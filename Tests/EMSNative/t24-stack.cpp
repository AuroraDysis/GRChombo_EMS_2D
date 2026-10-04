// Tests-only linked observer. Absent T24_STACK_PREFIX: no action at all.
#include <execinfo.h>
#include <dlfcn.h>
#include <fcntl.h>
#include <omp.h>
#include <pthread.h>
#include <signal.h>
#include <unistd.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

namespace
{
struct Thread { pthread_t id; int fd; };
std::vector<Thread> threads;

void dump(int)
{
    const auto self = pthread_self();
    for (const auto &t : threads)
        if (pthread_equal(self, t.id))
        {
            const char marker[] = "\nT24_SIGNAL_STACK_BEGIN\n";
            write(t.fd, marker, sizeof(marker)-1);
            void *frames[128];
            const int n = backtrace(frames, 128);
            backtrace_symbols_fd(frames, n, t.fd);
            const char end[] = "T24_SIGNAL_STACK_END\n";
            write(t.fd, end, sizeof(end)-1);
            break;
        }
}

void request(int)
{
    for (const auto &t : threads) pthread_kill(t.id, SIGUSR2);
}

// Diagnostic fallback for a sandbox that refuses debugserver launch. It uses
// pre-opened per-thread files and pre-warms the native unwinder. Signal-time
// unwinding remains best effort, not a substitute for a ptrace debugger.
struct Install
{
    Install()
    {
        const char *prefix = std::getenv("T24_STACK_PREFIX");
        if (!prefix) return;
        threads.resize(omp_get_max_threads());
        int failure = 0;
        #pragma omp parallel reduction(|:failure)
        {
            const int i = omp_get_thread_num();
            char path[4096];
            std::snprintf(path, sizeof(path), "%s.thread-%d.txt", prefix, i);
            const int fd = open(path, O_WRONLY|O_CREAT|O_APPEND, 0600);
            threads[i] = {pthread_self(), fd};
            if (fd < 0) failure = 1;
            void *frames[128];
            const int n = backtrace(frames, 128);
            if (fd >= 0)
            {
                Dl_info image{};
                dladdr(reinterpret_cast<void *>(&dump), &image);
                dprintf(fd, "T24_STACK_READY thread=%d image=%s load_address=%p\n",
                        i, image.dli_fname, image.dli_fbase);
                backtrace_symbols_fd(frames, n, fd);
            }
        }
        if (failure)
        {
            std::fprintf(stderr, "T24 stack observer could not open its files\n");
            std::exit(70);
        }
        struct sigaction action{};
        sigemptyset(&action.sa_mask);
        action.sa_flags = SA_RESTART;
        action.sa_handler = dump;
        sigaction(SIGUSR2, &action, nullptr);
        action.sa_handler = request;
        sigaction(SIGUSR1, &action, nullptr);
    }
} install;
}
