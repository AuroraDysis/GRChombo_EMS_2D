#include <signal.h>
#include <unistd.h>
int main()
{
    raise(SIGUSR1);
    sleep(1);
    return 0;
}
