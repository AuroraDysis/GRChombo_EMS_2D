#include "SetupFunctions.hpp"
#include "RHUnion.hpp"

// Run in an empty directory. Also works collectively in an MPI build.
int main(int argc, char **argv)
{
    mainSetup(argc, argv);
    const bool invalid = argc > 1 && std::string(argv[1]) == "invalid";
    if (procID() == 0)
    {
        std::ofstream("rh_surf_0.dat")
            << "# saved geometry\n1 0 2.25 found\n2 0 2.5 close\n3 0 99 found\n";
        std::ofstream("rh_f0.dat")
            << "# shape\n1 1 1 1 1\n2 0.8 0.9 1.1 "
            << (invalid ? "-1" : "1.2") << "\n3 9 9 9 9\n";
        std::ofstream("rh_surf_2.dat") << "1 2 4.5 dead\n";
        std::ofstream("rh_f2.dat") << "1 1 1 1 1\n";
    }
    RHUnion rh;
    rh.setup(3, {5, 6, 7}, {10, 11, 12}, {4, 4, 4}, {0, 0, 0},
             {1, 1, 1}, {0, 0, 0}, {1, 1, 1}, {0, 70, 0}, 2.);
    auto require = [](bool ok) { if (!ok) MayDay::Error("RH restart regression"); };
    const auto &s = rh.m_surfaces[0];
    require(s.m_centre[0] == 2.5 && s.m_f[s.m_NG] == .8 &&
            s.m_f[s.m_NG + 3] == 1.2 && !s.m_dead);
    require(s.m_f[1] == .8 && s.m_f[0] == .9 &&
            s.m_f[6] == 1.2 && s.m_f[7] == 1.1);
    require(rh.m_surfaces[1].m_centre[0] == 11 &&
            rh.m_surfaces[1].m_f[2] == 6 &&
            rh.m_surfaces[1].m_state == RHSurf::SolverState::DORMANT);
    require(rh.m_surfaces[2].m_dead && rh.m_surfaces[2].m_centre[0] == 4.5);
    if (procID() == 0)
        for (const auto *name : {"rh_surf_0.dat", "rh_f0.dat"})
        {
            std::ifstream file(name);
            std::string line;
            int rows = 0;
            while (std::getline(file, line))
            {
                std::istringstream row(line);
                double t;
                if (row >> t) { require(t <= 2.); ++rows; }
            }
            require(rows == 2);
        }
    pout() << "RH_RESTART_PASS\n";
    mainFinalize();
}
