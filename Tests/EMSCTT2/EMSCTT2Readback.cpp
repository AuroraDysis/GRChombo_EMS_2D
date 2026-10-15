#include "EMSBH_trumpet_read.hpp"
#include <iomanip>
#include <iostream>
int main(int argc, char **argv)
{
    try
    {
        if (argc != 3 && argc != 4)
            throw std::runtime_error(
                "usage: readback profile companion-or-plain");
        EMSTrumpetSolution_read profile;
        profile.main(argv[1]);
        if (argc == 4 && std::string(argv[3]) == "--check")
        {
            EMSCTTPulse_read parser;
            std::string reason;
            if (!parser.check_file(argv[2], profile, reason))
            {
                std::cout << reason << '\n';
                return 2;
            }
            return 0;
        }
        const auto c = profile.get_coupling_parameters();
        CouplingFunction::params_t coupling{c[0], c[1], c[2], c[3]};
        EMSBH_params_t p{};
        p.bh_mass = 1.;
        p.data_path = argv[1];
        p.ctt_data_path = std::string(argv[2]) == "plain" ? "" : argv[2];
        if (argc == 4)
        {
            const std::string mode = argv[3];
            if (mode == "--binary")
                p.binary = true;
            else if (mode == "--boosted")
                p.boosted = true;
            else if (mode == "--rapidity")
                p.rapidity = .1;
            else if (mode == "--mass")
                p.bh_mass = 2.;
            else if (mode == "--separation")
                p.separation = 1.;
            else if (mode == "--geometric-lapse")
                p.use_geometric_initial_lapse = true;
            else if (mode == "--maximal-lapse")
                p.use_maximal_initial_lapse = true;
            else if (mode == "--coupling")
                coupling.f2 += .1;
            else if (mode != "--domain" && mode != "--ccz4")
                throw std::runtime_error("unknown test mode");
        }
        EMSBH_trumpet_read reader(p, coupling, 1., 1., 0);
        reader.compute_1d_solution();
        if (argc == 4 && std::string(argv[3]) == "--domain")
        {
            EMSCTTPulse_read parser;
            std::string reason;
            if (!parser.check_file(argv[2], profile, reason))
                throw std::runtime_error(reason);
            for (double r : {0., 1e-22})
            {
                bool rejected = false;
                try
                {
                    parser.evaluate(r, 0., 0.);
                }
                catch (const std::runtime_error &)
                {
                    rejected = true;
                }
                if (!rejected)
                    throw std::runtime_error("unrepresented radius accepted");
            }
            for (int panel : {-1, 23})
            {
                bool rejected = false;
                try
                {
                    parser.evaluate(1., 0., 0., panel);
                }
                catch (const std::runtime_error &)
                {
                    rejected = true;
                }
                if (!rejected)
                    throw std::runtime_error("invalid panel accepted");
            }
            std::cout << "domain controls passed\n";
            return 0;
        }
        double x, y, z;
        int panel;
        std::cout << std::setprecision(17);
        while (std::cin >> x >> y >> z >> panel)
        {
            if (argc==4 && std::string(argv[3])=="--ccz4")
            {
                if (z!=0. || panel!=0) throw std::runtime_error("CCZ4 fixture requires cartoon plane and automatic panels");
                const auto v=reader.compute_single_bh_vars(x,y,1.,0.,0.);
                const double fields[]={v.chi,v.h[0][0],v.h[0][1],v.h[1][1],v.hww,v.K,
                    v.A[0][0],v.A[0][1],v.A[1][1],v.Aww,v.phi,v.Pi,v.Ex,v.Ey,v.Ez,
                    v.Bx,v.By,v.Bz,v.lapse,v.shift[0],v.shift[1]};
                for (double value:fields) std::cout<<value<<' ';
                std::cout<<'\n';continue;
            }
            const auto a =
                reader.compute_single_ems_adm_vars(x, y, z, 1., 0., 0., panel);
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    std::cout << a.gamma[i][j] << ' ';
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    std::cout << a.K[i][j] << ' ';
            std::cout << a.phi << ' ' << a.Pi << ' ';
            for (int i = 0; i < 3; ++i)
                std::cout << a.E[i] << ' ';
            for (int i = 0; i < 3; ++i)
                std::cout << a.B[i] << ' ';
            for (int i = 0; i < 3; ++i)
                std::cout << a.shift[i] << ' ';
            std::cout << '\n';
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
