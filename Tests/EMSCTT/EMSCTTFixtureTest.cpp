#include "EMSBH_trumpet_read.hpp"
#include <chrono>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

static void require(bool ok, const std::string &message)
{
    if (!ok)
        throw std::runtime_error(message);
}
static std::string read(const std::string &path)
{
    std::ifstream in(path, std::ios::binary);
    require(bool(in), "cannot read " + path);
    return std::string(std::istreambuf_iterator<char>(in), {});
}
static std::string replace(std::string text, const std::string &from,
                           const std::string &to)
{
    const auto pos = text.find(from);
    require(pos != std::string::npos, "mutation target missing: " + from);
    text.replace(pos, from.size(), to);
    return text;
}
static std::vector<std::string> split(const std::string &line,
                                      char delimiter = '\t')
{
    std::istringstream in(line);
    std::string word;
    std::vector<std::string> out;
    while (std::getline(in, word, delimiter))
        out.push_back(word);
    return out;
}
static std::vector<double> ccz4(const CCZ4CartoonVars::VarsWithGauge<double> &v)
{
    return {v.chi,      v.h[0][0], v.h[0][1],  v.h[1][1],  v.hww,   v.K,
            v.A[0][0],  v.A[0][1], v.A[1][1],  v.Aww,      v.Theta, v.Gamma[0],
            v.Gamma[1], v.lapse,   v.shift[0], v.shift[1], v.B[0],  v.B[1],
            v.phi,      v.Pi,      v.Ex,       v.Ey,       v.Ez,    v.Bx,
            v.By,       v.Bz,      v.Xi,       v.Lambda};
}
static void sha_controls()
{
    const std::pair<std::string, std::string> vectors[] = {
        {"",
         "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"},
        {"abc",
         "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"},
        {"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq",
         "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1"},
        {std::string(1, char(0xbd)),
         "68325720aabd7c82f30f554b313d0570c95accbb7dc4b5aae11204c08ffe732b"},
        {std::string(55, '\0'),
         "02779466cdec163811d078815c633f21901413081449002f24aa3e80f0b88ef7"},
        {std::string(56, '\0'),
         "d4817aa5497628e7c77e6b606107042bbba3130888c5f47a375e6179be789fbb"},
        {std::string(57, '\0'),
         "65a16cb7861335d5ace3c60718b5052e44660726da4cd13bb745381b235a1785"},
        {std::string(64, '\0'),
         "f5a5fd42d16a20302798ef6ed309979b43003d2320d9f0e8ea9831a92759fb4b"},
        {std::string(1000000, 'a'),
         "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0"}};
    for (const auto &v : vectors)
        require(SHA256::digest(v.first) == v.second,
                "SHA256 vector length " + std::to_string(v.first.size()));
    std::cout << "SHA256: 9 known-answer vectors passed\n";
}
static void parser_controls(const std::string &source, const std::string &path,
                            const EMSTrumpetSolution_read &profile,
                            const EMSCTTSolution_read::binding_t &binding)
{
    // Small zero-correction file retaining the authenticated layout, with 3x3
    // blocks. This is a parser control, never a solved binary fixture.
    std::istringstream in(source);
    std::ostringstream synthetic;
    std::string line;
    while (std::getline(in, line) && line.rfind("SHA256 ", 0) != 0)
    {
        if (line.rfind("PANEL ", 0) == 0)
        {
            auto w = split(line, ' ');
            w[10] = w[11] = "3";
            for (std::size_t i = 0; i < w.size(); ++i)
                synthetic << (i ? " " : "") << w[i];
            synthetic << '\n';
        }
        else if (line.rfind("FIELD ", 0) == 0)
        {
            auto w = split(line, ' ');
            synthetic << "FIELD " << w[1] << " 9\n";
            for (int i = 0; i < std::stoi(w[2]); ++i)
                require(bool(std::getline(in, line)), "source truncated");
            for (int i = 0; i < 9; ++i)
                synthetic << "0\n";
        }
        else
            synthetic << line << '\n';
    }
    const auto base = synthetic.str();
    auto write = [&](const std::string &s, bool rehash = true)
    {
        std::ofstream out(path, std::ios::binary);
        out << s;
        if (rehash)
            out << "SHA256 " << SHA256::digest(s) << '\n';
        require(bool(out), "cannot write synthetic control");
    };
    EMSCTTSolution_read reader;
    std::string reason;
    write(base);
    require(reader.check_file(path, profile, binding, reason),
            "valid synthetic: " + reason);
    require(reader.evaluate(0, 3, .2).logpsi == 0, "synthetic correction");
    const std::pair<std::string, std::pair<std::string, std::string>>
        mutants[] = {
            {"magic", {"EMSCTT 1", "EMSCTT 2"}},
            {"convention", {"ems-native-draft-v1", "unknown"}},
            {"schema", {"schema=EMSCTT/1", "schema=EMSCTT/2"}},
            {"missing key", {"# A_left=-0.0026627145861407193\n", ""}},
            {"duplicate key", {"# degree=28\n", "# degree=28\n# degree=28\n"}},
            {"unknown key", {"# degree=28\n", "# degree=28\n# x.unknown=1\n"}},
            {"block shape", {"FIELD 1 9", "FIELD 1 10"}},
            {"dimensions", {"sphere 3 3 5", "sphere 2 3 5"}},
            {"nonfinite coefficient", {"FIELD 1 9\n0", "FIELD 1 9\nNaN"}},
            {"nonfinite metadata", {"# degree=28", "# degree=Inf"}},
            {"integer control", {"# degree=28", "# degree=28.0"}},
            {"map topology", {"PANEL 1 log 1.0e-6", "PANEL 1 log 2.0e-6"}},
            {"END", {"END\n", "END \n"}},
            {"trailing data", {"END\n", "END\nextra\n"}},
            {"profile hash",
             {profile.get_source_sha256(), std::string(64, '0')}},
            {"parameter mismatch", {"# mass_left=1.0", "# mass_left=2.0"}},
            {"coupling mismatch", {"# ems_f2=-20", "# ems_f2=-19"}}};
    for (const auto &m : mutants)
    {
        write(replace(base, m.second.first, m.second.second));
        require(!reader.check_file(path, profile, binding, reason),
                "accepted " + m.first);
        require(!reason.empty(), "empty failure reason");
        std::cout << "parser " << m.first << ": rejected (" << reason << ")\n";
    }
    write(base.substr(0, base.find("FIELD 2")));
    require(!reader.check_file(path, profile, binding, reason),
            "accepted truncated payload");
    std::cout << "parser truncated payload: rejected (" << reason << ")\n";
    write(base + "SHA256 " + std::string(64, '0') + "\n", false);
    require(!reader.check_file(path, profile, binding, reason) &&
                reason == "companion checksum mismatch",
            "checksum control");
    write(base + "SHA256 " + SHA256::digest(base), false);
    require(!reader.check_file(path, profile, binding, reason),
            "accepted missing newline");
    // Every run binding is independently rejected, including the right hole.
    for (int field = 0; field < 10; ++field)
    {
        auto bad = binding;
        if (field < 2)
            bad.mass[field] += .01;
        else if (field < 4)
            bad.center[field - 2] += .01;
        else if (field < 6)
            bad.rapidity[field - 4] += .01;
        else
            bad.coupling[field - 6] += .01;
        write(base);
        require(!reader.check_file(path, profile, bad, reason),
                "accepted changed run binding");
    }
    write(base);
    require(reader.check_file(path, profile, binding, reason),
            "synthetic reload");
    for (double y : {0., 1e-22})
    {
        bool rejected = false;
        try
        {
            reader.evaluate(-16, y, 0);
        }
        catch (const std::domain_error &)
        {
            rejected = true;
        }
        require(rejected, "accepted unrepresented cylinder radius");
    }
    std::filesystem::remove(path);
    std::cout << "parser: 30 rejection controls passed; inner-domain controls "
                 "passed\n";
}
struct maximum_t
{
    double error = 0;
    std::string point = "all exact", field = "all exact";
};
static void parity(const std::string &dir, const EMSBH_trumpet_read &reader)
{
    int failures = 0;
    std::cout << "table,group,max_scaled_error,point,field\n";
    for (const std::string table : {"physical", "ccz4"})
    {
        std::ifstream in(dir + "/" + table + ".tsv");
        require(bool(in), "missing fixture " + table);
        std::string line;
        std::getline(in, line);
        const auto names = split(line);
        std::map<std::string, maximum_t> maxima;
        int count = 0;
        std::uint64_t fingerprint = 14695981039346656037ULL;
        while (std::getline(in, line))
        {
            auto row = split(line);
            require(row.size() == names.size(), "fixture row shape");
            const double x = std::stod(row[4]), y = std::stod(row[5]),
                         z = std::stod(row[6]);
            const int panel = std::stoi(row[3]);
            std::vector<double> actual;
            if (table == "physical")
            {
                const auto a = reader.compute_binary_ems_adm_vars(
                    x, y, 1., 32., .20273255, z, panel);
                for (const auto &tensor : {a.gamma, a.K})
                    for (int i = 0; i < 3; ++i)
                        for (int j = i; j < 3; ++j)
                            actual.push_back(tensor[i][j]);
                actual.push_back(a.phi);
                actual.push_back(a.Pi);
                for (const auto &v : {a.E, a.B, a.shift})
                    for (int i = 0; i < 3; ++i)
                        actual.push_back(v[i]);
                const auto q = reader.ctt_solution().evaluate(x, y, z, panel);
                actual.push_back(q.logpsi);
                for (int i = 0; i < 3; ++i)
                    for (int j = i; j < 3; ++j)
                        actual.push_back(q.C[i][j]);
            }
            else
                actual = ccz4(reader.compute_binary_bh_vars(x, y, 1., 32.,
                                                            .20273255, panel));
            require(actual.size() + 8 == row.size(), "fixture result shape");
            for (std::size_t k = 0; k < actual.size(); ++k)
            {
                std::uint64_t bits;
                std::memcpy(&bits, &actual[k], sizeof(bits));
                fingerprint = (fingerprint ^ bits) * 1099511628211ULL;
                const auto &name = names[k + 8];
                const double expected = std::stod(row[k + 8]);
                const double e =
                    std::abs(actual[k] - expected) / (1 + std::abs(expected));
                require(std::isfinite(e), "nonfinite parity");
                const std::string group = table == "physical"
                                              ? (k < 6     ? "gamma"
                                                 : k < 12  ? "Kij"
                                                 : k < 14  ? "scalar"
                                                 : k < 17  ? "E"
                                                 : k < 20  ? "B"
                                                 : k < 23  ? "shift"
                                                 : k == 23 ? "logpsi"
                                                           : "C")
                                              : (k < 5    ? "metric"
                                                 : k < 10 ? "curvature"
                                                 : k < 18 ? "gauge"
                                                 : k < 20 ? "scalar"
                                                 : k < 23 ? "E"
                                                 : k < 26 ? "B"
                                                          : "auxiliary");
                auto &m = maxima[group];
                if (e > m.error)
                    m = {e, row[0], name};
                if (e > 5e-13)
                {
                    ++failures;
                    std::cerr << table << " point=" << row[0] << " " << name
                              << " scaled_error=" << e << '\n';
                }
            }
            ++count;
        }
        require(count == 42, "fixture must contain 42 points");
        std::cout << table << " all-output fingerprint=" << std::hex
                  << fingerprint << std::dec << '\n';
        // Direct physical metric and the shared production CCZ4 conversion.
        if (table == "ccz4")
            require(fingerprint == 0xb3b32ee0d0a729dcULL,
                    "changed CTT CCZ4 bits");
        for (const auto &kv : maxima)
            std::cout << table << ',' << kv.first << ',' << kv.second.error
                      << ',' << kv.second.point << ',' << kv.second.field
                      << '\n';
    }
    require(failures == 0, "fixture parity failed");
}
int main(int argc, char **argv)
{
    try
    {
        require(argc == 5, "usage: EMSCTTFixtureTest profile companion "
                           "fixture_directory scratch_directory");
        const auto start = std::chrono::steady_clock::now();
        std::cout << std::setprecision(17);
        sha_controls();
        require(SHA256::digest(read(argv[1])) ==
                    "6f0820a576620f1f7230c701131af56a2312b130e31d5de2a752c32cef"
                    "e6d24f",
                "profile manifest hash");
        const auto source = read(argv[2]);
        require(SHA256::digest(source) == "ecba66b51c0bdb9a470afa50194dcf23bd0d"
                                          "7163a1af1e4b4b1fc5408b7e11b1",
                "companion manifest hash");
        require(SHA256::digest(read(std::string(argv[3]) + "/physical.tsv")) ==
                    "a0daf9a27596b96d29ddc302db193fe6919cb3d7b6ec5cee493b4a1538"
                    "5d1e90",
                "physical manifest hash");
        require(SHA256::digest(read(std::string(argv[3]) + "/ccz4.tsv")) ==
                    "bca759d4c2a6e7d8b3301a8b90f1b6b30de2c1431e0e1c254cb4146503"
                    "ef95c6",
                "ccz4 manifest hash");
        EMSBH_params_t params{};
        params.star_centre = {0, 0};
        params.bh_mass = 1;
        params.binary = params.boosted = true;
        params.separation = 32;
        params.rapidity = .20273255;
        params.data_path = argv[1];
        params.ctt_data_path = argv[2];
        CouplingFunction::params_t coupling{12.566370614359172, 0, 0, -20};
        EMSBH_trumpet_read reader(params, coupling, 1., 1., 0);
        reader.compute_1d_solution();
        const EMSCTTSolution_read::binding_t binding{
            {1, 1},
            {-16, 16},
            {.20273255, -.20273255},
            reader.m_1d_sol.get_coupling_parameters()};
        parser_controls(source, std::string(argv[4]) + "/synthetic.ctt",
                        reader.m_1d_sol, binding);
        // A translated point smaller than an ulp of the hole centre must still
        // carry its axial direction. Never form -16 + dx before evaluating.
        const auto a = reader.ctt_solution().evaluate(1e-18, 1e-18, 0, 0, -16);
        const auto b = reader.ctt_solution().evaluate(-1e-18, 1e-18, 0, 0, -16);
        require(a.C[0][1] != b.C[0][1], "lost hole-relative coordinate");
        parity(argv[3], reader);
        std::cout << "PASS seconds="
                  << std::chrono::duration<double>(
                         std::chrono::steady_clock::now() - start)
                         .count()
                  << '\n';
    }
    catch (const std::exception &e)
    {
        std::cerr << "FAIL: " << e.what() << '\n';
        return 1;
    }
}
