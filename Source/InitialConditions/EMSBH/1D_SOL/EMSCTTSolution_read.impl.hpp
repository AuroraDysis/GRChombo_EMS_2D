#if !defined(EMSCTTSOLUTION_READ_HPP_)
#error "Include EMSCTTSolution_read.hpp instead"
#endif
#ifndef EMSCTTSOLUTION_READ_IMPL_HPP_
#define EMSCTTSOLUTION_READ_IMPL_HPP_

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <tuple>

inline bool EMSCTTSolution_read::number(const std::string &text, double &value)
{
    char *end = nullptr;
    value = std::strtod(text.c_str(), &end);
    if (end == text.c_str() || !std::isfinite(value))
        return false;
    while (*end == ' ' || *end == '\t')
        ++end;
    return !*end && text.find('\0') == std::string::npos;
}
inline bool EMSCTTSolution_read::integer(const std::string &text, int &value,
                                         int lo, int hi)
{
    std::istringstream in(text);
    std::string token, extra;
    if (!(in >> token) || (in >> extra))
        return false;
    if (token[0] == '+')
        token.erase(0, 1);
    if (token.empty() || !std::all_of(token.begin(), token.end(), [](char c)
                                      { return c >= '0' && c <= '9'; }))
        return false;
    token.erase(0, token.find_first_not_of('0'));
    for (int i = lo; i <= hi; ++i)
        if (token == std::to_string(i))
        {
            value = i;
            return true;
        }
    return false;
}
inline bool EMSCTTSolution_read::hash_string(const std::string &text)
{
    return text.size() == 64 &&
           std::all_of(
               text.begin(), text.end(), [](char c)
               { return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f'); });
}
inline bool EMSCTTSolution_read::check_file(
    const std::string &a_path, const EMSTrumpetSolution_read &a_profile,
    const binding_t &a_binding, std::string &a_reason)
{
    m_panels.clear();
    m_values.clear();
    m_source_sha256.clear();
    a_reason.clear();
    const auto fail = [&](const std::string &message)
    {
        a_reason = message;
        return false;
    };
    std::ifstream file(a_path, std::ios::binary);
    if (!file)
        return fail("cannot open companion file");
    const std::string bytes((std::istreambuf_iterator<char>(file)), {});
    if (file.bad())
        return fail("cannot read companion file");
    if (bytes.empty() || bytes.back() != '\n')
        return fail("missing final newline");
    const auto trailer = bytes.rfind('\n', bytes.size() - 2);
    if (trailer == std::string::npos)
        return fail("missing checksum");
    const auto checksum = bytes.substr(trailer + 1, bytes.size() - trailer - 2);
    if (checksum.size() != 71 || checksum.substr(0, 7) != "SHA256 " ||
        !hash_string(checksum.substr(7)))
        return fail("malformed checksum");
    if (SHA256::digest(std::string_view(bytes.data(), trailer + 1)) !=
        checksum.substr(7))
        return fail("companion checksum mismatch");
    std::istringstream in(bytes.substr(0, trailer + 1));
    std::string line;
    if (!std::getline(in, line) || line != "EMSCTT 1")
        return fail("unknown companion magic");
    const std::map<std::string, std::string> required_text = {
        {"schema", "EMSCTT/1"},
        {"convention", "ems-native-draft-v1"},
        {"construction", "nonlinear-cylinder-dtn-B0"},
        {"end_branch", "B_A=0"},
        {"coefficients", "Chebyshev-Ti(t)-Tj(v)-full-a0-radial-fast"},
        {"fields", "logpsi,Cxx,Cxr_over_s,Cww,Crr_minus_Cww_over_s2"},
        {"tensor", "covariant-C=det(gbar)^(-1/3)*LW-lowered"},
        {"end_continuation", "resolved-log-panels-no-padding"}};
    const std::vector<std::string> numeric_keys = {
        "mass_left",        "mass_right",      "center_left",    "center_right",
        "rapidity_left",    "rapidity_right",  "ems_alpha",      "ems_f0",
        "ems_f1",           "ems_f2",          "psi_C_left",     "psi_C_right",
        "A_left",           "A_right",         "degree",         "cutoff",
        "newton_tolerance", "newton_residual", "tail_tolerance", "tail_degree",
        "tail_L",           "end_min_radius",  "panel_count"};
    const std::vector<std::string> hash_keys = {
        "profile_sha256", "solver_sha256", "end_solver_sha256",
        "checkpoint_sha256"};
    std::map<std::string, std::string> meta;
    while (std::getline(in, line) && line.rfind("# ", 0) == 0)
    {
        const auto pos = line.find('=', 2);
        if (pos == std::string::npos)
            return fail("malformed metadata");
        const auto key = line.substr(2, pos - 2), value = line.substr(pos + 1);
        if (!required_text.count(key) &&
            std::find(numeric_keys.begin(), numeric_keys.end(), key) ==
                numeric_keys.end() &&
            std::find(hash_keys.begin(), hash_keys.end(), key) ==
                hash_keys.end())
            return fail("unknown key " + key);
        if (!meta.emplace(key, value).second)
            return fail("duplicate key " + key);
    }
    if (meta.size() !=
        required_text.size() + numeric_keys.size() + hash_keys.size())
        return fail("missing required metadata");
    for (const auto &kv : required_text)
        if (meta.at(kv.first) != kv.second)
            return fail("unsupported " + kv.first);
    std::map<std::string, double> values;
    for (const auto &key : numeric_keys)
        if (!number(meta.at(key), values[key]))
            return fail("nonfinite " + key);
    for (const auto &key : hash_keys)
        if (!hash_string(meta.at(key)))
            return fail("malformed hash " + key);
    if (a_profile.get_source_sha256().empty() ||
        meta.at("profile_sha256") != a_profile.get_source_sha256())
        return fail("profile hash mismatch");
    const char *coupling_keys[] = {"ems_alpha", "ems_f0", "ems_f1", "ems_f2"};
    for (int i = 0; i < 4; ++i)
        if (values.at(coupling_keys[i]) != a_binding.coupling[i] ||
            values.at(coupling_keys[i]) !=
                a_profile.get_coupling_parameters()[i])
            return fail("coupling parameter mismatch");
    for (int side = 0; side < 2; ++side)
    {
        const std::string suffix = side == 0 ? "left" : "right";
        if (values.at("mass_" + suffix) != a_binding.mass[side] ||
            values.at("center_" + suffix) != a_binding.center[side] ||
            values.at("rapidity_" + suffix) != a_binding.rapidity[side])
            return fail("binary parameter mismatch");
    }
    const auto v = [&](const char *key) { return values.at(key); };
    if (v("mass_left") != 1 || v("mass_right") != 1)
        return fail("unsupported masses");
    if (v("center_left") != -v("center_right") || !(v("center_right") > 4))
        return fail("invalid centres");
    if (v("rapidity_left") != -v("rapidity_right"))
        return fail("invalid rapidities");
    for (const auto key :
         {"psi_C_left", "psi_C_right", "newton_tolerance", "tail_tolerance"})
        if (!(v(key) > 0))
            return fail("invalid controls");
    if (!(0 <= v("newton_residual") &&
          v("newton_residual") <= v("newton_tolerance")))
        return fail("unconverged source");
    if (!(0 < v("end_min_radius") && v("end_min_radius") < v("cutoff") &&
          v("cutoff") < .01))
        return fail("invalid end range");
    int count, unused;
    if (!integer(meta.at("panel_count"), count, 1, 128) ||
        !integer(meta.at("degree"), unused, 2, 128) ||
        !integer(meta.at("tail_degree"), unused, 2, 128) ||
        !integer(meta.at("tail_L"), unused, 1, 256))
        return fail("invalid integer controls");
    std::vector<panel_t> panels;
    for (int index = 1; index <= count; ++index)
    {
        panel_t p;
        std::istringstream row(line);
        std::vector<std::string> w{std::istream_iterator<std::string>(row), {}};
        int id, nf;
        if (w.size() != 13 || w[0] != "PANEL" ||
            !integer(w[1], id, index, index))
            return fail("malformed panel");
        p.kind = w[2];
        p.boundary = w[9];
        if ((p.kind != "log" && p.kind != "bridge" && p.kind != "inverse" &&
             p.kind != "end") ||
            (p.boundary != "sphere" && p.boundary != "plane"))
            return fail("unknown map");
        double *maps[] = {&p.lo,      &p.hi,    &p.center,
                          &p.stretch, &p.mu_lo, &p.mu_hi};
        for (int j = 0; j < 6; ++j)
            if (!number(w[j + 3], *maps[j]))
                return fail("nonfinite map");
        if (!(p.hi > p.lo && p.stretch > 0 && -1 <= p.mu_lo &&
              p.mu_lo < p.mu_hi && p.mu_hi <= 1))
            return fail("invalid map");
        if (p.kind == "inverse"
                ? !(p.lo == 0 && p.center == 0 && p.stretch == 1)
                : !(p.lo > 0))
            return fail("invalid radial map");
        if (p.kind == "end" && !(p.lo >= v("end_min_radius") * (1 - 1e-13) &&
                                 p.hi <= v("cutoff") * (1 + 1e-13)))
            return fail("invalid end map");
        if (!integer(w[10], p.nr, 3, 193) || !integer(w[11], p.na, 3, 193) ||
            !integer(w[12], nf, 5, 5))
            return fail("invalid block dimensions");
        for (int f = 0; f < 5; ++f)
        {
            if (!std::getline(in, line) ||
                line != "FIELD " + std::to_string(f + 1) + " " +
                            std::to_string(p.nr * p.na))
                return fail("malformed field block");
            p.coeff[f].resize(p.nr * p.na);
            for (auto &a : p.coeff[f])
            {
                if (!std::getline(in, line))
                    return fail("truncated coefficient block");
                if (!number(line, a))
                    return fail("nonfinite coefficient");
            }
        }
        panels.push_back(std::move(p));
        if (!std::getline(in, line))
            return fail("truncated companion");
    }
    if (line != "END" || std::getline(in, line))
        return fail("invalid companion ending");
    // The same v1 finite topology and resolved-end coverage as the Julia
    // reader.
    const double cut = v("cutoff"), center = v("center_right"),
                 stretch = std::cosh(v("rapidity_left")), outer = 4 * center;
    if (cut != 1e-4 && cut != 1e-6 && cut != 1e-8)
        return fail("unsupported cutoff layout");
    if (!std::isfinite(outer) || !std::isfinite(stretch))
        return fail("nonfinite binary map");
    std::vector<double> radii;
    if (cut == 1e-8)
        radii.push_back(1e-8);
    if (cut <= 1e-6)
        radii.push_back(1e-6);
    for (double r : {1e-4, .01, 1., center / 2})
        radii.push_back(r);
    using map_t = std::tuple<std::string, double, double, double, double,
                             double, double, std::string>;
    std::vector<map_t> expected, actual;
    for (double c : {-center, center})
    {
        const double mu =
            -stretch * c /
            std::sqrt((stretch * c) * (stretch * c) + outer * outer);
        const double edges[] = {-1, mu, 1};
        for (int sector = 0; sector < 2; ++sector)
        {
            for (std::size_t j = 0; j + 1 < radii.size(); ++j)
                expected.emplace_back("log", radii[j], radii[j + 1], c, stretch,
                                      edges[sector], edges[sector + 1],
                                      "sphere");
            expected.emplace_back("bridge", center / 2, outer, c, stretch,
                                  edges[sector], edges[sector + 1],
                                  (c < 0) == (sector == 0) ? "sphere"
                                                           : "plane");
        }
    }
    for (int sector = 0; sector < 2; ++sector)
        expected.emplace_back("inverse", 0, 1 / outer, 0, 1, sector - 1, sector,
                              "sphere");
    for (const auto &p : panels)
        if (p.kind != "end")
            actual.emplace_back(p.kind, p.lo, p.hi, p.center, p.stretch,
                                p.mu_lo, p.mu_hi, p.boundary);
        else if ((p.center != -center && p.center != center) ||
                 p.boundary != "sphere")
            return fail("unbound end panel");
    if (actual != expected)
        return fail("finite maps do not match binary binding");
    for (double c : {-center, center})
    {
        std::vector<const panel_t *> ends;
        for (const auto &p : panels)
            if (p.kind == "end" && p.center == c)
                ends.push_back(&p);
        if (ends.empty())
            return fail("missing resolved end");
        std::sort(ends.begin(), ends.end(),
                  [](const panel_t *a, const panel_t *b)
                  { return a->lo < b->lo; });
        for (const auto p : ends)
            if (p->stretch != stretch || p->mu_lo != -1 || p->mu_hi != 1)
                return fail("invalid end coverage");
        if (std::abs(ends.front()->lo / v("end_min_radius") - 1) > 1e-13 ||
            std::abs(ends.back()->hi / cut - 1) > 1e-13)
            return fail("incomplete end coverage");
        for (std::size_t j = 0; j + 1 < ends.size(); ++j)
            if (ends[j]->hi != ends[j + 1]->lo)
                return fail("end gap or overlap");
    }
    m_panels = std::move(panels);
    m_values = std::move(values);
    m_source_sha256 = SHA256::digest(bytes);
    return true;
}
inline bool EMSCTTSolution_read::matches(double mass, double separation,
                                         double rapidity) const
{
    return !m_panels.empty() && mass == m_values.at("mass_left") &&
           -separation / 2 == m_values.at("center_left") &&
           separation / 2 == m_values.at("center_right") &&
           rapidity == m_values.at("rapidity_left");
}
inline double EMSCTTSolution_read::clenshaw(const double *a, int n, double x)
{
    double b = 0, d = 0;
    for (int k = n - 1; k >= 1; --k)
    {
        const double next = 2 * x * b - d + a[k];
        d = b;
        b = next;
    }
    return a[0] + x * b - d;
}
inline EMSCTTSolution_read::correction_t
EMSCTTSolution_read::evaluate(double x, double y, double z, int panel,
                              double origin) const
{
    if (m_panels.empty() || panel < 0 || panel > int(m_panels.size()))
        throw std::invalid_argument(
            "EMSCTT: unloaded companion or invalid panel");
    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z) ||
        !std::isfinite(origin))
        throw std::domain_error("EMSCTT: nonfinite coordinates");
    const double eps = 8 * std::numeric_limits<double>::epsilon();
    // Reject the limiting cylinder and radii below the represented inner end;
    // no constant padding or clamping, as in the trumpet reader's inner rule.
    for (const char *key : {"center_left", "center_right"})
    {
        const double X = std::cosh(m_values.at("rapidity_left")) *
                         (x + (origin - m_values.at(key)));
        if (std::sqrt(X * X + y * y + z * z) < m_values.at("end_min_radius"))
            throw std::domain_error("EMSCTT: below minimum represented radius");
    }
    for (std::size_t k = 0; k < m_panels.size(); ++k)
    {
        if (panel && int(k + 1) != panel)
            continue;
        const auto &p = m_panels[k];
        const double X = p.stretch * (x + (origin - p.center));
        const double r = std::sqrt(X * X + y * y + z * z);
        if (!(r > 0))
            continue;
        const double mu = X / r;
        if (!panel && (mu < p.mu_lo - eps || mu > p.mu_hi + eps))
            continue;
        double hi = p.hi;
        if (p.kind == "bridge")
        {
            if (p.boundary == "plane")
                hi = -p.center * p.stretch / mu;
            else
            {
                const double a = 1 - mu * mu +
                                 mu * mu / (p.stretch * p.stretch),
                             b = p.center * mu / p.stretch;
                hi = (-b + std::sqrt(b * b +
                                     a * (p.hi * p.hi - p.center * p.center))) /
                     a;
            }
        }
        const double t =
            p.kind == "inverse"
                ? (2 / r - p.lo - p.hi) / (p.hi - p.lo)
                : (2 * std::log(r) - std::log(p.lo) - std::log(hi)) /
                      (std::log(hi) - std::log(p.lo));
        if (!std::isfinite(t) || (!panel && (t < -1 - eps || t > 1 + eps)))
            continue;
        const double v = (2 * mu - p.mu_lo - p.mu_hi) / (p.mu_hi - p.mu_lo);
        std::array<double, 5> a;
        std::array<double, 193> angular;
        for (int f = 0; f < 5; ++f)
        {
            for (int j = 0; j < p.na; ++j)
                angular[j] = clenshaw(p.coeff[f].data() + j * p.nr, p.nr, t);
            a[f] = clenshaw(angular.data(), p.na, v);
            if (!std::isfinite(a[f]))
                throw std::domain_error("EMSCTT: nonfinite correction");
        }
        correction_t q{a[0]};
        const double n[] = {y / r, z / r};
        q.C[0][0] = a[1];
        for (int i = 0; i < 2; ++i)
        {
            q.C[0][i + 1] = q.C[i + 1][0] = a[2] * n[i];
            for (int j = 0; j < 2; ++j)
                q.C[i + 1][j + 1] = (i == j ? a[3] : 0) + a[4] * n[i] * n[j];
        }
        return q;
    }
    throw std::domain_error("EMSCTT: outside represented companion domain");
}
#endif
