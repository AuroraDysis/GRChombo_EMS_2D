/* EMSCTT/2: file-only, unboosted mass-one initial data. */
#ifndef EMSCTTPULSE_READ_HPP_
#define EMSCTTPULSE_READ_HPP_
#include "EMSTrumpetSolution_read.hpp"
#include "SHA256.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iterator>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

class EMSCTTPulse_read
{
    struct panel_t
    {
        std::string kind;
        double lo, hi;
        int nr, na;
        std::array<std::vector<double>, 3> c;
    };
    std::vector<panel_t> m_panels;
    std::map<std::string, double> m_values;
    std::vector<double> m_phase;
    double m_loglo, m_loghi, m_center, m_left, m_right;
    static double number(const std::string &s)
    {
        size_t n;
        const double v = std::stod(s, &n);
        if (n != s.size() || !std::isfinite(v) ||
            s.find_first_not_of("+-0123456789.eE") != s.npos)
            throw std::runtime_error("invalid finite number");
        return v;
    }
    static int integer(const std::string &s, int lo, int hi)
    {
        if (s.empty() || s.find_first_not_of("0123456789") != s.npos)
            throw std::runtime_error("invalid integer");
        const double v = number(s);
        if (v < lo || v > hi)
            throw std::runtime_error("integer out of range");
        return static_cast<int>(v);
    }
    static bool hash(const std::string &s)
    {
        return s.size() == 64 &&
               s.find_first_not_of("0123456789abcdef") == s.npos;
    }
    static std::vector<std::string> tokens(const std::string &s)
    {
        std::istringstream in(s);
        return {std::istream_iterator<std::string>(in), {}};
    }
    // Full-a0 Chebyshev values and analytic first derivatives, including
    // endpoints.
    static void basis(double t, std::vector<double> &v, std::vector<double> &d)
    {
        v[0] = 1.;
        d[0] = 0.;
        if (v.size() > 1)
        {
            v[1] = t;
            d[1] = 1.;
        }
        for (size_t k = 2; k < v.size(); ++k)
        {
            v[k] = 2 * t * v[k - 1] - v[k - 2];
            d[k] = 2 * v[k - 1] + 2 * t * d[k - 1] - d[k - 2];
        }
    }
    std::array<double, 2> phase(double r) const
    {
        std::vector<double> v(m_phase.size()), d(v.size());
        basis((2 * std::log(r) - m_loglo - m_loghi) / (m_loghi - m_loglo), v,
              d);
        std::array<double, 2> out{};
        for (size_t k = 0; k < v.size(); ++k)
        {
            out[0] += m_phase[k] * v[k];
            out[1] += m_phase[k] * d[k];
        }
        out[1] *= 2 / (r * (m_loghi - m_loglo));
        return out;
    }
    static void transform(panel_t &p)
    {
        const double pi = std::acos(-1.);
        for (auto &c : p.c)
        {
            std::vector<double> tmp(c.size()), out(c.size());
            for (int j = 0; j < p.na; ++j)
                for (int k = 0; k < p.nr; ++k)
                {
                    double sum = 0.;
                    for (int i = 0; i < p.nr; ++i)
                        sum += c[i + p.nr * j] *
                               ((i == 0 || i == p.nr - 1) ? .5 : 1.) *
                               std::cos(k * pi * (p.nr - 1 - i) / (p.nr - 1));
                    tmp[k + p.nr * j] = sum * 2 / (p.nr - 1) *
                                        ((k == 0 || k == p.nr - 1) ? .5 : 1.);
                }
            for (int k = 0; k < p.na; ++k)
                for (int i = 0; i < p.nr; ++i)
                {
                    double sum = 0.;
                    for (int j = 0; j < p.na; ++j)
                        sum += tmp[i + p.nr * j] *
                               std::cos(k * pi * (2 * (p.na - 1 - j) + 1) /
                                        (2 * p.na));
                    out[i + p.nr * k] = sum * 2 / p.na * (k == 0 ? .5 : 1.);
                }
            c.swap(out);
        }
    }

  public:
    struct correction_t
    {
        double logpsi;
        std::array<std::array<double, 3>, 3> strain{};
    };
    bool check_file(const std::string &path,
                    const EMSTrumpetSolution_read &profile, std::string &reason)
    {
        m_panels.clear();
        m_values.clear();
        m_phase.clear();
        try
        {
            const auto require = [](bool ok, const char *why)
            {
                if (!ok)
                    throw std::runtime_error(why);
            };
            std::ifstream file(path, std::ios::binary);
            require(bool(file), "cannot open companion");
            const std::string bytes((std::istreambuf_iterator<char>(file)), {});
            require(!file.bad() && bytes.size() > 2 && bytes.back() == '\n',
                    "read/final newline");
            const auto end = bytes.rfind('\n', bytes.size() - 2);
            require(end != bytes.npos, "missing checksum");
            const auto trailer = bytes.substr(end + 1, bytes.size() - end - 2);
            require(trailer.size() == 71 && trailer.substr(0, 7) == "SHA256 " &&
                        hash(trailer.substr(7)),
                    "checksum syntax");
            require(SHA256::digest(std::string_view(bytes.data(), end + 1)) ==
                        trailer.substr(7),
                    "checksum mismatch");
            std::istringstream in(bytes.substr(0, end + 1));
            std::string line;
            const auto next = [&]()
            {
                require(bool(std::getline(in, line)), "truncated companion");
                return line;
            };
            require(next() == "EMSCTT 2", "unknown magic");
            const std::map<std::string, std::string> text = {
                {"schema", "EMSCTT/2-single-pulse"},
                {"geometry", "one-unboosted-trumpet"},
                {"convention", "ems-native-draft-v1"},
                {"construction", "nonlinear-cylinder-dtn-B0"},
                {"free_data", "trumpet-gamma-K-tau0-fixed-D0-B0-scalar-pulse"},
                {"pulse", "R-dphi=A*M*exp(-z^2/2)*scan-taper-Y20;xi=rstar+h"},
                {"momentum", "Pi=-n.dphi;dPi=-s.dphi-(s.dlnR+Ktheta).dphi"},
                {"fields", "nodal-logpsi,a,b;W=r/(1+r)*(a,y*b/r,z*b/r)"}};
            const std::set<std::string> numeric = {"mass",
                                                   "amplitude",
                                                   "radius_M",
                                                   "sigma_M",
                                                   "cutoff_M",
                                                   "end_min_radius_M",
                                                   "newton_tolerance",
                                                   "newton_residual",
                                                   "tail_degree",
                                                   "tail_L",
                                                   "panel_count"};
            const std::set<std::string> hashes = {
                "profile_sha256", "solver_sha256", "end_solver_sha256",
                "pulse_sha256", "checkpoint_sha256"};
            std::map<std::string, std::string> meta;
            next();
            while (line.rfind("# ", 0) == 0)
            {
                auto eq = line.find('=', 2);
                require(eq != line.npos, "metadata syntax");
                auto key = line.substr(2, eq - 2), value = line.substr(eq + 1);
                require(text.count(key) || numeric.count(key) ||
                            hashes.count(key),
                        "unknown metadata");
                require(meta.emplace(key, value).second, "duplicate metadata");
                next();
            }
            require(meta.size() == text.size() + numeric.size() + hashes.size(),
                    "missing metadata");
            for (const auto &kv : text)
                require(meta.at(kv.first) == kv.second,
                        "unsupported reconstruction");
            for (const auto &key : numeric)
                m_values[key] = number(meta.at(key));
            for (const auto &key : hashes)
                require(hash(meta.at(key)), "malformed provenance hash");
            require(meta.at("profile_sha256") == profile.get_source_sha256(),
                    "profile hash mismatch");
            const auto v = [&](const char *key) { return m_values.at(key); };
            require(v("mass") == 1 && v("radius_M") > 0 && v("sigma_M") > 0 &&
                        v("end_min_radius_M") > 0 &&
                        v("end_min_radius_M") < v("cutoff_M") &&
                        v("newton_tolerance") > 0 &&
                        v("newton_residual") >= 0 &&
                        v("newton_residual") <= v("newton_tolerance"),
                    "physical/solve controls");
            const int count = integer(meta.at("panel_count"), 2, 256);
            integer(meta.at("tail_degree"), 2, 192);
            integer(meta.at("tail_L"), 0, 192);
            auto row = tokens(line);
            require(row.size() == 7 && row[0] == "PHASE", "phase syntax");
            m_loglo = number(row[1]);
            m_loghi = number(row[2]);
            m_center = number(row[3]);
            m_left = number(row[4]);
            m_right = number(row[5]);
            int n = integer(row[6], 18, 194);
            require(m_loglo < m_loghi && std::exp(m_loglo) < m_left &&
                        m_left < m_center && m_center < m_right &&
                        m_right < std::exp(m_loghi),
                    "phase domain");
            for (int k = 0; k < n; ++k)
                m_phase.push_back(number(next()));
            require(std::abs(phase(m_center)[0]) < 1e-10 &&
                        std::abs(phase(m_left)[0] + 8 * v("sigma_M")) < 1e-9 &&
                        std::abs(phase(m_right)[0] - 8 * v("sigma_M")) < 1e-9,
                    "phase normalization");
            require(
                std::abs(m_center / profile.compute_radial_vars(m_center).X -
                         v("radius_M")) < 1e-10,
                "launch radius");
            for (int k = 0; k < count; ++k)
            {
                row = tokens(next());
                require(row.size() == 6 && row[0] == "PANEL", "panel syntax");
                panel_t p;
                p.kind = row[1];
                p.lo = number(row[2]);
                p.hi = number(row[3]);
                p.nr = integer(row[4], 3, 193);
                p.na = integer(row[5], 3, 193);
                require(p.hi > p.lo && ((p.kind == "log" && p.lo > 0) ||
                                        (p.kind == "inverse" && p.lo == 0)),
                        "panel map");
                require(next() == "VALUES " + std::to_string(3 * p.nr * p.na),
                        "values count");
                for (auto &c : p.c)
                    for (int j = 0; j < p.nr * p.na; ++j)
                        c.push_back(number(next()));
                m_panels.push_back(std::move(p));
            }
            require(next() == "END" && !std::getline(in, line),
                    "trailing/missing data");
            require(m_panels.front().lo == v("end_min_radius_M") &&
                        m_panels.back().kind == "inverse",
                    "end topology");
            bool cutoff = false;
            for (size_t j = 0; j + 1 < m_panels.size(); ++j)
            {
                auto &p = m_panels[j];
                require(p.kind == "log", "finite topology");
                cutoff |= p.lo == v("cutoff_M");
                if (j + 2 < m_panels.size())
                    require(p.hi == m_panels[j + 1].lo, "radial gap/overlap");
                else
                    require(std::abs(p.hi * m_panels.back().hi - 1) <=
                                16 * std::numeric_limits<double>::epsilon(),
                            "inverse gap/overlap");
            }
            require(cutoff, "missing cutoff face");
            for (auto &p : m_panels)
                transform(p);
            reason.clear();
            return true;
        }
        catch (const std::exception &e)
        {
            reason = e.what();
            m_panels.clear();
            m_values.clear();
            m_phase.clear();
            return false;
        }
    }
    correction_t evaluate(double x, double y, double z, int panel = 0) const
    {
        const double r = std::sqrt(x * x + y * y + z * z), mu = x / r;
        if (!std::isfinite(r) || !(r >= m_values.at("end_min_radius_M")))
            throw std::runtime_error("EMSCTT 2: outside represented end");
        int selected = panel - 1;
        if (!panel)
            for (size_t k = 0; k < m_panels.size(); ++k)
            {
                const auto &p = m_panels[k];
                const double t =
                    p.kind == "log"
                        ? (2 * std::log(r) - std::log(p.lo) - std::log(p.hi)) /
                              (std::log(p.hi) - std::log(p.lo))
                        : (2 / r - p.lo - p.hi) / (p.hi - p.lo);
                if (std::abs(t) <=
                    1 + 32 * std::numeric_limits<double>::epsilon())
                {
                    selected = k;
                    break;
                }
            }
        if (selected < 0 || selected >= int(m_panels.size()))
            throw std::runtime_error("EMSCTT 2: invalid panel/domain");
        const auto &p = m_panels[selected];
        const double t =
            p.kind == "log"
                ? (2 * std::log(r) - std::log(p.lo) - std::log(p.hi)) /
                      (std::log(p.hi) - std::log(p.lo))
                : (2 / r - p.lo - p.hi) / (p.hi - p.lo);
        const double tr = p.kind == "log"
                              ? 2 / (r * (std::log(p.hi) - std::log(p.lo)))
                              : -2 / (r * r * (p.hi - p.lo));
        std::vector<double> vr(p.nr), dr(p.nr), va(p.na), da(p.na);
        basis(t, vr, dr);
        basis(mu, va, da);
        double q[3][3] = {};
        for (int f = 0; f < 3; ++f)
            for (int j = 0; j < p.na; ++j)
                for (int i = 0; i < p.nr; ++i)
                {
                    const double c = p.c[f][i + p.nr * j];
                    q[f][0] += c * vr[i] * va[j];
                    q[f][1] += c * dr[i] * va[j] * tr;
                    q[f][2] += c * vr[i] * da[j];
                }
        const double pos[] = {x, y, z};
        double dW[3][3] = {}, div = 0.;
        for (int j = 0; j < 3; ++j)
        {
            const double rj = pos[j] / r,
                         muj = ((j == 0 ? 1. : 0.) - mu * rj) / r;
            const double aj = q[1][1] * rj + q[1][2] * muj,
                         bj = q[2][1] * rj + q[2][2] * muj;
            dW[0][j] = r / (1 + r) * aj + q[1][0] * rj / ((1 + r) * (1 + r));
            for (int i = 1; i < 3; ++i)
                dW[i][j] = pos[i] / (1 + r) * bj +
                           q[2][0] * ((i == j ? 1. : 0.) / (1 + r) -
                                      pos[i] * rj / ((1 + r) * (1 + r)));
            div += dW[j][j];
        }
        correction_t out;
        out.logpsi = q[0][0];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                out.strain[i][j] =
                    dW[i][j] + dW[j][i] - (i == j ? 2 * div / 3 : 0.);
        return out;
    }
    std::array<double, 2>
    perturbation(double r, double mu,
                 const EMSTrumpetSolution_read::ems_radial_vars_t &a) const
    {
        if (m_values.at("amplitude") == 0 || r <= m_left || r >= m_right)
            return {0., 0.};
        const auto xi = phase(r);
        const double sigma = m_values.at("sigma_M"), z = xi[0] / sigma,
                     az = std::abs(z);
        if (az >= 8)
            return {0., 0.};
        double w = 1., dw = 0.;
        if (az > 6)
        {
            const double t = (az - 6) / 2, l = std::exp(-1 / t),
                         h = std::exp(-1 / (1 - t));
            w = h / (l + h);
            dw = -w * (1 - w) * (1 / (t * t) + 1 / ((1 - t) * (1 - t))) *
                 std::copysign(.5, z);
        }
        const double g = std::exp(-z * z / 2),
                     Y = std::sqrt(5 / (16 * std::acos(-1.))) *
                         (3 * mu * mu - 1);
        const double scale = m_values.at("amplitude") * Y / (r / a.X),
                     phi = scale * g * w;
        return {phi, -scale * a.X * g * (dw - z * w) * xi[1] / sigma +
                         a.K_amplitude * phi};
    }
};
#endif
