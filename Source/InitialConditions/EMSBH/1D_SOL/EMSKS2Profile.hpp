#ifndef EMSKS2_PROFILE_HPP_
#define EMSKS2_PROFILE_HPP_

#include "SHA256.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <initializer_list>
#include <iterator>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace EMSKS2Detail
{
// Six source Taylor orders suffice to verify the first five endpoint
// coefficients after removing the cylinder's double root.
struct Jet
{
    std::array<double, 7> c{};
    Jet(double x = 0) { c[0] = x; }
    static Jet variable(double x) { Jet a(x); a.c[1] = 1; return a; }
};
inline Jet operator+(Jet a, Jet b)
{ for (int k = 0; k < 7; ++k) a.c[k] += b.c[k]; return a; }
inline Jet operator-(Jet a, Jet b)
{ for (int k = 0; k < 7; ++k) a.c[k] -= b.c[k]; return a; }
inline Jet operator*(Jet a, Jet b)
{
    Jet out;
    for (int k = 0; k < 7; ++k)
        for (int j = 0; j <= k; ++j) out.c[k] += a.c[j] * b.c[k - j];
    return out;
}
inline Jet inverse(Jet a)
{
    Jet out(1 / a.c[0]);
    for (int k = 1; k < 7; ++k)
    {
        for (int j = 1; j <= k; ++j) out.c[k] -= a.c[j] * out.c[k - j];
        out.c[k] /= a.c[0];
    }
    return out;
}
inline Jet operator/(Jet a, Jet b) { return a * inverse(b); }
inline Jet exponential(Jet a)
{
    Jet out(std::exp(a.c[0]));
    for (int k = 1; k < 7; ++k)
    {
        for (int j = 1; j <= k; ++j)
            out.c[k] += j * a.c[j] * out.c[k - j];
        out.c[k] /= k;
    }
    return out;
}
inline Jet square_root(Jet a)
{
    if (!(a.c[0] > 0)) throw std::runtime_error("nonpositive cylinder factor");
    Jet out(std::sqrt(a.c[0]));
    for (int k = 1; k < 7; ++k)
    {
        out.c[k] = a.c[k];
        for (int j = 1; j < k; ++j) out.c[k] -= out.c[j] * out.c[k - j];
        out.c[k] /= 2 * out.c[0];
    }
    return out;
}
inline Jet chebyshev(const std::vector<double> &c, Jet z)
{
    Jet b, d;
    for (int k = static_cast<int>(c.size()) - 1; k >= 1; --k)
    { Jet next = Jet(c[k]) + Jet(2) * z * b - d; d = b; b = next; }
    return Jet(c[0]) + z * b - d;
}
} // namespace EMSKS2Detail

// EMSKS 2 uses full-a0 Chebyshev coefficients throughout. This reader has no
// Chombo dependency so the serialized data can be checked before grid setup.
class EMSKS2Profile
{
  public:
    struct Sample
    {
        double rho, r, N, H, sigma, g, lambda, w, alpha, alpha_r;
        double S, S_r, k_r, k_t, K, phi, phi_r, Pi, E_rho, rho_a;
        double chi, lambda_rho, alpha_rho, beta_rho, q_star, taper;
    };
    struct Cartesian
    {
        Sample s;
        std::array<double, 9> h{}, A{};
        std::array<double, 3> Gamma{}, beta{}, E{};
        double hww, Aww;
    };

    Cartesian object(double x, double y, double z = 0) const
    {
        const double rho = std::sqrt(x * x + y * y + z * z);
        const auto s = sample(rho);
        Cartesian out{s};
        const double n[3] = {x / rho, y / rho, z / rho};
        const double T = std::pow(s.r / rho, 2), L = s.lambda * T;
        const double m = 1 / std::cbrt(s.lambda);
        const double G = 2 * s.lambda_rho /
                (3 * std::pow(s.lambda, 5. / 3.)) +
                2 * (s.lambda - 1) /
                (rho * std::pow(s.lambda, 2. / 3.));
        for (int i = 0; i < 3; ++i)
        {
            out.Gamma[i] = G * n[i];
            out.beta[i] = s.beta_rho * n[i];
            out.E[i] = s.E_rho * n[i];
            for (int j = 0; j < 3; ++j)
            {
                const double delta = i == j ? 1 : 0;
                const double gamma = T * delta + (L - T) * n[i] * n[j];
                const double Kij = T * s.k_t * delta +
                                    (L * s.k_r - T * s.k_t) * n[i] * n[j];
                out.h[3 * i + j] = s.chi * gamma;
                out.A[3 * i + j] = s.chi * (Kij - s.K * gamma / 3);
            }
        }
        out.hww = m;
        out.Aww = s.chi * T * (s.k_t - s.K / 3);
        return out;
    }

    void load(const std::string &path)
    {
        std::ifstream file(path, std::ios::binary);
        if (!file) throw std::runtime_error("cannot open KS2 file");
        const std::string bytes((std::istreambuf_iterator<char>(file)), {});
        if (file.bad() || bytes.empty() || bytes.back() != '\n')
            throw std::runtime_error("truncated KS2 file");
        const auto hashpos = bytes.rfind("SHA256 ");
        if (hashpos == std::string::npos || hashpos == 0 ||
            bytes.size() != hashpos + 72 || bytes[hashpos - 1] != '\n')
            throw std::runtime_error("unexpected trailing KS2 data");
        const std::string digest = bytes.substr(hashpos + 7, 64);
        if (!std::all_of(digest.begin(), digest.end(), [](char c) {
                return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f');
            }) || bytes[hashpos + 71] != '\n')
            throw std::runtime_error("bad KS2 hash syntax");
        if (SHA256::digest(std::string_view(bytes.data(), hashpos)) != digest)
            throw std::runtime_error("KS2 hash mismatch");
        fingerprint_ = SHA256::digest(bytes);
        std::vector<std::string> lines;
        std::size_t begin = 0;
        for (std::size_t end = bytes.find('\n'); end != std::string::npos;
             begin = end + 1, end = bytes.find('\n', begin))
            lines.push_back(bytes.substr(begin, end - begin));
        if (lines.empty() || lines[0] != "EMSKS 2")
            throw std::runtime_error("unknown KS version");
        std::size_t pos = 1;
        const std::map<std::string, std::string> text = {
            {"schema", "EMSKS/2"},
            {"convention", "ems-native-draft-v1"},
            {"slice", "KS-critical-maximal-S8-v1"},
            {"time_normalization", "delta(infinity)=0"},
            {"coefficient_basis", "Chebyshev-Tk-full-a0"},
            {"exterior_map", "identity-or-mass-sinh-v1"},
            {"interior_map", "affine-u-v1"},
            {"blend", "S8-(r-ra)/(rm-ra)-g-lambda-v1"},
            {"chart", "Z-anchored-antiderivative-of-Zprime-Chebyshev-v1"},
            {"gauge", "prescribed-Killing-shift-relative-harmonic-v1"},
            {"ko_taper", "z2+(1-z2)S8(z)-v1"}};
        const std::vector<std::string> numeric = {
            "ems_alpha", "ems_f0", "ems_f1", "ems_f2", "M", "Q", "e",
            "r_h", "r_exc", "phi_inf", "exterior_lo", "exterior_hi",
            "exterior_h", "exterior_beta", "interior_lo", "interior_hi",
            "exterior_degree", "interior_degree", "r_0", "D", "nu",
            "r_a", "r_m", "chart_degree", "kappa_alpha", "b_0", "q_limit"};
        std::string previous;
        while (pos < lines.size() && lines[pos].rfind("# ", 0) == 0)
        {
            const auto eq = lines[pos].find('=', 2);
            if (eq == std::string::npos) throw std::runtime_error("bad KS2 header");
            const auto key = lines[pos].substr(2, eq - 2);
            if (!previous.empty() && key <= previous)
                throw std::runtime_error("unsorted or duplicate KS2 metadata");
            previous = key;
            const auto value = lines[pos++].substr(eq + 1);
            if (auto it = text.find(key); it != text.end())
            {
                if (it->second != value)
                    throw std::runtime_error("unsupported " + key);
            }
            else if (std::find(numeric.begin(), numeric.end(), key) != numeric.end())
                values_[key] = number(value);
            else throw std::runtime_error("unknown KS2 metadata " + key);
            keys_.push_back(key);
        }
        if (keys_.size() != text.size() + numeric.size() ||
            values_.size() != numeric.size())
            throw std::runtime_error("missing KS2 metadata");
        auto degree = [&](const char *key, int lo, int hi) {
            const double v = get(key);
            if (v != std::floor(v) || v < lo || v > hi)
                throw std::runtime_error(std::string("invalid ") + key);
            return static_cast<int>(v);
        };
        const int ne = degree("exterior_degree", 2, 384) + 1;
        const int ni = degree("interior_degree", 2, 384) + 1;
        const int nc = degree("chart_degree", 24, 128) + 1;
        const double rh = get("r_h"), r0 = get("r_0"), ra = get("r_a"),
                     rm = get("r_m"), rexc = get("r_exc");
        if (!(get("M") > 0 && rexc > 0 && rexc < r0 && r0 < ra &&
              ra < rm && rm < rh && r0 - rexc >= .02 * rh))
            throw std::runtime_error("invalid KS2 domains");
        if (!(near(rm, .9 * rh, 1e-13 * rh) &&
              near(ra, r0 + (rm - r0) / 20, 1e-13 * rh)))
            throw std::runtime_error("invalid KS2 join");
        if (!(get("D") > 0 && get("nu") > 1 &&
              get("nu") < std::sqrt(2.) + 1e-12))
            throw std::runtime_error("invalid cylinder");
        if (!(get("exterior_lo") == 0 && get("exterior_hi") == 1 &&
              get("interior_lo") == 1 && get("interior_hi") > 1 &&
              near(rh / get("interior_hi"), rexc, 1e-13 * rh)))
            throw std::runtime_error("invalid element domains");
        if (!near(get("Q"), get("e") * rh,
                  1e-13 * std::max(1., std::abs(get("Q")))))
            throw std::runtime_error("charge mismatch");
        if (!(get("ems_alpha") > 0 && get("ems_f1") == 0 &&
              get("ems_f2") <= 0 && get("exterior_h") >= 0 &&
              get("exterior_beta") >= 0 &&
              (get("exterior_beta") == 0 || get("exterior_h") > 0)))
            throw std::runtime_error("invalid coupling or exterior map");
        auto block = [&](const std::string &label, int n,
                         std::initializer_list<const char *> names) {
            if (pos >= lines.size() || lines[pos++] != label + " " + std::to_string(n))
                throw std::runtime_error("malformed " + label + " block");
            std::vector<std::vector<double>> rows;
            for (const char *name : names)
            {
                if (pos >= lines.size() || lines[pos++] !=
                        std::string(name) + " " + std::to_string(n))
                    throw std::runtime_error(std::string("malformed ") + name + " block");
                std::vector<double> row;
                for (int k = 0; k < n; ++k)
                {
                    if (pos >= lines.size()) throw std::runtime_error("truncated coefficients");
                    row.push_back(number(lines[pos++]));
                }
                rows.push_back(std::move(row));
            }
            return rows;
        };
        exterior_ = block("EXTERIOR", ne, {"mu", "phi", "p", "delta", "V"});
        interior_ = block("INTERIOR", ni, {"mu", "phi", "p", "delta", "V"});
        if (pos >= lines.size() || lines[pos++] != "ENDPOINT 16")
            throw std::runtime_error("missing endpoint block");
        for (int k = 0; k < 16; ++k)
        {
            if (pos >= lines.size()) throw std::runtime_error("truncated endpoint");
            endpoint_[k] = number(lines[pos++]);
        }
        chart_ = block("CHART", nc + 1, {"CORE", "JOIN"});
        chartprime_ = block("CHART_DERIV", nc, {"CORE", "JOIN"});
        if (pos + 2 != lines.size() || lines[pos] != "END" ||
            lines[pos + 1].rfind("SHA256 ", 0) != 0)
            throw std::runtime_error("unexpected trailing KS2 data");
        validate();
        rho_a_ = chart(ra);
    }

    bool check_file(const std::string &path, std::string &reason)
    {
        try { EMSKS2Profile tmp; tmp.load(path); *this = std::move(tmp); reason.clear(); return true; }
        catch (const std::exception &e) { reason = e.what(); return false; }
    }
    double get(const std::string &key) const { return values_.at(key); }
    const std::string &fingerprint() const { return fingerprint_; }
    double rho_a() const { return rho_a_; }

    std::array<double, 5> source(double r, bool interior, int derivative = 0) const
    {
        const double u = get("r_h") / r;
        double y = u, yu = 1, yuu = 0;
        if (!interior && get("exterior_h") != 0)
        {
            const double h = get("exterior_h"), b = get("exterior_beta");
            const double t = (1 + h) * u / (h + u);
            const double tu = h * (1 + h) / ((h + u) * (h + u));
            const double tuu = -2 * tu / (h + u);
            if (b == 0) { y = t; yu = tu; yuu = tuu; }
            else
            {
                const double sh = std::sinh(b), v = (2 * t - 1) * sh;
                const double root = std::sqrt(1 + v * v);
                y = .5 + std::asinh(v) / (2 * b);
                const double yt = sh / (b * root);
                const double ytt = -2 * sh * sh * v / (b * root * root * root);
                yu = yt * tu;
                yuu = ytt * tu * tu + yt * tuu;
            }
        }
        const auto &rows = interior ? interior_ : exterior_;
        const double lo = interior ? get("interior_lo") : 0.;
        const double hi = interior ? get("interior_hi") : 1.;
        std::array<double, 5> out{};
        for (int j = 0; j < 5; ++j)
        {
            const double z = (2 * y - lo - hi) / (hi - lo);
            const auto p = cheb(rows[j], z);
            if (derivative == 0) out[j] = p[0];
            else if (derivative == 1) out[j] = p[1] * 2 * yu / (hi - lo);
            else out[j] = p[2] * std::pow(2 * yu / (hi - lo), 2) +
                          p[1] * 2 * yuu / (hi - lo);
        }
        return out;
    }

    Sample sample(double rho) const
    {
        if (!(rho > 0 && std::isfinite(rho)))
            throw std::runtime_error("KS2 puncture or nonfinite radius");
        const double r = areal(rho), r0 = get("r_0"), ra = get("r_a"),
                     rm = get("r_m"), rh = get("r_h"), x = r - r0;
        const auto d = source(r, r < rh);
        auto dr = source(r, r < rh, 1);
        for (double &v : dr) v *= -rh / (r * r);
        const double N = 1 - 2 * rh * d[0] / r, H = 2 - N;
        const double Hr = 2 * rh * (dr[0] - d[0] / r) / r;
        const double sigma = std::exp(-d[3]);
        const double A = N * sigma * sigma;
        const double Ar = sigma * sigma * (-Hr - 2 * N * dr[3]);
        double g, gp, lambda, lp, w, alpha, ap, S, Sp;
        if (r <= ra)
        {
            const auto a = endpoint(r - r0);
            g = r * a[0] / x;
            gp = g * (1 / r + a[1] / a[0] - 1 / x);
            lambda = 1; lp = 0; w = 0;
            alpha = sigma * x / (r * a[0]);
            ap = alpha * (-dr[3] + 1 / x - 1 / r - a[1] / a[0]);
            S = get("D") / (r * r); Sp = -2 * S / r;
        }
        else if (r < rm)
        {
            const double z = (r - ra) / (rm - ra);
            w = smooth(z);
            const double wp = std::pow(z * (1 - z), 8) *
                              218790. / (rm - ra); // 1/B(9,9)
            const auto a = endpoint(x);
            const double invht = r * a[0] / x;
            const double invhtp = invht * (1 / r + a[1] / a[0] - 1 / x);
            g = (1 - w) * invht + w;
            gp = (1 - w) * invhtp + wp * (1 - invht);
            lambda = std::exp(w * std::log(H));
            lp = lambda * (wp * std::log(H) + w * Hr / H);
            alpha = sigma / (g * std::sqrt(lambda));
            ap = alpha * (-dr[3] - gp / g - lp / (2 * lambda));
            const double S2 = alpha * alpha - A;
            if (!(S2 >= 0)) throw std::runtime_error("imaginary KS2 slice S");
            S = std::sqrt(S2); Sp = (alpha * ap - Ar / 2) / S;
        }
        else
        {
            g = 1; gp = 0; lambda = H; lp = Hr; w = 1;
            alpha = sigma / std::sqrt(H);
            ap = alpha * (-dr[3] - Hr / (2 * H));
            S = sigma * (H - 1) / std::sqrt(H);
            Sp = (alpha * ap - Ar / 2) / S;
        }
        const double drdrho = r / (g * rho);
        const double F = coupling(d[1]);
        const double q = get("Q") / std::sqrt(8 * pi);
        const double E = q * std::sqrt(lambda) / (F * r * rho);
        const double beta = rho * S / (r * std::sqrt(lambda));
        const double z = rho / rho_a_;
        const double taper = z >= 1 ? 1 : z * z + (1 - z * z) * smooth(z);
        return {rho, r, N, H, sigma, g, lambda, w, alpha, ap,
                S, Sp, Sp / sigma, S / (sigma * r),
                (Sp + 2 * S / r) / sigma, d[1], dr[1],
                (S / sigma) * dr[1], E, rho_a_,
                rho * rho / (r * r * std::cbrt(lambda)), lp * drdrho,
                ap * drdrho, beta, (S / sigma) * ap, taper};
    }

  private:
    static constexpr double pi = 3.141592653589793238462643383279502884;
    std::map<std::string, double> values_;
    std::vector<std::string> keys_;
    std::vector<std::vector<double>> exterior_, interior_, chart_, chartprime_;
    std::array<double, 16> endpoint_{};
    double rho_a_ = 0;
    std::string fingerprint_;
    static bool near(double a, double b, double tolerance)
    { return std::abs(a - b) < tolerance; }
    static double number(const std::string &s)
    {
        char *end = nullptr;
        const double v = std::strtod(s.c_str(), &end);
        if (end == s.c_str() || *end || !std::isfinite(v))
            throw std::runtime_error("malformed coefficient or metadata");
        return v;
    }
    static std::array<double, 3> cheb(const std::vector<double> &c, double z)
    {
        double b0 = 0, b1 = 0, b2 = 0, d0 = 0, d1 = 0, d2 = 0;
        double e0 = 0, e1 = 0, e2 = 0;
        for (int k = static_cast<int>(c.size()) - 1; k >= 1; --k)
        {
            b2 = b1; b1 = b0; b0 = c[k] + 2 * z * b1 - b2;
            e2 = e1; e1 = e0; e0 = 4 * d1 + 2 * z * e1 - e2;
            d2 = d1; d1 = d0; d0 = 2 * b1 + 2 * z * d1 - d2;
        }
        return {c[0] + z * b0 - b1, b0 + z * d0 - d1,
                2 * d0 + z * e0 - e1};
    }
    static double smooth(double z)
    {
        if (z <= 0) return 0;
        if (z >= 1) return 1;
        if (z > .5) return 1 - smooth(1 - z);
        double sum = 0;
        for (int k = 9; k <= 17; ++k)
        {
            int binom = 1;
            for (int j = 1; j <= k; ++j) binom = binom * (18 - j) / j;
            sum += binom * std::pow(z, k) * std::pow(1 - z, 17 - k);
        }
        return sum;
    }
    std::array<double, 2> endpoint(double x) const
    {
        double a = endpoint_[15], ap = 0;
        for (int k = 14; k >= 0; --k)
        { ap = ap * x + a; a = a * x + endpoint_[k]; }
        return {a, ap};
    }
    double series(const std::vector<double> &c, double r, double lo, double hi) const
    { return cheb(c, (2 * r - lo - hi) / (hi - lo))[0]; }
    double chart(double r) const
    {
        if (r >= get("r_m")) return r;
        const bool core = r <= get("r_a");
        const double Z = series(chart_[core ? 0 : 1], r,
                                core ? get("r_0") : get("r_a"),
                                core ? get("r_a") : get("r_m"));
        return get("r_m") * std::exp(std::log((r - get("r_0")) /
                       (get("r_m") - get("r_0"))) / get("nu") + Z);
    }
    double areal(double rho) const
    {
        if (rho >= get("r_m")) return rho;
        double lo = get("r_0"), hi = get("r_m");
        const double target = std::log(rho / hi);
        for (int k = 0; k < 100; ++k)
        {
            const double mid = (lo + hi) / 2;
            if (mid == lo) return hi;
            const bool core = mid <= get("r_a");
            const double Z = series(chart_[core ? 0 : 1], mid,
                                    core ? get("r_0") : get("r_a"),
                                    core ? get("r_a") : get("r_m"));
            const double log_rho = std::log((mid - get("r_0")) /
                       (get("r_m") - get("r_0"))) / get("nu") + Z;
            if (log_rho > target) hi = mid; else lo = mid;
        }
        return (lo + hi) / 2;
    }
    double coupling(double phi) const
    { return std::exp(-2 * get("ems_alpha") * (get("ems_f0") +
             get("ems_f1") * phi + get("ems_f2") * phi * phi)); }
    void validate() const
    {
        const double rh = get("r_h"), r0 = get("r_0"), ra = get("r_a"),
                     rm = get("r_m"), nu = get("nu"), D = get("D");
        if (!near(source(1e100, false)[0] * rh, get("M"), 1e-12 * get("M")) ||
            !near(source(1e100, false)[1], get("phi_inf"), 1e-12) ||
            !near(source(1e100, false)[3], 0, 1e-12))
            throw std::runtime_error("mass, scalar infinity or time normalization mismatch");
        const auto d = source(r0, true);
        const double N = 1 - 2 * rh * d[0] / r0;
        const double A = N * std::exp(-2 * d[3]);
        const double F = coupling(d[1]);
        const auto critical_B = [&](double r) {
            const auto v = source(r, r < rh), dv = source(r, r < rh, 1);
            const double n = 1 - 2 * rh * v[0] / r;
            const double phir = -rh * dv[1] / (r * r);
            return 1 + n * (3 + 8 * pi * r * r * phir * phir) -
                   get("Q") * get("Q") / (coupling(v[1]) * r * r);
        };
        double lo = get("r_exc"), hi = rh;
        if (!(critical_B(lo) < 0 && critical_B(hi) > 0))
            throw std::runtime_error("critical cylinder not bracketed");
        for (int k = 0; k < 100; ++k)
        {
            const double mid = (lo + hi) / 2;
            if (mid == lo || mid == hi) break;
            if (critical_B(mid) > 0) hi = mid; else lo = mid;
        }
        const double root = (lo + hi) / 2;
        if (!near(root, r0, 1e-12 * std::max(1., r0)) ||
            !near(D, std::sqrt(-r0 * r0 * r0 * r0 * A), 1e-12) ||
            !near(nu * nu, 2 - get("Q") * get("Q") / (F * r0 * r0), 1e-12))
            throw std::runtime_error("critical cylinder mismatch");
        if (!near(get("kappa_alpha") * get("M"), 1, 1e-12) ||
            !near(get("b_0"), D / (r0 * r0 * r0), 1e-12) ||
            !near(get("q_limit"), get("b_0") * nu, 1e-12))
            throw std::runtime_error("gauge metadata mismatch");
        if (!near(endpoint_[0], 1 / nu, 1e-10))
            throw std::runtime_error("endpoint/source mismatch");
        {
            using namespace EMSKS2Detail;
            const Jet r = Jet::variable(r0), u = Jet(rh) / r;
            const Jet z = (Jet(2) * u - Jet(get("interior_lo") +
                          get("interior_hi"))) /
                          Jet(get("interior_hi") - get("interior_lo"));
            const Jet mu = chebyshev(interior_[0], z);
            const Jet delta = chebyshev(interior_[3], z);
            const Jet N = Jet(1) - Jet(2) * u * mu;
            const Jet rr = r * r;
            const Jet U = Jet(-1) * rr * rr * N *
                          exponential(Jet(-2) * delta);
            Jet factor;
            for (int k = 0; k < 5; ++k) factor.c[k] = -U.c[k + 2];
            const Jet source_endpoint = r * exponential(Jet(-1) * delta) /
                                        square_root(factor);
            for (int k = 0; k < 5; ++k)
                if (!near(source_endpoint.c[k], endpoint_[k], 1e-10))
                    throw std::runtime_error("endpoint/source mismatch");
        }
        // The stored endpoint series must reconstruct the source's controlled
        // double-root factor throughout the resolved core and joining layer.
        for (int k = 1; k <= 64; ++k)
        {
            const double r = r0 + (rm - r0) * k / 64.;
            const auto d = source(r, true);
            const double N = 1 - 2 * rh * d[0] / r;
            const double sigma = std::exp(-d[3]);
            const double source_alpha = std::sqrt(N * sigma * sigma +
                                                  D * D / std::pow(r, 4));
            const double represented_alpha = sigma * (r - r0) /
                                               (r * endpoint(r - r0)[0]);
            if (!(source_alpha > 0 &&
                  std::abs(source_alpha - represented_alpha) /
                      source_alpha < 1e-10))
                throw std::runtime_error("endpoint/source mismatch");
        }
        auto integrated = [&](const std::vector<double> &c, double lo,
                              double hi, double zhi) {
            std::vector<double> b(c.size() + 1, 0);
            const double scale = (hi - lo) / 2;
            b[1] += scale * c[0];
            if (c.size() >= 2) b[2] += scale * c[1] / 4;
            for (std::size_t k = 2; k < c.size(); ++k)
            { b[k + 1] += scale * c[k] / (2 * (k + 1));
              b[k - 1] -= scale * c[k] / (2 * (k - 1)); }
            double sum = 0; for (std::size_t k = 1; k < b.size(); ++k) sum += b[k];
            b[0] = zhi - sum;
            return b;
        };
        const auto j = integrated(chartprime_[1], ra, rm, 0);
        const auto c = integrated(chartprime_[0], r0, ra,
                                  series(j, ra, ra, rm));
        for (int i = 0; i < 2; ++i)
            for (std::size_t k = 0; k < chart_[i].size(); ++k)
                if (!near(chart_[i][k], (i == 0 ? c : j)[k], 1e-11))
                    throw std::runtime_error("inconsistent chart antiderivative");
        for (int side = 0; side < 2; ++side)
        {
            const double lo = side ? ra : r0, hi = side ? rm : ra;
            for (int k = 0; k < 9; ++k)
            {
                const double r = lo + (hi - lo) * k / 8;
                const double zp = series(chartprime_[side], r, lo, hi);
                double expected;
                if (r == r0) expected = endpoint_[1];
                else if (r <= ra)
                {
                    const auto a = endpoint(r - r0);
                    expected = a[0] / (r - r0) - 1 / (nu * (r - r0));
                }
                else
                {
                    const auto d = source(r, true);
                    const double N = 1 - 2 * rh * d[0] / r;
                    const double sigma = std::exp(-d[3]);
                    const double a = endpoint(r - r0)[0];
                    const double z = (r - ra) / (rm - ra);
                    const double g = (1 - smooth(z)) * r * a / (r - r0) + smooth(z);
                    expected = g / r - 1 / (nu * (r - r0));
                    (void)N; (void)sigma;
                }
                if (!near(zp, expected, 1e-9))
                    throw std::runtime_error("chart derivative/source mismatch");
            }
        }
        for (int k = 1; k <= 4 * static_cast<int>(chartprime_[0].size()); ++k)
        {
            const double r = r0 + (rm - r0) * k /
                                     (4 * static_cast<double>(chartprime_[0].size()));
            const double rho = chart(r);
            if (!(rho > 0 && rho <= r + 1e-12 * rm))
                throw std::runtime_error("nonmonotone or expanding KS2 chart");
        }
    }
};

#endif
