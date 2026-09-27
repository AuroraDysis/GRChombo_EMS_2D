/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(EMSTRUMPETSOLUTION_READ_HPP_)
#error "This file should only be included through EMSTrumpetSolution_read.hpp"
#endif
#ifndef EMSTRUMPETSOLUTION_READ_IMPL_HPP_
#define EMSTRUMPETSOLUTION_READ_IMPL_HPP_

#include "MayDay.H"
#include "parstream.H"
#include "SHA256.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <limits>
#include <sstream>

inline bool
EMSTrumpetSolution_read::parse_finite_double(const std::string &a_text,
                                             double &a_value)
{
    char *end = nullptr;
    a_value = std::strtod(a_text.c_str(), &end);
    if (end == a_text.c_str() || !std::isfinite(a_value))
        return false;
    while (*end == ' ' || *end == '\t')
        ++end;
    return !*end;
}

inline bool EMSTrumpetSolution_read::parse_degree(const std::string &a_text,
                                                  int &a_value)
{
    std::istringstream in(a_text);
    std::string token, extra;
    if (!(in >> token) || (in >> extra) || token.empty())
        return false;
    if (token[0] == '+')
        token.erase(0, 1);
    if (token.empty() || !std::all_of(token.begin(), token.end(), [](char c)
                                      { return c >= '0' && c <= '9'; }))
        return false;
    // Compare strings to avoid integer overflow and nondecimal input.
    token.erase(0, token.find_first_not_of('0'));
    for (int allowed : {64, 96, 128, 192})
        if (token == std::to_string(allowed))
        {
            a_value = allowed;
            return true;
        }
    return false;
}

inline std::vector<double>
EMSTrumpetSolution_read::compute_derivative_coefficients(
    const std::vector<double> &a_coefficients)
{
    const int n = static_cast<int>(a_coefficients.size()) - 1;
    std::vector<double> b(a_coefficients.size(), 0.0);
    for (int k = n - 1; k >= 1; --k)
        b[k] = 2.0 * (k + 1) * a_coefficients[k + 1] +
               (k + 2 <= n ? b[k + 2] : 0.0);
    if (n > 0)
        b[0] = a_coefficients[1] + (n > 1 ? b[2] / 2.0 : 0.0);
    for (double &v : b)
        v *= 2.0; // x = 2s-1
    return b;
}

inline double EMSTrumpetSolution_read::evaluate_clenshaw(
    const std::vector<double> &a_coefficients, double a_s)
{
    const double x = 2.0 * a_s - 1.0;
    double b = 0.0, d = 0.0;
    for (int k = static_cast<int>(a_coefficients.size()) - 1; k >= 1; --k)
    {
        const double next = 2.0 * x * b - d + a_coefficients[k];
        d = b;
        b = next;
    }
    return x * b - d + a_coefficients[0];
}

inline double EMSTrumpetSolution_read::compute_coupling(double a_phi) const
{
    return std::exp(-2.0 * get_metadata_value("ems_alpha") *
                    (get_metadata_value("ems_f0") +
                     get_metadata_value("ems_f1") * a_phi +
                     get_metadata_value("ems_f2") * a_phi * a_phi));
}

inline void EMSTrumpetSolution_read::main(const std::string &a_data_path)
{
    std::string reason;
    if (!check_file(a_data_path, reason))
        MayDay::Error(("EMSTRUMPET 1: " + reason).c_str());
    pout() << "Read EMSTRUMPET 1 profile " << a_data_path << std::endl;
}

inline bool EMSTrumpetSolution_read::check_file(const std::string &a_data_path,
                                                std::string &a_reason)
{
    m_data_path = a_data_path;
    a_reason.clear();
    return read_from_file(a_reason);
}

inline bool EMSTrumpetSolution_read::read_from_file(std::string &a_reason)
{
    const auto fail = [&](const std::string &a_message)
    {
        a_reason = a_message;
        return false;
    };
    m_source_sha256.clear();
    std::ifstream file(m_data_path, std::ios::binary);
    if (!file)
        return fail("cannot open trumpet file");
    // Hash exactly the bytes parsed, not a second read of a mutable path.
    const std::string bytes((std::istreambuf_iterator<char>(file)), {});
    if (file.bad())
        return fail("cannot read trumpet file");
    std::istringstream in(bytes);
    std::string line;
    if (!std::getline(in, line))
        return fail("empty trumpet file");
    if (line != "EMSTRUMPET 1")
        return fail("unknown trumpet schema");
    const std::map<std::string, std::string> required_text = {
        {"schema", "EMSTRUMPET/1"},
        {"convention", "ems-native-draft-v1"},
        {"slice", "stationary-maximal-critical-trumpet-v1"},
        {"coordinate", "R=ell*s^(1/nu)/(1-s)"},
        {"coefficients", "Chebyshev-Tk(2s-1)-full-a0"}};
    const std::vector<std::string> numeric_keys = {
        "M",     "ell",       "r_C",    "nu",      "C",
        "Q",     "q_native",  "Q_s",    "phi_inf", "degree",
        "r_h",   "R_h",       "T_H",    "A_H",     "Phi_e",
        "phi_h", "ems_alpha", "ems_f0", "ems_f1",  "ems_f2"};
    std::map<std::string, std::string> meta;
    while (std::getline(in, line) && line.rfind("# ", 0) == 0)
    {
        const auto pos = line.find('=', 2);
        if (pos == std::string::npos)
            return fail("malformed metadata");
        const std::string key = line.substr(2, pos - 2);
        const std::string entry = line.substr(pos + 1);
        const bool known = required_text.count(key) ||
                           std::find(numeric_keys.begin(), numeric_keys.end(),
                                     key) != numeric_keys.end();
        if (!known && (key.rfind("x.", 0) != 0 || key.size() <= 2))
            return fail("unknown key " + key);
        if (key.rfind("x.", 0) == 0)
        {
            const std::string suffix = key.substr(2);
            if (required_text.count(suffix) ||
                std::find(numeric_keys.begin(), numeric_keys.end(), suffix) !=
                    numeric_keys.end())
                return fail("extension shadows required key " + key);
        }
        if (!meta.emplace(key, entry).second)
            return fail("duplicate key " + key);
        char *end = nullptr;
        const double maybe = std::strtod(entry.c_str(), &end);
        if (end != entry.c_str())
            while (*end == ' ' || *end == '\t')
                ++end;
        if (end != entry.c_str() && *end == '\0' && !std::isfinite(maybe))
            return fail("nonfinite metadata " + key);
    }
    for (const auto &kv : required_text)
    {
        if (!meta.count(kv.first))
            return fail("missing required metadata");
        if (meta.at(kv.first) != kv.second)
            return fail("unsupported " + kv.first);
    }
    for (const auto &key : numeric_keys)
    {
        if (!meta.count(key))
            return fail("missing required metadata");
        double number;
        if (!parse_finite_double(meta.at(key), number))
            return fail("nonfinite " + key);
        m_metadata_values[key] = number;
    }
    int polynomial_degree;
    if (!parse_degree(meta.at("degree"), polynomial_degree))
        return fail("invalid degree");
    for (int j = 0; j < 3; ++j)
    {
        const std::string label(1, "YDS"[j]);
        if (line != label + " " + std::to_string(polynomial_degree + 1))
            return fail("malformed " + label + " block");
        auto &coefficients = m_coefficients[j][0];
        coefficients.clear();
        coefficients.reserve(polynomial_degree + 1);
        for (int k = 0; k <= polynomial_degree; ++k)
        {
            if (!std::getline(in, line))
                return fail("truncated " + label + " block");
            std::istringstream row(line);
            std::string token, extra;
            if (!(row >> token) || (row >> extra))
                return fail("malformed coefficient");
            double coefficient;
            if (!parse_finite_double(token, coefficient))
                return fail("nonfinite coefficient");
            coefficients.push_back(coefficient);
        }
        for (int d = 1; d <= 3; ++d)
            m_coefficients[j][d] =
                compute_derivative_coefficients(m_coefficients[j][d - 1]);
        if (j != 2 && !std::getline(in, line))
            return fail("malformed " + std::string(1, "DS"[j]) + " block");
    }
    if (!std::getline(in, line) || line != "END" || std::getline(in, line))
        return fail("invalid file ending");
    if (get_metadata_value("M") != 1.0 || get_metadata_value("ell") <= 0 ||
        get_metadata_value("r_C") <= 0 ||
        get_metadata_value("r_h") <= get_metadata_value("r_C") ||
        get_metadata_value("C") >= 0 || get_metadata_value("R_h") <= 0)
        return fail("invalid trumpet geometry");
    constexpr double pi = 3.141592653589793238462643383279502884;
    if (std::abs(get_metadata_value("q_native") -
                 get_metadata_value("Q") / std::sqrt(8.0 * pi)) >= 2e-15)
        return fail("inconsistent native charge");
    const double ell = get_metadata_value("ell"), nu_ = get_nu();
    const std::array<double, 6> checks = {
        get_chebyshev_value(0, 0, 0) - get_metadata_value("r_C"),
        get_chebyshev_value(0, 1, 0) - ell,
        get_chebyshev_value(0, 1, 1) - (ell / nu_ - 1),
        get_chebyshev_value(2, 1, 0) - get_metadata_value("Q_s") / ell,
        get_chebyshev_value(1, 1, 0) -
            4 * pi * std::pow(get_metadata_value("Q_s") / ell, 2),
        nu_ * nu_ - (2 - std::pow(get_metadata_value("Q"), 2) /
                             (compute_coupling(get_metadata_value("phi_inf") +
                                               get_chebyshev_value(2, 0, 0)) *
                              std::pow(get_metadata_value("r_C"), 2)))};
    for (double check : checks)
        if (!(std::abs(check) < 1e-10))
            return fail("failed trumpet endpoint identity");
    double horizon_s;
    if (!invert_compactified_coordinate(get_metadata_value("R_h"), horizon_s))
        return fail("compactification inversion failed");
    if (!(std::abs(compute_radial_vars(get_metadata_value("R_h")).r -
                   get_metadata_value("r_h")) < 1e-10))
        return fail("inconsistent horizon radius");
    m_source_sha256 = SHA256::digest(bytes);
    return true;
}

inline double EMSTrumpetSolution_read::get_chebyshev_value(int a_field,
                                                           double a_s,
                                                           int a_order) const
{
    return evaluate_clenshaw(m_coefficients[a_field][a_order], a_s);
}

inline std::array<double, 4>
EMSTrumpetSolution_read::get_coupling_parameters() const
{
    return {get_metadata_value("ems_alpha"), get_metadata_value("ems_f0"),
            get_metadata_value("ems_f1"), get_metadata_value("ems_f2")};
}

inline bool
EMSTrumpetSolution_read::invert_compactified_coordinate(double a_R,
                                                        double &a_s) const
{
    if (!(a_R > 0) || !std::isfinite(a_R))
        return false;
    const double x = std::log(a_R / get_metadata_value("ell")), nu_ = get_nu();
    double lo = 0, hi = 1,
           u = x <= 0 ? std::exp(nu_ * x) / 2 : std::exp(-x) / 2;
    for (int i = 0; i < 100; ++i)
    {
        const double f = x <= 0 ? std::log(u) / nu_ - std::log1p(-u) - x
                                : std::log1p(-u) / nu_ - std::log(u) - x;
        if (std::abs(f) < 100 * std::numeric_limits<double>::epsilon())
        {
            a_s = x <= 0 ? u : 1 - u;
            return true;
        }
        if (x <= 0 ? f > 0 : f < 0)
            hi = u;
        else
            lo = u;
        const double slope =
            x <= 0 ? 1 / (nu_ * u) + 1 / (1 - u) : -1 / (nu_ * (1 - u)) - 1 / u;
        const double trial = u - f / slope;
        u = lo < trial && trial < hi ? trial : (lo + hi) / 2;
    }
    return false;
}

inline double
EMSTrumpetSolution_read::get_compactified_coordinate(double a_R) const
{
    double s = 0;
    if (!invert_compactified_coordinate(a_R, s))
        MayDay::Error("EMSTRUMPET compactification inversion failed");
    return s;
}

inline EMSTrumpetSolution_read::ems_radial_vars_t
EMSTrumpetSolution_read::compute_radial_vars(double a_R) const
{
    const double s = get_compactified_coordinate(a_R), t = 1 - s;
    const double Y = get_chebyshev_value(0, s),
                 Y1 = get_chebyshev_value(0, s, 1),
                 Y2 = get_chebyshev_value(0, s, 2);
    const double D = get_chebyshev_value(1, s),
                 D1 = get_chebyshev_value(1, s, 1);
    const double S = get_chebyshev_value(2, s),
                 S1 = get_chebyshev_value(2, s, 1);
    const double ell = get_metadata_value("ell"), nu_ = get_nu(),
                 C = get_metadata_value("C");
    const double G = Y + t * Y1, ks = 1 + (nu_ - 1) * s;
    const double r = Y / t, delta = t * t * D, phi = get_phi_inf() + t * S;
    const double sp = std::pow(s, 1 / nu_);
    const double X = ell * sp / Y;
    const double alpha = std::exp(-delta) * nu_ * s * G / (Y * ks);
    const double beta = -C * ell * sp * t * t / (Y * Y * Y);
    const double k = C * std::exp(delta) * t * t * t / (Y * Y * Y);
    const double phis = -S + t * S1;
    const double Pi =
        -C * std::exp(delta) * std::pow(t, 4) * phis / (Y * Y * G);
    const double ER =
        get_metadata_value("q_native") / (compute_coupling(phi) * r * a_R);
    const double sR = nu_ * s * t / (a_R * ks);
    const double rR = G / (t * t) * sR;
    const double deltas = -2 * t * D + t * t * D1;
    const double PR = Y / (ell * sp) * (Y1 / Y - 1 / (nu_ * s)) * sR;
    const double alphaR =
        alpha * (-deltas + 1 / s + t * Y2 / G - Y1 / Y - (nu_ - 1) / ks) * sR;
    const double betaR = beta * (1 / a_R - 3 * rR / r);
    return {s,  r,  delta, phi, X,  alpha,  beta, k,
            Pi, ER, G,     ks,  PR, alphaR, betaR};
}

#endif /* EMSTRUMPETSOLUTION_READ_IMPL_HPP_ */
