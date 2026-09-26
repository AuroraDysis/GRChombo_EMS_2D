/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "EMSBH_trumpet_read.hpp"
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct fixture_table_t
{
    std::string name;
    double max = 0;
    std::string worst;
    int rows = 0, columns = 0;
    std::uint64_t fingerprint = 14695981039346656037ULL;
};
static std::vector<std::string> split_tsv_row(const std::string &a_row)
{
    std::vector<std::string> fields;
    std::istringstream in(a_row);
    std::string field;
    while (std::getline(in, field, '\t'))
        fields.push_back(field);
    return fields;
}
static std::vector<double> parse_tsv_values(const std::string &a_row)
{
    std::vector<double> parsed_values;
    for (const auto &field : split_tsv_row(a_row))
        parsed_values.push_back(std::stod(field));
    return parsed_values;
}
static std::vector<double>
ccz4_components(const CCZ4CartoonVars::VarsWithGauge<double> &a_vars)
{
    return {a_vars.chi,      a_vars.h[0][0], a_vars.h[0][1],  a_vars.h[1][1],
            a_vars.hww,      a_vars.K,       a_vars.A[0][0],  a_vars.A[0][1],
            a_vars.A[1][1],  a_vars.Aww,     a_vars.Theta,    a_vars.Gamma[0],
            a_vars.Gamma[1], a_vars.lapse,   a_vars.shift[0], a_vars.shift[1],
            a_vars.B[0],     a_vars.B[1],    a_vars.phi,      a_vars.Pi,
            a_vars.Ex,       a_vars.Ey,      a_vars.Ez,       a_vars.Bx,
            a_vars.By,       a_vars.Bz,      a_vars.Xi,       a_vars.Lambda};
}
static const std::vector<std::string> ccz4_component_names = {
    "chi",    "h11",    "h12", "h22",   "hww",    "K",      "A11",
    "A12",    "A22",    "Aww", "Theta", "Gamma1", "Gamma2", "lapse",
    "shift1", "shift2", "B1",  "B2",    "phi",    "Pi",     "Ex",
    "Ey",     "Ez",     "Bx",  "By",    "Bz",     "Xi",     "Lambda"};
static double fixture_tolerance(const std::string &a_table,
                                const std::string &a_column)
{
    if (a_table == "binary_ccz4.tsv")
        return 5e-13;
    if (a_table == "jets.tsv")
    {
        if (a_column.find("_d1") != std::string::npos)
            return 1e-10;
        if (a_column.find("_d2") != std::string::npos ||
            a_column.find("_d3") != std::string::npos)
            return 1e-8;
    }
    if (a_table == "radial.tsv" &&
        (a_column == "PR" || a_column == "alphaR" ||
         a_column == "betaR_deriv" || a_column == "Pi"))
        return 1e-10;
    if (a_table == "objects.tsv" && (a_column[0] == 'K' || a_column == "Pi"))
        return 1e-10;
    if ((a_table == "single_ccz4.tsv" || a_table == "binary_ccz4.tsv") &&
        (a_column == "K" || a_column == "A11" || a_column == "A12" ||
         a_column == "A22" || a_column == "Aww" || a_column == "Pi"))
        return 1e-10;
    return 1e-12;
}
static std::vector<double>
compute_fixture_values(const std::string &a_table,
                       const EMSBH_trumpet_read &a_reader,
                       const std::vector<double> &a_values)
{
    if (a_table == "jets.tsv")
    {
        std::vector<double> computed_values = {a_values[0]};
        for (int j = 0; j < 3; ++j)
            for (int d = 0; d < 4; ++d)
                computed_values.push_back(
                    a_reader.m_1d_sol.get_chebyshev_value(j, a_values[0], d));
        return computed_values;
    }
    if (a_table == "radial.tsv")
    {
        const auto radial_vars =
            a_reader.m_1d_sol.compute_radial_vars(a_values[0]);
        return {a_values[0],
                radial_vars.s,
                radial_vars.r,
                radial_vars.delta,
                radial_vars.phi,
                radial_vars.X,
                radial_vars.lapse,
                radial_vars.shift_R,
                radial_vars.K_amplitude,
                radial_vars.Pi,
                radial_vars.E_R,
                radial_vars.dP_dR,
                radial_vars.dlapse_dR,
                radial_vars.dshift_dR};
    }
    if (a_table == "objects.tsv")
    {
        const auto adm_vars = a_reader.compute_ems_adm_vars(
            a_values[3], a_values[4], a_values[5], a_values[0], 0, a_values[1]);
        std::vector<double> computed_values(a_values.begin(),
                                            a_values.begin() + 6);
        for (int j = 0; j < 3; ++j)
            for (int i = 0; i < 3; ++i)
                computed_values.push_back(adm_vars.gamma[i][j]);
        for (int j = 0; j < 3; ++j)
            for (int i = 0; i < 3; ++i)
                computed_values.push_back(adm_vars.K[i][j]);
        for (int i = 0; i < 3; ++i)
            computed_values.push_back(adm_vars.shift[i]);
        computed_values.push_back(adm_vars.phi);
        computed_values.push_back(adm_vars.Pi);
        for (int i = 0; i < 3; ++i)
            computed_values.push_back(adm_vars.E[i]);
        for (int i = 0; i < 3; ++i)
            computed_values.push_back(adm_vars.B[i]);
        computed_values.push_back(adm_vars.lapse);
        computed_values.push_back(adm_vars.chi);
        computed_values.push_back(adm_vars.boost_metric_det);
        computed_values.push_back(adm_vars.boost_normal_factor);
        return computed_values;
    }
    const bool is_binary = a_table == "binary_ccz4.tsv";
    const auto ccz4_vars =
        is_binary ? a_reader.compute_binary_bh_vars(a_values[2], a_values[3], 1,
                                                    a_values[0], a_values[1])
                  : a_reader.compute_single_bh_vars(
                        a_values[3], a_values[4], a_values[0], 0, a_values[1]);
    std::vector<double> computed_values(a_values.begin(),
                                        a_values.begin() + (is_binary ? 4 : 5));
    auto tail = ccz4_components(ccz4_vars);
    computed_values.insert(computed_values.end(), tail.begin(), tail.end());
    return computed_values;
}
static fixture_table_t check_fixture_table(const std::string &a_fixture_dir,
                                           const std::string &a_fixture_name,
                                           const EMSBH_trumpet_read &a_reader)
{
    std::ifstream in(a_fixture_dir + "/" + a_fixture_name);
    if (!in)
        throw std::runtime_error("missing table " + a_fixture_name);
    std::string line;
    std::getline(in, line);
    auto names = split_tsv_row(line);
    if (a_fixture_name == "single_ccz4.tsv" ||
        a_fixture_name == "binary_ccz4.tsv")
    {
        names.resize(a_fixture_name == "single_ccz4.tsv" ? 5 : 4);
        names.insert(names.end(), ccz4_component_names.begin(),
                     ccz4_component_names.end());
    }
    fixture_table_t t{a_fixture_name, 0.0, "", 0, 0};
    while (std::getline(in, line))
    {
        if (line.empty())
            throw std::runtime_error("empty table row");
        auto expected = parse_tsv_values(line),
             actual =
                 compute_fixture_values(a_fixture_name, a_reader, expected);
        if (actual.size() != expected.size() || names.size() != actual.size())
            throw std::runtime_error("table width mismatch: " + a_fixture_name);
        for (size_t i = 0; i < expected.size(); ++i)
        {
            std::uint64_t bits;
            std::memcpy(&bits, &actual[i], sizeof(bits));
            t.fingerprint = (t.fingerprint ^ bits) * 1099511628211ULL;
            const double error =
                std::abs(actual[i] - expected[i]) / (1 + std::abs(expected[i]));
            if (!std::isfinite(error))
                throw std::runtime_error("nonfinite result: " + a_fixture_name +
                                         "/" + names[i]);
            if (error > t.max)
            {
                t.max = error;
                t.worst = names[i];
            }
            if (error > fixture_tolerance(a_fixture_name, names[i]))
            {
                std::ostringstream msg;
                msg << std::setprecision(17) << a_fixture_name << " row "
                    << (t.rows + 1) << " column " << names[i]
                    << " expected=" << expected[i] << " actual=" << actual[i]
                    << " normalized error=" << error;
                throw std::runtime_error(msg.str());
            }
            ++t.columns;
        }
        ++t.rows;
    }
    return t;
}
int main(int argc, char **argv)
{
    try
    {
        const std::string dir = argc > 1 ? argv[1] : "fixtures";
        EMSBH_params_t params{};
        params.star_centre = {0, 0};
        params.data_path = dir + "/reference.trumpet";
        CouplingFunction::params_t coupling{12.566370614359172, 0, 0, -20};
        EMSBH_trumpet_read reader(params, coupling, 1.0, 1.0, 0);
        reader.compute_1d_solution();
        std::cout << std::scientific << std::setprecision(4);
        for (const char *name : {"jets.tsv", "radial.tsv", "objects.tsv",
                                 "single_ccz4.tsv", "binary_ccz4.tsv"})
        {
            const fixture_table_t fixture_result =
                check_fixture_table(dir, name, reader);
            std::cout << name << ": rows=" << fixture_result.rows
                      << " columns=" << fixture_result.columns
                      << " max=" << fixture_result.max << " at "
                      << (fixture_result.worst.empty() ? "all exact"
                                                       : fixture_result.worst)
                      << "\n";
            if (std::string(name) == "objects.tsv" ||
                std::string(name) == "single_ccz4.tsv")
            {
                const std::uint64_t baseline =
                    std::string(name) == "objects.tsv" ? 0xe1107d452abc4e4bULL
                                                       : 0xc30cf22126f14d61ULL;
                if (fixture_result.fingerprint != baseline)
                    throw std::runtime_error(std::string(name) +
                                             " changed single-object bits");
                std::cout << name << " bit fingerprint=" << std::hex
                          << fixture_result.fingerprint << std::dec << '\n';
            }
        }
        const std::vector<std::pair<std::string, std::string>> malformed_files =
            {{"wrong-magic.trumpet", "unknown trumpet schema"},
             {"unknown-convention.trumpet", "unsupported convention"},
             {"missing-key.trumpet", "missing required metadata"},
             {"duplicate-key.trumpet", "duplicate key nu"},
             {"nan.trumpet", "nonfinite coefficient"},
             {"truncated-block.trumpet", "malformed S block"},
             {"wrong-charge.trumpet", "inconsistent native charge"},
             {"wrong-endpoint.trumpet", "failed trumpet endpoint identity"},
             {"wrong-nu.trumpet", "failed trumpet endpoint identity"}};
        for (const auto &malformed_file : malformed_files)
        {
            std::string reason;
            if (reader.m_1d_sol.check_file(dir + "/" + malformed_file.first,
                                           reason))
                throw std::runtime_error(malformed_file.first +
                                         " was accepted");
            if (reason.find(malformed_file.second) == std::string::npos)
                throw std::runtime_error(malformed_file.first + " expected '" +
                                         malformed_file.second + "', got '" +
                                         reason + "'");
            std::cout << malformed_file.first << ": rejected (" << reason
                      << ")\n";
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "FAIL: " << e.what() << '\n';
        return 1;
    }
}
