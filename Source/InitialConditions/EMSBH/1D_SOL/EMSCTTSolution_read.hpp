/* GRChombo - EMSCTT 1 file-only companion consumer. */
#ifndef EMSCTTSOLUTION_READ_HPP_
#define EMSCTTSOLUTION_READ_HPP_

#include "EMSTrumpetSolution_read.hpp"
#include "SHA256.hpp"
#include <array>
#include <map>
#include <string>
#include <vector>

class EMSCTTSolution_read
{
  public:
    struct binding_t
    {
        std::array<double, 2> mass, center, rapidity;
        std::array<double, 4> coupling;
    };
    struct correction_t
    {
        double logpsi;
        std::array<std::array<double, 3>, 3> C{};
    };
    bool check_file(const std::string &a_path,
                    const EMSTrumpetSolution_read &a_profile,
                    const binding_t &a_binding, std::string &a_reason);
    //! x is an offset from origin, translated BEFORE maps/Clenshaw. A positive
    //! panel selects a one-sided extension for fixture/interface checks only.
    correction_t evaluate(double x, double y, double z, int panel = 0,
                          double origin = 0) const;
    bool matches(double mass, double separation, double rapidity) const;
    const std::string &source_sha256() const { return m_source_sha256; }

  private:
    struct panel_t
    {
        std::string kind, boundary;
        double lo, hi, center, stretch, mu_lo, mu_hi;
        int nr, na;
        std::array<std::vector<double>, 5> coeff;
    };
    std::vector<panel_t> m_panels;
    std::map<std::string, double> m_values;
    std::string m_source_sha256;
    static bool number(const std::string &text, double &value);
    static bool integer(const std::string &text, int &value, int lo, int hi);
    static bool hash_string(const std::string &text);
    static double clenshaw(const double *a, int n, double x);
};
#include "EMSCTTSolution_read.impl.hpp"
#endif
