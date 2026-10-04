#ifndef T24_OBSERVER_HPP
#define T24_OBSERVER_HPP
#include "DisjointBoxLayout.H"
#include "LayoutIterator.H"
#include <cstdlib>
#include <fstream>
#include <iomanip>

// Private build only. Read the global layout; never touch an evolved FAB.
inline void t24_record_layout(const DisjointBoxLayout &grids, int level,
                              double time, double dx)
{
    const char *path = std::getenv("T24_LAYOUT_FILE");
    if (!path) return;
    std::ofstream out(path, std::ios::app);
    out << std::setprecision(17);
    for (LayoutIterator it = grids.layoutIterator(); it.ok(); ++it)
    {
        const Box &b = grids[it()];
        out << "postRegrid," << level << ',' << time << ',' << dx << ','
            << it().intCode() << ',' << grids.procID(it()) << ','
            << b.smallEnd(0) << ',' << b.smallEnd(1) << ','
            << b.bigEnd(0) << ',' << b.bigEnd(1) << '\n';
    }
    if (!out) MayDay::Error("T24 layout observer write failed");
}
#endif
