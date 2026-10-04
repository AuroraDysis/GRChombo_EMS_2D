#include "IntVectSet.H"
#include "BoxIterator.H"
#include <fstream>
#include <iostream>

// Probe the grow/minus-own-box/iteration portion of point-transfer ghost sets
// on stored layouts. This does not cover clipping or neighbour subtraction.
// No numerical state is read, and no production method is patched.
int main(int argc, char **argv)
{
    if (argc != 2) return 64;
    std::ifstream in(argv[1]);
    int level, x0, y0, x1, y1;
    unsigned boxes = 0;
    while (in >> level >> x0 >> y0 >> x1 >> y1)
    {
        const Box box(IntVect(D_DECL(x0,y0,0)), IntVect(D_DECL(x1,y1,0)));
        std::cout << "T24_SPARSE_BOX " << level << ' ' << box << std::endl;
        IntVectSet ghosts(grow(box,3));
        ghosts -= box;
        unsigned count = 0;
        for (IVSIterator it(ghosts); it.ok(); ++it) ++count;
        if (count != grow(box,3).numPts()-box.numPts()) return 65;
        ++boxes;
    }
    std::cout << "T24_SPARSE_COMPLETE boxes=" << boxes << std::endl;
    return 0;
}
