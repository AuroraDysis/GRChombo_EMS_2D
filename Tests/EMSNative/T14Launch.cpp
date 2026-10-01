// Frozen production main and level implementation; only the opt-in Tests
// factory selects dense initialization tag storage at levels >=13.
#include "T14DenseTags.hpp"
#define DefaultLevelFactory T14LevelFactory
#include "../../Examples/EMS/Main_EMSBH2DBH.cpp"
#undef DefaultLevelFactory
