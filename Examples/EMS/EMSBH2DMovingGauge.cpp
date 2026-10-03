#include "EMSBH2DLevel.hpp"
#include "BoxLoops.hpp"
#include "ComputePack.hpp"
#include "TraceARemovalCartoon.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "MovingPunctureGauge.hpp"
#include "CCZ4Cartoon.hpp"
#include "SetValue.hpp"
#include "BoxIterator.H"
#include <typeinfo>
#define EMS_RHS_METHOD specificEvalRHSMoving
#define EMS_RHS_GAUGE MovingPunctureGauge
#include "EMSBH2DRHS.impl.hpp"
