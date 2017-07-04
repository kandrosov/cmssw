#include "RecoPixelVertexing/PixelTriplets/interface/IntermediateHitTriplets.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include "RecoPixelVertexing/PixelLowPtUtilities/interface/LowPtClusterShapeSeedComparitor.h"

#include <vector>

namespace RecoPixelVertexing_PixelTriplets {
  struct dictionary {
    IntermediateHitTriplets iht;
    edm::Wrapper<IntermediateHitTriplets> wiht;
    LowPtClusterShapeSeedComparitor::ExtendedResult css_debug;
    LowPtClusterShapeSeedComparitor::ExtendedResultCollection css_debug_vec;
    edm::Wrapper<LowPtClusterShapeSeedComparitor::ExtendedResultCollection> css_debug_wrap;
  };
}
