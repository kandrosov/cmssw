#ifndef _TripletFilter_h_
#define _TripletFilter_h_

#include <vector>

#include "DataFormats/GeometryVector/interface/LocalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "RecoPixelVertexing/PixelLowPtUtilities/interface/ClusterShapeHitFilter.h"
#include "RecoPixelVertexing/PixelLowPtUtilities/interface/HitInfo.h"

namespace edm { class EventSetup; }
class TrackingRecHit;
class TrackerTopology;
class SiPixelClusterShapeCache;

class TripletFilter {
public:
    TripletFilter(const edm::EventSetup& es);
    template<typename Vector> // LocalVector or GlobalVector
    bool checkTrack(const std::vector<const TrackingRecHit*>& recHits, const std::vector<Vector>& dirs,
                    const TrackerTopology* tTopo, const SiPixelClusterShapeCache& clusterShapeCache)
    {
        bool ok = true;

        auto dir = dirs.begin();
        for(const auto& recHit : recHits) {
            const SiPixelRecHit* pixelRecHit = dynamic_cast<const SiPixelRecHit*>(recHit);
            if(!pixelRecHit->isValid()) {
                ok = false;
                break;
            }

            if(!theFilter->isCompatible(*pixelRecHit, *dir, clusterShapeCache)) {
                LogTrace("MinBiasTracking") << "  [TripletFilter] clusShape problem"
                                            << HitInfo::getInfo(*recHit,tTopo);
                ok = false;
                break;
            }

            ++dir;
        }

        return ok;
    }

private:
    const cluster_shape::ClusterShapeHitFilter* theFilter;
};

#endif

