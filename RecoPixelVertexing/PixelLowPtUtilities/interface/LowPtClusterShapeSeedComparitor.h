#ifndef _LowPtClusterShapeSeedComparitor_h_
#define _LowPtClusterShapeSeedComparitor_h_

#include "RecoTracker/TkSeedingLayers/interface/SeedComparitor.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelClusterShapeCache.h"
#include "RecoPixelVertexing/PixelLowPtUtilities/interface/ClusterShapeHitFilter.h"

class TrackerTopology;

namespace edm { class ParameterSet; class EventSetup; }

class LowPtClusterShapeSeedComparitor : public SeedComparitor {
public:

    struct ClusterMatchResult {
        size_t cluster_address;
        DetId detId;
        GlobalVector globalPos, globalDir;
        bool hasPixelRecHit = false, hasValidPixellRecHit = false;
        ClusterShapeHitFilter::ExtendedResult filterResult;
    };

    struct ExtendedResult {
        bool compatible = false, hasGlobalDir = false;
        std::vector<ClusterMatchResult> clusterMatchResults;
    };

    LowPtClusterShapeSeedComparitor(const edm::ParameterSet& ps, edm::ConsumesCollector& iC);
    virtual void init(const edm::Event& e, const edm::EventSetup& es);
    virtual bool compatible(const SeedingHitSet& hits) const override { return compatibleEx(hits, nullptr); }
    virtual bool compatible(const TrajectoryStateOnSurface&,
                            SeedingHitSet::ConstRecHitPointer) const override { return true; }
    virtual bool compatible(const SeedingHitSet&, const GlobalTrajectoryParameters&,
                            const FastHelix&) const override { return true; }

    bool compatibleEx(const SeedingHitSet& hits, ExtendedResult* outEx = nullptr) const;

private:
    /// something
    edm::ESHandle<ClusterShapeHitFilter> theShapeFilter;
    edm::ESHandle<TrackerTopology> theTTopo;
    edm::EDGetTokenT<SiPixelClusterShapeCache> thePixelClusterShapeCacheToken;
    edm::Handle<SiPixelClusterShapeCache> thePixelClusterShapeCache;
    std::string theShapeFilterLabel_;
};

#endif

