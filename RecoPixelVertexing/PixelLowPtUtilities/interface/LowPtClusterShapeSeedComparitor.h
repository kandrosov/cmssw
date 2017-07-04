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
        DetId detId;
        const SiPixelRecHit* recHit;
        GlobalPoint globalPos;
        GlobalVector globalDir;
        bool hasPixelRecHit = false, hasValidPixelRecHit = false;
        cluster_shape::ClusterShapeHitFilter::ExtendedResult filterResult;
    };

    struct ExtendedResult {
        bool compatible = false, hasGlobalDir = false;
        std::vector<ClusterMatchResult> clusterMatchResults;
    };

    using ExtendedResultCollection = std::vector<ExtendedResult>;

    LowPtClusterShapeSeedComparitor(const edm::ParameterSet& ps, edm::ConsumesCollector& iC);
    virtual void init(const edm::Event& e, const edm::EventSetup& es) override;
    virtual bool compatible(const SeedingHitSet& hits) const override { return compatibleEx(hits, nullptr); }
    virtual bool compatible(const TrajectoryStateOnSurface&,
                            SeedingHitSet::ConstRecHitPointer) const override { return true; }
    virtual bool compatible(const SeedingHitSet&, const GlobalTrajectoryParameters&,
                            const FastHelix&) const override { return true; }

    bool compatibleEx(const SeedingHitSet& hits, ExtendedResult* outEx = nullptr) const;

private:
    edm::ESHandle<cluster_shape::ClusterShapeHitFilter> theShapeFilter;
    edm::ESHandle<TrackerTopology> theTTopo;
    edm::EDGetTokenT<SiPixelClusterShapeCache> thePixelClusterShapeCacheToken;
    edm::Handle<SiPixelClusterShapeCache> thePixelClusterShapeCache;
    std::string theShapeFilterLabel_;
};

#endif

