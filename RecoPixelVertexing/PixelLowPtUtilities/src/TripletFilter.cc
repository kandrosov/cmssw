#include "RecoPixelVertexing/PixelLowPtUtilities/interface/TripletFilter.h"

#include "RecoPixelVertexing/PixelLowPtUtilities/interface/ClusterShapeHitFilter.h"
#include "RecoPixelVertexing/PixelLowPtUtilities/interface/HitInfo.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"

#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

#include "RecoTracker/Record/interface/CkfComponentsRecord.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"

TripletFilter::TripletFilter(const edm::EventSetup& es)
{
    // Get cluster shape hit filter
    edm::ESHandle<cluster_shape::ClusterShapeHitFilter> shape;
    es.get<CkfComponentsRecord>().get("ClusterShapeHitFilter", shape);
    theFilter = shape.product();
}
