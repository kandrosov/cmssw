#include "RecoPixelVertexing/PixelLowPtUtilities/interface/ClusterShapeHitFilterESProducer.h"

#include "RecoPixelVertexing/PixelLowPtUtilities/interface/ClusterShapeHitFilter.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "RecoTracker/Record/interface/CkfComponentsRecord.h"
#include "RecoLocalTracker/SiStripClusterizer/interface/ClusterChargeCut.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

ClusterShapeHitFilterESProducer::ClusterShapeHitFilterESProducer(const edm::ParameterSet& cfg) :
    minGoodPixelCharge_(0), minGoodStripCharge_((clusterChargeCut(cfg))),
    cutOnPixelCharge_(false), cutOnStripCharge_(minGoodStripCharge_ > 0),
    cutOnPixelShape_(cfg.exists("doPixelShapeCut") ? cfg.getParameter<bool>("doPixelShapeCut") : true),
    cutOnStripShape_(cfg.exists("doStripShapeCut") ? cfg.getParameter<bool>("doStripShapeCut") : true)
{
    use_PixelShapeFile = cfg.exists("PixelShapeFile")
            ? cfg.getParameter<std::string>("PixelShapeFile")
            : "RecoPixelVertexing/PixelLowPtUtilities/data/pixelShape.par";

    const std::string componentName = cfg.getParameter<std::string>("ComponentName");
    edm::LogInfo("ClusterShapeHitFilterESProducer") << " with name: " << componentName;
    setWhatProduced(this, componentName);
}


ClusterShapeHitFilterESProducer::ReturnType ClusterShapeHitFilterESProducer::produce(const Record &iRecord)
{
    // get all from SiStripLorentzAngle (why not!)

    // Retrieve magnetic field
    edm::ESHandle<MagneticField> field;
    iRecord.getRecord<TkStripCPERecord>().getRecord<IdealMagneticFieldRecord>().get(field);

    // Retrieve geometry
    edm::ESHandle<TrackerGeometry> geo;
    iRecord.getRecord<TkStripCPERecord>().getRecord<TrackerDigiGeometryRecord>().get(geo);

    // Retrieve pixel Lorentz
    edm::ESHandle<SiPixelLorentzAngle> pixel;
    iRecord.getRecord<TkPixelCPERecord>().getRecord<SiPixelLorentzAngleRcd>().get(pixel);

    // Retrieve strip Lorentz
    edm::ESHandle<SiStripLorentzAngle> strip;
    iRecord.getRecord<TkStripCPERecord>().getRecord<SiStripLorentzAngleDepRcd>().get(strip);

    // Produce the filter using the plugin factory
    ClusterShapeHitFilterESProducer::ReturnType aFilter(new cluster_shape::ClusterShapeHitFilter(
            geo.product(), field.product(), pixel.product(), strip.product(), &use_PixelShapeFile));
    aFilter->setShapeCuts(cutOnPixelShape_, cutOnStripShape_);
    aFilter->setChargeCuts(cutOnPixelCharge_, minGoodPixelCharge_, cutOnStripCharge_, minGoodStripCharge_);
    return aFilter;
}
