#include "RecoPixelVertexing/PixelLowPtUtilities/interface/ClusterShapeHitFilter.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"

#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"

#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"

#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"
#include "Geometry/CommonTopologies/interface/RectangularStripTopology.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "MagneticField/Engine/interface/MagneticField.h"

#include "CondFormats/SiPixelObjects/interface/SiPixelLorentzAngle.h"

#include "CondFormats/SiStripObjects/interface/SiStripLorentzAngle.h"
#include "DataFormats/SiStripCluster/interface/SiStripClusterTools.h"


#include <fstream>
#include<cassert>

namespace cluster_shape {

ClusterShapeHitFilter::ClusterShapeHitFilter(const TrackerGeometry* theTracker_,
                                             const MagneticField* theMagneticField_,
                                             const SiPixelLorentzAngle* theSiPixelLorentzAngle_,
                                             const SiStripLorentzAngle* theSiStripLorentzAngle_,
                                             const std::string* use_PixelShapeFile_) :
    theTracker(theTracker_),
    theMagneticField(theMagneticField_),
    theSiPixelLorentzAngle(theSiPixelLorentzAngle_),
    theSiStripLorentzAngle(theSiStripLorentzAngle_),
    pixelLimits(*use_PixelShapeFile_),
    stripLimits("RecoPixelVertexing/PixelLowPtUtilities/data/stripShape.par"),
    cutOnPixelCharge_(false), cutOnStripCharge_(false),
    cutOnPixelShape_(false), cutOnStripShape_(false)
{
    fillPixelData();
}

void PixelLimitsCollection::Load(const std::string& shapeFile)
{
    edm::FileInPath fileInPath(shapeFile);
    std::ifstream inFile(fileInPath.fullPath());

    while(!inFile.eof()) {
        int part,dx,dy;

        inFile >> part; // 0 or 1
        inFile >> dx;   // 0 to 10
        inFile >> dy;   // 0 to 15 ...

        const PixelKeys key(part,dx,dy);
        auto & pl = limits[key];

        for(int b = 0; b<2 ; b++) // branch
            for(int d = 0; d<2 ; d++) // direction
                for(int k = 0; k<2 ; k++) // lower and upper
                    inFile >> pl.data[b][d][k];

        double f;
        int d;

        inFile >> f; // density
        inFile >> d; // points
        inFile >> f; // density
        inFile >> d; // points
    }

    LogTrace("MinBiasTracking|ClusterShapeHitFilter") << " pixel-cluster-shape filter loaded";
}

void StripLimitsCollection::Load(const std::string& shapeFile)
{
    // Load strip
    edm::FileInPath fileInPath(shapeFile);
    std::ifstream inFile(fileInPath.fullPath());

    while(!inFile.eof()) {
        int dx;
        inFile >> dx;

        StripKeys key(dx);
        auto & sl = limits[key];

        for(int b = 0; b<2 ; b++) // branch
            for(int k = 0; k<2 ; k++) // lower and upper
                inFile >> sl.data[b][k];
    }

    LogTrace("MinBiasTracking|ClusterShapeHitFilter") << " strip-cluster-width filter loaded";
}

void ClusterShapeHitFilter::fillPixelData()
{
    //barrel
    for (auto det : theTracker->detsPXB()) {
        // better not to fail..
        const PixelGeomDetUnit* pixelDet = dynamic_cast<const PixelGeomDetUnit*>(det);
        assert(pixelDet);
        PixelData & pd = pixelData[pixelDet->geographicalId()];
        pd.det = pixelDet;
        pd.part=0;
        pd.cotangent=getCotangent(pixelDet);
        pd.drift=getDrift(pixelDet);
    }

    //endcap
    for (auto det : theTracker->detsPXF()) {
        // better not to fail..
        const PixelGeomDetUnit* pixelDet = dynamic_cast<const PixelGeomDetUnit*>(det);
        assert(pixelDet);
        PixelData& pd = pixelData[pixelDet->geographicalId()];
        pd.det = pixelDet;
        pd.part=1;
        pd.cotangent=getCotangent(pixelDet);
        pd.drift=getDrift(pixelDet);
    }
}

std::pair<float,float> ClusterShapeHitFilter::getCotangent(const PixelGeomDetUnit * pixelDet) const
{
    std::pair<float,float> cotangent;

    cotangent.first  = pixelDet->surface().bounds().thickness() /
                     pixelDet->specificTopology().pitch().first;
    cotangent.second = pixelDet->surface().bounds().thickness() /
                     pixelDet->specificTopology().pitch().second;

    return cotangent;
}

float ClusterShapeHitFilter::getCotangent(const StripGeomDetUnit* stripDet, const LocalPoint& pos) const
{
    // FIXME may be problematic in case of RadialStriptolopgy
    return stripDet->surface().bounds().thickness() /
         stripDet->specificTopology().localPitch(pos);
}

std::pair<float,float> ClusterShapeHitFilter::getDrift(const PixelGeomDetUnit* pixelDet) const
{
    LocalVector lBfield = pixelDet->surface().toLocal(
                          theMagneticField->inTesla(pixelDet->surface().position()));
    double theTanLorentzAnglePerTesla = theSiPixelLorentzAngle->getLorentzAngle(
                                        pixelDet->geographicalId().rawId());
    std::pair<float,float> dir;
    dir.first = - theTanLorentzAnglePerTesla * lBfield.y();
    dir.second = theTanLorentzAnglePerTesla * lBfield.x();
    return dir;
}

float ClusterShapeHitFilter::getDrift(const StripGeomDetUnit* stripDet) const
{
    LocalVector lBfield = stripDet->surface().toLocal(theMagneticField->inTesla(
                          stripDet->surface().position()));
    double theTanLorentzAnglePerTesla = theSiStripLorentzAngle->getLorentzAngle(
                                        stripDet->geographicalId().rawId());
    float dir = theTanLorentzAnglePerTesla * lBfield.y();
    return dir;
}

bool ClusterShapeHitFilter::isNormalOriented(const GeomDetUnit * geomDet) const
{
    if(geomDet->type().isBarrel()) { // barrel
        float perp0 = geomDet->toGlobal( Local3DPoint(0.,0.,0.) ).perp();
        float perp1 = geomDet->toGlobal( Local3DPoint(0.,0.,1.) ).perp();
        return (perp1 > perp0);
    } else { // endcap
        float rot = geomDet->toGlobal( LocalVector (0.,0.,1.) ).z();
        float pos = geomDet->toGlobal( Local3DPoint(0.,0.,0.) ).z();
        return (rot * pos > 0);
    }
}

bool ClusterShapeHitFilter::getSizes(const SiPixelRecHit& recHit, const LocalVector& ldir,
                                     const SiPixelClusterShapeCache& clusterShapeCache, int& part,
                                     ClusterData::ArrayType& meas, std::pair<float,float>& pred,
                                     const PixelData* ipd, ExtendedResult* outEx) const
{
    // Get detector
    const PixelData & pd = getpd(recHit, ipd);

    // Get shape information
    const SiPixelClusterShapeData& data = clusterShapeCache.get(recHit.cluster(), pd.det);
    bool usable = (data.isStraight() && data.isComplete());

    // Usable?
    //if(usable)
    {
        part = pd.part;

        // Predicted size
        pred.first  = ldir.x() / ldir.z();
        pred.second = ldir.y() / ldir.z();

        SiPixelClusterShapeData::Range sizeRange = data.size();
        if(sizeRange.first->second < 0)
          pred.second = - pred.second;

        meas.clear();
        assert(meas.capacity() >= std::distance(sizeRange.first, sizeRange.second));
        for(auto s=sizeRange.first; s != sizeRange.second; ++s) {
          meas.push_back_unchecked(*s);
        }
        if(sizeRange.first->second < 0) {
          for(auto& s: meas)
            s.second = -s.second;
        }

        // Take out drift
        std::pair<float,float> const & drift = pd.drift;
        pred.first  += drift.first;
        pred.second += drift.second;

        // Apply cotangent
        std::pair<float,float> const & cotangent = pd.cotangent;
        pred.first  *= cotangent.first;
        pred.second *= cotangent.second;
    }

    if(outEx) {
        outEx->part = pd.part;
        outEx->drift = {{ pd.drift.first, pd.drift.second }};
        outEx->cotangent = {{ pd.cotangent.first, pd.cotangent.second }};
        outEx->isStraight = data.isStraight();
        outEx->isComplete = data.isComplete();
        outEx->hasBigPixelsOnlyInside = data.hasBigPixelsOnlyInside();
        outEx->usable = usable;
        const auto sizes = data.size();
        for(auto s = sizes.first; s != sizes.second; ++s) {
            outEx->clusterSizes.push_back({{ s->first, s->second }});
        }
        outEx->pred = {{ pred.first, pred.second }};
    }

    // Usable?
    return usable;
}

bool ClusterShapeHitFilter::isCompatibleImpl(const SiPixelRecHit& recHit, const LocalVector& ldir,
                                             const SiPixelClusterShapeCache& clusterShapeCache, const PixelData* ipd,
                                             ExtendedResult* outEx) const
{
    // Get detector
    if (cutOnPixelCharge_ && (!checkClusterCharge(recHit.geographicalId(), *(recHit.cluster()), ldir))) return false;
    if (!cutOnPixelShape_) return true;

    const PixelData & pd = getpd(recHit, ipd);
    int part;
    ClusterData::ArrayType meas;
    std::pair<float,float> pred;

    if(getSizes(recHit, ldir, clusterShapeCache, part,meas, pred, &pd, outEx)) {
        if(outEx) {
            outEx->hasPredInsideLimits = false;
            outEx->hasInvalidKey = false;
            for(const auto& m : meas) {
                const PixelKeys key(part, m.first, m.second);
                LimitResult limitResult;
                limitResult.limitKey = key;
                limitResult.limitKeyValid = key.isValid();
                if(key.isValid()) {
                    limitResult.predInsideLimits = pixelLimits[key].isInside(pred);
                    if(limitResult.predInsideLimits)
                        outEx->hasPredInsideLimits = true;
                    for(unsigned i = 0; i < 2; ++i) {
                        const auto limit = pixelLimits[key].data[i];
                        limitResult.limits[i].low = {{ limit[0][0], limit[1][0] }};
                        limitResult.limits[i].high = {{ limit[0][1], limit[1][1] }};
                    }
                } else {
                    outEx->hasInvalidKey = true;
                }
                outEx->limits.push_back(limitResult);
            }
        }

        for(const auto& m: meas) {
            PixelKeys key(part, m.first, m.second);
            if (!key.isValid()) return true; // FIXME original logic
            if (pixelLimits[key].isInside(pred)) return true;
        }
        // none of the choices worked
        return false;
    }
    // not usable
    return true;
}


bool ClusterShapeHitFilter::isCompatible(const SiPixelRecHit& recHit, const LocalVector& ldir,
                                         const SiPixelClusterShapeCache& clusterShapeCache,
                                         const PixelData* ipd, ExtendedResult* outEx) const
{
    const bool compatible = isCompatibleImpl(recHit, ldir, clusterShapeCache, ipd, outEx);
    if(outEx) {
        outEx->detId = recHit.geographicalId();
        outEx->compatible = compatible;
        outEx->localPos = recHit.localPosition();
        outEx->localDir = ldir;
    }
    return compatible;
}

bool ClusterShapeHitFilter::isCompatible(const SiPixelRecHit & recHit, const GlobalVector & gdir,
                                         const SiPixelClusterShapeCache& clusterShapeCache,
                                         const PixelData* ipd, ExtendedResult* outEx) const
{
    // Get detector
    const PixelData & pd = getpd(recHit,ipd);
    LocalVector ldir =pd.det->toLocal(gdir);
    return isCompatible(recHit, ldir, clusterShapeCache, &pd, outEx);
}

bool ClusterShapeHitFilter::getSizes(DetId id, const SiStripCluster& cluster, const LocalPoint &lpos,
                                     const LocalVector & ldir, int & meas, float & pred) const
{
    // Get detector
    const StripGeomDetUnit* stripDet = dynamic_cast<const StripGeomDetUnit*> (theTracker->idToDet(id));
    // Measured width
    meas = cluster.amplitudes().size();

    // Usable?
    int fs = cluster.firstStrip();
    int ns = stripDet->specificTopology().nstrips();
    // bool usable = (fs > 1 && fs + meas - 1 < ns);
    bool usable = (fs >= 1 && fs + meas - 1 <= ns);

    // Usable?
    //if(usable)
    {
        // Predicted width
        pred = ldir.x() / ldir.z();

        // Take out drift
        float drift = getDrift(stripDet);
        pred += drift;

        // Apply cotangent
        pred *= getCotangent(stripDet,lpos);
    }

    return usable;
}   

bool ClusterShapeHitFilter::isCompatible(DetId detId, const SiStripCluster& cluster, const LocalPoint& lpos,
                                         const LocalVector& ldir) const
{
    int meas;
    float pred;

    if (cutOnStripCharge_ && (!checkClusterCharge(detId, cluster, ldir))) return false;
    if (!cutOnStripShape_) return true;

    if(getSizes(detId, cluster, lpos, ldir, meas, pred)) {
        StripKeys key(meas);
        if (key.isValid())
            return stripLimits[key].isInside(pred);
    }

    // Not usable or no limits
    return true;
}

bool ClusterShapeHitFilter::isCompatible(DetId detId, const SiStripCluster& cluster, const GlobalPoint &gpos,
                                         const GlobalVector & gdir) const
{
    const GeomDet *det = theTracker->idToDet(detId);
    LocalVector ldir = det->toLocal(gdir);
    LocalPoint  lpos = det->toLocal(gpos);
    // now here we do the transformation
    lpos -= ldir * lpos.z()/ldir.z();
    return isCompatible(detId, cluster, lpos, ldir);
}

bool ClusterShapeHitFilter::isCompatible(DetId detId, const SiStripCluster & cluster, const GlobalVector & gdir) const
{
    return isCompatible(detId, cluster, theTracker->idToDet(detId)->toLocal(gdir));
}

bool ClusterShapeHitFilter::checkClusterCharge(DetId detId, const SiStripCluster& cluster,
                                               const LocalVector & ldir) const
{
    return siStripClusterTools::chargePerCM(detId, cluster, ldir) >  minGoodStripCharge_;
}

bool ClusterShapeHitFilter::checkClusterCharge(DetId detId, const SiPixelCluster& cluster,
                                               const LocalVector & ldir) const
{
    return siStripClusterTools::chargePerCM(detId, cluster, ldir) >  minGoodPixelCharge_;
}

} // namespace cluster_shape

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Utilities/interface/typelookup.h"
TYPELOOKUP_DATA_REG(cluster_shape::ClusterShapeHitFilter);

