#ifndef _ClusterShapeHitFilter_h_
#define _ClusterShapeHitFilter_h_

#include "TrackingTools/TrajectoryFiltering/interface/TrajectoryFilter.h"

#include "DataFormats/GeometryVector/interface/LocalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelClusterShapeCache.h"

#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"

#include "RecoPixelVertexing/PixelLowPtUtilities/interface/ClusterData.h"

#include <utility>
#include <unordered_map>
#include <cstring>

namespace edm { class EventSetup; }

class TrackerGeometry;
class MagneticField;
class SiPixelLorentzAngle;
class SiStripLorentzAngle;
class PixelGeomDetUnit;
class StripGeomDetUnit;

namespace cluster_shape {
class PixelKeys {    
public:
    static constexpr int offset_endcap_dy = 5;
    static constexpr int offset_endcap_dx = 10;
    static constexpr int N_endcap = 55;
    static constexpr int N_barrel = 137;
    static constexpr int N = N_barrel + N_endcap;

    using KeyType = unsigned char;

    PixelKeys(int part, int dx, int dy) :
        key(part==0 ? barrelPacking(dx,dy) : endcapPacking(dx,dy)) {}
  
    static KeyType endcapPacking(int dx, int dy)
    {
        if(dx < 0 || dy < 0 || dx > 10 || dy > 4) return N;
        return N_barrel + dx * offset_endcap_dy + dy;  // max 11*5 = 55
    }

    static KeyType barrelPacking(int dx, int dy)
    {
        if(dx < 0 || dy < 0 || dx > 10 || dy > 15) return N;
        if(dx < 8) return dx * 16 + dy;  // max 8*16=128
        if(dy > 2) return N;
        return 128 + (dx - 8) * 3 + dy; // max = 128+9 = 137
    }

    bool isValid() const { return key < N; }
    operator unsigned int() const { return key; }
    bool operator<(const PixelKeys& right) const { return key < right.key; }

private:
    KeyType key;
};

class StripKeys {
public:
    static constexpr int N = 40;

    using KeyType = unsigned char;

    StripKeys(int width) : key(width > 0 ? width - 1 : N) {}
    bool isValid() const { return key < N; }
    operator unsigned int() const { return key; }
    bool operator<(const StripKeys& right) const { return key < right.key; }

private:
    KeyType key;  // max 40;
};

struct PixelLimits {
    PixelLimits()
    {
        // init to make sure inside is true;
        auto limit = data[0];
        limit[0][0] = -10e10;
        limit[0][1] =  10e10;
        limit[1][0] = -10e10;
        limit[1][1] =  10e10;
        limit = data[1];
        limit[0][0] = -10e10;
        limit[0][1] =  10e10;
        limit[1][0] = -10e10;
        limit[1][1] =  10e10;
    }

    bool isInside( const std::pair<float,float> & pred) const
    {
        auto limit = data[0];
        bool one = pred.first > limit[0][0] && pred.first < limit[0][1]
                && pred.second > limit[1][0] && pred.second < limit[1][1];
        limit = data[1];
        bool two = pred.first > limit[0][0] && pred.first < limit[0][1]
                && pred.second > limit[1][0] && pred.second < limit[1][1];
        return one || two;
    }

    float data[2][2][2];
};

struct StripLimits {
    StripLimits()
    {
        data[0][0] = -10e10;
        data[0][1] =  10e10;
        data[1][0] = -10e10;
        data[1][1] =  10e10;
    }

    bool isInside(float pred) const
    {
        auto limit = data[0];
        bool one = pred > limit[0] && pred < limit[1];
        limit = data[1];
        bool two = pred > limit[0] && pred < limit[1];
        return one || two;
    }

    float data[2][2];
};

class PixelLimitsCollection {
public:
    explicit PixelLimitsCollection(const std::string& shapeFile) { Load(shapeFile); }
    const PixelLimits& operator[](size_t n) const { return limits.at(n); }
    size_t size() const { return limits.size(); }

private:
    void Load(const std::string& shapeFile);
    std::array<PixelLimits, PixelKeys::N + 1> limits; // [2][2][2]
};

class StripLimitsCollection {
public:
    explicit StripLimitsCollection(const std::string& shapeFile) { Load(shapeFile); }
    const StripLimits& operator[](size_t n) const { return limits.at(n); }
    size_t size() const { return limits.size(); }

private:
    void Load(const std::string& shapeFile);
    std::array<StripLimits, StripKeys::N + 1> limits; // [2][2]
};

class ClusterShapeHitFilter {
public:
    struct PixelData {
        const PixelGeomDetUnit * det;
        unsigned int part;
        std::pair<float, float> drift;
        std::pair<float, float> cotangent;
    };

    template<typename T>
    using Quantity2D = std::array<T, 2>;

    struct Limit2D {
        Quantity2D<float> low, high;
    };

    struct LimitResult {
        unsigned limitKey;
        bool limitKeyValid;
        bool predInsideLimits;
        std::array<Limit2D, 2> limits;
    };

    struct ExtendedResult {
        bool compatible = false;
        LocalPoint localPos;
        LocalVector localDir;
        DetId detId;
        bool isStraight = false, isComplete = false, hasBigPixelsOnlyInside = false, usable = false;
        unsigned part = 0;
        Quantity2D<float> drift;
        Quantity2D<float> cotangent;
        Quantity2D<float> pred;
        std::vector<Quantity2D<int>> clusterSizes;
        std::vector<LimitResult> limits;
        bool hasPredInsideLimits = false, hasInvalidKey = false;
    };

    ClusterShapeHitFilter(const TrackerGeometry* theTracker_,
                          const MagneticField* theMagneticField_,
                          const SiPixelLorentzAngle* theSiPixelLorentzAngle_,
                          const SiStripLorentzAngle* theSiStripLorentzAngle_,
                          const std::string* use_PixelShapeFile_);
 
    void setShapeCuts(bool cutOnPixelShape, bool cutOnStripShape)
    {
        cutOnPixelShape_ = cutOnPixelShape;
        cutOnStripShape_ = cutOnStripShape;
    }

    void setChargeCuts(bool cutOnPixelCharge, float minGoodPixelCharge, bool cutOnStripCharge, float minGoodStripCharge)
    {
        cutOnPixelCharge_ = cutOnPixelCharge;
        minGoodPixelCharge_= minGoodPixelCharge;
        cutOnStripCharge_ = cutOnStripCharge;
        minGoodStripCharge_= minGoodStripCharge;
    }

    bool getSizes(const SiPixelRecHit& recHit, const LocalVector& ldir,
                  const SiPixelClusterShapeCache& clusterShapeCache, int& part, ClusterData::ArrayType& meas,
                  std::pair<float,float> & predr, const PixelData* pd = nullptr, ExtendedResult* outEx = nullptr) const;
    bool isCompatible(const SiPixelRecHit& recHit, const LocalVector& ldir,
                      const SiPixelClusterShapeCache& clusterShapeCache, const PixelData* pd = nullptr,
                      ExtendedResult* outEx = nullptr) const;
    bool isCompatible(const SiPixelRecHit& recHit, const GlobalVector& gdir,
                      const SiPixelClusterShapeCache& clusterShapeCache, const PixelData* pd = nullptr,
                      ExtendedResult* outEx = nullptr) const;


    bool getSizes(DetId detId, const SiStripCluster& cluster, const LocalPoint& lpos, const LocalVector& ldir,
                  int& meas, float& pred) const;
    bool getSizes(const SiStripRecHit2D& recHit, const LocalPoint& lpos, const LocalVector& ldir, int& meas,
                  float& pred) const
    {
        return getSizes(recHit.geographicalId(), recHit.stripCluster(), lpos, ldir, meas, pred);
    }

    bool isCompatible(DetId detId, const SiStripCluster& cluster, const LocalPoint& lpos,
                      const LocalVector& ldir) const;
    bool isCompatible(DetId detId, const SiStripCluster& cluster, const LocalVector& ldir) const
    {
        return isCompatible(detId, cluster, LocalPoint(0,0,0), ldir);
    }

    bool isCompatible(DetId detId, const SiStripCluster& cluster, const GlobalPoint& gpos,
                      const GlobalVector& gdir) const;
    bool isCompatible(DetId detId, const SiStripCluster& cluster, const GlobalVector& gdir) const;
    bool isCompatible(const SiStripRecHit2D& recHit, const LocalPoint& lpos, const LocalVector& ldir) const
    {
        return isCompatible(recHit.geographicalId(), recHit.stripCluster(), lpos, ldir);
    }
    bool isCompatible(const SiStripRecHit2D& recHit, const LocalVector& ldir) const
    {
        return isCompatible(recHit.geographicalId(), recHit.stripCluster(), ldir);
    }
    bool isCompatible(const SiStripRecHit2D& recHit, const GlobalPoint& gpos, const GlobalVector& gdir) const
    {
        return isCompatible(recHit.geographicalId(), recHit.stripCluster(), gpos, gdir);
    }
    bool isCompatible(const SiStripRecHit2D& recHit, const GlobalVector& gdir) const
    {
        return isCompatible(recHit.geographicalId(), recHit.stripCluster(), gdir);
    }

private:
    const PixelData & getpd(const SiPixelRecHit& recHit, const PixelData* pd = nullptr) const
    {
        if (pd) return *pd;
        // Get detector
        DetId id = recHit.geographicalId();
        auto p = pixelData.find(id);
        return (*p).second;
    }

    void fillPixelData();

    std::pair<float,float> getCotangent(const PixelGeomDetUnit* pixelDet) const;
    float getCotangent(const StripGeomDetUnit* stripDet, const LocalPoint& p = LocalPoint(0,0,0)) const;

    std::pair<float,float> getDrift(const PixelGeomDetUnit* pixelDet) const;
    float getDrift(const StripGeomDetUnit* stripDet) const;

    bool isNormalOriented(const GeomDetUnit* geomDet) const;

    bool isCompatibleImpl(const SiPixelRecHit& recHit, const LocalVector& ldir,
                          const SiPixelClusterShapeCache& clusterShapeCache, const PixelData* pd,
                          ExtendedResult* outEx) const;


    const TrackerGeometry* theTracker;
    const MagneticField* theMagneticField;

    const SiPixelLorentzAngle* theSiPixelLorentzAngle;
    const SiStripLorentzAngle* theSiStripLorentzAngle;

    std::unordered_map<unsigned int, PixelData> pixelData;

    PixelLimitsCollection pixelLimits;
    StripLimitsCollection stripLimits;

    float theAngle[6];
    bool cutOnPixelCharge_, cutOnStripCharge_;
    float minGoodPixelCharge_, minGoodStripCharge_;
    bool cutOnPixelShape_, cutOnStripShape_;
    bool checkClusterCharge(DetId detId, const SiStripCluster& cluster, const LocalVector& ldir) const;
    bool checkClusterCharge(DetId detId, const SiPixelCluster& cluster, const LocalVector& ldir) const;
};

} // namespace cluster_shape

#endif
