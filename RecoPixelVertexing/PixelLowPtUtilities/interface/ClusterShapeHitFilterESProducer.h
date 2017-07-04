#ifndef RecoPixelVertexing_PixelLowPtUtilities_ClusterShapeHitFilterESProducer_H
#define RecoPixelVertexing_PixelLowPtUtilities_ClusterShapeHitFilterESProducer_H


// -*- C++ -*-
//
// Package:    ClusterShapeHitFilterESProducer
// Class:      ClusterShapeHitFilterESProducer
// 
/**\class ClusterShapeHitFilterESProducer ClusterShapeHitFilterESProducer.h
 * TrackingTools/ClusterShapeHitFilterESProducer/src/ClusterShapeHitFilterESProducer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jean-Roch Vlimant
//         Created:  Fri Sep 28 18:07:52 CEST 2007
//
//

#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "RecoPixelVertexing/PixelLowPtUtilities/interface/ClusterShapeHitFilter.h"

class ClusterShapeHitFilterESProducer : public edm::ESProducer {
public:
    using ReturnType = std::unique_ptr<cluster_shape::ClusterShapeHitFilter>;
    using Record = TrajectoryFilter::Record;

    ClusterShapeHitFilterESProducer(const edm::ParameterSet&);
    ReturnType produce(const Record &);

private:
    std::string use_PixelShapeFile;
    float minGoodPixelCharge_, minGoodStripCharge_;
    bool cutOnPixelCharge_, cutOnStripCharge_;
    bool cutOnPixelShape_, cutOnStripShape_;
};

#endif
