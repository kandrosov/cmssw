#ifndef CaloTowerCreator_CaloTowerFromL1TCreatorForTauHLT_h
#define CaloTowerCreator_CaloTowerFromL1TCreatorForTauHLT_h

/** \class CaloTowerFromL1TCreatorForTauHLT
 *
 * Framework module that produces a collection
 * of calo towers in the region of interest for Tau HLT reconnstruction,
 * depending on tau type trigger:
 *                   Tau1 - take location of 1st L1 Tau
 *                   Tau2 - take location of 2nd L1 Tau; if does not exists,
 *                          take location of 1st Calo Tower
 *                   ETau - take L1 Tau candidate which is not collinear
 *                          to HLT (or L1) electron candidate.
 *
 * \author A. Nikitenko. IC.   based on L. Lista and J. Mans
 *
 */

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/L1Trigger/interface/Tau.h"

class CaloTowerFromL1TCreatorForTauHLT : public edm::global::EDProducer<> {
public:
    CaloTowerFromL1TCreatorForTauHLT(const edm::ParameterSet&);
    virtual void produce(edm::StreamID, edm::Event& evt, const edm::EventSetup&) const override;
    static void fillDescriptions(edm::ConfigurationDescriptions& desc);

private:
    const edm::EDGetTokenT<CaloTowerCollection> mtowers_token;
    const edm::EDGetTokenT<l1t::TauBxCollection> mTauTrigger_token;
    const int mBX;
    const double mCone, mConeSquare;
    const double mEtThreshold;
    const double mEThreshold;
    const int selectedTauId;
};

#endif
