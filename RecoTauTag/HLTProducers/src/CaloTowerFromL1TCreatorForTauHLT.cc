// makes CaloTowerCandidates from CaloTowers
// original author: L.Lista INFN, modifyed by: F.Ratnikov UMd
// Author for regionality A. Nikitenko
// Modified by S. Gennai

#include "DataFormats/RecoCandidate/interface/RecoCaloTowerCandidate.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoTauTag/HLTProducers/interface/CaloTowerFromL1TCreatorForTauHLT.h"
#include "Math/GenVector/VectorUtil.h"

CaloTowerFromL1TCreatorForTauHLT::CaloTowerFromL1TCreatorForTauHLT(const edm::ParameterSet& p) :
    mtowers_token(consumes<CaloTowerCollection>(p.getParameter<edm::InputTag>("towers"))),
    mTauTrigger_token(consumes<l1t::TauBxCollection>(p.getParameter<edm::InputTag>("TauTrigger"))),
    mBX(p.getParameter<int>("BX")),
    mCone(p.getParameter<double>("UseTowersInCone")),
    mConeSquare(mCone * mCone),
    mEtThreshold(p.getParameter<double>("minimumEt")),
    mEThreshold(p.getParameter<double>("minimumE")),
    selectedTauId(p.getParameter<int>("TauId")),
    legacyMode(selectedTauId > 0)
{
    if(legacyMode) {
        produces<CaloTowerCollection>();
        return;
    }

    for(size_t tauId = 0; tauId < NumberOfTauIndexes; ++tauId)
        produces<CaloTowerCollection>(GetCollectionName(tauId));
}

void CaloTowerFromL1TCreatorForTauHLT::produce(edm::StreamID sid, edm::Event& evt, const edm::EventSetup& stp) const
{
    TowerCollectionList tauTowers;
    splitTowerCollection(evt, tauTowers);

    if(legacyMode) {
        if(size_t(selectedTauId) >= tauTowers.size()) {
            evt.put(TowerCollectionPtr( new CaloTowerCollection ));
        } else {
            auto iter = tauTowers.begin();
            std::advance(iter, selectedTauId);
            evt.put(TowerCollectionPtr(std::move(*iter)));
        }
        return;
    }

    size_t tauId = 0;
    for(auto& cands : tauTowers)
        evt.put(std::move(cands), GetCollectionName(tauId++));

    for(; tauId < NumberOfTauIndexes; ++tauId)
        evt.put(TowerCollectionPtr( new CaloTowerCollection ), GetCollectionName(tauId));
}

void CaloTowerFromL1TCreatorForTauHLT::splitTowerCollection(const edm::Event& evt, TowerCollectionList& tauTowers) const
{
    edm::Handle<CaloTowerCollection> h_caloTowers;
    evt.getByToken(mtowers_token, h_caloTowers);
    if(!h_caloTowers.isValid()) {
        edm::LogWarning("MissingProduct") << "Calo tower collection not found.";
        return;
    }
    const auto& caloTowers = *h_caloTowers;

    edm::Handle<l1t::TauBxCollection> h_jetsgen;
    evt.getByToken(mTauTrigger_token, h_jetsgen);
    if(!h_jetsgen.isValid()) {
        edm::LogWarning("MissingProduct") << "L1Upgrade jet bx collection not found.";
        return;
    }
    const auto& jetsgen = *h_jetsgen;

    for (auto l1Jet = jetsgen.begin(mBX); l1Jet != jetsgen.end(mBX); ++l1Jet) {
        TowerCollectionPtr towers( new CaloTowerCollection );
        for(const auto& tower : caloTowers) {
            if(tower.et() < mEtThreshold || tower.energy() < mEThreshold) continue;
            const double deltaR2  = ROOT::Math::VectorUtil::DeltaR2(l1Jet->polarP4(), tower.polarP4());
            if(deltaR2 < mConeSquare)
                towers->push_back(tower);
        }
        tauTowers.push_back(std::move(towers));
    }
}

std::string CaloTowerFromL1TCreatorForTauHLT::GetCollectionName(size_t tauId)
{
    static const std::string NamePrefix = "Tau";
    std::ostringstream s_name;
    s_name << NamePrefix << tauId;
    return s_name.str();
}

void CaloTowerFromL1TCreatorForTauHLT::fillDescriptions( edm::ConfigurationDescriptions & desc )
{
    edm::ParameterSetDescription aDesc;

    aDesc.add<edm::InputTag>("TauTrigger", edm::InputTag("caloStage2Digis"))->setComment("L1 Tau collection for seeding");
    aDesc.add<edm::InputTag>("towers", edm::InputTag("towerMaker"))->setComment("Input tower collection");
    aDesc.add<double>("UseTowersInCone", 0.8)->setComment("Radius of cone around seed");
    aDesc.add<double>("minimumE", 0.8)->setComment("Minimum tower energy");
    aDesc.add<double>("minimumEt", 0.5)->setComment("Minimum tower ET");
    aDesc.add<int>("BX", 0)->setComment("Set bunch crossing; 0 = in time, -1 = previous, 1 = following");

    // Legacy
    aDesc.add<int>("TauId", -1)->setComment("Procude output only for the one selected L1 tau");
    aDesc.addUntracked<int>("verbose", 0)->setComment("Verbosity level");

    desc.add("CaloTowerFromL1TCreatorForTauHLT", aDesc);
    desc.setComment("Produce tower collection around L1 particle seed.");
}
