#include "RecoTauTag/HLTProducers/interface/L2TauJetsMerger.h"
#include "Math/GenVector/VectorUtil.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

L2TauJetsMerger::L2TauJetsMerger(const edm::ParameterSet& iConfig) :
    jetSrc( iConfig.getParameter<vtag>("JetSrc") ),
    mEt_Min( iConfig.getParameter<double>("EtMin") )
{
    for(const auto& source : jetSrc)
        jetSrc_token.push_back(consumes<reco::CaloJetCollection>(source));
    produces<reco::CaloJetCollection>();
}

// Getting the Collections of L2ReconstructedJets from L1Seeds and removing the collinear jets.
void L2TauJetsMerger::produce(edm::StreamID iSId, edm::Event& iEvent, const edm::EventSetup& iES) const
{
    static constexpr int TauPdgId = 15;
    static constexpr double DeltaR_cut = 0.1;
    static constexpr double DeltaR2_cut = DeltaR_cut * DeltaR_cut;
    static const auto ptOrderer = [](const reco::CaloJet& jet1 , const reco::CaloJet& jet2) {
        return jet1.pt() > jet2.pt();
    };
    using CaloJetSet = std::multiset<reco::CaloJet, decltype(ptOrderer)>;

    CaloJetSet allJets(ptOrderer);
    for(const auto& token : jetSrc_token) {
        edm::Handle<reco::CaloJetCollection> h_tauJets;
        iEvent.getByToken(token, h_tauJets);
        for(const auto& tau : *h_tauJets) {
            if(tau.et() <= mEt_Min) continue;
            reco::CaloJet jet(tau);
            jet.setPdgId(TauPdgId);
            allJets.insert(std::move(jet));
        }
    }

    // Removing the collinear jets.
    std::unique_ptr<reco::CaloJetCollection> l2Jets(new reco::CaloJetCollection);
    for(const auto& jet : allJets) {
        bool pass_deltaR = true;
        for(auto l2Jet_iter = l2Jets->begin(); pass_deltaR && l2Jet_iter != l2Jets->end(); ++l2Jet_iter) {
            const double deltaR2 = ROOT::Math::VectorUtil::DeltaR2(jet.polarP4(), l2Jet_iter->polarP4());
            pass_deltaR = deltaR2 > DeltaR2_cut;
        }
        if(pass_deltaR)
            l2Jets->push_back(jet);
    }

    iEvent.put(std::move(l2Jets));
}

void L2TauJetsMerger::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    edm::ParameterSetDescription desc;
    std::vector<edm::InputTag> inputTags;
    inputTags.push_back( edm::InputTag("hltAkIsoTau1Regional") );
    inputTags.push_back( edm::InputTag("hltAkIsoTau2Regional") );
    inputTags.push_back( edm::InputTag("hltAkIsoTau3Regional") );
    inputTags.push_back( edm::InputTag("hltAkIsoTau4Regional") );
    desc.add< std::vector<edm::InputTag> >("JetSrc",inputTags)->setComment("CaloJet collections to merge");
    desc.add<double>("EtMin",20.0)->setComment("Minimal ET of jet to merge");
    descriptions.setComment("Merges CaloJet collections removing duplicates");
    descriptions.add("L2TauJetsMerger",desc);
}
