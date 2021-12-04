/*
 * \class TrackFilterForTauReg
 *
 * Select pixel tracks and calo jets around L1 tau for the regional tau reconstruction
 *
 */

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"

class TrackFilterForTauReg : public edm::stream::EDProducer<> {
public:
  explicit TrackFilterForTauReg(const edm::ParameterSet& cfg)
      : l1TausToken_(consumes<trigger::TriggerFilterObjectWithRefs>(cfg.getParameter<edm::InputTag>("L1Taus"))),
        pixelTracksToken_(consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("pixelTracks"))),
        caloJetsToken_(consumes<edm::View<reco::CaloJet>>(cfg.getParameter<edm::InputTag>("caloJets"))),
        maxDeltaR2_pixelTrack_(std::pow(cfg.getParameter<double>("maxDeltaR_pixelTrack"), 2)),
        maxDeltaR2_caloJet_(std::pow(cfg.getParameter<double>("maxDeltaR_caloJet"), 2)),
        minPt_caloJet_(cfg.getParameter<double>("minPt_caloJet"))
  {
      produces<reco::TrackCollection>("PixelTracks");
      produces<reco::CaloJetCollection>("CaloJets");
  }
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("L1Taus", edm::InputTag("hltL1sTauVeryBigOR"))->setComment("L1 taus");
    desc.add<edm::InputTag>("pixelTracks", edm::InputTag("hltPixelTracks"))->setComment("pixel tracks");
    desc.add<edm::InputTag>("caloJets", edm::InputTag("hltAK4CaloJetsPFEt5"))->setComment("calo jets");
    desc.add<double>("maxDeltaR_pixelTrack", 0.8)->setComment("cone size around the L1 tau seed");
    desc.add<double>("maxDeltaR_caloJet", 0.4)->setComment("cone size around the L1 tau seed");
    desc.add<double>("minPt_caloJet", 20)->setComment("min pt for calo jets");
    descriptions.addWithDefaultLabel(desc);
  }

private:
    void produce(edm::Event& event, const edm::EventSetup&) override
    {
        l1t::TauVectorRef l1Taus;
        const auto& l1TriggeredTaus = event.get(l1TausToken_);
        l1TriggeredTaus.getObjects(trigger::TriggerL1Tau, l1Taus);

        const auto& tracks = event.get(pixelTracksToken_);
        auto selectedTracks = std::make_unique<reco::TrackCollection>();
        for(const auto& track : tracks) {
            for(const auto& tau : l1Taus) {
                if(reco::deltaR2(track.momentum(), tau->polarP4()) < maxDeltaR2_pixelTrack_) {
                    selectedTracks->emplace_back(track);
                    break;
                }
            }
        }
        event.put(std::move(selectedTracks), "PixelTracks");

        const auto& caloJets = event.get(caloJetsToken_);
        auto selectedCaloJets = std::make_unique<reco::CaloJetCollection>();
        for(const auto& caloJet : caloJets) {
            for(const auto& tau : l1Taus) {
                if(caloJet.pt() > minPt_caloJet_
                        && reco::deltaR2(caloJet.polarP4(), tau->polarP4()) < maxDeltaR2_caloJet_) {
                    selectedCaloJets->emplace_back(caloJet);
                    break;
                }
            }
        }
        event.put(std::move(selectedCaloJets), "CaloJets");
    }

private:
  const edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> l1TausToken_;
  const edm::EDGetTokenT<reco::TrackCollection> pixelTracksToken_;
  const edm::EDGetTokenT<edm::View<reco::CaloJet>> caloJetsToken_;
  const double maxDeltaR2_pixelTrack_, maxDeltaR2_caloJet_, minPt_caloJet_;
};

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TrackFilterForTauReg);
