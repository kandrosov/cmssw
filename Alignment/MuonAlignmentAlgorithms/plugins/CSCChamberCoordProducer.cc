#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"

#include "TFile.h"
#include "TTree.h"

class CSCChamberCoordProducer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit CSCChamberCoordProducer(const edm::ParameterSet& ps)
      : inputFile_(ps.getUntrackedParameter<std::string>("inputFile")),
        outputFile_(ps.getUntrackedParameter<std::string>("outputFile")),
        treeName_(ps.getUntrackedParameter<std::string>("treeName", "csc_layer_ttree")),
        cscGeometryToken_(esConsumes<CSCGeometry, MuonGeometryRecord>()),
        globalGeometryToken_(esConsumes<GlobalTrackingGeometry, GlobalTrackingGeometryRecord>()) {}

  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void beginJob() override;
  void endJob() override;

private:
  std::string inputFile_;
  std::string outputFile_;
  std::string treeName_;

  edm::ESGetToken<CSCGeometry, MuonGeometryRecord> cscGeometryToken_;
  edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> globalGeometryToken_;

  TFile* inputFileHandle_ = nullptr;
  TFile* outputFileHandle_ = nullptr;
  TTree* inputTree_ = nullptr;
  TTree* outputTree_ = nullptr;

  Int_t endcap_ = 0;
  Int_t station_ = 0;
  Int_t ring_ = 0;
  Int_t chamber_ = 0;
  Float_t hit_x_[6]{};
  Float_t hit_y_[6]{};
  Float_t res_x_[6]{};
  Float_t res_y_[6]{};
  Float_t hit_chamber_x_[6]{};
  Float_t hit_chamber_y_[6]{};
};

void CSCChamberCoordProducer::beginJob() {
  inputFileHandle_ = TFile::Open(inputFile_.c_str());
  if (!inputFileHandle_ || inputFileHandle_->IsZombie()) {
    throw cms::Exception("CSCChamberCoordProducer") << "Unable to open input file: " << inputFile_;
  }

  inputTree_ = dynamic_cast<TTree*>(inputFileHandle_->Get(treeName_.c_str()));
  if (!inputTree_) {
    throw cms::Exception("CSCChamberCoordProducer") << "Unable to find tree '" << treeName_ << "' in " << inputFile_;
  }

  outputFileHandle_ = TFile::Open(outputFile_.c_str(), "RECREATE");
  if (!outputFileHandle_ || outputFileHandle_->IsZombie()) {
    throw cms::Exception("CSCChamberCoordProducer") << "Unable to create output file: " << outputFile_;
  }

  inputTree_->SetBranchAddress("endcap", &endcap_);
  inputTree_->SetBranchAddress("station", &station_);
  inputTree_->SetBranchAddress("ring", &ring_);
  inputTree_->SetBranchAddress("chamber", &chamber_);
  inputTree_->SetBranchAddress("hit_x", hit_x_);
  inputTree_->SetBranchAddress("hit_y", hit_y_);
  inputTree_->SetBranchAddress("res_x", res_x_);
  inputTree_->SetBranchAddress("res_y", res_y_);

  outputTree_ = inputTree_->CloneTree(0);
  if (!outputTree_) {
    throw cms::Exception("CSCChamberCoordProducer") << "Failed to clone input tree.";
  }
  inputTree_->CopyAddresses(outputTree_);
  outputTree_->Branch("hit_chamber_x", hit_chamber_x_, "hit_chamber_x[6]/F");
  outputTree_->Branch("hit_chamber_y", hit_chamber_y_, "hit_chamber_y[6]/F");
  outputTree_->SetDirectory(outputFileHandle_);
}

void CSCChamberCoordProducer::analyze(const edm::Event&, const edm::EventSetup& iSetup) {
  const CSCGeometry* cscGeometry = &iSetup.getData(cscGeometryToken_);
  const GlobalTrackingGeometry* globalGeometry = &iSetup.getData(globalGeometryToken_);

  for (Long64_t i = 0; i < inputTree_->GetEntries(); ++i) {
    inputTree_->GetEntry(i);

    for (int layer = 0; layer < 6; ++layer) {
      hit_chamber_x_[layer] = -999.;
      hit_chamber_y_[layer] = -999.;

      if (hit_x_[layer] < -900.) {
        continue;
      }

      CSCDetId layerId(endcap_, station_, ring_, chamber_, layer + 1);
      CSCDetId chamberId(endcap_, station_, ring_, chamber_, 0);

      const auto* layerDet = cscGeometry->idToDet(layerId);
      const auto* chamberDet = globalGeometry->idToDet(chamberId);
      if (!layerDet || !chamberDet) {
        continue;
      }

      const auto globalPoint = layerDet->toGlobal(Geom::LocalPoint(hit_x_[layer], hit_y_[layer], 0.));
      const auto chamberLocalPoint = chamberDet->surface().toLocal(globalPoint);
      hit_chamber_x_[layer] = chamberLocalPoint.x();
      hit_chamber_y_[layer] = chamberLocalPoint.y();
    }

    outputTree_->Fill();
  }
}

void CSCChamberCoordProducer::endJob() {
  outputFileHandle_->cd();
  outputTree_->Write();
  outputFileHandle_->Close();
  inputFileHandle_->Close();
}

DEFINE_FWK_MODULE(CSCChamberCoordProducer);
