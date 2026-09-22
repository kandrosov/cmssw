#include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentAlgorithmBase.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonTopologies/interface/GeomDet.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"

#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/DTSuperLayerId.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "TrackingTools/PatternTools/interface/Trajectory.h"

#include "Alignment/CommonAlignment/interface/Alignable.h"
#include "Alignment/CommonAlignment/interface/AlignableDetUnit.h"
#include "Alignment/CommonAlignment/interface/AlignableNavigator.h"
#include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentParameterStore.h"
#include "Alignment/MuonAlignment/interface/AlignableDTSuperLayer.h"
#include "Alignment/MuonAlignment/interface/AlignableDTChamber.h"
#include "Alignment/MuonAlignment/interface/AlignableDTStation.h"
#include "Alignment/MuonAlignment/interface/AlignableDTWheel.h"
#include "Alignment/MuonAlignment/interface/AlignableCSCChamber.h"
#include "Alignment/MuonAlignment/interface/AlignableCSCRing.h"
#include "Alignment/MuonAlignment/interface/AlignableCSCStation.h"
#include "Alignment/MuonAlignment/interface/AlignableMuon.h"

#include "Alignment/MuonAlignmentAlgorithms/interface/MuonResidualsFromTrack.h"
#include "Alignment/MuonAlignmentAlgorithms/interface/MuonResiduals6DOFFitter.h"
#include "Alignment/MuonAlignmentAlgorithms/interface/MuonResiduals5DOFFitter.h"
#include "Alignment/MuonAlignmentAlgorithms/interface/MuonResiduals6DOFrphiFitter.h"
#include "Alignment/MuonAlignmentAlgorithms/interface/MuonResidualsTwoBin.h"
#include "Alignment/MuonAlignmentAlgorithms/interface/MuonResidualsBigStructureFitter.h"

#include "Alignment/MuonAlignmentAlgorithms/interface/DTTTree.h"
#include "Alignment/MuonAlignmentAlgorithms/interface/CSCTTree.h"
#include "Alignment/MuonAlignmentAlgorithms/interface/FlatOccupancy.h"

#include "TrackingTools/Records/interface/DetIdAssociatorRecord.h"

#include "TFile.h"
#include "TTree.h"
#include "TStopwatch.h"

#include "Alignment/MuonAlignmentAlgorithms/interface/Tracer.hpp"
#include "CondFormats/Alignment/interface/AlignTransform.h"
#include "Alignment/CommonAlignment/interface/Utilities.h"

#include <map>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <iterator> 
#include <chrono>
#include <memory>
#include <random>

class MuonBigStructureAlignment : public AlignmentAlgorithmBase {
public:
  MuonBigStructureAlignment(const edm::ParameterSet& cfg, edm::ConsumesCollector& iC);
  ~MuonBigStructureAlignment() override;

  void initialize(const edm::EventSetup& iSetup,
                  AlignableTracker* alignableTracker,
                  AlignableMuon* alignableMuon,
                  AlignableExtras* extras,
                  AlignmentParameterStore* alignmentParameterStore) override;

  void startNewLoop() override {};
  void run(const edm::EventSetup& iSetup, const EventInfo& eventInfo) override {};
  void terminate(const edm::EventSetup& iSetup) override;

  FlatOccupancy myOccupancyMap;

private:
  bool numeric(std::string s);
  int number(std::string s);
  std::string chamberPrettyNameFromId(unsigned int idx);

  void readTmpFiles();
  void writeTmpFiles();

  void doGlobalAlignment();

  void fiducialCuts();
  void eraseNotSelectedResiduals();

  // tokens
  const edm::ESGetToken<CSCGeometry, MuonGeometryRecord> m_cscGeometryToken;
  const edm::ESGetToken<DTGeometry, MuonGeometryRecord> m_dtGeometryToken; // token added for accessing dt geometry; include DTGeometry.h
  const edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> m_globTackingToken;
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> m_MagFieldToken;
  const edm::ESGetToken<Propagator, TrackingComponentsRecord> m_propToken;
  const MuonResidualsFromTrack::BuilderToken m_builderToken;

  // configutarion paramenters:
  edm::InputTag m_muonCollectionTag;
  std::vector<std::string> m_reference;
  double m_minTrackPt;
  double m_maxTrackPt;
  double m_minTrackP;
  double m_maxTrackP;
  double m_maxDxy;
  int m_minTrackerHits;
  double m_maxTrackerRedChi2;
  bool m_allowTIDTEC;
  int m_minNCrossedChambers;
  int m_minDT13Hits;
  int m_minDT2Hits;
  int m_minCSCHits;
  std::string m_writeTemporaryFile;
  std::vector<std::string> m_readTemporaryFiles;
  bool m_doAlignment;
  int m_strategy;
  std::string m_residualsModel;
  int m_minAlignmentHits;
  bool m_twoBin;
  bool m_combineME11;
  bool m_weightAlignment;
  std::string m_reportFileName;
  double m_maxResSlopeY;
  bool m_createNtuple;
  double m_peakNSigma;
  int m_BFieldCorrection;
  bool m_doDT;
  bool m_doCSC;
  std::string m_useResiduals;

  // global alignment parameters
  bool m_doGlobalAlignment = true;
  std::string m_CSCEndcaps;
  std::string m_CSCRings;
  std::string m_CSCStations;
  std::string m_DTStations;
  std::string m_DTWheels;
  double m_min_pt = 0.0;
  double m_max_pt = 10e4;
  double m_errdef = 0.5;
  DTGeometry const* m_dtGeometry;
  CSCGeometry const* m_cscGeometry;

  //Layer Plots
  bool m_createLayerNtuple_DT;
  bool m_createLayerNtuple_CSC;

  // utility objects
  AlignableNavigator* m_alignableNavigator;
  AlignmentParameterStore* m_alignmentParameterStore;
  align::Alignables m_alignables;
  std::map<Alignable*, Alignable*> m_me11map;
  std::map<Alignable*, MuonResidualsTwoBin*> m_fitters;
  std::vector<unsigned int> m_indexes;
  std::map<unsigned int, MuonResidualsTwoBin*> m_fitterOrder;

  // counters
  long m_counter_events;
  long m_counter_tracks;
  long m_counter_trackmomentum;
  long m_counter_trackdxy;
  long m_counter_trackerhits;
  long m_counter_trackerchi2;
  long m_counter_trackertidtec;
  long m_counter_minchambers;
  long m_counter_totchambers;
  long m_counter_station123;
  long m_counter_station123valid;
  long m_counter_station123dt13hits;
  long m_counter_station123dt2hits;
  long m_counter_station123aligning;
  long m_counter_station4;
  long m_counter_station4valid;
  long m_counter_station4hits;
  long m_counter_station4aligning;
  long m_counter_csc;
  long m_counter_cscvalid;
  long m_counter_cschits;
  long m_counter_cscaligning;
  long m_counter_resslopey;

  // Ntuples for debug
  void bookNtuple();
  void bookNtupleLayers_DT();
  void bookNtupleLayers_CSC();
  TTree* m_ttree;
  TTree* m_ttree_DT_layers;
  TTree* m_ttree_CSC_layers;
  MuonResidualsFitter::MuonAlignmentTreeRow m_tree_row;
  DTLayerData layerData_DT;
  CSCLayerData layerData_CSC;

  bool m_debug;
};

MuonBigStructureAlignment::MuonBigStructureAlignment(const edm::ParameterSet& cfg, edm::ConsumesCollector& iC)
    : AlignmentAlgorithmBase(cfg, iC),
      m_cscGeometryToken(iC.esConsumes<edm::Transition::BeginRun>()),
      m_dtGeometryToken(iC.esConsumes<edm::Transition::BeginRun>()), // initialize dtGeometry token
      m_globTackingToken(iC.esConsumes()),
      m_MagFieldToken(iC.esConsumes()),
      m_propToken(iC.esConsumes(edm::ESInputTag("", "SteppingHelixPropagatorAny"))),
      m_builderToken(iC.esConsumes(MuonResidualsFromTrack::builderESInputTag())),
      m_muonCollectionTag(cfg.getParameter<edm::InputTag>("muonCollectionTag")),
      m_reference(cfg.getParameter<std::vector<std::string>>("reference")),
      m_minTrackPt(cfg.getParameter<double>("minTrackPt")),
      m_maxTrackPt(cfg.getParameter<double>("maxTrackPt")),
      m_minTrackP(cfg.getParameter<double>("minTrackP")),
      m_maxTrackP(cfg.getParameter<double>("maxTrackP")),
      m_maxDxy(cfg.getParameter<double>("maxDxy")),
      m_minTrackerHits(cfg.getParameter<int>("minTrackerHits")),
      m_maxTrackerRedChi2(cfg.getParameter<double>("maxTrackerRedChi2")),
      m_allowTIDTEC(cfg.getParameter<bool>("allowTIDTEC")),
      m_minNCrossedChambers(cfg.getParameter<int>("minNCrossedChambers")),
      m_minDT13Hits(cfg.getParameter<int>("minDT13Hits")),
      m_minDT2Hits(cfg.getParameter<int>("minDT2Hits")),
      m_minCSCHits(cfg.getParameter<int>("minCSCHits")),
      m_writeTemporaryFile(cfg.getParameter<std::string>("writeTemporaryFile")),
      m_readTemporaryFiles(cfg.getParameter<std::vector<std::string> >("readTemporaryFiles")),
      m_doAlignment(cfg.getParameter<bool>("doAlignment")),
      m_strategy(cfg.getParameter<int>("strategy")),
      m_residualsModel(cfg.getParameter<std::string>("residualsModel")),
      m_minAlignmentHits(cfg.getParameter<int>("minAlignmentHits")),
      m_twoBin(cfg.getParameter<bool>("twoBin")),
      m_combineME11(cfg.getParameter<bool>("combineME11")),
      m_weightAlignment(cfg.getParameter<bool>("weightAlignment")),
      m_reportFileName(cfg.getParameter<std::string>("reportFileName")),
      m_maxResSlopeY(cfg.getParameter<double>("maxResSlopeY")),
      m_createNtuple(cfg.getParameter<bool>("createNtuple")),
      m_peakNSigma(cfg.getParameter<double>("peakNSigma")),
      m_BFieldCorrection(cfg.getParameter<int>("bFieldCorrection")),
      m_doDT(cfg.getParameter<bool>("doDT")),
      m_doCSC(cfg.getParameter<bool>("doCSC")),
      m_useResiduals(cfg.getParameter<std::string>("useResiduals")),
      m_CSCEndcaps(cfg.getParameter<std::string>("cscEndcaps")),
      m_CSCRings(cfg.getParameter<std::string>("cscRings")),
      m_CSCStations(cfg.getParameter<std::string>("cscStations")),
      m_DTStations(cfg.getParameter<std::string>("dtStations")),
      m_DTWheels(cfg.getParameter<std::string>("dtWheels")),
      m_min_pt(cfg.getParameter<double>("minPt")),
      m_max_pt(cfg.getParameter<double>("maxPt")),
      m_errdef(cfg.getParameter<double>("errDef")),
      m_createLayerNtuple_DT(cfg.getParameter<bool>("createLayerNtupleDT")),
      m_createLayerNtuple_CSC(cfg.getParameter<bool>("createLayerNtupleCSC")){
  // alignment requires a TFile to provide plots to check the fit output
  // just filling the residuals lists does not
  // but we don't want to wait until the end of the job to find out that the TFile is missing
  if (m_doAlignment || m_createNtuple) {
    edm::Service<TFileService> fs;
    TFile& tfile = fs->file();
    tfile.ls();
  }

  m_ttree = nullptr;
  if (m_createNtuple)
    bookNtuple();
  //Layer Ntuples
  m_ttree_DT_layers = nullptr;   //By default the TTree is nullptr
  m_ttree_CSC_layers = nullptr;  //By default the TTree is nullptr
  if (m_createLayerNtuple_DT) {
    layerData_DT.doFill = true;
    bookNtupleLayers_DT();
  }
  if (m_createLayerNtuple_CSC) {
    layerData_CSC.doFill = true;
    bookNtupleLayers_CSC();
  }

  m_counter_events = 0;
  m_counter_tracks = 0;
  m_counter_trackmomentum = 0;
  m_counter_trackdxy = 0;
  m_counter_trackerhits = 0;
  m_counter_trackerchi2 = 0;
  m_counter_trackertidtec = 0;
  m_counter_minchambers = 0;
  m_counter_totchambers = 0;
  m_counter_station123 = 0;
  m_counter_station123valid = 0;
  m_counter_station123dt13hits = 0;
  m_counter_station123dt2hits = 0;
  m_counter_station123aligning = 0;
  m_counter_station4 = 0;
  m_counter_station4valid = 0;
  m_counter_station4hits = 0;
  m_counter_station4aligning = 0;
  m_counter_csc = 0;
  m_counter_cscvalid = 0;
  m_counter_cschits = 0;
  m_counter_cscaligning = 0;
  m_counter_resslopey = 0;

  //myOccupancyMap.LoadWeigths("/afs/cern.ch/work/l/lpernie/MuonAlign/WD/FixOccupancy/CMSSW_8_0_17/src/LayerPlots/Output_Occupancy/Output_weights.root"); //Place a root file containing the weight map for the 250 DTs as a function of X and Y, if you want flat occupancy (and decomment the weight part in MuonResidualsXDOFFitter.cc).
  myOccupancyMap.LoadWeigths(
      "");  //Place a root file containing the weight map for the 250 DTs as a function of X and Y, if you want flat occupancy (and decomment the weight part in MuonResidualsXDOFFitter.cc).
  m_debug = false;
}

MuonBigStructureAlignment::~MuonBigStructureAlignment() { delete m_alignableNavigator; }

void MuonBigStructureAlignment::bookNtuple() {
  edm::Service<TFileService> fs;
  m_ttree = fs->make<TTree>("mual_ttree", "mual_ttree");
  m_ttree->Branch("is_plus", &m_tree_row.is_plus, "is_plus/O");
  m_ttree->Branch("is_dt", &m_tree_row.is_dt, "is_dt/O");
  m_ttree->Branch("station", &m_tree_row.station, "station/b");
  m_ttree->Branch("ring_wheel", &m_tree_row.ring_wheel, "ring_wheel/B");
  m_ttree->Branch("sector", &m_tree_row.sector, "sector/b");
  m_ttree->Branch("res_x", &m_tree_row.res_x, "res_x/F");
  m_ttree->Branch("res_y", &m_tree_row.res_y, "res_y/F");
  m_ttree->Branch("res_slope_x", &m_tree_row.res_slope_x, "res_slope_x/F");
  m_ttree->Branch("res_slope_y", &m_tree_row.res_slope_y, "res_slope_y/F");
  m_ttree->Branch("pos_x", &m_tree_row.pos_x, "pos_x/F");
  m_ttree->Branch("pos_y", &m_tree_row.pos_y, "pos_y/F");
  m_ttree->Branch("angle_x", &m_tree_row.angle_x, "angle_x/F");
  m_ttree->Branch("angle_y", &m_tree_row.angle_y, "angle_y/F");
  m_ttree->Branch("pz", &m_tree_row.pz, "pz/F");
  m_ttree->Branch("pt", &m_tree_row.pt, "pt/F");
  m_ttree->Branch("q", &m_tree_row.q, "q/B");
  m_ttree->Branch("OccuWeight", &m_tree_row.OccuWeight, "OccuWeight/F");
  m_ttree->Branch("select", &m_tree_row.select, "select/O");
}

void MuonBigStructureAlignment::bookNtupleLayers_DT() {
  edm::Service<TFileService> fs;
  m_ttree_DT_layers = fs->make<TTree>("dt_layer_ttree", "dt_layer_ttree");
  m_ttree_DT_layers->Branch("charge", &(layerData_DT.charge), "charge/i");
  m_ttree_DT_layers->Branch("nEvent", &(layerData_DT.nEvent), "nEvent/i");
  m_ttree_DT_layers->Branch("nlayers", &(layerData_DT.nlayers), "nlayers/i");
  m_ttree_DT_layers->Branch("nDT", &(layerData_DT.nDT), "nDT/i");
  m_ttree_DT_layers->Branch("nCSC", &(layerData_DT.nDT), "nCSC/i");
  m_ttree_DT_layers->Branch("nTracker", &(layerData_DT.nTracker), "nTracker/i");

  m_ttree_DT_layers->Branch("wheel", &(layerData_DT.wheel), "wheel/b");
  m_ttree_DT_layers->Branch("station", &(layerData_DT.station), "station/b");
  m_ttree_DT_layers->Branch("sector", &(layerData_DT.sector), "sector/b");

  m_ttree_DT_layers->Branch("pt", &(layerData_DT.pt), "pt/F");
  m_ttree_DT_layers->Branch("pz", &(layerData_DT.pz), "pz/F");
  m_ttree_DT_layers->Branch("eta", &(layerData_DT.eta), "eta/F");
  m_ttree_DT_layers->Branch("phi", &(layerData_DT.phi), "phi/F");

  m_ttree_DT_layers->Branch("hit_x", &(layerData_DT.v_hitx), "hit_x[8]/F");
  m_ttree_DT_layers->Branch("hit_y", &(layerData_DT.v_hity), "hit_y[4]/F");
  m_ttree_DT_layers->Branch("track_x", &(layerData_DT.v_trackx), "track_x[8]/F");
  m_ttree_DT_layers->Branch("track_y", &(layerData_DT.v_tracky), "track_y[4]/F");
  m_ttree_DT_layers->Branch("track_y_x_layer", &(layerData_DT.v_tracky_x_layer), "track_y_x_layer[4]/F");
}

void MuonBigStructureAlignment::bookNtupleLayers_CSC() {
  edm::Service<TFileService> fs;
  m_ttree_CSC_layers = fs->make<TTree>("csc_layer_ttree", "csc_layer_ttree");
  m_ttree_CSC_layers->Branch("charge", &(layerData_CSC.charge), "charge/i");
  m_ttree_CSC_layers->Branch("nEvent", &(layerData_CSC.nEvent), "nEvent/i");
  m_ttree_CSC_layers->Branch("nlayers", &(layerData_CSC.nlayers), "nlayers/i");

  m_ttree_CSC_layers->Branch("nDT", &(layerData_CSC.nDT), "nDT/i");
  m_ttree_CSC_layers->Branch("nCSC", &(layerData_CSC.nCSC), "nCSC/i");
  m_ttree_CSC_layers->Branch("nTracker", &(layerData_CSC.nTracker), "nTracker/i");

  m_ttree_CSC_layers->Branch("endcap", &(layerData_CSC.endcap), "endcap/b");
  m_ttree_CSC_layers->Branch("station", &(layerData_CSC.station), "station/b");
  m_ttree_CSC_layers->Branch("ring", &(layerData_CSC.ring), "ring/b");
  m_ttree_CSC_layers->Branch("chamber", &(layerData_CSC.chamber), "chamber/b");

  m_ttree_CSC_layers->Branch("pt", &(layerData_CSC.pt), "pt/F");
  m_ttree_CSC_layers->Branch("pz", &(layerData_CSC.pz), "pz/F");
  m_ttree_CSC_layers->Branch("eta", &(layerData_CSC.eta), "eta/F");
  m_ttree_CSC_layers->Branch("phi", &(layerData_CSC.phi), "phi/F");

  m_ttree_CSC_layers->Branch("hit_x", &(layerData_CSC.v_hitx), "hit_x[6]/F");
  m_ttree_CSC_layers->Branch("hit_y", &(layerData_CSC.v_hity), "hit_y[6]/F");
  m_ttree_CSC_layers->Branch("res_x", &(layerData_CSC.v_resx), "res_x[6]/F");
  m_ttree_CSC_layers->Branch("res_y", &(layerData_CSC.v_resy), "res_y[6]/F");
}

bool MuonBigStructureAlignment::numeric(std::string s) { return s.length() == 1 && std::isdigit(s[0]); }

int MuonBigStructureAlignment::number(std::string s) {
  if (!numeric(s))
    assert(false);
  return atoi(s.c_str());
}

void MuonBigStructureAlignment::initialize(const edm::EventSetup& iSetup,
                                            AlignableTracker* alignableTracker,
                                            AlignableMuon* alignableMuon,
                                            AlignableExtras* extras,
                                            AlignmentParameterStore* alignmentParameterStore) {
  if (alignableMuon == nullptr)
    throw cms::Exception("MuonBigStructureAlignment") << "doMuon must be set to True" << std::endl;

  m_alignableNavigator = new AlignableNavigator(alignableMuon);
  m_alignmentParameterStore = alignmentParameterStore;
  m_alignables = m_alignmentParameterStore->alignables();

  int residualsModel;
  if (m_residualsModel == std::string("pureGaussian"))
    residualsModel = MuonResidualsFitter::kPureGaussian;
  else if (m_residualsModel == std::string("pureGaussian2D"))
    residualsModel = MuonResidualsFitter::kPureGaussian2D;
  else if (m_residualsModel == std::string("powerLawTails"))
    residualsModel = MuonResidualsFitter::kPowerLawTails;
  else if (m_residualsModel == std::string("ROOTVoigt"))
    residualsModel = MuonResidualsFitter::kROOTVoigt;
  else if (m_residualsModel == std::string("GaussPowerTails"))
    residualsModel = MuonResidualsFitter::kGaussPowerTails;
  else
    throw cms::Exception("MuonBigStructureAlignment")
        << "unrecognized residualsModel: \"" << m_residualsModel << "\"" << std::endl;

  int useResiduals;
  if (m_useResiduals == std::string("1111"))
    useResiduals = MuonResidualsFitter::k1111;
  else if (m_useResiduals == std::string("1110"))
    useResiduals = MuonResidualsFitter::k1110;
  else if (m_useResiduals == std::string("1100"))
    useResiduals = MuonResidualsFitter::k1100;
  else if (m_useResiduals == std::string("1000"))
    useResiduals = MuonResidualsFitter::k1000;
  else if (m_useResiduals == std::string("1010"))
    useResiduals = MuonResidualsFitter::k1010;
  else if (m_useResiduals == std::string("0010"))
    useResiduals = MuonResidualsFitter::k0010;
  else
    throw cms::Exception("MuonBigStructureAlignment")
        << "unrecognized useResiduals: \"" << m_useResiduals << "\"" << std::endl;

  const CSCGeometry* cscGeometry = &iSetup.getData(m_cscGeometryToken);
  m_dtGeometry = &iSetup.getData(m_dtGeometryToken);
  m_cscGeometry = cscGeometry;

  // set up the MuonResidualsFitters (which also collect residuals for fitting)
  m_me11map.clear();
  m_fitters.clear();
  m_indexes.clear();
  m_fitterOrder.clear();

  for (const auto& ali : m_alignables) {
    bool made_fitter = false;

    // fitters for DT
    if (ali->alignableObjectId() == align::AlignableDTChamber) {
      DTChamberId id(ali->geomDetId().rawId());

      if (id.station() == 4) {
        m_fitters[ali] = new MuonResidualsTwoBin(
            m_twoBin,
            new MuonResiduals5DOFFitter(residualsModel, m_minAlignmentHits, useResiduals, m_weightAlignment),
            new MuonResiduals5DOFFitter(residualsModel, m_minAlignmentHits, useResiduals, m_weightAlignment));
        made_fitter = true;
      } else {
        m_fitters[ali] = new MuonResidualsTwoBin(
            m_twoBin,
            new MuonResiduals6DOFFitter(residualsModel, m_minAlignmentHits, useResiduals, m_weightAlignment),
            new MuonResiduals6DOFFitter(residualsModel, m_minAlignmentHits, useResiduals, m_weightAlignment));
        made_fitter = true;
      }
    }

    // fitters for CSC
    else if (ali->alignableObjectId() == align::AlignableCSCChamber) {
      auto thisali = ali;
      CSCDetId id(ali->geomDetId().rawId());

      // take care of ME1/1a
      if (m_combineME11 && id.station() == 1 && id.ring() == 4) {
        CSCDetId pairid(id.endcap(), 1, 1, id.chamber());

        for (const auto& ali2 : m_alignables) {
          if (ali2->alignableObjectId() == align::AlignableCSCChamber && ali2->geomDetId().rawId() == pairid.rawId()) {
            thisali = ali2;
            break;
          }
        }
        m_me11map[ali] = thisali;  // points from each ME1/4 chamber to the corresponding ME1/1 chamber
      }

      if (thisali == ali)  // don't make fitters for ME1/4; they get taken care of in ME1/1
      {
        m_fitters[ali] = new MuonResidualsTwoBin(
            m_twoBin,
            new MuonResiduals6DOFrphiFitter(
                residualsModel, m_minAlignmentHits, useResiduals, cscGeometry, m_weightAlignment),
            new MuonResiduals6DOFrphiFitter(
                residualsModel, m_minAlignmentHits, useResiduals, cscGeometry, m_weightAlignment));
        made_fitter = true;
      }
    }

    else
      throw cms::Exception("MuonBigStructureAlignment")
          << "only DTChambers and CSCChambers can be aligned with this module" << std::endl;

    if (made_fitter) {
      m_fitters[ali]->setStrategy(m_strategy);

      int index = ali->geomDetId().rawId();
      m_indexes.push_back(index);
      m_fitterOrder[index] = m_fitters[ali];
    }
  }  // end loop over chambers chosen for alignment

  // cannonical order of fitters in the file
  std::sort(m_indexes.begin(), m_indexes.end());

  // de-weight all chambers but the reference
  const auto& all_DT_chambers = alignableMuon->DTChambers();
  const auto& all_CSC_chambers = alignableMuon->CSCChambers();

  align::Alignables reference;
  alignmentParameterStore->setAlignmentPositionError(all_DT_chambers, 100000000., 0.);
  alignmentParameterStore->setAlignmentPositionError(all_CSC_chambers, 100000000., 0.);
  alignmentParameterStore->setAlignmentPositionError(reference, 0., 0.);
}

void MuonBigStructureAlignment::terminate(const edm::EventSetup& iSetup) {
  bool m_debug = false;

  // one-time print-out
  std::cout << "Counters:" << std::endl
            << "COUNT{ events: " << m_counter_events << " }" << std::endl
            << "COUNT{  tracks: " << m_counter_tracks << " }" << std::endl
            << "COUNT{   trackppt: " << m_counter_trackmomentum << " }" << std::endl
            << "COUNT{    trackdxy: " << m_counter_trackdxy << " }" << std::endl
            << "COUNT{     trackerhits: " << m_counter_trackerhits << " }" << std::endl
            << "COUNT{      trackerchi2: " << m_counter_trackerchi2 << " }" << std::endl
            << "COUNT{       trackertidtec: " << m_counter_trackertidtec << " }" << std::endl
            << "COUNT{        minnchambers: " << m_counter_minchambers << " }" << std::endl
            << "COUNT{          totchambers: " << m_counter_totchambers << " }" << std::endl
            << "COUNT{            station123:             " << m_counter_station123 << " }" << std::endl
            << "COUNT{             station123valid:       " << m_counter_station123valid << " }" << std::endl
            << "COUNT{              station123dt13hits:   " << m_counter_station123dt13hits << " }" << std::endl
            << "COUNT{               station123dt2hits:   " << m_counter_station123dt2hits << " }" << std::endl
            << "COUNT{                station123aligning: " << m_counter_station123aligning << " }" << std::endl
            << "COUNT{                 resslopey: " << m_counter_resslopey << " }" << std::endl
            << "COUNT{            station4:            " << m_counter_station4 << " }" << std::endl
            << "COUNT{             station4valid:      " << m_counter_station4valid << " }" << std::endl
            << "COUNT{              station4hits:      " << m_counter_station4hits << " }" << std::endl
            << "COUNT{               station4aligning: " << m_counter_station4aligning << " }" << std::endl
            << "COUNT{            csc:            " << m_counter_csc << " }" << std::endl
            << "COUNT{             cscvalid:      " << m_counter_cscvalid << " }" << std::endl
            << "COUNT{              cschits:      " << m_counter_cschits << " }" << std::endl
            << "COUNT{               cscaligning: " << m_counter_cscaligning << " }" << std::endl
            << "That's all!" << std::endl;

  TStopwatch stop_watch;

  // collect temporary files
  if (!m_readTemporaryFiles.empty()) {
    stop_watch.Start();
    readTmpFiles();
    if (m_debug)
      std::cout << "readTmpFiles took " << stop_watch.CpuTime() << " sec" << std::endl;
    stop_watch.Stop();
  }

  if (m_doAlignment) {
    stop_watch.Start();
    fiducialCuts();
    if (m_debug)
      std::cout << "fiducialCuts took " << stop_watch.CpuTime() << " sec" << std::endl;
    stop_watch.Stop();
  }

  if (m_doAlignment) {
    stop_watch.Start();
    eraseNotSelectedResiduals();
    if (m_debug)
      std::cout << "eraseNotSelectedResiduals took " << stop_watch.CpuTime() << " sec" << std::endl;
    stop_watch.Stop();
  }

  if (m_doGlobalAlignment) {
    stop_watch.Start();
    doGlobalAlignment();
    if (m_debug)
      std::cout << "doGlobalAlignment took " << stop_watch.CpuTime() << " sec" << std::endl;
    stop_watch.Stop();
  }

  // write out the pseudontuples for a later job to collect
  if (m_writeTemporaryFile != std::string(""))
    writeTmpFiles();
  if (m_debug)
    std::cout << "end: MuonBigStructureAlignment::terminate()" << std::endl;
}

void MuonBigStructureAlignment::readTmpFiles() {
  for (std::vector<std::string>::const_iterator fileName = m_readTemporaryFiles.begin();
       fileName != m_readTemporaryFiles.end();
       ++fileName) {
    FILE* file;
    int size;
    file = fopen((*fileName).c_str(), "r");
    if (file == nullptr)
      throw cms::Exception("MuonBigStructureAlignment")
          << "file \"" << *fileName << "\" can't be opened (doesn't exist?)" << std::endl;

    fread(&size, sizeof(int), 1, file);
    if (int(m_indexes.size()) != size)
      throw cms::Exception("MuonBigStructureAlignment")
          << "file \"" << *fileName << "\" has " << size << " fitters, but this job has " << m_indexes.size()
          << " fitters (probably corresponds to the wrong alignment job)" << std::endl;

    int i = 0;
    for (std::vector<unsigned int>::const_iterator index = m_indexes.begin(); index != m_indexes.end(); ++index, ++i) {
      MuonResidualsTwoBin* fitter = m_fitterOrder[*index];
      unsigned int index_toread;
      fread(&index_toread, sizeof(unsigned int), 1, file);
      if (*index != index_toread)
        throw cms::Exception("MuonBigStructureAlignment")
            << "file \"" << *fileName << "\" has index " << index_toread << " at position " << i
            << ", but this job is expecting " << *index << " (probably corresponds to the wrong alignment job)"
            << std::endl;
      fitter->read(file, i);
    }

    fclose(file);
  }
}

void MuonBigStructureAlignment::doGlobalAlignment()
{
  using namespace std::chrono;
  // configure GPR fitter
  std::vector<std::string> options{ m_DTWheels, m_DTStations, m_CSCEndcaps, m_CSCRings, m_CSCStations };
  MuonResidualsBigStructureFitter big_structure_fitter(m_dtGeometry, m_cscGeometry, options, m_min_pt, m_max_pt, m_errdef);

  if (m_doCSC && !m_doDT && m_CSCEndcaps == "11")
  {
    Tracer::instance() << "fitting CSC endcaps separately!\n";

    big_structure_fitter.SetOption(2, "10");

    auto start = high_resolution_clock::now();
    big_structure_fitter.CopyData(m_fitters);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    std::cout << "Copying took " << duration.count() << " ms" << "\n";
    std::cout << "Track count " << big_structure_fitter.TrackCount() << "\n";
    Tracer::instance() << "Copying took " << duration.count() << " ms\n";

    Tracer::instance() << "Current track count: " << big_structure_fitter.TrackCount() << "\n";
    Tracer::instance() << "Number of chambers selected for alignment: " << big_structure_fitter.NumberOfChambers() << "\n";

    big_structure_fitter.CalcStats();

    std::ofstream file;
    file.open("gpr_params.csv");
    
    start = high_resolution_clock::now();
    bool gpr_fit_done = big_structure_fitter.Fit();
    stop = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(stop - start);

    if (gpr_fit_done)
    {
      Tracer::instance() << "Global Alignment of CSC Endcap 1 done!" << "\n";
      Tracer::instance() << "Global Alignment of CSC Endcap 1 took " << duration.count() << " ms" << "\n";
      Tracer::instance() << "Global Alignment of CSC Endcap 1 parameters: ";
      for (size_t i = 0; i < 6; ++i)
      {
        Tracer::instance() << big_structure_fitter.GetParamValue(i) << " ";
        file << big_structure_fitter.GetParamValue(i) << ",";
      }

      file << "\n";
      Tracer::instance() << "\n";
    }

    big_structure_fitter.ReleaseData();
    Tracer::instance() << "Pointers to data have been released\n";
    Tracer::instance() << "Current track count: " << big_structure_fitter.TrackCount() << "\n";
    Tracer::instance() << "Number of chambers selected for alignment: " << big_structure_fitter.NumberOfChambers() << "\n";

    big_structure_fitter.SetOption(2, "01");

    start = high_resolution_clock::now();
    big_structure_fitter.CopyData(m_fitters);
    stop = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(stop - start);
    std::cout << "Copying took " << duration.count() << " ms" << "\n";
    Tracer::instance() << "Copying took " << duration.count() << " ms\n";
    Tracer::instance() << "Current track count: " << big_structure_fitter.TrackCount() << "\n";
    Tracer::instance() << "Number of chambers selected for alignment: " << big_structure_fitter.NumberOfChambers() << "\n";

    big_structure_fitter.CalcStats();

    start = high_resolution_clock::now();
    gpr_fit_done = big_structure_fitter.Fit();
    stop = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(stop - start);

    if (gpr_fit_done)
    {
      Tracer::instance() << "Global Alignment of CSC Endcap 2 done!" << "\n";
      Tracer::instance() << "Global Alignment of CSC Endcap 2 took " << duration.count() << " ms" << "\n";
      Tracer::instance() << "Global Alignment of CSC Endcap 2 parameters: ";
      for (size_t i = 0; i < 6; ++i)
      {
        Tracer::instance() << big_structure_fitter.GetParamValue(i) << " ";
        file << big_structure_fitter.GetParamValue(i) << ",";
      }

      file << "\n";
      Tracer::instance() << "\n";
    }
    
    file.close();
  }

  if (m_doDT && !m_doCSC) 
  {
    Tracer::instance() << "Fitting DT barrel\n";
    auto start = high_resolution_clock::now();
    big_structure_fitter.CopyData(m_fitters);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    std::cout << "Copying took " << duration.count() << " ms" << "\n";
    Tracer::instance() << "Copying took " << duration.count() << " ms\n";

    Tracer::instance() << "Current track count: " << big_structure_fitter.TrackCount() << "\n";
    Tracer::instance() << "Number of chambers selected for alignment: " << big_structure_fitter.NumberOfChambers() << "\n";

    big_structure_fitter.CalcStats();

    std::ofstream file;
    file.open("gpr_params.csv");
    
    start = high_resolution_clock::now();
    bool gpr_fit_done = big_structure_fitter.Fit();
    stop = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(stop - start);

    if (gpr_fit_done)
    {
      Tracer::instance() << "Global Alignment of DT barrel done!" << "\n";
      Tracer::instance() << "Global Alignment of DT barrel took " << duration.count() << " ms" << "\n";
      Tracer::instance() << "Global Alignment of DT barrel parameters: ";
      for (size_t i = 0; i < 6; ++i)
      {
        Tracer::instance() << big_structure_fitter.GetParamValue(i) << " ";
        file << big_structure_fitter.GetParamValue(i) << ",";
      }

      file << "\n";
      Tracer::instance() << "\n";
    }
    file.close();
  }
}

void MuonBigStructureAlignment::writeTmpFiles() {
  FILE* file;
  file = fopen(m_writeTemporaryFile.c_str(), "w");
  int size = m_indexes.size();          // m_indexes is a vector of int. Each int is the Raw ID of the chamber.
  fwrite(&size, sizeof(int), 1, file);  //WRITE the size of m_indexes

  int i = 0;
  for (std::vector<unsigned int>::const_iterator index = m_indexes.begin(); index != m_indexes.end(); ++index, ++i) {
    MuonResidualsTwoBin* fitter = m_fitterOrder[*index];
    unsigned int index_towrite = *index;
    fwrite(&index_towrite, sizeof(unsigned int), 1, file);  // WRITE the RAW_ID of the chamber (one by one)
    fitter->write(
        file,
        i);  //All the other informations are written inside this method (you just need the file and the index to the element to write)
  }

  fclose(file);
}

void MuonBigStructureAlignment::fiducialCuts() {
  for (std::vector<unsigned int>::const_iterator index = m_indexes.begin(); index != m_indexes.end(); ++index) {
    if (m_debug)
      std::cout << "applying fiducial cuts in " << chamberPrettyNameFromId(*index) << std::endl;
    MuonResidualsTwoBin* fitter = m_fitterOrder[*index];
    fitter->fiducialCuts(*index);
  }
}

void MuonBigStructureAlignment::eraseNotSelectedResiduals() {
  for (std::vector<unsigned int>::const_iterator index = m_indexes.begin(); index != m_indexes.end(); ++index) {
    if (m_debug)
      std::cout << "erasing in " << chamberPrettyNameFromId(*index) << std::endl;
    MuonResidualsTwoBin* fitter = m_fitterOrder[*index];
    fitter->eraseNotSelectedResiduals();
  }
}

std::string MuonBigStructureAlignment::chamberPrettyNameFromId(unsigned int idx) {
  DetId id(idx);
  char cname[40];
  if (id.subdetId() == MuonSubdetId::DT) {
    DTChamberId chamberId(id.rawId());
    sprintf(cname, "MB%+d/%d/%02d", chamberId.wheel(), chamberId.station(), chamberId.sector());
  } else if (id.subdetId() == MuonSubdetId::CSC) {
    CSCDetId chamberId(id.rawId());
    sprintf(cname,
            "ME%s%d/%d/%02d",
            (chamberId.endcap() == 1 ? "+" : "-"),
            chamberId.station(),
            chamberId.ring(),
            chamberId.chamber());
  }
  return std::string(cname);
}

#include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentAlgorithmPluginFactory.h"
DEFINE_EDM_PLUGIN(AlignmentAlgorithmPluginFactory, MuonBigStructureAlignment, "MuonBigStructureAlignment");
