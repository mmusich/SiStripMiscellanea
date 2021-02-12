// -*- C++ -*-
//
// Package:    SiStripMiscellanea/SiStripApvGainInspector
// Class:      SiStripApvGainInspector
//
/**\class SiStripApvGainInspector SiStripApvGainInspector.cc SiStripMiscellanea/SiStripApvGainInspector/plugins/SiStripApvGainInspector.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marco Musich
//         Created:  Tue, 05 Jun 2018 15:46:15 GMT
//
//

// system include files
#include <memory>

// user include files
#include "CalibFormats/SiStripObjects/interface/SiStripDetCabling.h"
#include "CalibFormats/SiStripObjects/interface/SiStripGain.h"
#include "CalibFormats/SiStripObjects/interface/SiStripQuality.h"
#include "CalibTracker/Records/interface/SiStripDetCablingRcd.h"
#include "CalibTracker/Records/interface/SiStripGainRcd.h"
#include "CalibTracker/Records/interface/SiStripQualityRcd.h"
#include "CondFormats/SiStripObjects/interface/SiStripApvGain.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "CommonTools/TrackerMap/interface/TrackerMap.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CondTools/SiStrip/plugins/SiStripMiscalibrateHelper.cc"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

#include <cmath> /* log */

// ROOT includes
#include "TStyle.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2S.h"
#include "TProfile.h"
#include "TF1.h"

// user includes
#include "CalibTracker/SiStripChannelGain/interface/APVGainStruct.h"

class SiStripApvGainInspector : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit SiStripApvGainInspector(const edm::ParameterSet&);
  ~SiStripApvGainInspector() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  void checkBookAPVColls(const edm::EventSetup& es);
  void checkAndRetrieveTopology(const edm::EventSetup& setup);
  bool isGoodLandauFit(double* FitResults);
  void getPeakOfLandau(TH1* InputHisto, double* FitResults, double LowRange = 50, double HighRange = 5400);
  void storeOnTree(TFileService* tfs);
  void makeNicePlotStyle(TH1F* plot);

  // ----------member data ---------------------------

  TFileService* tfs;

  // map the APV ids to the charge plots
  std::map<std::pair<unsigned char, uint32_t>, TH1F*> histoMap_;

  edm::ESHandle<TrackerGeometry> tkGeom_;
  const TrackerGeometry* bareTkGeomPtr_;  // ugly hack to fill APV colls only once, but checks
  const TrackerTopology* tTopo_;

  int NStripAPVs;
  int NPixelDets;

  unsigned int GOOD;
  unsigned int BAD;
  unsigned int MASKED;

  std::vector<std::shared_ptr<stAPVGain> > APVsCollOrdered;
  std::unordered_map<unsigned int, std::shared_ptr<stAPVGain> > APVsColl;

  const TH2F* Charge_Vs_Index;
  TFile* fin;
  const std::string filename_;
  double minNrEntries;
  std::vector<unsigned int> wantedmods;

  std::unique_ptr<TrackerMap> ratio_map;
  std::unique_ptr<TrackerMap> old_payload_map;
  std::unique_ptr<TrackerMap> new_payload_map;
  std::unique_ptr<TrackerMap> mpv_map;
  std::unique_ptr<TrackerMap> mpv_err_map;
  std::unique_ptr<TrackerMap> entries_map;
  std::unique_ptr<TrackerMap> fitChi2_map;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
SiStripApvGainInspector::SiStripApvGainInspector(const edm::ParameterSet& iConfig)
    : bareTkGeomPtr_(nullptr),
      tTopo_(nullptr),
      GOOD(0),
      BAD(0),
      filename_(iConfig.getUntrackedParameter<std::string>("inputFile")),
      minNrEntries(iConfig.getUntrackedParameter<double>("minNrEntries", 20)),
      wantedmods(iConfig.getUntrackedParameter<std::vector<unsigned int> >("selectedModules")) {
  sort(wantedmods.begin(), wantedmods.end());

  edm::LogInfo("SelectedModules") << "Selected module list";
  for (std::vector<unsigned int>::const_iterator mod = wantedmods.begin(); mod != wantedmods.end(); mod++) {
    edm::LogVerbatim("SelectedModules") << *mod;
  }

  //now do what ever initialization is needed
  fin = TFile::Open(filename_.c_str(), "READ");
  Charge_Vs_Index = (TH2F*)fin->Get("DQMData/Run 999999/AlCaReco/Run summary/SiStripGainsAAG/Charge_Vs_Index_AagBunch");

  ratio_map = std::make_unique<TrackerMap>("ratio");
  ratio_map->setTitle("Average by module of the G2 Gain payload ratio (new/old)");
  ratio_map->setPalette(1);

  new_payload_map = std::make_unique<TrackerMap>("new_payload");
  new_payload_map->setTitle("Tracker Map of Updated G2 Gain payload averaged by module");
  new_payload_map->setPalette(1);

  old_payload_map = std::make_unique<TrackerMap>("old_payload");
  old_payload_map->setTitle("Tracker Map of Starting G2 Gain Payload averaged by module");
  old_payload_map->setPalette(1);

  // fit quality maps

  mpv_map = std::make_unique<TrackerMap>("MPV");
  mpv_map->setTitle("Landau Fit MPV average value per module [ADC counts/mm]");
  mpv_map->setPalette(1);

  mpv_err_map = std::make_unique<TrackerMap>("MPVerr");
  mpv_err_map->setTitle("Landau Fit MPV average error per module [ADC counts/mm]");
  mpv_err_map->setPalette(1);

  entries_map = std::make_unique<TrackerMap>("Entries");
  entries_map->setTitle("log_{10}(entries) average per module");
  entries_map->setPalette(1);

  fitChi2_map = std::make_unique<TrackerMap>("FitChi2");
  fitChi2_map->setTitle("log_{10}(Fit #chi^{2}/ndf) average per module");
  fitChi2_map->setPalette(1);
}

SiStripApvGainInspector::~SiStripApvGainInspector() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  fin->Close();
  delete fin;
}

//
// member functions
//

// ------------ method called for each event  ------------
void SiStripApvGainInspector::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  this->checkBookAPVColls(iSetup);  // check whether APV colls are booked and do so if not yet done
  this->checkAndRetrieveTopology(iSetup);

  edm::ESHandle<SiStripGain> gainHandle;
  iSetup.get<SiStripGainRcd>().get(gainHandle);
  if (!gainHandle.isValid()) {
    edm::LogError("SiStripApvGainInspector") << "gainHandle is not valid\n";
    exit(0);
  }

  edm::ESHandle<SiStripQuality> SiStripQuality_;
  iSetup.get<SiStripQualityRcd>().get(SiStripQuality_);

  for (unsigned int a = 0; a < APVsCollOrdered.size(); a++) {
    std::shared_ptr<stAPVGain> APV = APVsCollOrdered[a];

    if (APV->SubDet == PixelSubdetector::PixelBarrel || APV->SubDet == PixelSubdetector::PixelEndcap)
      continue;

    APV->isMasked = SiStripQuality_->IsApvBad(APV->DetId, APV->APVId);

    if (gainHandle->getNumberOfTags() != 2) {
      edm::LogError("SiStripApvGainInspector") << "NUMBER OF GAIN TAG IS EXPECTED TO BE 2\n";
      fflush(stdout);
      exit(0);
    };
    float newPreviousGain = gainHandle->getApvGain(APV->APVId, gainHandle->getRange(APV->DetId, 1), 1);
    if (APV->PreviousGain != 1 and newPreviousGain != APV->PreviousGain)
      edm::LogWarning("SiStripApvGainInspector") << "WARNING: ParticleGain in the global tag changed\n";
    APV->PreviousGain = newPreviousGain;

    float newPreviousGainTick = gainHandle->getApvGain(APV->APVId, gainHandle->getRange(APV->DetId, 0), 0);
    if (APV->PreviousGainTick != 1 and newPreviousGainTick != APV->PreviousGainTick) {
      edm::LogWarning("SiStripApvGainInspector")
          << "WARNING: TickMarkGain in the global tag changed\n"
          << std::endl
          << " APV->SubDet: " << APV->SubDet << " APV->APVId:" << APV->APVId << std::endl
          << " APV->PreviousGainTick: " << APV->PreviousGainTick << " newPreviousGainTick: " << newPreviousGainTick
          << std::endl;
    }
    APV->PreviousGainTick = newPreviousGainTick;
  }

  unsigned int I = 0;
  TH1F* Proj = nullptr;
  double FitResults[6];
  double MPVmean = 300;

  if (Charge_Vs_Index == nullptr) {
    edm::LogError("SiStripGainsPCLHarvester") << "Harvesting: could not find input histogram " << std::endl;
    return;
  }

  printf("Progressing Bar              :0%%       20%%       40%%       60%%       80%%       100%%\n");
  printf("Fitting Charge Distribution  :");
  int TreeStep = APVsColl.size() / 50;

  for (auto it = APVsColl.begin(); it != APVsColl.end(); it++, I++) {
    if (I % TreeStep == 0) {
      printf(".");
      fflush(stdout);
    }
    std::shared_ptr<stAPVGain> APV = it->second;
    if (APV->Bin < 0)
      APV->Bin = Charge_Vs_Index->GetXaxis()->FindBin(APV->Index);

    Proj = (TH1F*)(Charge_Vs_Index->ProjectionY(
        "", Charge_Vs_Index->GetXaxis()->FindBin(APV->Index), Charge_Vs_Index->GetXaxis()->FindBin(APV->Index), "e"));
    if (!Proj)
      continue;

    getPeakOfLandau(Proj, FitResults);
    APV->FitMPV = FitResults[0];
    APV->FitMPVErr = FitResults[1];
    APV->FitWidth = FitResults[2];
    APV->FitWidthErr = FitResults[3];
    APV->FitChi2 = FitResults[4];
    APV->FitNorm = FitResults[5];
    APV->NEntries = Proj->GetEntries();

    for (const auto& mod : wantedmods) {
      if (mod == APV->DetId) {
        edm::LogInfo("ModuleFound") << " module " << mod << " found! Storing... " << std::endl;
        histoMap_[std::make_pair(APV->APVId, APV->DetId)] = (TH1F*)Proj->Clone(Form("hClone_%s", Proj->GetName()));
      }
    }

    if (isGoodLandauFit(FitResults)) {
      APV->Gain = APV->FitMPV / MPVmean;
      if (APV->SubDet > 2)
        GOOD++;
    } else {
      APV->Gain = APV->PreviousGain;
      if (APV->SubDet > 2)
        BAD++;
    }
    if (APV->Gain <= 0)
      APV->Gain = 1;

    delete Proj;
  }
  printf("\n");
}

//********************************************************************************//
// ------------ method called once each job just before starting event loop  ------------
void SiStripApvGainInspector::checkBookAPVColls(const edm::EventSetup& es) {
  es.get<TrackerDigiGeometryRecord>().get(tkGeom_);
  const TrackerGeometry* newBareTkGeomPtr = &(*tkGeom_);
  if (newBareTkGeomPtr == bareTkGeomPtr_)
    return;  // already filled APVColls, nothing changed

  if (!bareTkGeomPtr_) {  // pointer not yet set: called the first time => fill the APVColls
    auto const& Det = newBareTkGeomPtr->dets();

    unsigned int Index = 0;

    for (unsigned int i = 0; i < Det.size(); i++) {
      DetId Detid = Det[i]->geographicalId();
      int SubDet = Detid.subdetId();

      if (SubDet == StripSubdetector::TIB || SubDet == StripSubdetector::TID || SubDet == StripSubdetector::TOB ||
          SubDet == StripSubdetector::TEC) {
        auto DetUnit = dynamic_cast<const StripGeomDetUnit*>(Det[i]);
        if (!DetUnit)
          continue;

        const StripTopology& Topo = DetUnit->specificTopology();
        unsigned int NAPV = Topo.nstrips() / 128;

        for (unsigned int j = 0; j < NAPV; j++) {
          auto APV = std::make_shared<stAPVGain>();
          APV->Index = Index;
          APV->Bin = -1;
          APV->DetId = Detid.rawId();
          APV->APVId = j;
          APV->SubDet = SubDet;
          APV->FitMPV = -1;
          APV->FitMPVErr = -1;
          APV->FitWidth = -1;
          APV->FitWidthErr = -1;
          APV->FitChi2 = -1;
          APV->FitNorm = -1;
          APV->Gain = -1;
          APV->PreviousGain = 1;
          APV->PreviousGainTick = 1;
          APV->x = DetUnit->position().basicVector().x();
          APV->y = DetUnit->position().basicVector().y();
          APV->z = DetUnit->position().basicVector().z();
          APV->Eta = DetUnit->position().basicVector().eta();
          APV->Phi = DetUnit->position().basicVector().phi();
          APV->R = DetUnit->position().basicVector().transverse();
          APV->Thickness = DetUnit->surface().bounds().thickness();
          APV->NEntries = 0;
          APV->isMasked = false;

          APVsCollOrdered.push_back(APV);
          APVsColl[(APV->DetId << 4) | APV->APVId] = APV;
          Index++;
          NStripAPVs++;
        }  // loop on APVs
      }    // if is Strips
    }      // loop on dets

    for (unsigned int i = 0; i < Det.size();
         i++) {  //Make two loop such that the Pixel information is added at the end --> make transition simpler
      DetId Detid = Det[i]->geographicalId();
      int SubDet = Detid.subdetId();
      if (SubDet == PixelSubdetector::PixelBarrel || SubDet == PixelSubdetector::PixelEndcap) {
        auto DetUnit = dynamic_cast<const PixelGeomDetUnit*>(Det[i]);
        if (!DetUnit)
          continue;

        const PixelTopology& Topo = DetUnit->specificTopology();
        unsigned int NROCRow = Topo.nrows() / (80.);
        unsigned int NROCCol = Topo.ncolumns() / (52.);

        for (unsigned int j = 0; j < NROCRow; j++) {
          for (unsigned int i = 0; i < NROCCol; i++) {
            auto APV = std::make_shared<stAPVGain>();
            APV->Index = Index;
            APV->Bin = -1;
            APV->DetId = Detid.rawId();
            APV->APVId = (j << 3 | i);
            APV->SubDet = SubDet;
            APV->FitMPV = -1;
            APV->FitMPVErr = -1;
            APV->FitWidth = -1;
            APV->FitWidthErr = -1;
            APV->FitChi2 = -1;
            APV->Gain = -1;
            APV->PreviousGain = 1;
            APV->PreviousGainTick = 1;
            APV->x = DetUnit->position().basicVector().x();
            APV->y = DetUnit->position().basicVector().y();
            APV->z = DetUnit->position().basicVector().z();
            APV->Eta = DetUnit->position().basicVector().eta();
            APV->Phi = DetUnit->position().basicVector().phi();
            APV->R = DetUnit->position().basicVector().transverse();
            APV->Thickness = DetUnit->surface().bounds().thickness();
            APV->isMasked = false;  //SiPixelQuality_->IsModuleBad(Detid.rawId());
            APV->NEntries = 0;

            APVsCollOrdered.push_back(APV);
            APVsColl[(APV->DetId << 4) | APV->APVId] = APV;
            Index++;
            NPixelDets++;

          }  // loop on ROC cols
        }    // loop on ROC rows
      }      // if Pixel
    }        // loop on Dets
  }          //if (!bareTkGeomPtr_) ...
  bareTkGeomPtr_ = newBareTkGeomPtr;
}

void SiStripApvGainInspector::storeOnTree(TFileService* tfs) {
  unsigned int tree_Index;
  unsigned int tree_Bin;
  unsigned int tree_DetId;
  unsigned char tree_APVId;
  unsigned char tree_SubDet;
  float tree_x;
  float tree_y;
  float tree_z;
  float tree_Eta;
  float tree_R;
  float tree_Phi;
  float tree_Thickness;
  float tree_FitMPV;
  float tree_FitMPVErr;
  float tree_FitWidth;
  float tree_FitWidthErr;
  float tree_FitChi2NDF;
  float tree_FitNorm;
  double tree_Gain;
  double tree_PrevGain;
  double tree_PrevGainTick;
  double tree_NEntries;
  bool tree_isMasked;

  TTree* MyTree;
  MyTree = tfs->make<TTree>("APVGain", "APVGain");
  MyTree->Branch("Index", &tree_Index, "Index/i");
  MyTree->Branch("Bin", &tree_Bin, "Bin/i");
  MyTree->Branch("DetId", &tree_DetId, "DetId/i");
  MyTree->Branch("APVId", &tree_APVId, "APVId/b");
  MyTree->Branch("SubDet", &tree_SubDet, "SubDet/b");
  MyTree->Branch("x", &tree_x, "x/F");
  MyTree->Branch("y", &tree_y, "y/F");
  MyTree->Branch("z", &tree_z, "z/F");
  MyTree->Branch("Eta", &tree_Eta, "Eta/F");
  MyTree->Branch("R", &tree_R, "R/F");
  MyTree->Branch("Phi", &tree_Phi, "Phi/F");
  MyTree->Branch("Thickness", &tree_Thickness, "Thickness/F");
  MyTree->Branch("FitMPV", &tree_FitMPV, "FitMPV/F");
  MyTree->Branch("FitMPVErr", &tree_FitMPVErr, "FitMPVErr/F");
  MyTree->Branch("FitWidth", &tree_FitWidth, "FitWidth/F");
  MyTree->Branch("FitWidthErr", &tree_FitWidthErr, "FitWidthErr/F");
  MyTree->Branch("FitChi2NDF", &tree_FitChi2NDF, "FitChi2NDF/F");
  MyTree->Branch("FitNorm", &tree_FitNorm, "FitNorm/F");
  MyTree->Branch("Gain", &tree_Gain, "Gain/D");
  MyTree->Branch("PrevGain", &tree_PrevGain, "PrevGain/D");
  MyTree->Branch("PrevGainTick", &tree_PrevGainTick, "PrevGainTick/D");
  MyTree->Branch("NEntries", &tree_NEntries, "NEntries/D");
  MyTree->Branch("isMasked", &tree_isMasked, "isMasked/O");

  uint32_t cachedId(0);
  SiStripMiscalibrate::Entry gain_ratio;
  SiStripMiscalibrate::Entry o_gain;
  SiStripMiscalibrate::Entry n_gain;
  SiStripMiscalibrate::Entry mpv;
  SiStripMiscalibrate::Entry mpv_err;
  SiStripMiscalibrate::Entry entries;
  SiStripMiscalibrate::Entry fitChi2;

  for (unsigned int a = 0; a < APVsCollOrdered.size(); a++) {
    std::shared_ptr<stAPVGain> APV = APVsCollOrdered[a];
    if (APV == nullptr)
      continue;
    //     printf(      "%i | %i | PreviousGain = %7.5f NewGain = %7.5f (#clusters=%8.0f)\n", APV->DetId,APV->APVId,APV->PreviousGain,APV->Gain, APV->NEntries);
    //fprintf(Gains,"%i | %i | PreviousGain = %7.5f(tick) x %7.5f(particle) NewGain (particle) = %7.5f (#clusters=%8.0f)\n", APV->DetId,APV->APVId,APV->PreviousGainTick, APV->PreviousGain,APV->Gain, APV->NEntries);

    // do not fill the Pixel
    if (APV->SubDet == PixelSubdetector::PixelBarrel || APV->SubDet == PixelSubdetector::PixelEndcap)
      continue;

    tree_Index = APV->Index;
    tree_Bin = Charge_Vs_Index->GetXaxis()->FindBin(APV->Index);
    tree_DetId = APV->DetId;
    tree_APVId = APV->APVId;
    tree_SubDet = APV->SubDet;
    tree_x = APV->x;
    tree_y = APV->y;
    tree_z = APV->z;
    tree_Eta = APV->Eta;
    tree_R = APV->R;
    tree_Phi = APV->Phi;
    tree_Thickness = APV->Thickness;
    tree_FitMPV = APV->FitMPV;
    tree_FitMPVErr = APV->FitMPVErr;
    tree_FitWidth = APV->FitWidth;
    tree_FitWidthErr = APV->FitWidthErr;
    tree_FitChi2NDF = APV->FitChi2;
    tree_FitNorm = APV->FitNorm;
    tree_Gain = APV->Gain;
    tree_PrevGain = APV->PreviousGain;
    tree_PrevGainTick = APV->PreviousGainTick;
    tree_NEntries = APV->NEntries;
    tree_isMasked = APV->isMasked;

    // flush the counters
    if (cachedId != 0 && tree_DetId != cachedId) {
      //ratio_map->fill(cachedId,gain_ratio.mean());
      ratio_map->fill(cachedId, o_gain.mean() / n_gain.mean());
      old_payload_map->fill(cachedId, o_gain.mean());
      new_payload_map->fill(cachedId, n_gain.mean());

      if (entries.mean() > 0) {
        mpv_map->fill(cachedId, mpv.mean());
        mpv_err_map->fill(cachedId, mpv_err.mean());
        entries_map->fill(cachedId, log10(entries.mean()));
        if (fitChi2.mean() > 0) {
          fitChi2_map->fill(cachedId, log10(fitChi2.mean()));
        } else {
          fitChi2_map->fill(cachedId, -1);
        }
      }

      gain_ratio.reset();
      o_gain.reset();
      n_gain.reset();

      mpv.reset();
      mpv_err.reset();
      entries.reset();
      fitChi2.reset();
    }

    cachedId = tree_DetId;
    gain_ratio.add(tree_PrevGain / tree_Gain);
    o_gain.add(tree_PrevGain);
    n_gain.add(tree_Gain);
    mpv.add(tree_FitMPV);
    mpv_err.add(tree_FitMPVErr);
    entries.add(tree_NEntries);
    fitChi2.add(tree_FitChi2NDF);

    if (tree_DetId == 402673324) {
      printf("%i | %i : %f --> %f  (%f)\n", tree_DetId, tree_APVId, tree_PrevGain, tree_Gain, tree_NEntries);
    }

    MyTree->Fill();
  }
}

//********************************************************************************//
void SiStripApvGainInspector::checkAndRetrieveTopology(const edm::EventSetup& setup) {
  if (!tTopo_) {
    edm::ESHandle<TrackerTopology> TopoHandle;
    setup.get<TrackerTopologyRcd>().get(TopoHandle);
    tTopo_ = TopoHandle.product();
  }
}

//********************************************************************************//
void SiStripApvGainInspector::getPeakOfLandau(TH1* InputHisto, double* FitResults, double LowRange, double HighRange) {
  FitResults[0] = -0.5;  //MPV
  FitResults[1] = 0;     //MPV error
  FitResults[2] = -0.5;  //Width
  FitResults[3] = 0;     //Width error
  FitResults[4] = -0.5;  //Fit Chi2/NDF
  FitResults[5] = 0;     //Normalization

  if (InputHisto->GetEntries() < minNrEntries)
    return;

  // perform fit with standard landau
  TF1 MyLandau("MyLandau", "landau", LowRange, HighRange);
  MyLandau.SetParameter(1, 300);
  InputHisto->Fit(&MyLandau, "QR WW");

  // MPV is parameter 1 (0=constant, 1=MPV, 2=Sigma)
  FitResults[0] = MyLandau.GetParameter(1);                     //MPV
  FitResults[1] = MyLandau.GetParError(1);                      //MPV error
  FitResults[2] = MyLandau.GetParameter(2);                     //Width
  FitResults[3] = MyLandau.GetParError(2);                      //Width error
  FitResults[4] = MyLandau.GetChisquare() / MyLandau.GetNDF();  //Fit Chi2/NDF
  FitResults[5] = MyLandau.GetParameter(0);
}

//********************************************************************************//
bool SiStripApvGainInspector::isGoodLandauFit(double* FitResults) {
  if (FitResults[0] <= 0)
    return false;
  //   if(FitResults[1] > MaxMPVError   )return false;
  //   if(FitResults[4] > MaxChi2OverNDF)return false;
  return true;
}

/*--------------------------------------------------------------------*/
void SiStripApvGainInspector::makeNicePlotStyle(TH1F* plot)
/*--------------------------------------------------------------------*/
{
  plot->GetXaxis()->CenterTitle(true);
  plot->GetYaxis()->CenterTitle(true);
  plot->GetXaxis()->SetTitleFont(42);
  plot->GetYaxis()->SetTitleFont(42);
  plot->GetXaxis()->SetTitleSize(0.05);
  plot->GetYaxis()->SetTitleSize(0.05);
  plot->GetXaxis()->SetTitleOffset(0.9);
  plot->GetYaxis()->SetTitleOffset(1.3);
  plot->GetXaxis()->SetLabelFont(42);
  plot->GetYaxis()->SetLabelFont(42);
  plot->GetYaxis()->SetLabelSize(.05);
  plot->GetXaxis()->SetLabelSize(.05);
}

// ------------ method called once each job just before starting event loop  ------------
void SiStripApvGainInspector::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void SiStripApvGainInspector::endJob() {
  edm::LogVerbatim("SelectedModules") << "Selected APVs:" << histoMap_.size() << std::endl;
  for (const auto& plot : histoMap_) {
    TCanvas* c1 = new TCanvas(Form("c1_%i_%i", plot.first.second, plot.first.first),
                              Form("c1_%i_%i", plot.first.second, plot.first.first),
                              800,
                              600);
    // Define common things for the different fits

    gStyle->SetOptFit(1011);
    c1->Clear();

    c1->SetLeftMargin(0.15);
    c1->SetRightMargin(0.10);
    plot.second->SetTitle(Form("Cluster Charge (%i,%i)", plot.first.second, plot.first.first));
    plot.second->GetXaxis()->SetTitle("Normalized Cluster Charge [ADC counts/mm]");
    plot.second->GetYaxis()->SetTitle("On-track clusters");
    plot.second->GetXaxis()->SetRangeUser(0., 1000.);

    this->makeNicePlotStyle(plot.second);
    plot.second->Draw();
    edm::LogVerbatim("SelectedModules") << " DetId: " << plot.first.second << " (" << plot.first.first << ")"
                                        << std::endl;
    ;

    c1->Print(Form("c1_%i_%i.png", plot.first.second, plot.first.first));
    c1->Print(Form("c1_%i_%i.pdf", plot.first.second, plot.first.first));
  }

  tfs = edm::Service<TFileService>().operator->();
  storeOnTree(tfs);

  auto range = SiStripMiscalibrate::getTruncatedRange(ratio_map.get());

  ratio_map->save(true, range.first, range.second, "G2_gain_ratio_map.pdf");
  ratio_map->save(true, range.first, range.second, "G2_gain_ratio_map.png");

  range = SiStripMiscalibrate::getTruncatedRange(old_payload_map.get());

  old_payload_map->save(true, range.first, range.second, "starting_G2_gain_payload_map.pdf");
  old_payload_map->save(true, range.first, range.second, "starting_G2_gain_payload_map.png");

  range = SiStripMiscalibrate::getTruncatedRange(new_payload_map.get());

  new_payload_map->save(true, range.first, range.second, "new_G2_gain_payload_map.pdf");
  new_payload_map->save(true, range.first, range.second, "new_G2_gain_payload_map.png");

  mpv_map->save(true, 250, 350., "mpv_map.pdf");
  mpv_map->save(true, 250, 350., "mpv_map.png");

  mpv_err_map->save(true, 0., 3., "mpv_err_map.pdf");
  mpv_err_map->save(true, 0., 3., "mpv_err_map.png");

  entries_map->save(true, 0, 0, "entries_map.pdf");
  entries_map->save(true, 0, 0, "entries_map.png");

  fitChi2_map->save(true, 0., 0., "fitChi2_map.pdf");
  fitChi2_map->save(true, 0., 0., "fitChi2_map.png");
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void SiStripApvGainInspector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SiStripApvGainInspector);
