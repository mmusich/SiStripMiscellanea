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
#include "CondFormats/SiStripObjects/interface/SiStripApvGain.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
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

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

// ROOT includes
#include "TFile.h"
#include "TH1F.h"
#include "TH2S.h"
#include "TProfile.h"
#include "TF1.h"

// user includes
#include "CalibTracker/SiStripChannelGain/interface/APVGainStruct.h"

class SiStripApvGainInspector : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit SiStripApvGainInspector(const edm::ParameterSet&);
      ~SiStripApvGainInspector();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      void checkBookAPVColls(const edm::EventSetup& es);
      void checkAndRetrieveTopology(const edm::EventSetup& setup);

      // ----------member data ---------------------------

      edm::ESHandle<TrackerGeometry> tkGeom_;
      const TrackerGeometry *bareTkGeomPtr_;  // ugly hack to fill APV colls only once, but checks
      const TrackerTopology* tTopo_;

      int NStripAPVs;
      int NPixelDets;

      std::vector<std::shared_ptr<stAPVGain> > APVsCollOrdered;
      std::unordered_map<unsigned int, std::shared_ptr<stAPVGain> > APVsColl; 	

      const TH2F* Charge_Vs_Index;
      TFile *fin;
      const std::string filename_;
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
SiStripApvGainInspector::SiStripApvGainInspector(const edm::ParameterSet& iConfig):
  filename_(iConfig.getUntrackedParameter<std::string> ("inputFile"))
{
   //now do what ever initialization is needed
  fin = TFile::Open(filename_.c_str(),"READ");
  Charge_Vs_Index = (TH2F*)fin->Get("DQMData/Run 999999/AlCaReco/Run summary/SiStripGainsAAG/Charge_Vs_Index_AagBunch");

}


SiStripApvGainInspector::~SiStripApvGainInspector()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  fin->Close();
  delete fin;
}


//
// member functions
//

// ------------ method called for each event  ------------
void
SiStripApvGainInspector::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   this->checkBookAPVColls(iSetup); // check whether APV colls are booked and do so if not yet done
   this->checkAndRetrieveTopology(iSetup);

   unsigned int I=0;
   TH1F* Proj = nullptr;
   double FitResults[6];
   double MPVmean = 300;
   
   if ( Charge_Vs_Index==nullptr ) {
     edm::LogError("SiStripGainsPCLHarvester") << "Harvesting: could not find input histogram "<< std::endl;
     return;
   }
    
  printf("Progressing Bar              :0%%       20%%       40%%       60%%       80%%       100%%\n");
  printf("Fitting Charge Distribution  :");
  int TreeStep = APVsColl.size()/50;
  
  for(auto it = APVsColl.begin();it!=APVsColl.end();it++,I++){
    
    if(I%TreeStep==0){printf(".");fflush(stdout);}
    std::shared_ptr<stAPVGain> APV = it->second;
    if(APV->Bin<0) APV->Bin = Charge_Vs_Index->GetXaxis()->FindBin(APV->Index);
        
    Proj = (TH1F*)(Charge_Vs_Index->ProjectionY("",Charge_Vs_Index->GetXaxis()->FindBin(APV->Index),Charge_Vs_Index->GetXaxis()->FindBin(APV->Index),"e"));
    if(!Proj)continue;

    
  }
   
}


//********************************************************************************//
// ------------ method called once each job just before starting event loop  ------------
void
SiStripApvGainInspector::checkBookAPVColls(const edm::EventSetup& es){

  es.get<TrackerDigiGeometryRecord>().get( tkGeom_ );
  const TrackerGeometry *newBareTkGeomPtr = &(*tkGeom_);
  if (newBareTkGeomPtr == bareTkGeomPtr_) return; // already filled APVColls, nothing changed

  if (!bareTkGeomPtr_) { // pointer not yet set: called the first time => fill the APVColls
    auto const & Det = newBareTkGeomPtr->dets();
    
    unsigned int Index=0;

    for(unsigned int i=0;i<Det.size();i++){
      
      DetId  Detid  = Det[i]->geographicalId(); 
      int    SubDet = Detid.subdetId();
      
      if( SubDet == StripSubdetector::TIB ||  SubDet == StripSubdetector::TID ||
	  SubDet == StripSubdetector::TOB ||  SubDet == StripSubdetector::TEC  ){
	
	auto DetUnit     = dynamic_cast<const StripGeomDetUnit*> (Det[i]);
	if(!DetUnit)continue;
	
	const StripTopology& Topo     = DetUnit->specificTopology();	
	unsigned int         NAPV     = Topo.nstrips()/128;
	
	for(unsigned int j=0;j<NAPV;j++){
	  auto APV = std::make_shared<stAPVGain>();
	  APV->Index         = Index;
	  APV->Bin           = -1;
	  APV->DetId         = Detid.rawId();
	  APV->APVId         = j;
	  APV->SubDet        = SubDet;
	  APV->FitMPV        = -1;
	  APV->FitMPVErr     = -1;
	  APV->FitWidth      = -1;
	  APV->FitWidthErr   = -1;
	  APV->FitChi2       = -1;
	  APV->FitNorm       = -1;
	  APV->Gain          = -1;
	  APV->PreviousGain  = 1;
	  APV->PreviousGainTick  = 1;
	  APV->x             = DetUnit->position().basicVector().x();
	  APV->y             = DetUnit->position().basicVector().y();
	  APV->z             = DetUnit->position().basicVector().z();
	  APV->Eta           = DetUnit->position().basicVector().eta();
	  APV->Phi           = DetUnit->position().basicVector().phi();
	  APV->R             = DetUnit->position().basicVector().transverse();
	  APV->Thickness     = DetUnit->surface().bounds().thickness();
	  APV->NEntries	   = 0;
	  APV->isMasked      = false;
	  
	  APVsCollOrdered.push_back(APV);
	  APVsColl[(APV->DetId<<4) | APV->APVId] = APV;
	  Index++;
	  NStripAPVs++;
	} // loop on APVs
      } // if is Strips
    } // loop on dets
    
    for(unsigned int i=0;i<Det.size();i++){  //Make two loop such that the Pixel information is added at the end --> make transition simpler
      DetId  Detid  = Det[i]->geographicalId();
      int    SubDet = Detid.subdetId();
      if( SubDet == PixelSubdetector::PixelBarrel || SubDet == PixelSubdetector::PixelEndcap ){
	auto DetUnit     = dynamic_cast<const PixelGeomDetUnit*> (Det[i]);
	if(!DetUnit) continue;
	
	const PixelTopology& Topo     = DetUnit->specificTopology();
	unsigned int         NROCRow  = Topo.nrows()/(80.);
	unsigned int         NROCCol  = Topo.ncolumns()/(52.);
	
	for(unsigned int j=0;j<NROCRow;j++){
	  for(unsigned int i=0;i<NROCCol;i++){
	    auto APV = std::make_shared<stAPVGain>();
	    APV->Index         = Index;
	    APV->Bin           = -1;
	    APV->DetId         = Detid.rawId();
	    APV->APVId         = (j<<3 | i);
	    APV->SubDet        = SubDet;
	    APV->FitMPV        = -1;
	    APV->FitMPVErr     = -1;
	    APV->FitWidth      = -1;
	    APV->FitWidthErr   = -1;
	    APV->FitChi2       = -1;
	    APV->Gain          = -1;
	    APV->PreviousGain  = 1;
	    APV->PreviousGainTick = 1;
	    APV->x             = DetUnit->position().basicVector().x();
	    APV->y             = DetUnit->position().basicVector().y();
	    APV->z             = DetUnit->position().basicVector().z();
	    APV->Eta           = DetUnit->position().basicVector().eta();
	    APV->Phi           = DetUnit->position().basicVector().phi();
	    APV->R             = DetUnit->position().basicVector().transverse();
	    APV->Thickness     = DetUnit->surface().bounds().thickness();
	    APV->isMasked      = false; //SiPixelQuality_->IsModuleBad(Detid.rawId());
	    APV->NEntries      = 0;
	    
	    APVsCollOrdered.push_back(APV);
	    APVsColl[(APV->DetId<<4) | APV->APVId] = APV;
	    Index++;
	    NPixelDets++;

	  } // loop on ROC cols
	} // loop on ROC rows
      } // if Pixel
    } // loop on Dets  
  }  //if (!bareTkGeomPtr_) ... 
  bareTkGeomPtr_ = newBareTkGeomPtr;
}

void SiStripApvGainInspector::checkAndRetrieveTopology(const edm::EventSetup& setup) {
  if( !tTopo_ ) {
    edm::ESHandle<TrackerTopology> TopoHandle;
    setup.get<TrackerTopologyRcd>().get( TopoHandle );
    tTopo_ = TopoHandle.product();
  }
}



// ------------ method called once each job just before starting event loop  ------------
void
SiStripApvGainInspector::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
SiStripApvGainInspector::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SiStripApvGainInspector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
