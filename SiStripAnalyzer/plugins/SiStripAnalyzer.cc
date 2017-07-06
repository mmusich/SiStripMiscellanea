// -*- C++ -*-
//
// Package:    CalibTracker/SiStripAnalyzer
// Class:      SiStripAnalyzer
// 
/**\class SiStripAnalyzer SiStripAnalyzer.cc CalibTracker/SiStripAnalyzer/plugins/SiStripAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marco Musich
//         Created:  Thu, 06 Jul 2017 08:19:57 GMT
//
//


// system include files
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Scalers/interface/LumiScalers.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TH1I.h"
#include "TProfile.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class SiStripAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit SiStripAnalyzer(const edm::ParameterSet&);
      ~SiStripAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      edm::EDGetTokenT< LumiScalersCollection > scalerToken_; 
      // ----------member data ---------------------------
      int maxLS,maxBx;
      float maxPU,maxLumi;

      edm::Service<TFileService> fs;
      TProfile* p_instlumi_per_bx;   
      TProfile* p_PU_per_bx;        
      TProfile* p_lumiVsLS_;
      TH1I*     h_events_per_bx;   
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
SiStripAnalyzer::SiStripAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   usesResource("TFileService");
   scalerToken_ = consumes< LumiScalersCollection >(iConfig.getParameter<edm::InputTag>("lumiScalers"));
   maxLS=-999; 
   maxBx=-999;
   maxPU=-1.;
   maxLumi=-1.;
}


SiStripAnalyzer::~SiStripAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
SiStripAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   int bx = iEvent.bunchCrossing();
   int ls = iEvent.id().luminosityBlock();

   // Luminosity informations
   edm::Handle< LumiScalersCollection > lumiScalers;
   float instLumi_=-1.; float PU_=-1.;
   iEvent.getByToken(scalerToken_, lumiScalers); 
   if(lumiScalers.isValid()){
     if (lumiScalers->begin() != lumiScalers->end()) {
       instLumi_ = lumiScalers->begin()->instantLumi();
       PU_       = lumiScalers->begin()->pileup();
     }
   } else {
     edm::LogInfo("SiStripAnalyzer") 
       << "LumiScalers collection not found in the event; will write dummy values";
   }
   
   if(bx>maxBx)  maxBx = bx;
   if(ls>maxLS)  maxLS = ls;
   if(PU_>maxPU) maxPU = PU_;
   if(instLumi_>maxLumi) maxLumi = instLumi_;

   p_lumiVsLS_->Fill(ls,instLumi_);
   p_instlumi_per_bx->Fill(bx,instLumi_);
   p_PU_per_bx->Fill(bx,PU_);   
   h_events_per_bx->Fill(bx);

   //std::cout<<"bx:"<<bx<<" ls:"<<ls<<" instLumi:"<<instLumi_<<" PU:"<<PU_<<std::endl;
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
SiStripAnalyzer::beginJob()
{
  TH1F::SetDefaultSumw2(kTRUE);
  p_instlumi_per_bx = fs->make<TProfile>("instlumi_per_bx","instantaneous lumi per bx",3500,-0.5,3495.);
  p_PU_per_bx       = fs->make<TProfile>("PU_per_bx","PU per bx",3500,-0.5,3495.);
  h_events_per_bx   = fs->make<TH1I>("events_per_bx","events per bx",3500,-0.5,3495.);
  p_lumiVsLS_       = fs->make<TProfile>("lumiVsLS","scal lumi vs LS;LS;scal inst lumi E30 [Hz cm^{-2}]",2500,0,2500); 	      

}

// ------------ method called once each job just after ending the event loop  ------------
void 
SiStripAnalyzer::endJob() 
{
  std::cout<<"bx:"<<maxBx<<" ls:"<<maxLS<<" instLumi:"<<maxLumi<<" PU:"<<maxPU<<std::endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SiStripAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SiStripAnalyzer);
