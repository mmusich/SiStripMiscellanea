// -*- C++ -*-
//
// Package:    SiStripMiscellanea/SiStripResolutionAnalyzer
// Class:      SiStripResolutionAnalyzer
// 
/**\class SiStripResolutionAnalyzer SiStripResolutionAnalyzer.cc SiStripMiscellanea/SiStripResolutionAnalyzer/plugins/SiStripResolutionAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marco Musich
//         Created:  Wed, 12 Jul 2017 12:34:49 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"


#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "TrackingTools/PatternTools/interface/TrajectoryMeasurement.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryStateCombiner.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"

#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "CommonTools/UtilAlgos/interface/DetIdSelector.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "TH1F.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class SiStripResolutionAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit SiStripResolutionAnalyzer(const edm::ParameterSet&);
      ~SiStripResolutionAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      TH1F* m_ptrk;
      TH1F* m_etatrk;
      TH1F* m_nhits;
      TH1F* m_chi2Prob;

      edm::EDGetTokenT<TrajTrackAssociationCollection> m_ttacollToken;
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
SiStripResolutionAnalyzer::SiStripResolutionAnalyzer(const edm::ParameterSet& iConfig) :
  m_ttacollToken(consumes<TrajTrackAssociationCollection>(iConfig.getParameter<edm::InputTag>("trajTrackAssoCollection")))
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> tfserv;

   m_ptrk     = tfserv->make<TH1F>("trkmomentum","Refitted Track  momentum",100,0.,200.);
   m_etatrk   = tfserv->make<TH1F>("trketa","Refitted Track pseudorapidity",100,-4.,4.);
   m_nhits    = tfserv->make<TH1F>("trknhits","n. hits",30,0,30);
   m_chi2Prob = tfserv->make<TH1F>("chi2Prob","chi^{2} probability",50,0.,1.);
}


SiStripResolutionAnalyzer::~SiStripResolutionAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
SiStripResolutionAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  // loop on trajectories and plot TSOS local coordinate
  
  TrajectoryStateCombiner combiner;
  
  // Trajectory Handle

  Handle<TrajTrackAssociationCollection> ttac;
  iEvent.getByToken(m_ttacollToken,ttac);
  
  for(TrajTrackAssociationCollection::const_iterator pair=ttac->begin();pair!=ttac->end();++pair) {
    
    const edm::Ref<std::vector<Trajectory> > & traj = pair->key;
    const reco::TrackRef & trk = pair->val;
    const std::vector<TrajectoryMeasurement> & tmcoll = traj->measurements();
        
    double ProbChi2 = ChiSquaredProbability((double)( traj->chiSquared() ),(double)( traj->ndof(false) ));

    if (traj->foundHits() < 6 ) continue;
    if (ProbChi2 < 0.001 ) continue;
    if (trk->p()<3.) continue;

    m_ptrk->Fill(trk->p());
    m_etatrk->Fill(trk->eta());
    m_nhits->Fill(traj->foundHits());
    m_chi2Prob->Fill(ProbChi2);

    for(std::vector<TrajectoryMeasurement>::const_iterator measurement = tmcoll.begin() ; measurement!= tmcoll.end() ; ++measurement) {
      
      if(!measurement->updatedState().isValid()) continue;
      
      const TrajectoryStateOnSurface tsos = measurement->updatedState();
      const TrajectoryStateOnSurface unbiased = combiner(measurement->forwardPredictedState(), measurement->backwardPredictedState());

      TransientTrackingRecHit::ConstRecHitPointer hit = measurement->recHit();
      
    }
  }   
}


// ------------ method called once each job just before starting event loop  ------------
void 
SiStripResolutionAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SiStripResolutionAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SiStripResolutionAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SiStripResolutionAnalyzer);
