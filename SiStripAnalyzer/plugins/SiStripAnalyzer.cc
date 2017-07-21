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
#include <map>
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
#include "FWCore/Framework/interface/ESHandle.h"

#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "CondFormats/SiStripObjects/interface/SiStripLorentzAngle.h"
#include "CondFormats/DataRecord/interface/SiStripLorentzAngleRcd.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"


#include "DataFormats/GeometrySurface/interface/TrapezoidalPlaneBounds.h"
#include "DataFormats/GeometrySurface/interface/RectangularPlaneBounds.h"

#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "CalibTracker/SiStripCommon/interface/ShallowGainCalibration.h"


#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/SiStripDetId/interface/SiStripSubStructure.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/TrackReco/interface/DeDxHit.h"
#include "DataFormats/TrackReco/interface/TrackDeDxHits.h"

#include "CalibFormats/SiStripObjects/interface/SiStripGain.h" 
#include "CalibTracker/Records/interface/SiStripGainRcd.h"  


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
      double thickness(DetId id);
      bool   IsFarFromBorder(TrajectoryStateOnSurface* trajState, const uint32_t detid, const edm::EventSetup* iSetup);

      edm::EDGetTokenT< LumiScalersCollection > scalerToken_; 
      edm::EDGetTokenT< edm::View<reco::Track> > tracks_token_;
      edm::EDGetTokenT< TrajTrackAssociationCollection > association_token_;

      // ----------member data ---------------------------

      double       MinTrackMomentum;
      double       MaxTrackMomentum;
      double       MinTrackEta;
      double       MaxTrackEta;
      unsigned int MaxNrStrips;
      unsigned int MinTrackHits;
      double       MaxTrackChiOverNdf;
      int          MaxTrackingIteration;
      bool AllowSaturation;

      int maxLS,maxBx;
      float maxPU,maxLumi;

      const TrackerGeometry* m_tracker;
      std::map<DetId,double> m_thicknessMap;
      
      edm::Service<TFileService> fs;
      TProfile* p_instlumi_per_bx;   
      TProfile* p_PU_per_bx;        
      TProfile* p_lumiVsLS_;
      
      TProfile* p_nClust_per_bx;        
      TProfile* p_nClust_per_LS;
      TProfile* p_nClust_vs_lumi;

      TProfile* p_nUsedClust_per_bx;        
      TProfile* p_nUsedClust_per_LS;
      TProfile* p_nUsedClust_vs_lumi;


      TH1I*     h_events_per_bx; 
      TH1F*     h_CChargeOverPath;
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

   tracks_token_       =  consumes< edm::View<reco::Track> > (iConfig.getParameter<edm::InputTag>("Tracks"));
   association_token_  =  consumes< TrajTrackAssociationCollection >(iConfig.getParameter<edm::InputTag>("Tracks"));

   MinTrackMomentum        = iConfig.getUntrackedParameter<double>  ("minTrackMomentum"     ,  3.0);
   MaxTrackMomentum        = iConfig.getUntrackedParameter<double>  ("maxTrackMomentum"     ,  99999.0);
   MinTrackEta             = iConfig.getUntrackedParameter<double>  ("minTrackEta"          , -5.0);
   MaxTrackEta             = iConfig.getUntrackedParameter<double>  ("maxTrackEta"          ,  5.0);
   MaxNrStrips             = iConfig.getUntrackedParameter<unsigned>("maxNrStrips"          ,  2);
   MinTrackHits            = iConfig.getUntrackedParameter<unsigned>("MinTrackHits"         ,  8);
   MaxTrackChiOverNdf      = iConfig.getUntrackedParameter<double>  ("MaxTrackChiOverNdf"   ,  3);
   MaxTrackingIteration    = iConfig.getUntrackedParameter<int>     ("MaxTrackingIteration" ,  7);
   AllowSaturation         = iConfig.getUntrackedParameter<bool>    ("AllowSaturation"      , false);

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
  
  
  unsigned int nStripClus     = 0;
  unsigned int nUsedStripClus = 0;
  
  //std::cout<<"bx:"<<bx<<" ls:"<<ls<<" instLumi:"<<instLumi_<<" PU:"<<PU_<<std::endl;
  
  edm::ESHandle<TrackerGeometry> theTrackerGeometry;         iSetup.get<TrackerDigiGeometryRecord>().get( theTrackerGeometry );  
  m_tracker=&(* theTrackerGeometry );
  edm::ESHandle<SiStripGain> gainHandle;                     iSetup.get<SiStripGainRcd>().get(gainHandle);
  edm::Handle<edm::View<reco::Track> > tracks;	              iEvent.getByToken(tracks_token_, tracks);	  
  edm::Handle<TrajTrackAssociationCollection> associations;  iEvent.getByToken(association_token_, associations);
  
  for( TrajTrackAssociationCollection::const_iterator association = associations->begin(); association != associations->end(); association++) {
    const Trajectory*  traj  = association->key.get();
    const reco::Track* track = association->val.get();
    
    double	   trackchi2        = track->chi2();             
    double	   trackndof        = track->ndof();             
    double	   trackchi2ndof    = track->chi2()/track->ndof();
    float	   trackcharge      = track->charge();           
    float	   trackmomentum    = track->p();                
    float	   trackpt          = track->pt();               
    float	   trackpterr       = track->ptError();          
    unsigned int  trackhitsvalid   = track->numberOfValidHits();
    unsigned int  trackhitslost    = track->numberOfLostHits(); 
    double	   tracktheta       = track->theta();            
    double	   trackthetaerr    = track->thetaError();       
    double	   trackphi         = track->phi();              
    double	   trackphierr      = track->phiError();         
    double	   tracketa         = track->eta();              
    double	   tracketaerr      = track->etaError();         
    double	   trackdxy         = track->dxy();              
    double	   trackdxyerr      = track->dxyError();         
    double	   trackdsz         = track->dsz();              
    double	   trackdszerr      = track->dszError();         
    double	   trackqoverp      = track->qoverp();           
    double	   trackqoverperr   = track->qoverpError();      
    double	   trackvx          = track->vx();               
    double	   trackvy          = track->vy();               
    double	   trackvz          = track->vz();               
    int           trackalgo        = (int)track->algo();
    
    std::vector<TrajectoryMeasurement> measurements = traj->measurements();
    for(std::vector<TrajectoryMeasurement>::const_iterator measurement_it = measurements.begin(); measurement_it!=measurements.end(); measurement_it++){
      TrajectoryStateOnSurface trajState = measurement_it->updatedState();
      if( !trajState.isValid() ) continue;     
      
      const TrackingRecHit*         hit                 = (*measurement_it->recHit()).hit();
      const SiStripRecHit1D*        sistripsimple1dhit  = dynamic_cast<const SiStripRecHit1D*>(hit);
      const SiStripRecHit2D*        sistripsimplehit    = dynamic_cast<const SiStripRecHit2D*>(hit);
      const SiStripMatchedRecHit2D* sistripmatchedhit   = dynamic_cast<const SiStripMatchedRecHit2D*>(hit);
      const SiPixelRecHit*          sipixelhit          = dynamic_cast<const SiPixelRecHit*>(hit);
      
      const SiPixelCluster*   PixelCluster = NULL;
      const SiStripCluster*   StripCluster = NULL;
      uint32_t                DetId = 0;
      
      for(unsigned int h=0;h<2;h++){
	if(!sistripmatchedhit && h==1){
	  continue;
	}else if(sistripmatchedhit  && h==0){
	  StripCluster = &sistripmatchedhit->monoCluster();
	  DetId = sistripmatchedhit->monoId();
	}else if(sistripmatchedhit  && h==1){
	  StripCluster = &sistripmatchedhit->stereoCluster();;
	  DetId = sistripmatchedhit->stereoId();
	}else if(sistripsimplehit){
	  StripCluster = (sistripsimplehit->cluster()).get();
	  DetId = sistripsimplehit->geographicalId().rawId();
	}else if(sistripsimple1dhit){
	  StripCluster = (sistripsimple1dhit->cluster()).get();
	  DetId = sistripsimple1dhit->geographicalId().rawId();
	}else if(sipixelhit){
	  PixelCluster = (sipixelhit->cluster()).get();
	  DetId = sipixelhit->geographicalId().rawId();
	}else{
	  continue;
	}
	
	LocalVector             trackDirection = trajState.localDirection();
	double                  cosine         = trackDirection.z()/trackDirection.mag();
	bool                    Saturation     = false;
	bool                    Overlapping    = false;
	unsigned int            Charge         = 0;
	double                  Path           = (10.0*thickness(DetId))/fabs(cosine);
	double                  PrevGain       = -1;
	double                  PrevGainTick   = -1;
	int                     FirstStrip     = 0;
	unsigned int            NStrips        = 0;
	
	std::vector<unsigned char> amplitude;
	
	if(StripCluster){
	  
	  nStripClus++;
	  
	  const auto &Ampls          = StripCluster->amplitudes();
	  FirstStrip                 = StripCluster->firstStrip();
	  NStrips                    = Ampls.size();
	  int APVId                  = FirstStrip/128;
	  
	  /*
	    if(gainHandle.isValid()){ 
	    PrevGain     =  gainHandle->getApvGain(APVId,gainHandle->getRange(DetId, 1),1); 
	    PrevGainTick =  gainHandle->getApvGain(APVId,gainHandle->getRange(DetId, 0),1);           
	    }
	  */
	  
	  for(unsigned int a=0;a<Ampls.size();a++){               
	    Charge+=Ampls[a];
	    if(Ampls[a] >=254)Saturation =true;
	    amplitude.push_back( Ampls[a] );
	  }
	  
	  if(FirstStrip==0                                  )Overlapping=true;
	  if(FirstStrip==128                                )Overlapping=true;
	  if(FirstStrip==256                                )Overlapping=true;
	  if(FirstStrip==384                                )Overlapping=true;
	  if(FirstStrip==512                                )Overlapping=true;
	  if(FirstStrip==640                                )Overlapping=true;
	  
	  if(FirstStrip<=127 && FirstStrip+Ampls.size()>127)Overlapping=true;
	  if(FirstStrip<=255 && FirstStrip+Ampls.size()>255)Overlapping=true;
	  if(FirstStrip<=383 && FirstStrip+Ampls.size()>383)Overlapping=true;
	  if(FirstStrip<=511 && FirstStrip+Ampls.size()>511)Overlapping=true;
	  if(FirstStrip<=639 && FirstStrip+Ampls.size()>639)Overlapping=true;
	  
	  if(FirstStrip+Ampls.size()==127                   )Overlapping=true;
	  if(FirstStrip+Ampls.size()==255                   )Overlapping=true;
	  if(FirstStrip+Ampls.size()==383                   )Overlapping=true;
	  if(FirstStrip+Ampls.size()==511                   )Overlapping=true;
	  if(FirstStrip+Ampls.size()==639                   )Overlapping=true;
	  if(FirstStrip+Ampls.size()==767                   )Overlapping=true;
	  
	  bool farfromedge = IsFarFromBorder(&trajState,DetId, &iSetup);
	  
	  // starts here the selection
	  
	  if(tracketa        < MinTrackEta          )continue;
	  if(tracketa        > MaxTrackEta          )continue;
	  if(trackmomentum   < MinTrackMomentum     )continue;
	  if(trackmomentum   > MaxTrackMomentum     )continue;
	  if(trackhitsvalid  < MinTrackHits         )continue;
	  if(trackchi2ndof   > MaxTrackChiOverNdf   )continue;
	  if(trackalgo       > MaxTrackingIteration )continue;
	  
	  if(farfromedge     == false               )continue;
	  if(Overlapping     == true                )continue;
	  if(Saturation      && !AllowSaturation    )continue;
	  if(NStrips         >  MaxNrStrips         )continue;
	  
	  double ChargeOverPath = (double)Charge / Path ;
	  
	  nUsedStripClus++;
	  h_CChargeOverPath->Fill(ChargeOverPath);
	  
	} else if(PixelCluster){
	  
	  const auto&             Ampls          = PixelCluster->pixelADC();
	  int                     FirstRow       = PixelCluster->minPixelRow();
	  int                     FirstCol       = PixelCluster->minPixelCol();
	  FirstStrip                             = ((FirstRow/80)<<3 | (FirstCol/52)) * 128; //Hack to save the APVId
	  NStrips                                = 0;
	  Saturation                             = false;
	  Overlapping                            = false;
	  
	  for(unsigned int a=0;a<Ampls.size();a++){
	    Charge+=Ampls[a];
	    if(Ampls[a] >=254)Saturation =true;
	  } // loop on amplitudes
	} // if it's pixel   	 
      } // h-index
    } // loop on TM
  } // loop on tracks
  
  // fill histograms
  p_nClust_per_bx->Fill(bx,nStripClus); 
  p_nClust_per_LS->Fill(ls,nStripClus); 
  p_nClust_vs_lumi->Fill(instLumi_,nStripClus);

  p_nUsedClust_per_bx->Fill(bx,nUsedStripClus); 
  p_nUsedClust_per_LS->Fill(ls,nUsedStripClus); 
  p_nUsedClust_vs_lumi->Fill(instLumi_,nUsedStripClus);

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
  h_CChargeOverPath = fs->make<TH1F>("clusterChargeOverPath","cluster charge over path;cluster charge / path;# clusters",100,0.,1000.);


  p_nClust_per_bx  = fs->make<TProfile>("nStripClustersVsBx","n. of strip clusters vs bx;bx id;n. strip clusters",3500,-0.5,3495.);     
  p_nClust_per_LS = fs->make<TProfile>("nStripClustersVsLS","n. of strip clusters vs LS;LS; n. strip clusters",2500,0,2500);    
  p_nClust_vs_lumi = fs->make<TProfile>("nStripClustersVsLumi","n. of strip clusters vs inst. lumi;scal inst lumi E30 [Hz cm^{-2}]; n. strip clusters",100,0.,15000.);    
                       
  p_nUsedClust_per_bx  = fs->make<TProfile>("nUsedStripClustersVsBx","n. of used strip clusters vs bx;bx id;n. strip clusters",3500,-0.5,3495.);     		    
  p_nUsedClust_per_LS  = fs->make<TProfile>("nUsesStripClustersVsLS","n. of used strip clusters vs LS;LS; n. strip clusters",2500,0,2500);    			    
  p_nUsedClust_vs_lumi = fs->make<TProfile>("nUsedStripClustersVsLumi","n. of used strip clusters vs inst. lumi;scal inst lumi E30 [Hz cm^{-2}]; n. strip clusters",100,0.,15000.);   

}

// ------------ method called once each job just after ending the event loop  ------------
void 
SiStripAnalyzer::endJob() 
{
  std::cout<<"bx:"<<maxBx<<" ls:"<<maxLS<<" instLumi:"<<maxLumi<<" PU:"<<maxPU<<std::endl;
}

// ------------ method to get the detector thickness ------------
double SiStripAnalyzer::thickness(DetId id)
{
  std::map<DetId,double>::iterator th=m_thicknessMap.find(id);
  if(th!=m_thicknessMap.end())
    return (*th).second;
  else {
    double detThickness=1.;
    //compute thickness normalization
    const GeomDetUnit* it = m_tracker->idToDetUnit(DetId(id));
    bool isPixel = dynamic_cast<const PixelGeomDetUnit*>(it)!=0;
    bool isStrip = dynamic_cast<const StripGeomDetUnit*>(it)!=0;
    if (!isPixel && ! isStrip) {
      //FIXME throw exception
      edm::LogWarning("DeDxHitsProducer") << "\t\t this detID doesn't seem to belong to the Tracker";
      detThickness = 1.;
    }else{
      detThickness = it->surface().bounds().thickness();
    }
    
    m_thicknessMap[id]=detThickness;//computed value
    return detThickness;
  }
}

bool SiStripAnalyzer::IsFarFromBorder(TrajectoryStateOnSurface* trajState, const uint32_t detid, const edm::EventSetup* iSetup)
{ 
  edm::ESHandle<TrackerGeometry> tkGeom; iSetup->get<TrackerDigiGeometryRecord>().get( tkGeom );

  LocalPoint  HitLocalPos   = trajState->localPosition();
  LocalError  HitLocalError = trajState->localError().positionError() ;

  const GeomDetUnit* it = tkGeom->idToDetUnit(DetId(detid));
  if (dynamic_cast<const StripGeomDetUnit*>(it)==0 && dynamic_cast<const PixelGeomDetUnit*>(it)==0) {
     std::cout << "this detID doesn't seem to belong to the Tracker" << std::endl;
     return false;
  }

  const BoundPlane plane = it->surface();
  const TrapezoidalPlaneBounds* trapezoidalBounds( dynamic_cast<const TrapezoidalPlaneBounds*>(&(plane.bounds())));
  const RectangularPlaneBounds* rectangularBounds( dynamic_cast<const RectangularPlaneBounds*>(&(plane.bounds())));

  double DistFromBorder = 1.0;    
  double HalfLength     = it->surface().bounds().length() /2.0;

  if(trapezoidalBounds)
  {
      std::array<const float, 4> const & parameters = (*trapezoidalBounds).parameters();
     HalfLength     = parameters[3];
  }else if(rectangularBounds){
     HalfLength     = it->surface().bounds().length() /2.0;
  }else{return false;}

  if (fabs(HitLocalPos.y())+HitLocalError.yy() >= (HalfLength - DistFromBorder) ) return false;

  return true;
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
