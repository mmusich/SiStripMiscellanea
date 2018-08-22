// -*- C++ -*-
//
// Package:    CalibTracker/SiStripClusterAnalyzer
// Class:      SiStripClusterAnalyzer
// 
/**\class SiStripClusterAnalyzer SiStripClusterAnalyzer.cc CalibTracker/SiStripAnalyzer/plugins/SiStripClusterAnalyzer.cc

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

#include "CondTools/SiStrip/plugins/SiStripMiscalibrateHelper.cc"
#include "RecoLocalTracker/SiStripClusterizer/interface/SiStripClusterInfo.h"
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

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"

#include "TH1I.h"
#include "TProfile.h"

// ROOT includes
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"

// RooFit includes
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooFFTConvPdf.h"

#define HIPDEBUG true

namespace SiStripHelper {
  
  enum layer {TIBL1=1,
	      TIBL2=2,
	      TIBL3=3,
	      TIBL4=4,	       
	      TOBL1=5,
	      TOBL2=6,
	      TOBL3=7,
	      TOBL4=8,
	      TOBL5=9,
	      TOBL6=10,
	      TIDD1=11,
	      TIDD2=12,
	      TIDD3=13,
	      TECD1=14,
	      TECD2=15,
	      TECD3=16,
	      TECD4=17,
	      TECD5=18,
	      TECD6=19,
	      TECD7=20,
	      TECD8=21,
	      TECD9=22,
	      NUM_OF_TYPES=23,
  };

  TF1 * f2 = NULL;
  TF1 * f3 = NULL;
}

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class SiStripClusterAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit SiStripClusterAnalyzer(const edm::ParameterSet&);
      ~SiStripClusterAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      double thickness(DetId id);
      bool   IsFarFromBorder(TrajectoryStateOnSurface* trajState, const uint32_t detid, const edm::EventSetup* iSetup);
      int checkLayer( unsigned int iidd, const TrackerTopology* tTopo); 
      std::string getStringFromEnum(SiStripHelper::layer e);
      template <class OBJECT_TYPE> int GetIndex(const std::vector<OBJECT_TYPE*>& vec, const TString& name);
      std::string myreplace(const std::string &s,const std::string &toReplace,const std::string &replaceWith);
      void fitStoN(TH1F *hist);
      static Double_t function_sum(Double_t *x, Double_t *par);
      void makeNicePlotStyle(RooPlot* plot);

      edm::EDGetTokenT< LumiScalersCollection > scalerToken_; 
      edm::EDGetTokenT< edm::View<reco::Track> > tracks_token_;
      edm::EDGetTokenT< TrajTrackAssociationCollection > association_token_;
      edm::FileInPath gfp_;          /*!< File Path for the ideal geometry. */ 

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

      bool m_verbose_fit;   

      // cluster analysis
      bool   applyClusterQuality_;
      double sToNLowerLimit_;
      double sToNUpperLimit_;
      double widthLowerLimit_;
      double widthUpperLimit_;

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

      TProfile* p_clusterWitdh_vs_Layer;

      TH1F*     h_nClusters;
      TH1F*     h_nUsedClusters;
      TH1I*     h_events_per_bx; 
      TH1F*     h_CChargeOverPath;
      TH1I*     h_CCharge;
      TH1I*     h_SeedCharge;
      TH1I*     h_SingleStripCharge;  
      TH1F*     h_SeedNoise;
      TH1F*     h_SeedGain;

      TH1F*     h_StoN;
      TH1F*     h_Noise;
      TH1F*     h_StoNCorr;
      TH1F*     h_ClusterWidth;

      std::vector<TH1*> vTrackHistos_;

      // std::map<SiStripHelper::layer,TH1F*>  h_StoN_layer;
      // std::map<SiStripHelper::layer,TH1F*>  h_Noise_layer;
      // std::map<SiStripHelper::layer,TH1F*>  h_StoNCorr_layer;
      // std::map<SiStripHelper::layer,TH1F*>  h_CCOverP_layer;
      // std::map<SiStripHelper::layer,TH1F*>  h_ClusterWidth_layer;

      std::map<sistripsummary::TrackerRegion,TH1F*>  h_StoN_region;	    
      std::map<sistripsummary::TrackerRegion,TH1F*>  h_Noise_region;	    
      std::map<sistripsummary::TrackerRegion,TH1F*>  h_StoNCorr_region;    
      std::map<sistripsummary::TrackerRegion,TH1F*>  h_CCOverP_region;	    
      std::map<sistripsummary::TrackerRegion,TH1F*>  h_ClusterWidth_region;
      std::map<sistripsummary::TrackerRegion,TH1I*>  h_SeedCharge_region;
      std::map<sistripsummary::TrackerRegion,TH1I*>  h_SingleStripCharge_region;

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
SiStripClusterAnalyzer::SiStripClusterAnalyzer(const edm::ParameterSet& iConfig):
  gfp_(iConfig.getUntrackedParameter<edm::FileInPath>("geoFile",edm::FileInPath("CalibTracker/SiStripCommon/data/SiStripDetInfo.dat")))
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
   m_verbose_fit           = iConfig.getParameter<bool>("verbose_fit");

   // // cluster quality conditions
   // edm::ParameterSet cluster_condition = iConfig.getParameter<edm::ParameterSet>("ClusterConditions");
   // applyClusterQuality_ = cluster_condition.getParameter<bool>("On");
   // sToNLowerLimit_      = cluster_condition.getParameter<double>("minStoN");
   // sToNUpperLimit_      = cluster_condition.getParameter<double>("maxStoN");
   // widthLowerLimit_     = cluster_condition.getParameter<double>("minWidth");
   // widthUpperLimit_     = cluster_condition.getParameter<double>("maxWidth");


   maxLS=-999; 
   maxBx=-999;
   maxPU=-1.;
   maxLumi=-1.;
}


SiStripClusterAnalyzer::~SiStripClusterAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
SiStripClusterAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
    edm::LogInfo("SiStripClusterAnalyzer") 
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
  
  // Get the geometry
  edm::ESHandle<TrackerTopology> tTopo_handle;
  iSetup.get<TrackerTopologyRcd>().get(tTopo_handle);
  const TrackerTopology* tTopo = tTopo_handle.product();

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
    unsigned int   trackhitsvalid   = track->numberOfValidHits();
    unsigned int   trackhitslost    = track->numberOfLostHits(); 
    double	   tracktheta       = track->theta();            
    double	   trackthetaerr    = track->thetaError();       
    double	   trackphi         = track->phi();              
    double	   trackphierr      = track->phiError();         
    double	   tracketa         = track->eta();              
    double	   tracketaerr      = track->etaError();         
    double	   trackdxy         = track->dxy();              
    double	   trackdxyerr      = track->dxyError();         
    double	   trackdz          = track->dz();              
    double	   trackdzerr       = track->dzError();           
    double	   trackqoverp      = track->qoverp();           
    double	   trackqoverperr   = track->qoverpError();      
    double	   trackvx          = track->vx();               
    double	   trackvy          = track->vy();               
    double	   trackvz          = track->vz();               
    int            trackalgo        = (int)track->algo();
    
    // fill control plots
    static const int chi2index = this->GetIndex(vTrackHistos_,"h_chi2");
    vTrackHistos_[chi2index]->Fill(trackchi2);
    
    static const int ndofindex = this->GetIndex(vTrackHistos_,"h_trackndof");
    vTrackHistos_[ndofindex]->Fill(trackndof);
    
    static const int chargeindex = this->GetIndex(vTrackHistos_,"h_trackcharge");
    vTrackHistos_[chargeindex]->Fill(trackcharge);

    static const int thetaindex = this->GetIndex(vTrackHistos_,"h_tracktheta");
    vTrackHistos_[thetaindex]->Fill(tracktheta);

    static const int dxyindex = this->GetIndex(vTrackHistos_,"h_trackdxy");
    vTrackHistos_[dxyindex]->Fill(trackdxy);

    static const int dzindex = this->GetIndex(vTrackHistos_,"h_trackdz");
    vTrackHistos_[dzindex]->Fill(trackdz);

    static const int dxyerrindex = this->GetIndex(vTrackHistos_,"h_trackdxyerr");
    vTrackHistos_[dxyerrindex]->Fill(trackdxyerr);

    static const int dzerrindex = this->GetIndex(vTrackHistos_,"h_trackdzerr");
    vTrackHistos_[dzerrindex]->Fill(trackdzerr);
    
    static const int vxindex = this->GetIndex(vTrackHistos_,"h_trackvx");
    vTrackHistos_[vxindex]->Fill(trackvx);

    static const int vyindex = this->GetIndex(vTrackHistos_,"h_trackvy");
    vTrackHistos_[vyindex]->Fill(trackvy);

    static const int vzindex = this->GetIndex(vTrackHistos_,"h_trackvz");
    vTrackHistos_[vzindex]->Fill(trackvz);

    static const int thetaerrindex = this->GetIndex(vTrackHistos_,"h_trackthetaerr");
    vTrackHistos_[thetaerrindex]->Fill(trackthetaerr);

    static const int phierrindex = this->GetIndex(vTrackHistos_,"h_trackphierr");
    vTrackHistos_[phierrindex]->Fill(trackphierr);

    static const int kappaerrindex = this->GetIndex(vTrackHistos_,"h_curvatureerr");
    vTrackHistos_[kappaerrindex]->Fill(trackqoverperr);

    static const int algoindex = this->GetIndex(vTrackHistos_,"h_trackalgo");
    vTrackHistos_[algoindex]->Fill(trackalgo);
   
    static const int normchi2index = this->GetIndex(vTrackHistos_,"h_normchi2");
    vTrackHistos_[normchi2index]->Fill(trackchi2ndof);
    static const int ptindex = this->GetIndex(vTrackHistos_,"h_pt");
    if (trackpterr != 0.) {
      static const int ptResolutionindex = this->GetIndex(vTrackHistos_,"h_ptResolution");
      vTrackHistos_[ptResolutionindex]->Fill(trackpterr/trackpt);
    }
    vTrackHistos_[ptindex]->Fill(trackpt);

    static const int etaindex = this->GetIndex(vTrackHistos_,"h_tracketa");
    vTrackHistos_[etaindex]->Fill(tracketa);

    static const int etaerrindex = this->GetIndex(vTrackHistos_,"h_tracketaerr");
    vTrackHistos_[etaerrindex]->Fill(tracketaerr);

    static const int phiindex = this->GetIndex(vTrackHistos_,"h_trackphi");
    vTrackHistos_[phiindex]->Fill(trackphi);

    static const int numOfValidHitsindex = this->GetIndex(vTrackHistos_,"h_trackNumberOfValidHits");
    vTrackHistos_[numOfValidHitsindex]->Fill(trackhitsvalid);

    static const int numOfLostHitsindex = this->GetIndex(vTrackHistos_,"h_trackNumberOfLostHits");
    vTrackHistos_[numOfLostHitsindex]->Fill(trackhitslost);

    static const int kappaindex = this->GetIndex(vTrackHistos_,"h_curvature");
    vTrackHistos_[kappaindex]->Fill(trackqoverp);

    std::vector<TrajectoryMeasurement> measurements = traj->measurements();
    for(std::vector<TrajectoryMeasurement>::const_iterator measurement_it = measurements.begin(); measurement_it!=measurements.end(); measurement_it++){
      TrajectoryStateOnSurface trajState = measurement_it->updatedState();
      if( !trajState.isValid() ) continue;     
      
      const TrackingRecHit*         hit                 = (*measurement_it->recHit()).hit();
      const SiStripRecHit1D*        sistripsimple1dhit  = dynamic_cast<const SiStripRecHit1D*>(hit);
      const SiStripRecHit2D*        sistripsimplehit    = dynamic_cast<const SiStripRecHit2D*>(hit);
      const SiStripMatchedRecHit2D* sistripmatchedhit   = dynamic_cast<const SiStripMatchedRecHit2D*>(hit);
      const SiPixelRecHit*          sipixelhit          = dynamic_cast<const SiPixelRecHit*>(hit);
      
      const uint32_t& detid = hit->geographicalId().rawId();
      //      SiStripHelper::layers layer = static_cast<std::underlying_type_t<SiStripHelper::layers>>(checkLayer(detid,tTopo)); 
      auto layer  = checkLayer(detid,tTopo);
      auto regions = SiStripMiscalibrate::getRegionsFromDetId(tTopo,detid); 

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
	double                  cosRZ          = fabs(trackDirection.z())/trackDirection.mag();
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
	  
	  SiStripClusterInfo SiStripClusterInfo_(*StripCluster,iSetup,detid);
	  float    StoN     = SiStripClusterInfo_.signalOverNoise();
	  float    noise    = SiStripClusterInfo_.noiseRescaledByGain();
	  uint16_t charge   = SiStripClusterInfo_.charge();
	  uint16_t width    = SiStripClusterInfo_.width();
	  float position    = SiStripClusterInfo_.baryStrip();
	  
	  uint16_t seedcharge = SiStripClusterInfo_.maxCharge();
	  float seednoise     = SiStripClusterInfo_.stripNoisesRescaledByGain().at(SiStripClusterInfo_.maxIndex());
	  float seednoisepure = SiStripClusterInfo_.stripNoises().at(SiStripClusterInfo_.maxIndex());
	  float seedgain      = SiStripClusterInfo_.stripGains().at(SiStripClusterInfo_.maxIndex());
	  
	  h_CCharge->Fill(charge);
	  h_SeedCharge->Fill(seedcharge);
	  h_SeedNoise->Fill(seednoise);
	  h_SeedGain->Fill(seedgain);

	  h_StoN->Fill(StoN);
	  h_StoNCorr->Fill(StoN*cosRZ);
	  h_ClusterWidth->Fill(width);
	  h_Noise->Fill(noise);
	  
	  p_clusterWitdh_vs_Layer->Fill(layer,width);

	  for(unsigned int a=0;a<Ampls.size();a++){               
	    h_SingleStripCharge->Fill(Ampls[a]);
	  }

	  //if( (applyClusterQuality_) &&
	  //      (SiStripClusterInfo_.signalOverNoise() > sToNLowerLimit_ &&
	  // SiStripClusterInfo_.signalOverNoise() < sToNUpperLimit_ &&
	  // SiStripClusterInfo_.width() > widthLowerLimit_ &&
	  // SiStripClusterInfo_.width() < widthUpperLimit_) ) { 
	 
	  for (unsigned int j=0; j<regions.size();j++){

	    h_StoN_region[regions[j]]->Fill(StoN); 
	    h_StoNCorr_region[regions[j]]->Fill(StoN*cosRZ);
	    h_ClusterWidth_region[regions[j]]->Fill(width);
	    h_Noise_region[regions[j]]->Fill(noise);
	    h_SeedCharge_region[regions[j]]->Fill(seedcharge);
	    
	    for(unsigned int a=0;a<Ampls.size();a++){               
	      h_SingleStripCharge_region[regions[j]]->Fill(Ampls[a]);
	    }
	  }

	  // }	  
	  
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
	  
	  double ChargeOverPath = (double)Charge / Path ;
	  for (unsigned int j=0; j<regions.size();j++){
	    h_CCOverP_region[regions[j]]->Fill(ChargeOverPath);
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

  h_nClusters->Fill(nStripClus);	
  h_nUsedClusters->Fill(nUsedStripClus);

}


// ------------ method called once each job just before starting event loop  ------------
void 
SiStripClusterAnalyzer::beginJob()
{
  TH1F::SetDefaultSumw2(kTRUE);
  p_instlumi_per_bx = fs->make<TProfile>("instlumi_per_bx","instantaneous lumi per bx",3500,-0.5,3495.);
  p_PU_per_bx       = fs->make<TProfile>("PU_per_bx","PU per bx",3500,-0.5,3495.);
  h_events_per_bx   = fs->make<TH1I>("events_per_bx","events per bx",3500,-0.5,3495.);
  p_lumiVsLS_       = fs->make<TProfile>("lumiVsLS","scal lumi vs LS;LS;scal inst lumi E30 [Hz cm^{-2}]",2500,0,2500); 	
  h_CChargeOverPath = fs->make<TH1F>("clusterChargeOverPath","cluster charge over path;cluster charge / path;# clusters",150,0.,1500.);

  h_CCharge         = fs->make<TH1I>("clusterCharge","cluster charge;cluster charge [ADC counts]",100,0,500);
  h_SeedCharge      = fs->make<TH1I>("seedCharge","cluster seed charge;cluster seed charge [ADC counts]",256,0,256);
  h_SingleStripCharge = fs->make<TH1I>("signleStripCharge","single strip charge;single strip charge [ADC counts]",256,0,256);
  h_SeedNoise       = fs->make<TH1F>("seedNoise","cluster seed noise;cluster seed noise [ADC counts]",100,0.,10.);
  h_SeedGain        = fs->make<TH1F>("seedGain","cluster seed gain;cluster seed gain",100,0.,2.);

  h_StoN            = fs->make<TH1F>("clusterStoN","cluster Raw S/N;cluster raw S/N;# clusters",100,0.,100.);
  h_Noise           = fs->make<TH1F>("clusterNoise","cluster Noise;cluster Noise [ADC counts];# clusters",100,0.,10.);
  h_StoNCorr        = fs->make<TH1F>("clusterStoNCorr","cluster S/N corrected for path length;cluster S/N;# clusters",100,0.,100.);

  h_ClusterWidth    = fs->make<TH1F>("clusterWidth","cluster width;cluster witdh;# clusters",20,-0.5,19.5);

  p_clusterWitdh_vs_Layer = fs->make<TProfile>("clusterWidthPerLayer","Tracker Layer;Average Cluster width",22,0.5,22.5);
  for ( int fooInt = SiStripHelper::TIBL1; fooInt != SiStripHelper::NUM_OF_TYPES; fooInt++ ){
    SiStripHelper::layer layer = static_cast<SiStripHelper::layer>(fooInt);
    std::string s_layer = getStringFromEnum(layer);
    p_clusterWitdh_vs_Layer->GetXaxis()->SetBinLabel(fooInt,s_layer.c_str());
  }

  std::vector<std::string> parts = { "Tracker","TIB" ,"TIB_1" ,"TIB_2" ,"TIB_3" ,"TIB_4"  ,"TID"    ,"TIDP"   ,"TIDP_1" ,"TIDP_2" ,"TIDP_3" ,"TIDM"   ,"TIDM_1" ,"TIDM_2" ,"TIDM_3" ,"TOB"    ,"TOB_1"  ,"TOB_2"  ,"TOB_3"  ,"TOB_4"  ,"TOB_5"  ,"TOB_6"  ,"TEC"    ,"TECP"   ,"TECP_1" ,"TECP_2" ,"TECP_3" ,"TECP_4" ,"TECP_5" ,"TECP_6" ,"TECP_7" ,"TECP_8" ,"TECP_9" ,"TECM"   ,"TECM_1" ,"TECM_2" ,"TECM_3" ,"TECM_4" ,"TECM_5" ,"TECM_6" ,"TECM_7" ,"TECM_8" ,"TECM_9"}; 

  for(const auto &part : parts){
    std::string append = part; //myreplace(part,"_"," ");
    auto region = SiStripMiscalibrate::getRegionFromString(part);
   
    h_StoN_region[region]     = fs->make<TH1F>(Form("clusterStoN_%s",part.c_str()),Form("cluster Raw S/N (%s); %s cluster raw S/N;# clusters",append.c_str(),append.c_str()),100,0.,100.);
    h_Noise_region[region]    = fs->make<TH1F>(Form("clusterNoise_%s",part.c_str()),Form("cluster Noise (%s); %s cluster Noise [ADC counts];# clusters",append.c_str(),append.c_str()),100,0.,10.);
    h_StoNCorr_region[region] = fs->make<TH1F>(Form("clusterStoNCorr_%s",part.c_str()),Form("cluster S/N (%s) corrected for path length;%s cluster S/N;# clusters",append.c_str(),append.c_str()),100,0.,100.);
    h_CCOverP_region[region]  = fs->make<TH1F>(Form("clusterChargeOverPath_%s",part.c_str()),Form("cluster charge over path (%s);%s cluster charge / path;# clusters",append.c_str(),append.c_str()),150,0.,1500.);
    h_ClusterWidth_region[region] = fs->make<TH1F>(Form("clusterWidth_%s",part.c_str()),Form("cluster width (%s);%s cluster witdh;# clusters",append.c_str(),append.c_str()),20,-0.5,19.5);
    h_SeedCharge_region[region]   = fs->make<TH1I>(Form("seedCharge_%s",part.c_str()),Form("cluster seed charge (%s);%s cluster seed charge [ADC counts]",append.c_str(),append.c_str()),256,0,256);
    h_SingleStripCharge_region[region] = fs->make<TH1I>(Form("signleStripCharge_%s",part.c_str()),Form("single strip charge (%s);%s single strip charge [ADC counts]",append.c_str(),append.c_str()),256,0,256);
  
  }

  // for ( int fooInt = SiStripHelper::TIBL1; fooInt != SiStripHelper::NUM_OF_TYPES; fooInt++ ){
  //   SiStripHelper::layer layer = static_cast<SiStripHelper::layer>(fooInt);
  //   std::string s_layer = getStringFromEnum(layer);
  //   std::string append = myreplace(s_layer," ","_");

  //   h_StoN_layer[layer]     = fs->make<TH1F>(Form("clusterStoN_%s",append.c_str()),Form("cluster Raw S/N (%s); %s cluster raw S/N;# clusters",s_layer.c_str(),s_layer.c_str()),100,0.,100.);

  //   h_Noise_layer[layer]    = fs->make<TH1F>(Form("clusterNoise_%s",append.c_str()),Form("cluster Noise (%s); %s cluster Noise [ADC counts];# clusters",s_layer.c_str(),s_layer.c_str()),100,0.,10.);
  //   h_StoNCorr_layer[layer] = fs->make<TH1F>(Form("clusterStoNCorr_%s",append.c_str()),Form("cluster S/N (%s) corrected for path length;%s cluster S/N;# clusters",s_layer.c_str(),s_layer.c_str()),100,0.,100.);
  //   h_CCOverP_layer[layer]  = fs->make<TH1F>(Form("clusterChargeOverPath_%s",append.c_str()),Form("cluster charge over path (%s);%s cluster charge / path;# clusters",s_layer.c_str(),s_layer.c_str()),100,0.,1000.);

  //   h_ClusterWidth_layer[layer] = fs->make<TH1F>(Form("clusterWidth_%s",append.c_str()),Form("cluster width (%s);%s cluster witdh;# clusters",s_layer.c_str(),s_layer.c_str()),20,-0.5,19.5);

  // }  

  h_nClusters       = fs->make<TH1F>("nStripClusters","n. of Strip clusters;n. Strip clusters; events",250,0.,5000.);
  h_nUsedClusters   = fs->make<TH1F>("nUsedStripClusters","n. of selected Strip clusters;n. selected Strip clusters; events",200,0.,2000.);

  p_nClust_per_bx  = fs->make<TProfile>("nStripClustersVsBx","n. of Strip clusters vs bx;bx id;#LT n.clusters/event #GT",3500,-0.5,3495.);     
  p_nClust_per_LS  = fs->make<TProfile>("nStripClustersVsLS","n. of Strip clusters vs LS;LS; #LT n.clusters/event #GT",2500,0,2500);    
  p_nClust_vs_lumi = fs->make<TProfile>("nStripClustersVsLumi","n. of Strip clusters vs inst. lumi;scal inst lumi E30 [Hz cm^{-2}]; #LT n.clusters/event #GT",100,0.,15000.);        
  p_nUsedClust_per_bx  = fs->make<TProfile>("nUsedStripClustersVsBx","n. of used Strip clusters vs bx;bx id;#LT n.clusters/event #GT",3500,-0.5,3495.);     		    
  p_nUsedClust_per_LS  = fs->make<TProfile>("nUsesStripClustersVsLS","n. of used Strip clusters vs LS;LS; #LT n.clusters/event #GT",2500,0,2500);    			    
  p_nUsedClust_vs_lumi = fs->make<TProfile>("nUsedStripClustersVsLumi","n. of used Strip clusters vs inst. lumi;scal inst lumi E30 [Hz cm^{-2}]; #LT n.clusters/event #GT",100,0.,15000.);   

    TFileDirectory tfd = fs->mkdir("trackControl");

  vTrackHistos_.push_back(tfd.make<TH1F>("h_tracketa",
					 "Track #eta;#eta_{Track};Number of Tracks",
					 90,-3.,3.));

  vTrackHistos_.push_back(tfd.make<TH1F>("h_tracketaerr",
					 "Track #eta error;err(#eta); Number of Tracks",
					 100,0.,0.1));

  vTrackHistos_.push_back(tfd.make<TH1F>("h_trackphi",
					 "Track #phi;#phi_{Track};Number of Tracks",
					 90,-3.15,3.15));
  
  vTrackHistos_.push_back(tfd.make<TH1F>("h_trackphierr",
					 "Track #phi error;err(#phi) [rad]; Number of Tracks",
					 100,0.,0.1));
  
  vTrackHistos_.push_back(tfd.make<TH1F>("h_trackNumberOfValidHits",
					 "Track # of valid hits;# of valid hits _{Track};Number of Tracks",
					 40,0.,40.));

  vTrackHistos_.push_back(tfd.make<TH1F>("h_trackNumberOfLostHits",
					 "Track # of lost hits;# of lost hits _{Track};Number of Tracks",
					 10,0.,10.));
  
  vTrackHistos_.push_back(tfd.make<TH1F>("h_curvature",
					 "Curvature #kappa;#kappa_{Track} [GeV^{-1}];Number of Tracks",
					 100,-0.5,0.5));

  vTrackHistos_.push_back(tfd.make<TH1F>("h_curvatureerr",
					 "Track Curvature error;err(#kappa) [GeV^{-1}]",
					 100,0.,0.1));

  vTrackHistos_.push_back(tfd.make<TH1F>("h_chi2",
					 "#chi^{2};#chi^{2}_{Track};Number of Tracks",
					 500,-0.01,500.));

  vTrackHistos_.push_back(tfd.make<TH1F>("h_normchi2",
					 "#chi^{2}/ndof;#chi^{2}/ndof;Number of Tracks",
					 100,-0.01,10.));     

  vTrackHistos_.push_back(tfd.make<TH1F>("h_pt",
					 "p_{T}^{track};p_{T}^{track} [GeV];Number of Tracks",
					 250,0.,250));           

  vTrackHistos_.push_back(tfd.make<TH1F>("h_ptResolution",
					 "#delta_{p_{T}}/p_{T}^{track};#delta_{p_{T}}/p_{T}^{track};Number of Tracks",
					 100,0.,0.5)); 

  vTrackHistos_.push_back(tfd.make<TH1F>("h_trackndof",
					 "track # of DOF;track n. DOF;Number of Tracks",
					 30,0.5,30.5)); 
  
  vTrackHistos_.push_back(tfd.make<TH1F>("h_trackcharge",
					 "track charge;track charge;Number of Tracks",
					 3,-1.5,1.5));
  
  vTrackHistos_.push_back(tfd.make<TH1F>("h_tracktheta",
					 "track #theta;track #theta angle;Number of Tracks",
					 100,0.,-3.15));

  vTrackHistos_.push_back(tfd.make<TH1F>("h_trackthetaerr",
					 "track #theta error;err(#theta) [rad]; Number of Tracks",
					 100,0.,0.1));
     
  vTrackHistos_.push_back(tfd.make<TH1F>("h_trackdxy",
					 "Transverse Impact Parameter;d_{xy} [cm]; Number of Tracks",
					 200,-1.,1.));
  
  vTrackHistos_.push_back(tfd.make<TH1F>("h_trackdz",
					 "Longitudinal Impact Parameter;d_{z} [cm]; Number of Tracks",
					 200,-30.,30.));

  vTrackHistos_.push_back(tfd.make<TH1F>("h_trackvx",
					 "track x-coordinate reference point;track v_{x} [cm]; Number of Tracks",
					 200,-1.,1.));
  
  vTrackHistos_.push_back(tfd.make<TH1F>("h_trackvy",
					 "track y-coordinate reference point;track v_{y} [cm]; Number of Tracks",
					 200,-1.,1.));

  vTrackHistos_.push_back(tfd.make<TH1F>("h_trackvz",
					 "track z-coordinate reference point;track v_{z} [cm]; Number of Tracks",
					 200,-10.,10.));

  vTrackHistos_.push_back(tfd.make<TH1F>("h_trackdxyerr",
					 "Transverse Impact Parameter error;err(d_{xy}) [cm]; Number of Tracks",
					 100,0.,0.2));
  
  vTrackHistos_.push_back(tfd.make<TH1F>("h_trackdzerr",
					 "Longitudinal Impact Parameter error;err(d_{z}) [cm]; Number of Tracks",
					 100,0.,2.));
  
  vTrackHistos_.push_back(tfd.make<TH1F>("h_trackalgo",
					 "tracking algorithm;tracking algorithm; Number of Tracks",
					 15,-0.5,14.5));

}

// ------------ method called once each job just after ending the event loop  ------------
void 
SiStripClusterAnalyzer::endJob() 
{
  std::cout<<"bx:"<<maxBx<<" ls:"<<maxLS<<" instLumi:"<<maxLumi<<" PU:"<<maxPU<<std::endl;

  // for ( int fooInt = SiStripHelper::TIBL1; fooInt != SiStripHelper::NUM_OF_TYPES; fooInt++ ){
  //   SiStripHelper::layer layer = static_cast<SiStripHelper::layer>(fooInt);
  //   fitStoN( h_StoNCorr_layer[layer]);
  // }
}

// ------------ method to get the detector thickness ------------
//****************************************************************/
double SiStripClusterAnalyzer::thickness(DetId id)
//****************************************************************/
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

//****************************************************************/
bool SiStripClusterAnalyzer::IsFarFromBorder(TrajectoryStateOnSurface* trajState, const uint32_t detid, const edm::EventSetup* iSetup)
//****************************************************************/
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

//****************************************************************/
int SiStripClusterAnalyzer::checkLayer( unsigned int iidd, const TrackerTopology* tTopo) 
//****************************************************************/
{
  StripSubdetector strip=StripSubdetector(iidd);
  unsigned int subid=strip.subdetId();
  if (subid ==  StripSubdetector::TIB) { 
    
    return tTopo->tibLayer(iidd);
  }
  if (subid ==  StripSubdetector::TOB) { 
    
    return tTopo->tobLayer(iidd) + 4 ; 
  }
  if (subid ==  StripSubdetector::TID) { 
    
    return tTopo->tidWheel(iidd) + 10;
  }
  if (subid ==  StripSubdetector::TEC) { 
    
    return tTopo->tecWheel(iidd) + 13 ; 
  }
  return 0;
}

//****************************************************************/
template <class OBJECT_TYPE>  
int SiStripClusterAnalyzer::GetIndex(const std::vector<OBJECT_TYPE*> &vec, const TString &name)
//****************************************************************/
{
  int result = 0;
  for (typename std::vector<OBJECT_TYPE*>::const_iterator iter = vec.begin(), iterEnd = vec.end();
       iter != iterEnd; ++iter, ++result) {
    if (*iter && (*iter)->GetName() == name) return result;
  }
  edm::LogError("SiStripGainsValidator") << "@SUB=SiStripGainsValidator::GetIndex" << " could not find " << name;
  return -1;
}

// -------------- method to get the topology from the detID ------------------------------
//****************************************************************/
std::string SiStripClusterAnalyzer::getStringFromEnum(SiStripHelper::layer e)
//****************************************************************/
{
  switch(e)
    {
    case SiStripHelper::TIBL1:  return "TIB L1";
    case SiStripHelper::TIBL2:  return "TIB L2";
    case SiStripHelper::TIBL3:  return "TIB L3";
    case SiStripHelper::TIBL4:  return "TIB L4";
    case SiStripHelper::TOBL1:  return "TOB L1";
    case SiStripHelper::TOBL2:  return "TOB L2";
    case SiStripHelper::TOBL3:  return "TOB L3";
    case SiStripHelper::TOBL4:  return "TOB L4";
    case SiStripHelper::TOBL5:  return "TOB L5";
    case SiStripHelper::TOBL6:  return "TOB L6";
    case SiStripHelper::TIDD1:  return "TID Disk 1";
    case SiStripHelper::TIDD2:  return "TID Disk 2";
    case SiStripHelper::TIDD3:  return "TID Disk 3";
    case SiStripHelper::TECD1:  return "TEC Disk 1";
    case SiStripHelper::TECD2:  return "TEC Disk 2";
    case SiStripHelper::TECD3:  return "TEC Disk 3";  
    case SiStripHelper::TECD4:  return "TEC Disk 4";
    case SiStripHelper::TECD5:  return "TEC Disk 5";
    case SiStripHelper::TECD6:  return "TEC Disk 6";
    case SiStripHelper::TECD7:  return "TEC Disk 7";
    case SiStripHelper::TECD8:  return "TEC Disk 8";
    case SiStripHelper::TECD9:  return "TEC Disk 9";
    default: 
      edm::LogWarning("LogicError") << "Unknown partition: " <<  e;
      return "";
    }
}

//****************************************************************/
std::string SiStripClusterAnalyzer::myreplace(const std::string &s,const std::string &toReplace,const std::string &replaceWith)
//****************************************************************/
{
  std::string replacement=s;
  return(replacement.replace(replacement.find(toReplace),toReplace.length(),replaceWith));
}

//****************************************************************/
Double_t SiStripClusterAnalyzer::function_sum(Double_t *x, Double_t *par) 
//****************************************************************/
{
    const Double_t xx =x[0];
    return (1 - par[0]) * SiStripHelper::f2->Eval(xx) / par[1] + (par[0]) * SiStripHelper::f3->Eval(xx) / par[2];
    //return (par[0]) * f2->Eval(xx) + (1 - par[0]) * f3->Eval(xx);
}


//****************************************************************/
void SiStripClusterAnalyzer::fitStoN(TH1F *hist)
//****************************************************************/
{
  std::cout << "## Fitting TH1 histograms ##" << std::endl;
  if (!m_verbose_fit) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  }

  TCanvas *c1 = new TCanvas();
  // Define common things for the different fits
  
  c1->Clear();

  c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.10);

  Double_t xmin = 0;
  Double_t xmax = 100;
  RooRealVar StoN("StoN", "cluster S/N", xmin, xmax);
  RooPlot* frame = StoN.frame();
  RooDataHist datahist("datahist", "datahist", StoN, RooFit::Import(*hist));
  datahist.plotOn(frame);
  
  // Try Landau convolved with a gaussian
  RooRealVar meanG("#mu_{gauss}", "#mu_{gauss}", 0);
  RooRealVar sigmaG("#sigma_{gauss}", "#sigma_{gauss}", 1, 0.1, 20);
  RooRealVar meanL("MPV", "MPV", 30, 10, 50);
  RooRealVar sigmaL("#sigma_{landau}", "#sigma_{landau}", 0.1, 100);
  RooGaussian gaussian("gaussian", "gaussian", StoN, meanG, sigmaG);
  RooLandau landau("landau", "landau", StoN, meanL, sigmaL);
  // Set #bins to be used for FFT sampling to 10000
  StoN.setBins(10000, "cache");
  RooFFTConvPdf model("model", "landau (X) gauss", StoN, landau, gaussian) ;
  if (m_verbose_fit) {
    model.fitTo(datahist,RooFit::Range(12.,40.));
  } else {
    model.fitTo(datahist,RooFit::Range(12.,40.),RooFit::PrintLevel(-1));
  }
  model.plotOn(frame, RooFit::LineColor(kRed));
  // Get the maxima
  // lxg == landau (X) gauss
  TF1 *lxg = model.asTF(RooArgList(StoN));

  model.paramOn(frame,RooFit::Layout(0.55,0.9,0.8)) ;
  Double_t xmax_lxg = lxg->GetMaximumX();
  Double_t ymax_lxg = lxg->GetMaximum();
  Double_t x1_lxg = lxg->GetX(ymax_lxg / 2., xmin, xmax_lxg);
  Double_t x2_lxg = lxg->GetX(ymax_lxg / 2., xmax_lxg, xmax);
  if (HIPDEBUG) {
    std::cout << "From TF1, maximum at (x, y) = (" << xmax_lxg << " , " << ymax_lxg << ")" << std::endl;
    std::cout << "From TF1, FWHM at (x1, x2) = (" << x1_lxg << " , " << x2_lxg << ")" << std::endl;
  }
  
  makeNicePlotStyle(frame);

  // Try landau + a gaussian
  /*
    RooRealVar meanG2("meanG2", "meanG2", 30, 10, 50);
    RooRealVar sigmaG2("sigmaG2", "sigmaG2", 1, 1, 50);
    RooRealVar meanL2("meanL2", "meanL2", 30, 10, 50);
    RooRealVar sigmaL2("sigmaL2", "sigmaL2", 0.1, 100);
    RooGaussian gaussian2("gaussian2", "gaussian2", StoN, meanG2, sigmaG2);
    RooLandau landau2("landau2", "landau2", StoN, meanL2, sigmaL2);
    RooRealVar x("x", "x", 0.1, 0, 0.4);
    RooAddPdf model2("model2", "model2", gaussian2, landau2, x);
    if (m_verbose_fit) {
      model2.fitTo(datahist);
    } else {
      model2.fitTo(datahist, RooFit::PrintLevel(-1));
    }
    model2.plotOn(frame, RooFit::LineColor(kRed));
    // Get the maxima
    // lpg == landau + gauss
    f2 = landau2.asTF(RooArgList(StoN));
    f3 = gaussian2.asTF(RooArgList(StoN));
    TF1 *lpg = new TF1("Sum", function_sum, 0, 100, 3);
    lpg->SetParameter(0, x.getVal());
    lpg->SetParameter(1, f2->Integral(0, 100));
    lpg->SetParameter(2, f3->Integral(0, 100));
    Double_t xmax_lpg = lpg->GetMaximumX();
    Double_t ymax_lpg = lpg->GetMaximum();
    Double_t x1_lpg = lpg->GetX(ymax_lpg / 2., xmin, xmax_lpg);
    Double_t x2_lpg = lpg->GetX(ymax_lpg / 2., xmax_lpg, xmax);
    if (HIPDEBUG) {
      std::cout << "From TF1, maximum at (x, y) = (" << xmax_lpg << " , " << ymax_lpg << ")" << std::endl;
      std::cout << "From TF1, FWHM at (x1, x2) = (" << x1_lpg << " , " << x2_lpg << ")" << std::endl;
      }
  */    

  // Redraw data on top and print / store everything
  datahist.plotOn(frame);
  frame->GetYaxis()->SetTitle("n. of on-track clusters");
  TString histName = hist->GetName();
  frame->SetName("frame"+histName);
  frame->SetTitle(hist->GetTitle());
  frame->Draw();
  //    m_output->cd();
  //frame->Write();
  if (HIPDEBUG) {
    c1->Print("fit_debug"+histName+".pdf");
  }
  
  delete lxg;
  //delete lpg;
  delete c1;
}

/*--------------------------------------------------------------------*/
void SiStripClusterAnalyzer::makeNicePlotStyle(RooPlot* plot)
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


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SiStripClusterAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SiStripClusterAnalyzer);
