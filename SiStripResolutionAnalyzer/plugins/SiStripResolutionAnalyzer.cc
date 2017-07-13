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

#include "Alignment/OfflineValidation/interface/TrackerValidationVariables.h"
#include "Alignment/TrackerAlignment/interface/TrackerAlignableId.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "CommonTools/UtilAlgos/interface/DetIdSelector.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CondFormats/Alignment/interface/Definitions.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/CommonTopologies/interface/RadialStripTopology.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementError.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"
#include "DataFormats/GeometrySurface/interface/RectangularPlaneBounds.h"
#include "DataFormats/GeometrySurface/interface/TrapezoidalPlaneBounds.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajectoryMeasurement.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryStateCombiner.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"

#include "SiStripMiscellanea/SiStripResolutionAnalyzer/interface/SiStripResolutionHelper.h"

#include "TMath.h"
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
      void fillHitQuantities(const Trajectory* trajectory, std::vector<TrackerValidationVariables::AVHitStruct> & v_avhitout, bool isBiased);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      TH1F* m_ptrk;
      TH1F* m_etatrk;
      TH1F* m_nhits;
      TH1F* m_chi2Prob;

      edm::EDGetTokenT<std::vector<Trajectory>> trajCollectionToken_;
      edm::EDGetTokenT<std::vector<reco::Track>> tracksToken_;
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
  trajCollectionToken_ (consumes<std::vector<Trajectory> >(edm::InputTag(iConfig.getParameter<std::string>("trajectoryInput")))),
  tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("Tracks")))
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
  
  /*
    edm::ESHandle<MagneticField> magneticField;
    iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
  */  

  // retrieve tracker topology from geometry
  edm::ESHandle<TrackerTopology> tTopoHandle;
  iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);
  const TrackerTopology* tTopo = tTopoHandle.product();
   
  // get the tracks
  edm::Handle<reco::TrackCollection> tracksH;
  iEvent.getByToken(tracksToken_, tracksH);

  if(!tracksH.isValid()) return;
  auto const & tracks = *tracksH;
  auto ntrk = tracks.size();
  LogDebug("TrackerValidationVariables") << "Track collection size " << ntrk;
  
  edm::Handle<std::vector<Trajectory>> trajsH;
  iEvent.getByToken(trajCollectionToken_, trajsH);
  bool yesTraj = trajsH.isValid();
  std::vector<Trajectory> const * trajs = nullptr;
  if (yesTraj) trajs = &(*trajsH);
  if (yesTraj) assert (trajs->size()==tracks.size());

  Trajectory const * traj =  nullptr;
  for (unsigned int i=0; i<ntrk; ++i) {
    auto const & track = tracks[i];  
    if (yesTraj) traj = &(*trajs)[i];

    double ProbChi2 = ChiSquaredProbability((double)( traj->chiSquared() ),(double)( traj->ndof(false) ));

    if (traj->foundHits() < 6 ) continue;
    if (ProbChi2 < 0.001 ) continue;
    if (track.pt()<3.) continue;

    m_ptrk->Fill(track.pt());
    m_etatrk->Fill(track.eta());
    m_nhits->Fill(traj->foundHits());
    m_chi2Prob->Fill(ProbChi2);

    std::vector<TrackerValidationVariables::AVHitStruct> vhits_unbiased;
    std::vector<TrackerValidationVariables::AVHitStruct> vhits_biased;

    fillHitQuantities(traj,vhits_biased,true);
    fillHitQuantities(traj,vhits_unbiased,false);

    assert(vhits_biased.size()==vhits_unbiased.size());

    for (auto& it : vhits_biased) {
      auto id = DetId(it.rawDetId);
      auto isStrip = id.subdetId() > 2; 
      if (!isStrip) continue; 

      // auto packedTopo = SiStripResol::typeAndLayerFromDetId(id,tTopo);
      // std::cout<<" packedTopo =("<< packedTopo.second.first<<","<<packedTopo.second.first<<")"<<std::endl;

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

// ------------ method called once each job just after ending the event loop  ------------
void
SiStripResolutionAnalyzer::fillHitQuantities(const Trajectory* trajectory, std::vector<TrackerValidationVariables::AVHitStruct> & v_avhitout, bool isBiased)
{
  TrajectoryStateCombiner tsoscomb;

  const std::vector<TrajectoryMeasurement> &tmColl = trajectory->measurements();
  for(std::vector<TrajectoryMeasurement>::const_iterator itTraj = tmColl.begin();
      itTraj != tmColl.end();
      ++itTraj) {
    
    if (!itTraj->updatedState().isValid()) continue;

    TrajectoryStateOnSurface tsos;

    if(isBiased){
      tsos = itTraj->updatedState();
    } else {
      tsos = tsoscomb( itTraj->forwardPredictedState(), itTraj->backwardPredictedState());
    }

    if(!tsos.isValid()) continue;
    TransientTrackingRecHit::ConstRecHitPointer hit = itTraj->recHit();
    
    if(!hit->isValid() || hit->geographicalId().det() != DetId::Tracker) continue;
    
    TrackerValidationVariables::AVHitStruct hitStruct;
    const DetId& hit_detId = hit->geographicalId();
    unsigned int IntRawDetID = (hit_detId.rawId());	
    unsigned int IntSubDetID = (hit_detId.subdetId());
    
    if(IntSubDetID == 0) continue;
    
    //first calculate residuals in cartesian coordinates in the local module coordinate system
    
    LocalPoint lPHit = hit->localPosition();
    LocalPoint lPTrk = tsos.localPosition();
    LocalVector lVTrk = tsos.localDirection();

    hitStruct.localAlpha = atan2(lVTrk.x(), lVTrk.z()); // wrt. normal tg(alpha)=x/z
    hitStruct.localBeta  = atan2(lVTrk.y(), lVTrk.z()); // wrt. normal tg(beta)= y/z

    LocalError errHit = hit->localPositionError();
    // no need to add  APE to hitError anymore
    // AlgebraicROOTObject<2>::SymMatrix mat = asSMatrix<2>(hit->parametersError());
    // LocalError errHit = LocalError( mat(0,0),mat(0,1),mat(1,1) );
    LocalError errTrk = tsos.localError().positionError();
    
    //check for negative error values: track error can have negative value, if matrix inversion fails (very rare case)
    //hit error should always give positive values
    if(errHit.xx()<0. || errHit.yy()<0. || errTrk.xx()<0. || errTrk.yy()<0.){
      edm::LogError("SiStripResolutionAnalyzer") << "@SUB=SiStripResolutionAnalyzer::fillHitQuantities"
						  << "One of the squared error methods gives negative result"
						  << "\n\terrHit.xx()\terrHit.yy()\terrTrk.xx()\terrTrk.yy()"
						  << "\n\t" << errHit.xx()
						  << "\t" << errHit.yy()
						  << "\t" << errTrk.xx()
						  << "\t" << errTrk.yy();
      continue;
    }
    
    align::LocalVector res = lPTrk - lPHit;
    
    float resXErr = std::sqrt( errHit.xx() + errTrk.xx() );
    float resYErr = std::sqrt( errHit.yy() + errTrk.yy() );
    
    hitStruct.resX = res.x();
    hitStruct.resY = res.y();
    hitStruct.resErrX = resXErr;
    hitStruct.resErrY = resYErr;

    // hitStruct.localX = lPhit.x();
    // hitStruct.localY = lPhit.y();
    // EM: use predictions for local coordinates
    hitStruct.localX = lPTrk.x();
    hitStruct.localY = lPTrk.y();

    // now calculate residuals taking global orientation of modules and radial topology in TID/TEC into account
    float resXprime(999.F), resYprime(999.F);
    float resXatTrkY(999.F);
    float resXprimeErr(999.F), resYprimeErr(999.F);
    
    if(hit->detUnit()){ // is it a single physical module?
      const GeomDetUnit& detUnit = *(hit->detUnit());
      float uOrientation(-999.F), vOrientation(-999.F);
      float resXTopol(999.F), resYTopol(999.F);
      float resXatTrkYTopol(999.F);       

      const Surface& surface = hit->detUnit()->surface();
      const BoundPlane& boundplane = hit->detUnit()->surface();
      const Bounds& bound = boundplane.bounds();
      
      float length = 0;
      float width = 0;

      LocalPoint lPModule(0.,0.,0.), lUDirection(1.,0.,0.), lVDirection(0.,1.,0.);
      GlobalPoint gPModule    = surface.toGlobal(lPModule),
	          gUDirection = surface.toGlobal(lUDirection),
	          gVDirection = surface.toGlobal(lVDirection);
      
      if (IntSubDetID == PixelSubdetector::PixelBarrel ||
	  IntSubDetID == StripSubdetector::TIB || 
	  IntSubDetID == StripSubdetector::TOB) {
	
	uOrientation = deltaPhi(gUDirection.barePhi(),gPModule.barePhi()) >= 0. ? +1.F : -1.F;
	vOrientation = gVDirection.z() - gPModule.z() >= 0 ? +1.F : -1.F;
	resXTopol = res.x();
	resXatTrkYTopol = res.x();
	resYTopol = res.y();
	resXprimeErr = resXErr;
	resYprimeErr = resYErr;

	const RectangularPlaneBounds *rectangularBound = dynamic_cast<const RectangularPlaneBounds*>(&bound);
	if (rectangularBound!=NULL) {
	  hitStruct.inside = rectangularBound->inside(lPTrk);
	  length = rectangularBound->length();
	  width = rectangularBound->width();
	  hitStruct.localXnorm = 2*hitStruct.localX/width;
	  hitStruct.localYnorm = 2*hitStruct.localY/length;
	} else {
	  throw cms::Exception("Geometry Error")
	    << "[SiStripResolutionAnalyzer] Cannot cast bounds to RectangularPlaneBounds as expected for TPB, TIB and TOB";
	}

      } else if (IntSubDetID == PixelSubdetector::PixelEndcap) {
	
	uOrientation = gUDirection.perp() - gPModule.perp() >= 0 ? +1.F : -1.F;
	vOrientation = deltaPhi(gVDirection.barePhi(),gPModule.barePhi()) >= 0. ? +1.F : -1.F;
	resXTopol = res.x();
	resXatTrkYTopol = res.x();
	resYTopol = res.y();
	resXprimeErr = resXErr;
	resYprimeErr = resYErr;
	
	const RectangularPlaneBounds *rectangularBound = dynamic_cast<const RectangularPlaneBounds*>(&bound);
	if (rectangularBound!=NULL) {
	  hitStruct.inside = rectangularBound->inside(lPTrk);
	  length = rectangularBound->length();
	  width = rectangularBound->width();
	  hitStruct.localXnorm = 2*hitStruct.localX/width;
	  hitStruct.localYnorm = 2*hitStruct.localY/length;
	} else {
	  throw cms::Exception("Geometry Error")
	    << "[SiStripResolutionAnalyzer] Cannot cast bounds to RectangularPlaneBounds as expected for TPE";
	}

      } else if (IntSubDetID == StripSubdetector::TID ||
		 IntSubDetID == StripSubdetector::TEC) {
	
	uOrientation = deltaPhi(gUDirection.barePhi(),gPModule.barePhi()) >= 0. ? +1.F : -1.F;
	vOrientation = gVDirection.perp() - gPModule.perp() >= 0. ? +1.F : -1.F;
	
	if (!dynamic_cast<const RadialStripTopology*>(&detUnit.type().topology()))continue;
	const RadialStripTopology& topol = dynamic_cast<const RadialStripTopology&>(detUnit.type().topology());
	
	MeasurementPoint measHitPos = topol.measurementPosition(lPHit);
	MeasurementPoint measTrkPos = topol.measurementPosition(lPTrk);
	
	MeasurementError measHitErr = topol.measurementError(lPHit,errHit);
	MeasurementError measTrkErr = topol.measurementError(lPTrk,errTrk);
	
	if (measHitErr.uu()<0. ||
	    measHitErr.vv()<0. ||
	    measTrkErr.uu()<0. ||
	    measTrkErr.vv()<0.){
	  edm::LogError("SiStripResolutionAnalyzer") << "@SUB=SiStripResolutionAnalyzer::fillHitQuantities"
						      << "One of the squared error methods gives negative result"
						      << "\n\tmeasHitErr.uu()\tmeasHitErr.vv()\tmeasTrkErr.uu()\tmeasTrkErr.vv()"
						      << "\n\t" << measHitErr.uu()
						      << "\t" << measHitErr.vv()
						      << "\t" << measTrkErr.uu()
						      << "\t" << measTrkErr.vv();
	  continue;
	}
	  
	float localStripLengthHit = topol.localStripLength(lPHit);
	float localStripLengthTrk = topol.localStripLength(lPTrk);
	float phiHit = topol.stripAngle(measHitPos.x());
	float phiTrk = topol.stripAngle(measTrkPos.x());
	float r_0 = topol.originToIntersection();
	
	resXTopol = (phiTrk-phiHit)*r_0;
	//resXTopol = (tan(phiTrk)-tan(phiHit))*r_0;
        
	LocalPoint LocalHitPosCor= topol.localPosition(MeasurementPoint(measHitPos.x(), measTrkPos.y()));                       
	resXatTrkYTopol = lPTrk.x() - LocalHitPosCor.x();  

    
	//resYTopol = measTrkPos.y()*localStripLengthTrk - measHitPos.y()*localStripLengthHit;
	float cosPhiHit(cos(phiHit)), cosPhiTrk(cos(phiTrk)),
	      sinPhiHit(sin(phiHit)), sinPhiTrk(sin(phiTrk));
	float l_0 = r_0 - topol.detHeight()/2;
	resYTopol = measTrkPos.y()*localStripLengthTrk - measHitPos.y()*localStripLengthHit + l_0*(1/cosPhiTrk - 1/cosPhiHit);
	
	resXprimeErr = std::sqrt(measHitErr.uu()+measTrkErr.uu())*topol.angularWidth()*r_0;
	//resYprimeErr = std::sqrt(measHitErr.vv()*localStripLengthHit*localStripLengthHit + measTrkErr.vv()*localStripLengthTrk*localStripLengthTrk);
	float helpSummand = l_0*l_0*topol.angularWidth()*topol.angularWidth()*(sinPhiHit*sinPhiHit/pow(cosPhiHit,4)*measHitErr.uu()
									       + sinPhiTrk*sinPhiTrk/pow(cosPhiTrk,4)*measTrkErr.uu() );
	resYprimeErr = std::sqrt(measHitErr.vv()*localStripLengthHit*localStripLengthHit
				 + measTrkErr.vv()*localStripLengthTrk*localStripLengthTrk + helpSummand );  


	const TrapezoidalPlaneBounds *trapezoidalBound = dynamic_cast < const TrapezoidalPlaneBounds * >(& bound);
	if (trapezoidalBound!=NULL) {
	  hitStruct.inside = trapezoidalBound->inside(lPTrk);
	  length = trapezoidalBound->length();
	  width  = trapezoidalBound->width();
	  //float widthAtHalfLength = trapezoidalBound->widthAtHalfLength();

	  //	int yAxisOrientation=trapezoidalBound->yAxisOrientation(); 
	  // for trapezoidal shape modules, scale with as function of local y coordinate 
	  //	float widthAtlocalY=width-(1-yAxisOrientation*2*lPTrk.y()/length)*(width-widthAtHalfLength); 
	  //	hitStruct.localXnorm = 2*hitStruct.localX/widthAtlocalY;  
	  hitStruct.localXnorm = 2*hitStruct.localX/width; 
	  hitStruct.localYnorm = 2*hitStruct.localY/length;
	} else {
	  throw cms::Exception("Geometry Error")
	    << "[SiStripResolutionAnalyzer] Cannot cast bounds to TrapezoidalPlaneBounds as expected for TID and TEC";
	}

      } else {
	edm::LogWarning("SiStripResolutionAnalyzer") << "@SUB=SiStripResolutionAnalyzer::fillHitQuantities" 
						      << "No valid tracker subdetector " << IntSubDetID;
	continue;
      }  
      
      resXprime = resXTopol*uOrientation;
      resXatTrkY = resXatTrkYTopol;
      resYprime = resYTopol*vOrientation;
      
    } else { 
      // not a detUnit, so must be a virtual 2D-Module
      // FIXME: at present only for det units residuals are calculated and filled in the hitStruct
      // But in principle this method should also be useable for the gluedDets (2D modules in TIB, TID, TOB, TEC)
      // In this case, only orientation should be taken into account for primeResiduals, but not the radial topology
      // At present, default values (999.F) are given out
    }
    
    hitStruct.resXprime = resXprime;
    hitStruct.resXatTrkY = resXatTrkY;
    hitStruct.resYprime = resYprime;
    hitStruct.resXprimeErr = resXprimeErr;
    hitStruct.resYprimeErr = resYprimeErr;
    
    hitStruct.rawDetId = IntRawDetID;
    hitStruct.phi = tsos.globalDirection().phi();
    hitStruct.eta = tsos.globalDirection().eta();
    
    v_avhitout.push_back(hitStruct);
  } 
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
