// -*- C++ -*-
//
// Package:    SiStripMiscellanea/SiStripGainsValidator
// Class:      SiStripGainsValidator
//
/**\class SiStripGainsValidator SiStripGainsValidator.cc SiStripMiscellanea/SiStripGainsValidator/plugins/SiStripGainsValidator.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marco Musich
//         Created:  Fri, 18 Aug 2017 11:28:24 GMT
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
#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
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
#include "CalibFormats/SiStripObjects/interface/SiStripQuality.h"
#include "CalibTracker/Records/interface/SiStripGainRcd.h"

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"

#include "TF1.h"
#include "TH2.h"
#include "TH1I.h"
#include "TProfile.h"

namespace SiStripHelper {

  enum layer {
    TIBL1 = 1,
    TIBL2 = 2,
    TIBL3 = 3,
    TIBL4 = 4,
    TOBL1 = 5,
    TOBL2 = 6,
    TOBL3 = 7,
    TOBL4 = 8,
    TOBL5 = 9,
    TOBL6 = 10,
    TIDD1 = 11,
    TIDD2 = 12,
    TIDD3 = 13,
    TECD1 = 14,
    TECD2 = 15,
    TECD3 = 16,
    TECD4 = 17,
    TECD5 = 18,
    TECD6 = 19,
    TECD7 = 20,
    TECD8 = 21,
    TECD9 = 22,
    NUM_OF_TYPES = 23,
  };
}

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class SiStripGainsValidator : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit SiStripGainsValidator(const edm::ParameterSet&);
  ~SiStripGainsValidator() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  double thickness(DetId id);
  bool IsFarFromBorder(TrajectoryStateOnSurface* trajState, const uint32_t detid, const edm::EventSetup* iSetup);
  int checkLayer(unsigned int iidd, const TrackerTopology* tTopo);
  std::string getStringFromEnum(SiStripHelper::layer e);
  std::string myreplace(const std::string& s, const std::string& toReplace, const std::string& replaceWith);
  template <class OBJECT_TYPE>
  int GetIndex(const std::vector<OBJECT_TYPE*>& vec, const TString& name);
  void getPeakOfLandau(TH1* InputHisto, double* FitResults, double LowRange, double HighRange);

  edm::EDGetTokenT<edm::View<reco::Track> > tracks_token_;
  edm::EDGetTokenT<TrajTrackAssociationCollection> association_token_;

  SiStripClusterInfo siStripClusterInfo_;

  // ----------member data ---------------------------

  double MinTrackMomentum;
  double MaxTrackMomentum;
  double MinTrackEta;
  double MaxTrackEta;
  unsigned int MaxNrStrips;
  unsigned int MinTrackHits;
  double MaxTrackChiOverNdf;
  int MaxTrackingIteration;
  bool AllowSaturation;

  const TrackerGeometry* m_tracker;
  std::map<DetId, double> m_thicknessMap;

  edm::Service<TFileService> fs;

  TH1F* h_nClusters;
  TH1F* h_nUsedClusters;

  TH1F* h_APVGain;
  TH1F* h_SumStripGain;
  TH1F* h_diffAPVfromClusterGain;
  TH1F* h_diffAPVfromClusterGainNoOverlap;

  TH1F* h_chargeFromStrips;
  TH1F* h_diffChargeMethod;

  TH1F* h_CChargeOverPath;
  TH1F* h_CChargeOverPathNewG2;
  TH1F* h_CChargeOverPathNoG2;
  TH1F* h_CChargeOverPathNoG1;
  TH1F* h_CChargeOverPathNoG1G2;

  TH1F* h_StoN;
  TH1F* h_StoNCorr;

  TH1F* h_Noise;
  TH1F* h_chargeFromClusterInfo;
  TH1F* h_clusterwidth;
  TH1F* h_clusterposition;

  TProfile* p_ClusterChargeVsLS;
  TH2F* h2_ClusterChargeVsLS;

  std::vector<TH1*> vTrackHistos_;

  std::map<SiStripHelper::layer, TH1F*> h_StoN_layer;
  std::map<SiStripHelper::layer, TH1F*> h_StoNCorr_layer;

  std::map<SiStripHelper::layer, TH1F*> h_CChargeOverPath_layer;
  std::map<SiStripHelper::layer, TH1F*> h_CChargeOverPathNewG2_layer;
  std::map<SiStripHelper::layer, TH1F*> h_CChargeOverPathNoG2_layer;
  std::map<SiStripHelper::layer, TH1F*> h_CChargeOverPathNoG1_layer;
  std::map<SiStripHelper::layer, TH1F*> h_CChargeOverPathNoG1G2_layer;
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
SiStripGainsValidator::SiStripGainsValidator(const edm::ParameterSet& iConfig)
    : siStripClusterInfo_(consumesCollector()) {
  //now do what ever initialization is needed
  usesResource("TFileService");

  tracks_token_ = consumes<edm::View<reco::Track> >(iConfig.getParameter<edm::InputTag>("Tracks"));
  association_token_ = consumes<TrajTrackAssociationCollection>(iConfig.getParameter<edm::InputTag>("Tracks"));

  MinTrackMomentum = iConfig.getUntrackedParameter<double>("minTrackMomentum", 3.0);
  MaxTrackMomentum = iConfig.getUntrackedParameter<double>("maxTrackMomentum", 99999.0);
  MinTrackEta = iConfig.getUntrackedParameter<double>("minTrackEta", -5.0);
  MaxTrackEta = iConfig.getUntrackedParameter<double>("maxTrackEta", 5.0);
  MaxNrStrips = iConfig.getUntrackedParameter<unsigned>("maxNrStrips", 2);
  MinTrackHits = iConfig.getUntrackedParameter<unsigned>("MinTrackHits", 8);
  MaxTrackChiOverNdf = iConfig.getUntrackedParameter<double>("MaxTrackChiOverNdf", 3);
  MaxTrackingIteration = iConfig.getUntrackedParameter<int>("MaxTrackingIteration", 7);
  AllowSaturation = iConfig.getUntrackedParameter<bool>("AllowSaturation", false);
}

SiStripGainsValidator::~SiStripGainsValidator() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void SiStripGainsValidator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  // initialize SiStripClusterInfo
  siStripClusterInfo_.initEvent(iSetup);

  //  int bx = iEvent.bunchCrossing();
  int ls = iEvent.id().luminosityBlock();

  unsigned int nStripClus = 0;
  unsigned int nUsedStripClus = 0;

  // Get the geometry
  edm::ESHandle<TrackerTopology> tTopo_handle;
  iSetup.get<TrackerTopologyRcd>().get(tTopo_handle);
  const TrackerTopology* tTopo = tTopo_handle.product();

  edm::ESHandle<TrackerGeometry> theTrackerGeometry;
  iSetup.get<TrackerDigiGeometryRecord>().get(theTrackerGeometry);
  m_tracker = &(*theTrackerGeometry);
  edm::ESHandle<SiStripGain> gainHandle;
  iSetup.get<SiStripGainRcd>().get(gainHandle);

  edm::ESHandle<SiStripQuality> qualityHandle;
  iSetup.get<SiStripQualityRcd>().get("", qualityHandle);
  const SiStripQuality* stripQuality = qualityHandle.product();

  // gain to be validated
  edm::ESHandle<SiStripApvGain> g3Handle;
  iSetup.get<SiStripApvGain3Rcd>().get(g3Handle);

  edm::Handle<edm::View<reco::Track> > tracks;
  iEvent.getByToken(tracks_token_, tracks);
  edm::Handle<TrajTrackAssociationCollection> associations;
  iEvent.getByToken(association_token_, associations);

  for (TrajTrackAssociationCollection::const_iterator association = associations->begin();
       association != associations->end();
       association++) {
    const Trajectory* traj = association->key.get();
    const reco::Track* track = association->val.get();

    double trackchi2 = track->chi2();
    double trackndof = track->ndof();
    double trackchi2ndof = track->chi2() / track->ndof();
    float trackcharge = track->charge();
    float trackmomentum = track->p();
    float trackpt = track->pt();
    float trackpterr = track->ptError();
    unsigned int trackhitsvalid = track->numberOfValidHits();
    unsigned int trackhitslost = track->numberOfLostHits();
    double tracktheta = track->theta();
    double trackthetaerr = track->thetaError();
    double trackphi = track->phi();
    double trackphierr = track->phiError();
    double tracketa = track->eta();
    double tracketaerr = track->etaError();
    double trackdxy = track->dxy();
    double trackdxyerr = track->dxyError();
    double trackdz = track->dz();
    double trackdzerr = track->dzError();
    double trackqoverp = track->qoverp();
    double trackqoverperr = track->qoverpError();
    double trackvx = track->vx();
    double trackvy = track->vy();
    double trackvz = track->vz();
    int trackalgo = (int)track->algo();

    if (tracketa < MinTrackEta)
      continue;
    if (tracketa > MaxTrackEta)
      continue;
    if (trackmomentum < MinTrackMomentum)
      continue;
    if (trackmomentum > MaxTrackMomentum)
      continue;
    if (trackhitsvalid < MinTrackHits)
      continue;
    if (trackchi2ndof > MaxTrackChiOverNdf)
      continue;
    if (trackalgo > MaxTrackingIteration)
      continue;

    // fill control plots
    static const int chi2index = this->GetIndex(vTrackHistos_, "h_chi2");
    vTrackHistos_[chi2index]->Fill(trackchi2);

    static const int ndofindex = this->GetIndex(vTrackHistos_, "h_trackndof");
    vTrackHistos_[ndofindex]->Fill(trackndof);

    static const int chargeindex = this->GetIndex(vTrackHistos_, "h_trackcharge");
    vTrackHistos_[chargeindex]->Fill(trackcharge);

    static const int thetaindex = this->GetIndex(vTrackHistos_, "h_tracktheta");
    vTrackHistos_[thetaindex]->Fill(tracktheta);

    static const int dxyindex = this->GetIndex(vTrackHistos_, "h_trackdxy");
    vTrackHistos_[dxyindex]->Fill(trackdxy);

    static const int dzindex = this->GetIndex(vTrackHistos_, "h_trackdz");
    vTrackHistos_[dzindex]->Fill(trackdz);

    static const int dxyerrindex = this->GetIndex(vTrackHistos_, "h_trackdxyerr");
    vTrackHistos_[dxyerrindex]->Fill(trackdxyerr);

    static const int dzerrindex = this->GetIndex(vTrackHistos_, "h_trackdzerr");
    vTrackHistos_[dzerrindex]->Fill(trackdzerr);

    static const int vxindex = this->GetIndex(vTrackHistos_, "h_trackvx");
    vTrackHistos_[vxindex]->Fill(trackvx);

    static const int vyindex = this->GetIndex(vTrackHistos_, "h_trackvy");
    vTrackHistos_[vyindex]->Fill(trackvy);

    static const int vzindex = this->GetIndex(vTrackHistos_, "h_trackvz");
    vTrackHistos_[vzindex]->Fill(trackvz);

    static const int thetaerrindex = this->GetIndex(vTrackHistos_, "h_trackthetaerr");
    vTrackHistos_[thetaerrindex]->Fill(trackthetaerr);

    static const int phierrindex = this->GetIndex(vTrackHistos_, "h_trackphierr");
    vTrackHistos_[phierrindex]->Fill(trackphierr);

    static const int kappaerrindex = this->GetIndex(vTrackHistos_, "h_curvatureerr");
    vTrackHistos_[kappaerrindex]->Fill(trackqoverperr);

    static const int algoindex = this->GetIndex(vTrackHistos_, "h_trackalgo");
    vTrackHistos_[algoindex]->Fill(trackalgo);

    static const int normchi2index = this->GetIndex(vTrackHistos_, "h_normchi2");
    vTrackHistos_[normchi2index]->Fill(trackchi2ndof);
    static const int ptindex = this->GetIndex(vTrackHistos_, "h_pt");
    if (trackpterr != 0.) {
      static const int ptResolutionindex = this->GetIndex(vTrackHistos_, "h_ptResolution");
      vTrackHistos_[ptResolutionindex]->Fill(trackpterr / trackpt);
    }
    vTrackHistos_[ptindex]->Fill(trackpt);

    static const int etaindex = this->GetIndex(vTrackHistos_, "h_tracketa");
    vTrackHistos_[etaindex]->Fill(tracketa);

    static const int etaerrindex = this->GetIndex(vTrackHistos_, "h_tracketaerr");
    vTrackHistos_[etaerrindex]->Fill(tracketaerr);

    static const int phiindex = this->GetIndex(vTrackHistos_, "h_trackphi");
    vTrackHistos_[phiindex]->Fill(trackphi);

    static const int numOfValidHitsindex = this->GetIndex(vTrackHistos_, "h_trackNumberOfValidHits");
    vTrackHistos_[numOfValidHitsindex]->Fill(trackhitsvalid);

    static const int numOfLostHitsindex = this->GetIndex(vTrackHistos_, "h_trackNumberOfLostHits");
    vTrackHistos_[numOfLostHitsindex]->Fill(trackhitslost);

    static const int kappaindex = this->GetIndex(vTrackHistos_, "h_curvature");
    vTrackHistos_[kappaindex]->Fill(trackqoverp);

    std::vector<TrajectoryMeasurement> measurements = traj->measurements();
    for (std::vector<TrajectoryMeasurement>::const_iterator measurement_it = measurements.begin();
         measurement_it != measurements.end();
         measurement_it++) {
      TrajectoryStateOnSurface trajState = measurement_it->updatedState();
      if (!trajState.isValid())
        continue;

      const TrackingRecHit* hit = (*measurement_it->recHit()).hit();
      const SiStripRecHit1D* sistripsimple1dhit = dynamic_cast<const SiStripRecHit1D*>(hit);
      const SiStripRecHit2D* sistripsimplehit = dynamic_cast<const SiStripRecHit2D*>(hit);
      const SiStripMatchedRecHit2D* sistripmatchedhit = dynamic_cast<const SiStripMatchedRecHit2D*>(hit);
      const SiPixelRecHit* sipixelhit = dynamic_cast<const SiPixelRecHit*>(hit);

      const uint32_t& detid = hit->geographicalId().rawId();
      auto layer = checkLayer(detid, tTopo);
      SiStripHelper::layer myLayer = static_cast<SiStripHelper::layer>(layer);

      const SiPixelCluster* PixelCluster = nullptr;
      const SiStripCluster* StripCluster = nullptr;
      uint32_t DetId = 0;

      for (unsigned int h = 0; h < 2; h++) {
        if (!sistripmatchedhit && h == 1) {
          continue;
        } else if (sistripmatchedhit && h == 0) {
          StripCluster = &sistripmatchedhit->monoCluster();
          DetId = sistripmatchedhit->monoId();
        } else if (sistripmatchedhit && h == 1) {
          StripCluster = &sistripmatchedhit->stereoCluster();
          ;
          DetId = sistripmatchedhit->stereoId();
        } else if (sistripsimplehit) {
          StripCluster = (sistripsimplehit->cluster()).get();
          DetId = sistripsimplehit->geographicalId().rawId();
        } else if (sistripsimple1dhit) {
          StripCluster = (sistripsimple1dhit->cluster()).get();
          DetId = sistripsimple1dhit->geographicalId().rawId();
        } else if (sipixelhit) {
          PixelCluster = (sipixelhit->cluster()).get();
          DetId = sipixelhit->geographicalId().rawId();
        } else {
          continue;
        }

        LocalVector trackDirection = trajState.localDirection();
        double cosine = trackDirection.z() / trackDirection.mag();
        double cosRZ = fabs(trackDirection.z()) / trackDirection.mag();
        bool Saturation = false;
        bool Overlapping = false;
        unsigned int Charge = 0;
        double Path = (10.0 * thickness(DetId)) / fabs(cosine);
        double PrevGain = -1;
        double PrevGainTick = -1;
        double NewGain = -1;
        int FirstStrip = 0;
        unsigned int NStrips = 0;

        std::vector<unsigned char> amplitude;

        if (StripCluster) {
          unsigned int FirstAmplitude = 0;
          nStripClus++;

          const auto& Ampls = StripCluster->amplitudes();
          FirstStrip = StripCluster->firstStrip();
          NStrips = Ampls.size();
          int APVId = FirstStrip / 128;

          FirstAmplitude += NStrips;

          siStripClusterInfo_.setCluster(*StripCluster, detid);
          float StoN = siStripClusterInfo_.signalOverNoise();
          float noise = siStripClusterInfo_.noiseRescaledByGain();
          uint16_t charge = siStripClusterInfo_.charge();
          uint16_t width = siStripClusterInfo_.width();
          float position = siStripClusterInfo_.baryStrip();

          h_Noise->Fill(noise);
          h_chargeFromClusterInfo->Fill(charge);
          h_clusterwidth->Fill(width);
          h_clusterposition->Fill(position);

          h_StoN->Fill(StoN);
          h_StoNCorr->Fill(StoN * cosRZ);

          h_StoN_layer[myLayer]->Fill(StoN);
          h_StoNCorr_layer[myLayer]->Fill(StoN * cosRZ);

          if (gainHandle.isValid()) {
            PrevGain = gainHandle->getApvGain(APVId, gainHandle->getRange(DetId, 1), 1);
            PrevGainTick = gainHandle->getApvGain(APVId, gainHandle->getRange(DetId, 0), 1);
          }

          // take from the G3 handle the new gain
          if (g3Handle.isValid()) {
            NewGain = g3Handle->getApvGain(APVId, g3Handle->getRange(DetId));
          }

          //std::cout<<"Found in the ES:"<<gainHandle->getNumberOfTags()<<std::endl;
          //std::cout<<"PrevGain: "<<PrevGain<<" PrevGainTick: "<<PrevGainTick<< " NewGain:"<< NewGain <<std::endl;

          for (unsigned int a = 0; a < Ampls.size(); a++) {
            Charge += Ampls[a];
            if (Ampls[a] >= 254)
              Saturation = true;
            amplitude.push_back(Ampls[a]);
          }

          h_chargeFromStrips->Fill(Charge);
          //h_diffChargeMethod->Fill(Charge-charge));

          // Fill monitoring histograms
          int mCharge1 = 0;
          int mCharge2 = 0;
          int mCharge3 = 0;
          int mCharge4 = 0;
          for (unsigned int s = 0; s < NStrips; s++) {
            int StripCharge = Ampls[FirstAmplitude - NStrips + s];
            //int StripCharge = Ampls[s];
            if (StripCharge > 1024)
              StripCharge = 255;
            else if (StripCharge > 254)
              StripCharge = 254;
            mCharge1 += StripCharge;
            mCharge2 += StripCharge;
            mCharge3 += StripCharge;
            mCharge4 += StripCharge;
          }
          // Revome gains for monitoring
          mCharge2 *= PrevGain;                   // remove G2
          mCharge3 *= PrevGainTick;               // remove G1
          mCharge4 *= (PrevGain * PrevGainTick);  // remove G1 and G2
          mCharge1 *= (PrevGain / NewGain);       // remove old G2 and apply new G2

          // Getting raw charge with strip gain.

          double clustergain = 0;
          double cleanclustergain = 0;
          // SiStripClusterInfo.stripCharges() <==> SiStripCluster.amplitudes()
          for (size_t chidx = 0; chidx < siStripClusterInfo_.stripCharges().size(); ++chidx) {
            if (siStripClusterInfo_.stripCharges()[chidx] <= 0) {
              continue;
            }  // nonzero amplitude
            clustergain +=
                gainHandle->getStripGain(siStripClusterInfo_.firstStrip() + chidx, gainHandle->getRange(DetId));
            if (stripQuality->IsStripBad(stripQuality->getRange(DetId), siStripClusterInfo_.firstStrip() + chidx)) {
              continue;
            }
            cleanclustergain +=
                gainHandle->getStripGain(siStripClusterInfo_.firstStrip() + chidx, gainHandle->getRange(DetId));
          }
          clustergain /= double(siStripClusterInfo_.stripCharges().size());  // calculating average gain inside cluster
          cleanclustergain /= double(siStripClusterInfo_.stripCharges().size());

          //std::cout<< "effect of bad strips: "<<clustergain - cleanclustergain << std::endl;

          //if(clustergain!=(PrevGain*PrevGainTick)) std::cout<<"clustergain: "<<clustergain<<" | apvgain:"<< (PrevGain*PrevGainTick) <<std::endl;

          h_APVGain->Fill(PrevGain * PrevGainTick);
          h_SumStripGain->Fill(clustergain);
          h_diffAPVfromClusterGain->Fill(fabs(clustergain - (PrevGain * PrevGainTick)));

          if (FirstStrip == 0)
            Overlapping = true;
          if (FirstStrip == 128)
            Overlapping = true;
          if (FirstStrip == 256)
            Overlapping = true;
          if (FirstStrip == 384)
            Overlapping = true;
          if (FirstStrip == 512)
            Overlapping = true;
          if (FirstStrip == 640)
            Overlapping = true;

          if (FirstStrip <= 127 && FirstStrip + Ampls.size() > 127)
            Overlapping = true;
          if (FirstStrip <= 255 && FirstStrip + Ampls.size() > 255)
            Overlapping = true;
          if (FirstStrip <= 383 && FirstStrip + Ampls.size() > 383)
            Overlapping = true;
          if (FirstStrip <= 511 && FirstStrip + Ampls.size() > 511)
            Overlapping = true;
          if (FirstStrip <= 639 && FirstStrip + Ampls.size() > 639)
            Overlapping = true;

          if (FirstStrip + Ampls.size() == 127)
            Overlapping = true;
          if (FirstStrip + Ampls.size() == 255)
            Overlapping = true;
          if (FirstStrip + Ampls.size() == 383)
            Overlapping = true;
          if (FirstStrip + Ampls.size() == 511)
            Overlapping = true;
          if (FirstStrip + Ampls.size() == 639)
            Overlapping = true;
          if (FirstStrip + Ampls.size() == 767)
            Overlapping = true;

          bool farfromedge = IsFarFromBorder(&trajState, DetId, &iSetup);

          // starts here the cluster selection

          if (farfromedge == false)
            continue;
          if (Overlapping == true)
            continue;
          if (Saturation && !AllowSaturation)
            continue;
          if (NStrips > MaxNrStrips)
            continue;

          h_diffAPVfromClusterGainNoOverlap->Fill(fabs(clustergain - (PrevGain * PrevGainTick)));

          double ChargeOverPath = (double)Charge / Path;

          p_ClusterChargeVsLS->Fill(ls, ChargeOverPath);
          h2_ClusterChargeVsLS->Fill(ls, ChargeOverPath);

          h_CChargeOverPath->Fill(ChargeOverPath);
          h_CChargeOverPathNewG2->Fill(double(mCharge1) / Path);
          h_CChargeOverPathNoG2->Fill(double(mCharge2) / Path);
          h_CChargeOverPathNoG1->Fill(double(mCharge3) / Path);
          h_CChargeOverPathNoG1G2->Fill(double(mCharge4) / Path);

          h_CChargeOverPath_layer[myLayer]->Fill(ChargeOverPath);
          h_CChargeOverPathNewG2_layer[myLayer]->Fill(double(mCharge1) / Path);
          h_CChargeOverPathNoG2_layer[myLayer]->Fill(double(mCharge2) / Path);
          h_CChargeOverPathNoG1_layer[myLayer]->Fill(double(mCharge3) / Path);
          h_CChargeOverPathNoG1G2_layer[myLayer]->Fill(double(mCharge4) / Path);

          nUsedStripClus++;

        } else if (PixelCluster) {
          const auto& Ampls = PixelCluster->pixelADC();
          int FirstRow = PixelCluster->minPixelRow();
          int FirstCol = PixelCluster->minPixelCol();
          FirstStrip = ((FirstRow / 80) << 3 | (FirstCol / 52)) * 128;  //Hack to save the APVId
          NStrips = 0;
          Saturation = false;
          Overlapping = false;

          for (unsigned int a = 0; a < Ampls.size(); a++) {
            Charge += Ampls[a];
            if (Ampls[a] >= 254)
              Saturation = true;
          }  // loop on amplitudes

        }  // if it's pixel
      }    // h-index
    }      // loop on TM
  }        // loop on tracks

  h_nClusters->Fill(nStripClus);
  h_nUsedClusters->Fill(nUsedStripClus);
}

// ------------ method called once each job just before starting event loop  ------------
void SiStripGainsValidator::beginJob() {
  TH1F::SetDefaultSumw2(kTRUE);

  h_StoN = fs->make<TH1F>("clusterStoN", "cluster Raw S/N;cluster raw S/N;# clusters", 100, 0., 100.);
  h_StoNCorr =
      fs->make<TH1F>("clusterStoNCorr", "cluster S/N corrected for path lenght;cluster S/N;# clusters", 100, 0., 100.);

  h_Noise = fs->make<TH1F>("clusterNoise", "cluster noise;cluster noise;# clusters", 100, 0., 100.);
  h_chargeFromClusterInfo = fs->make<TH1F>("clusterCharge", "cluster charge;cluster charge;# clusters", 100, 0., 1000.);
  h_chargeFromStrips =
      fs->make<TH1F>("sumOfStripsCharge", "cluster charge from strips;cluster charge;# clusters", 100, 0., 1000.);

  h_diffChargeMethod =
      fs->make<TH1F>("diffChargeMethod",
                     "#Delta(cluster charge,#Sigma strips charge);#Delta(Q_{c},#Sigma_{s} Q_{s});clusters",
                     100,
                     0.,
                     100.);

  h_clusterwidth = fs->make<TH1F>("clusetWidth", "cluster width;cluster width;# clusters", 50, -0.5, 49.5);
  h_clusterposition =
      fs->make<TH1F>("clusterPosition", "cluster position; cluster position; # clusters", 100, 0., 1000.);

  h_APVGain = fs->make<TH1F>("APVGain", "average cluster gain (from APV);cluster gain;clusters", 100, 0., 2.);
  h_SumStripGain =
      fs->make<TH1F>("sumStripGain", "average sum of strip gains;average sums of strips gains;clusters", 100, 0., 2.);

  h_diffAPVfromClusterGain = fs->make<TH1F>(
      "diffClusterFromAPVgain", "#Delta(cluster gain,APV gain);#Delta(G_{cluster},G_{APV});clusters", 100, 0., 1.);

  h_diffAPVfromClusterGainNoOverlap =
      fs->make<TH1F>("diffClusterFromAPVgainNoOverlap",
                     "#Delta(cluster gain,APV gain);#Delta(G_{cluster},G_{APV});clusters",
                     100,
                     0.,
                     1.);

  h_CChargeOverPath = fs->make<TH1F>(
      "clusterChargeOverPath", "cluster charge over path;cluster charge / path;# clusters", 100, 0., 1000.);

  p_ClusterChargeVsLS = fs->make<TProfile>(
      "p_ClusterChargeVsLS", "cluster charge / path per LS;LS; #LT cluster charge / path #GT", 2500, 0, 2500);
  h2_ClusterChargeVsLS = fs->make<TH2F>(
      "h2_ClusterChargeVsLS", "cluster charge / path per LS;LS; cluster charge / path", 2500, 0, 2500, 100, 0., 1000.);

  h_CChargeOverPathNewG2 = fs->make<TH1F>("clusterChargeOverPathNewG2",
                                          "cluster charge over path (new G2);cluster charge / path;# clusters",
                                          100,
                                          0.,
                                          1000.);
  h_CChargeOverPathNoG2 = fs->make<TH1F>("clusterChargeOverPathNoG2",
                                         "cluster charge over path (G2 removed);cluster charge / path;# clusters",
                                         100,
                                         0.,
                                         1000.);
  h_CChargeOverPathNoG1 = fs->make<TH1F>("clusterChargeOverPathNoG1",
                                         "cluster charge over path (G1 removed);cluster charge / path;# clusters",
                                         100,
                                         0.,
                                         1000.);
  h_CChargeOverPathNoG1G2 =
      fs->make<TH1F>("clusterChargeOverPathNoG1G2",
                     "cluster charge over path (all gains removed);cluster charge / path;# clusters",
                     100,
                     0.,
                     1000.);

  TFileDirectory StoN = fs->mkdir("SignalToNoise");
  TFileDirectory ClusterCharge = fs->mkdir("ClusterCharge");

  for (int fooInt = SiStripHelper::TIBL1; fooInt != SiStripHelper::NUM_OF_TYPES; fooInt++) {
    SiStripHelper::layer layer = static_cast<SiStripHelper::layer>(fooInt);
    std::string s_layer = getStringFromEnum(layer);
    std::string append = myreplace(s_layer, " ", "_");

    h_StoN_layer[layer] =
        StoN.make<TH1F>(Form("clusterStoN_%s", append.c_str()),
                        Form("cluster Raw S/N (%s); %s cluster raw S/N;# clusters", s_layer.c_str(), s_layer.c_str()),
                        100,
                        0.,
                        100.);
    h_StoNCorr_layer[layer] = StoN.make<TH1F>(
        Form("clusterStoNCorr_%s", append.c_str()),
        Form("cluster S/N (%s) corrected for path lenght;%s cluster S/N;# clusters", s_layer.c_str(), s_layer.c_str()),
        100,
        0.,
        100.);
    h_CChargeOverPath_layer[layer] = ClusterCharge.make<TH1F>(
        Form("clusterChargeOverPath_%s", append.c_str()),
        Form("cluster charge over path (%s); %s cluster charge / path;# clusters", s_layer.c_str(), s_layer.c_str()),
        100,
        0.,
        1000.);
    h_CChargeOverPathNewG2_layer[layer] =
        ClusterCharge.make<TH1F>(Form("clusterChargeOverPathNewG2_%s", append.c_str()),
                                 Form("cluster charge over path (new G2) (%s); %s cluster charge / path;# clusters",
                                      s_layer.c_str(),
                                      s_layer.c_str()),
                                 100,
                                 0.,
                                 1000.);
    h_CChargeOverPathNoG2_layer[layer] =
        ClusterCharge.make<TH1F>(Form("clusterChargeOverPathNoG2_%s", append.c_str()),
                                 Form("cluster charge over path (G2 removed)(%s); %s cluster charge / path;# clusters",
                                      s_layer.c_str(),
                                      s_layer.c_str()),
                                 100,
                                 0.,
                                 1000.);
    h_CChargeOverPathNoG1_layer[layer] =
        ClusterCharge.make<TH1F>(Form("clusterChargeOverPathNoG1_%s", append.c_str()),
                                 Form("cluster charge over path (G1 removed)(%s); %s cluster charge / path;# clusters",
                                      s_layer.c_str(),
                                      s_layer.c_str()),
                                 100,
                                 0.,
                                 1000.);
    h_CChargeOverPathNoG1G2_layer[layer] = ClusterCharge.make<TH1F>(
        Form("clusterChargeOverPathNoG1G2_%s", append.c_str()),
        Form("cluster charge over path (all gains removed)(%s); %s cluster charge / path;# clusters",
             s_layer.c_str(),
             s_layer.c_str()),
        100,
        0.,
        1000.);
  }

  TFileDirectory tfd = fs->mkdir("trackControl");

  vTrackHistos_.push_back(tfd.make<TH1F>("h_tracketa", "Track #eta;#eta_{Track};Number of Tracks", 90, -3., 3.));

  vTrackHistos_.push_back(
      tfd.make<TH1F>("h_tracketaerr", "Track #eta error;err(#eta); Number of Tracks", 100, 0., 0.1));

  vTrackHistos_.push_back(tfd.make<TH1F>("h_trackphi", "Track #phi;#phi_{Track};Number of Tracks", 90, -3.15, 3.15));

  vTrackHistos_.push_back(
      tfd.make<TH1F>("h_trackphierr", "Track #phi error;err(#phi) [rad]; Number of Tracks", 100, 0., 0.1));

  vTrackHistos_.push_back(tfd.make<TH1F>(
      "h_trackNumberOfValidHits", "Track # of valid hits;# of valid hits _{Track};Number of Tracks", 40, 0., 40.));

  vTrackHistos_.push_back(tfd.make<TH1F>(
      "h_trackNumberOfLostHits", "Track # of lost hits;# of lost hits _{Track};Number of Tracks", 10, 0., 10.));

  vTrackHistos_.push_back(
      tfd.make<TH1F>("h_curvature", "Curvature #kappa;#kappa_{Track} [GeV^{-1}];Number of Tracks", 100, -0.5, 0.5));

  vTrackHistos_.push_back(
      tfd.make<TH1F>("h_curvatureerr", "Track Curvature error;err(#kappa) [GeV^{-1}]", 100, 0., 0.1));

  vTrackHistos_.push_back(tfd.make<TH1F>("h_chi2", "#chi^{2};#chi^{2}_{Track};Number of Tracks", 500, -0.01, 500.));

  vTrackHistos_.push_back(
      tfd.make<TH1F>("h_normchi2", "#chi^{2}/ndof;#chi^{2}/ndof;Number of Tracks", 100, -0.01, 10.));

  vTrackHistos_.push_back(tfd.make<TH1F>("h_pt", "p_{T}^{track};p_{T}^{track} [GeV];Number of Tracks", 250, 0., 250));

  vTrackHistos_.push_back(tfd.make<TH1F>(
      "h_ptResolution", "#delta_{p_{T}}/p_{T}^{track};#delta_{p_{T}}/p_{T}^{track};Number of Tracks", 100, 0., 0.5));

  vTrackHistos_.push_back(tfd.make<TH1F>("h_trackndof", "track # of DOF;track n. DOF;Number of Tracks", 30, 0.5, 30.5));

  vTrackHistos_.push_back(tfd.make<TH1F>("h_trackcharge", "track charge;track charge;Number of Tracks", 3, -1.5, 1.5));

  vTrackHistos_.push_back(
      tfd.make<TH1F>("h_tracktheta", "track #theta;track #theta angle;Number of Tracks", 100, 0., -3.15));

  vTrackHistos_.push_back(
      tfd.make<TH1F>("h_trackthetaerr", "track #theta error;err(#theta) [rad]; Number of Tracks", 100, 0., 0.1));

  vTrackHistos_.push_back(
      tfd.make<TH1F>("h_trackdxy", "Transverse Impact Parameter;d_{xy} [cm]; Number of Tracks", 200, -1., 1.));

  vTrackHistos_.push_back(
      tfd.make<TH1F>("h_trackdz", "Longitudinal Impact Parameter;d_{z} [cm]; Number of Tracks", 200, -30., 30.));

  vTrackHistos_.push_back(tfd.make<TH1F>(
      "h_trackvx", "track x-coordinate reference point;track v_{x} [cm]; Number of Tracks", 200, -1., 1.));

  vTrackHistos_.push_back(tfd.make<TH1F>(
      "h_trackvy", "track y-coordinate reference point;track v_{y} [cm]; Number of Tracks", 200, -1., 1.));

  vTrackHistos_.push_back(tfd.make<TH1F>(
      "h_trackvz", "track z-coordinate reference point;track v_{z} [cm]; Number of Tracks", 200, -10., 10.));

  vTrackHistos_.push_back(tfd.make<TH1F>(
      "h_trackdxyerr", "Transverse Impact Parameter error;err(d_{xy}) [cm]; Number of Tracks", 100, 0., 0.2));

  vTrackHistos_.push_back(tfd.make<TH1F>(
      "h_trackdzerr", "Longitudinal Impact Parameter error;err(d_{z}) [cm]; Number of Tracks", 100, 0., 2.));

  vTrackHistos_.push_back(
      tfd.make<TH1F>("h_trackalgo", "tracking algorithm;tracking algorithm; Number of Tracks", 15, -0.5, 14.5));

  h_nClusters = fs->make<TH1F>("nStripClusters", "n. of Strip clusters;n. Strip clusters; events", 5000, 0., 5000.);
  h_nUsedClusters = fs->make<TH1F>(
      "nUsedStripClusters", "n. of selected Strip clusters;n. selected Strip clusters; events", 2000, 0., 2000.);
}

// ------------ method to get the detector thickness ------------
//****************************************************************/
double SiStripGainsValidator::thickness(DetId id)
//****************************************************************/
{
  std::map<DetId, double>::iterator th = m_thicknessMap.find(id);
  if (th != m_thicknessMap.end())
    return (*th).second;
  else {
    double detThickness = 1.;
    //compute thickness normalization
    const GeomDetUnit* it = m_tracker->idToDetUnit(DetId(id));
    bool isPixel = dynamic_cast<const PixelGeomDetUnit*>(it) != nullptr;
    bool isStrip = dynamic_cast<const StripGeomDetUnit*>(it) != nullptr;
    if (!isPixel && !isStrip) {
      //FIXME throw exception
      edm::LogWarning("DeDxHitsProducer") << "\t\t this detID doesn't seem to belong to the Tracker";
      detThickness = 1.;
    } else {
      detThickness = it->surface().bounds().thickness();
    }

    m_thicknessMap[id] = detThickness;  //computed value
    return detThickness;
  }
}

//****************************************************************/
bool SiStripGainsValidator::IsFarFromBorder(TrajectoryStateOnSurface* trajState,
                                            const uint32_t detid,
                                            const edm::EventSetup* iSetup)
//****************************************************************/
{
  edm::ESHandle<TrackerGeometry> tkGeom;
  iSetup->get<TrackerDigiGeometryRecord>().get(tkGeom);

  LocalPoint HitLocalPos = trajState->localPosition();
  LocalError HitLocalError = trajState->localError().positionError();

  const GeomDetUnit* it = tkGeom->idToDetUnit(DetId(detid));
  if (dynamic_cast<const StripGeomDetUnit*>(it) == nullptr && dynamic_cast<const PixelGeomDetUnit*>(it) == nullptr) {
    std::cout << "this detID doesn't seem to belong to the Tracker" << std::endl;
    return false;
  }

  const BoundPlane plane = it->surface();
  const TrapezoidalPlaneBounds* trapezoidalBounds(dynamic_cast<const TrapezoidalPlaneBounds*>(&(plane.bounds())));
  const RectangularPlaneBounds* rectangularBounds(dynamic_cast<const RectangularPlaneBounds*>(&(plane.bounds())));

  double DistFromBorder = 1.0;
  double HalfLength = it->surface().bounds().length() / 2.0;

  if (trapezoidalBounds) {
    std::array<const float, 4> const& parameters = (*trapezoidalBounds).parameters();
    HalfLength = parameters[3];
  } else if (rectangularBounds) {
    HalfLength = it->surface().bounds().length() / 2.0;
  } else {
    return false;
  }

  if (fabs(HitLocalPos.y()) + HitLocalError.yy() >= (HalfLength - DistFromBorder))
    return false;

  return true;
}

//****************************************************************/
int SiStripGainsValidator::checkLayer(unsigned int iidd, const TrackerTopology* tTopo)
//****************************************************************/
{
  StripSubdetector strip = StripSubdetector(iidd);
  unsigned int subid = strip.subdetId();
  if (subid == StripSubdetector::TIB) {
    return tTopo->tibLayer(iidd);
  }
  if (subid == StripSubdetector::TOB) {
    return tTopo->tobLayer(iidd) + 4;
  }
  if (subid == StripSubdetector::TID) {
    return tTopo->tidWheel(iidd) + 10;
  }
  if (subid == StripSubdetector::TEC) {
    return tTopo->tecWheel(iidd) + 13;
  }
  return 0;
}

//****************************************************************/
template <class OBJECT_TYPE>
int SiStripGainsValidator::GetIndex(const std::vector<OBJECT_TYPE*>& vec, const TString& name)
//****************************************************************/
{
  int result = 0;
  for (typename std::vector<OBJECT_TYPE*>::const_iterator iter = vec.begin(), iterEnd = vec.end(); iter != iterEnd;
       ++iter, ++result) {
    if (*iter && (*iter)->GetName() == name)
      return result;
  }
  edm::LogError("SiStripGainsValidator") << "@SUB=SiStripGainsValidator::GetIndex"
                                         << " could not find " << name;
  return -1;
}

// -------------- method to get the topology from the detID ------------------------------
//****************************************************************/
std::string SiStripGainsValidator::getStringFromEnum(SiStripHelper::layer e)
//****************************************************************/
{
  switch (e) {
    case SiStripHelper::TIBL1:
      return "TIB L1";
    case SiStripHelper::TIBL2:
      return "TIB L2";
    case SiStripHelper::TIBL3:
      return "TIB L3";
    case SiStripHelper::TIBL4:
      return "TIB L4";
    case SiStripHelper::TOBL1:
      return "TOB L1";
    case SiStripHelper::TOBL2:
      return "TOB L2";
    case SiStripHelper::TOBL3:
      return "TOB L3";
    case SiStripHelper::TOBL4:
      return "TOB L4";
    case SiStripHelper::TOBL5:
      return "TOB L5";
    case SiStripHelper::TOBL6:
      return "TOB L6";
    case SiStripHelper::TIDD1:
      return "TID Disk 1";
    case SiStripHelper::TIDD2:
      return "TID Disk 2";
    case SiStripHelper::TIDD3:
      return "TID Disk 3";
    case SiStripHelper::TECD1:
      return "TEC Disk 1";
    case SiStripHelper::TECD2:
      return "TEC Disk 2";
    case SiStripHelper::TECD3:
      return "TEC Disk 3";
    case SiStripHelper::TECD4:
      return "TEC Disk 4";
    case SiStripHelper::TECD5:
      return "TEC Disk 5";
    case SiStripHelper::TECD6:
      return "TEC Disk 6";
    case SiStripHelper::TECD7:
      return "TEC Disk 7";
    case SiStripHelper::TECD8:
      return "TEC Disk 8";
    case SiStripHelper::TECD9:
      return "TEC Disk 9";
    default:
      edm::LogWarning("LogicError") << "Unknown partition: " << e;
      return "";
  }
}

//****************************************************************/
std::string SiStripGainsValidator::myreplace(const std::string& s,
                                             const std::string& toReplace,
                                             const std::string& replaceWith)
//****************************************************************/
{
  std::string replacement = s;
  return (replacement.replace(replacement.find(toReplace), toReplace.length(), replaceWith));
}

//****************************************************************/
void SiStripGainsValidator::getPeakOfLandau(TH1* InputHisto, double* FitResults, double LowRange, double HighRange)
//****************************************************************/
{
  FitResults[0] = -0.5;  //MPV
  FitResults[1] = 0;     //MPV error
  FitResults[2] = -0.5;  //Width
  FitResults[3] = 0;     //Width error
  FitResults[4] = -0.5;  //Fit Chi2/NDF
  FitResults[5] = 0;     //Normalization

  if (InputHisto->GetEntries() < 10.)
    return;

  // perform fit with standard landau
  TF1 MyLandau("MyLandau", "landau", LowRange, HighRange);
  MyLandau.SetParameter(1, 300);

  if (0 == InputHisto->Fit(&MyLandau, "0QR WW")) {
    if (InputHisto->GetFunction(MyLandau.GetName())) {  // Take care that it is later on drawn:
      InputHisto->GetFunction(MyLandau.GetName())->ResetBit(TF1::kNotDraw);
    }
  }

  //InputHisto->Fit(&MyLandau,"0QR WW");

  // MPV is parameter 1 (0=constant, 1=MPV, 2=Sigma)
  FitResults[0] = MyLandau.GetParameter(1);                     //MPV
  FitResults[1] = MyLandau.GetParError(1);                      //MPV error
  FitResults[2] = MyLandau.GetParameter(2);                     //Width
  FitResults[3] = MyLandau.GetParError(2);                      //Width error
  FitResults[4] = MyLandau.GetChisquare() / MyLandau.GetNDF();  //Fit Chi2/NDF
  FitResults[5] = MyLandau.GetParameter(0);
}

// ------------ method called once each job just after ending the event loop  ------------
void SiStripGainsValidator::endJob() {
  for (int fooInt = SiStripHelper::TIBL1; fooInt != SiStripHelper::NUM_OF_TYPES; fooInt++) {
    SiStripHelper::layer layer = static_cast<SiStripHelper::layer>(fooInt);
    std::string s_layer = getStringFromEnum(layer);

    double FitResults_Old[6], FitResults_New[6];
    getPeakOfLandau(h_CChargeOverPath_layer[layer], FitResults_Old, 250., 400.);
    getPeakOfLandau(h_CChargeOverPathNewG2_layer[layer], FitResults_New, 250., 400.);

    std::cout << "*********************************************************" << std::endl;
    std::cout << "Summary for: " << s_layer << std::endl;
    std::cout << "Old MPV: " << FitResults_Old[0] << " +/-" << FitResults_Old[1] << std::endl;
    std::cout << "New MPV: " << FitResults_New[0] << " +/-" << FitResults_New[1] << std::endl;
  }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void SiStripGainsValidator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SiStripGainsValidator);
