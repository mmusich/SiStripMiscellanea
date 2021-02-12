// -*- C++ -*-
//
// Package:    SiStripMiscellanea/SiStripDigisAnalyzer
// Class:      SiStripDigisAnalyzer
//
/**\class SiStripDigisAnalyzer SiStripDigisAnalyzer.cc SiStripMiscellanea/SiStripDigisAnalyzer/plugins/SiStripDigisAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marco Musich
//         Created:  Fri, 12 Feb 2021 09:26:21 GMT
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
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/SiStripDigi/interface/SiStripDigi.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "CommonTools/TrackerMap/interface/TrackerMap.h"

//
// class declaration
//

class SiStripDigisAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit SiStripDigisAnalyzer(const edm::ParameterSet&);
  ~SiStripDigisAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<edm::DetSetVector<SiStripDigi>> edmDetSetVector_SiStripDigi_Token_;
  std::unique_ptr<TrackerMap> tmap;
};

//
// constructors and destructor
//
SiStripDigisAnalyzer::SiStripDigisAnalyzer(const edm::ParameterSet& iConfig)
    : edmDetSetVector_SiStripDigi_Token_(
          consumes<edm::DetSetVector<SiStripDigi>>(iConfig.getParameter<edm::InputTag>("src"))) {
  tmap = std::make_unique<TrackerMap>("Strip");
  tmap->setTitle("Strip digis entries");
  tmap->setPalette(1);
}

SiStripDigisAnalyzer::~SiStripDigisAnalyzer() {}

//
// member functions
//

// ------------ method called for each event  ------------
void SiStripDigisAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  std::string digiProducer = "siStripDigis";
  edm::Handle<edm::DetSetVector<SiStripDigi>> stripDigis;
  iEvent.getByToken(edmDetSetVector_SiStripDigi_Token_, stripDigis);
  edm::DetSetVector<SiStripDigi>::const_iterator DSViter = stripDigis->begin();
  for (; DSViter != stripDigis->end(); DSViter++) {
    unsigned int id = DSViter->id;
    DetId detId(id);
    tmap->fill(detId.rawId(), DSViter->size());
  }
}

// ------------ method called once each job just before starting event loop  ------------
void SiStripDigisAnalyzer::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void SiStripDigisAnalyzer::endJob() {
  tmap->save(true, 0, 0, "StripDigiMap.pdf");
  tmap->save(true, 0, 0, "StripDigiMap.png");
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void SiStripDigisAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("simSiStripDigis", "ZeroSuppressed"));
  descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SiStripDigisAnalyzer);
