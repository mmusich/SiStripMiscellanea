#include "SiStripMiscellanea/SiStripResolutionAnalyzer/interface/SiStripResolutionHelper.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// -------------- method to get the topology from the detID ------------------------------
std::pair<int, std::pair<int, int> > SiStripResol::typeAndLayerFromDetId(const DetId& detId,
                                                                         const TrackerTopology* tTopo) {
  int layerNumber = 0;
  int sideNumber = 0;
  unsigned int subdetId = static_cast<unsigned int>(detId.subdetId());

  if (subdetId == StripSubdetector::TIB) {
    layerNumber = tTopo->tibLayer(detId.rawId());
    sideNumber = tTopo->tibSide(detId.rawId());
  } else if (subdetId == StripSubdetector::TOB) {
    layerNumber = tTopo->tobLayer(detId.rawId());
    sideNumber = tTopo->tobSide(detId.rawId());
  } else if (subdetId == StripSubdetector::TID) {
    layerNumber = tTopo->tidWheel(detId.rawId());
    sideNumber = tTopo->tidSide(detId.rawId());
  } else if (subdetId == StripSubdetector::TEC) {
    layerNumber = tTopo->tecWheel(detId.rawId());
    sideNumber = tTopo->tecSide(detId.rawId());
  } else if (subdetId == PixelSubdetector::PixelBarrel) {
    layerNumber = tTopo->pxbLayer(detId.rawId());
  } else if (subdetId == PixelSubdetector::PixelEndcap) {
    layerNumber = tTopo->pxfDisk(detId.rawId());
    sideNumber = tTopo->pxfSide(detId.rawId());
  } else
    edm::LogWarning("LogicError") << "Unknown subdetid: " << subdetId;

  return std::make_pair(subdetId, std::make_pair(layerNumber, sideNumber));
}

// -------------- method to get the topology from the detID ------------------------------
std::string SiStripResol::getStringFromEnum(SiStripResol::partitions e) {
  switch (e) {
    case partitions::TIB:
      return "TIB";
    case partitions::TIBL1:
      return "TIB L1";
    case partitions::TIBL2:
      return "TIB L2";
    case partitions::TIBL3:
      return "TIB L3";
    case partitions::TIBL4:
      return "TIB L4";
    case partitions::TOB:
      return "TOB";
    case partitions::TOBL1:
      return "TOB L1";
    case partitions::TOBL2:
      return "TOB L2";
    case partitions::TOBL3:
      return "TOB L3";
    case partitions::TOBL4:
      return "TOB L4";
    case partitions::TOBL5:
      return "TOB L5";
    case partitions::TOBL6:
      return "TOB L6";
    case partitions::TIDP:
      return "TIDplus";
    case partitions::TIDM:
      return "TIDminus";
    case partitions::TIDPD1:
      return "TID plus Disk 1";
    case partitions::TIDPD2:
      return "TID plus Disk 2";
    case partitions::TIDPD3:
      return "TID plus Disk 3";
    case partitions::TIDMD1:
      return "TID minus Disk 1";
    case partitions::TIDMD2:
      return "TID minus Disk 2";
    case partitions::TIDMD3:
      return "TID minus Disk 3";
    case partitions::TECP:
      return "TECplus";
    case partitions::TECM:
      return "TECminus";
    default:
      edm::LogWarning("LogicError") << "Unknown partition: " << e;
      return "";
  }
}
