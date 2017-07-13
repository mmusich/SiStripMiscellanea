#ifndef SISTRIPMISCELLANEA_SISTRIPRESOLUTIONANALYZER_HELPER_H
#define SISTRIPMISCELLANEA_SISTRIPRESOLUTIONANALYZER_HELPER_H

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"


#include <string>
#include <vector>
#include <utility>
#include <stdint.h>

namespace SiStripResol {
  
  enum partitions {TIB,
		   TIBL1,
		   TIBL2,
		   TIBL3,
		   TIBL4,
		   TOB,
		   TOBL1,
		   TOBL2,
		   TOBL3,
		   TOBL4,
		   TOBL5,
		   TOBL6,
		   TIDP,
		   TIDM,
		   TIDPD1,
		   TIDPD2,
		   TIDPD3,
		   TIDMD1,
		   TIDMD2,
		   TIDMD3,
		   TECP,
		   TECM,
		   NUM_OF_TYPES};
  
  std::pair<int,std::pair<int,int> > typeAndLayerFromDetId( const DetId& detId , const TrackerTopology* tTopo);
  std::string getStringFromEnum(partitions e);
};

#endif
