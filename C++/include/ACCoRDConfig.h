/* ====================================================================
 * Authors:      George Hagen, Rick Butler and Jeff Maddalon,
 * Organization: NASA/Langley Research Center
 * Website:      http://shemesh.larc.nasa.gov/fm/
 *
 * Copyright (c) 2011-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 * ====================================================================
 */

#ifndef ACCORDCONFIG_H_
#define ACCORDCONFIG_H_
#include "Real.h"

namespace larcfm {

/**
 * This class captures the default parameters for ACCORD.
 *
 */
class ACCoRDConfig {

private:

public:

  static const real minHorizExitSpeedLoS_default;   //  = Units.from("kn",100.0);
  static const real minVertExitSpeedLoS_default;
  static const real maxGs_default;                  //  = Units.from("kn",700.0);
  static const real minGs_default;                  // = Units.from("kn",150.0);       // must be greater than 0
  static const real maxVs_default;
  static real gsSearchLosDiscard;     // = Units.from("nm",1.5);
  static real vsDiscretization_default; //  = Units::from("fpm",10.0);

  static const real NMAC_D; // Defined in RTCA SC-147
  static const real NMAC_H; // Defined in RTCA SC-147

  static void setLosRepulsiveCrit(bool flag);
  static void setNewCesarAlgorithm(bool flag);

  static void setGsSearchLosDiscard(real val);

}; // ACCORDCONFIG
}

#endif /* ACCORDCONFIG_H_ */
