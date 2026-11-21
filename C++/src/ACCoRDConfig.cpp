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
 *
 * ACCoRDConfig 
 *
 */

#include "Real.h"
#include "ACCoRDConfig.h"
#include "Units.h"

namespace larcfm {

const real ACCoRDConfig::minHorizExitSpeedLoS_default = Units::from("kn",real(100.0));
const real ACCoRDConfig::minVertExitSpeedLoS_default = Units::from("fpm",real(1499.99999999));
const real ACCoRDConfig::maxGs_default = Units::from("kn",real(700.0));
const real ACCoRDConfig::minGs_default = Units::from("kn",real(150.0));       // must be greater than 0
const real ACCoRDConfig::maxVs_default = Units::from("fpm",real(5000.0));       // must be greater than 0
real ACCoRDConfig::gsSearchLosDiscard = Units::from("nm",real(1.5));
real ACCoRDConfig::vsDiscretization_default = Units::from("fpm",real(10.0));

const real ACCoRDConfig::NMAC_D = Units::from("ft",500); // Defined in RTCA SC-147
const real ACCoRDConfig::NMAC_H = Units::from("ft",100); // Defined in RTCA SC-147

// internal units right now
void ACCoRDConfig::setGsSearchLosDiscard(real val) {
  gsSearchLosDiscard = val;
}


}


