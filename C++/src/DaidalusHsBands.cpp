/*
 * Copyright (c) 2015-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */

#include "Real.h"
#include "DaidalusHsBands.h"
#include "Velocity.h"
#include "Position.h"
#include "Detection3D.h"
#include "DaidalusRealBands.h"
#include "IntervalSet.h"
#include "Tuple5.h"
#include "Kinematics.h"
#include "ProjectedKinematics.h"
#include <vector>

namespace larcfm {

DaidalusHsBands::DaidalusHsBands() {}

DaidalusHsBands::DaidalusHsBands(const DaidalusHsBands& b) : DaidalusRealBands(b) {}

bool DaidalusHsBands::do_recovery(const DaidalusParameters& parameters) const {
  return parameters.isEnabledRecoveryHorizontalSpeedBands();
}

real DaidalusHsBands::get_step(const DaidalusParameters& parameters) const {
  return parameters.getHorizontalSpeedStep();
}

real DaidalusHsBands::get_min(const DaidalusParameters& parameters) const {
  return parameters.getMinHorizontalSpeed();
}

real DaidalusHsBands::get_max(const DaidalusParameters& parameters) const {
  return parameters.getMaxHorizontalSpeed();
}

real DaidalusHsBands::get_min_rel(const DaidalusParameters& parameters) const {
  return parameters.getBelowRelativeHorizontalSpeed();
}

real DaidalusHsBands::get_max_rel(const DaidalusParameters& parameters) const{
  return parameters.getAboveRelativeHorizontalSpeed();
}

bool DaidalusHsBands::instantaneous_bands(const DaidalusParameters& parameters) const {
  return parameters.getHorizontalAcceleration() == real(0.0);
}

real DaidalusHsBands::own_val(const TrafficState& ownship) const {
  return ownship.velocityXYZ().gs();
}

real DaidalusHsBands::time_step(const DaidalusParameters& parameters, const TrafficState& ownship) const {
  // This function is never called when horizontal acceleration is zero
  return parameters.getHorizontalAcceleration() == real(0.0) ? real(0.0) : get_step(parameters)/parameters.getHorizontalAcceleration();
}

std::pair<Vect3, Vect3> DaidalusHsBands::trajectory(const DaidalusParameters& parameters, const TrafficState& ownship, real time, bool dir, int target_step, bool instantaneous) const {
  std::pair<Position,Velocity> posvel;
  if (time == real(0.0) && target_step == real(0.0)) {
    return std::pair<Vect3, Vect3>(ownship.get_s(),ownship.get_v());
  } else if (instantaneous) {
    real gs = ownship.velocityXYZ().gs()+(dir?1:-1)*target_step*get_step(parameters);
    posvel = std::pair<Position, Velocity>(ownship.positionXYZ(),ownship.velocityXYZ().mkGs(gs));
  } else {
    posvel = ProjectedKinematics::gsAccel(ownship.positionXYZ(),ownship.velocityXYZ(),time,
        (dir?1:-1)*parameters.getHorizontalAcceleration());
  }
  return std::pair<Vect3, Vect3>(ownship.pos_to_s(posvel.first),ownship.vel_to_v(posvel.first,posvel.second));
}


real DaidalusHsBands::max_delta_resolution(const DaidalusParameters& parameters) const {
  return parameters.getPersistencePreferredHorizontalSpeedResolution();
}

}

