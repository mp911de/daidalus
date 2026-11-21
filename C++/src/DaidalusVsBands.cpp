/*
 * Copyright (c) 2015-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */

#include "Real.h"
#include "DaidalusVsBands.h"
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

DaidalusVsBands::DaidalusVsBands() {}

DaidalusVsBands::DaidalusVsBands(const DaidalusVsBands& b) : DaidalusRealBands(b) {}

bool DaidalusVsBands::do_recovery(const DaidalusParameters& parameters) const {
  return parameters.isEnabledRecoveryVerticalSpeedBands();
}

real DaidalusVsBands::get_step(const DaidalusParameters& parameters) const {
  return parameters.getVerticalSpeedStep();
}

real DaidalusVsBands::get_min(const DaidalusParameters& parameters) const {
  return parameters.getMinVerticalSpeed();
}

real DaidalusVsBands::get_max(const DaidalusParameters& parameters) const {
  return parameters.getMaxVerticalSpeed();
}

real DaidalusVsBands::get_min_rel(const DaidalusParameters& parameters) const {
  return parameters.getBelowRelativeVerticalSpeed();
}

real DaidalusVsBands::get_max_rel(const DaidalusParameters& parameters) const {
  return parameters.getAboveRelativeVerticalSpeed();
}

void DaidalusVsBands::set_special_configuration(const DaidalusParameters& parameters, const SpecialBandFlags& special_flags) {
  if (special_flags.dta_status > 0) {
    set_min_max_rel(0,-1);
  }
}

bool DaidalusVsBands::instantaneous_bands(const DaidalusParameters& parameters) const {
  return parameters.getVerticalAcceleration() == real(0.0);
}

real DaidalusVsBands::own_val(const TrafficState& ownship) const {
  return ownship.velocityXYZ().vs();
}

real DaidalusVsBands::time_step(const DaidalusParameters& parameters, const TrafficState& ownship) const {
  // This function is never called when vertical acceleration is zero
  return parameters.getVerticalAcceleration() == real(0.0) ? real(0.0) : get_step(parameters)/parameters.getVerticalAcceleration();
}

std::pair<Vect3, Vect3> DaidalusVsBands::trajectory(const DaidalusParameters& parameters, const TrafficState& ownship, real time, bool dir, int target_step, bool instantaneous) const {
  std::pair<Position,Velocity> posvel;
  if (time == real(0.0) && target_step == real(0.0)) {
    return std::pair<Vect3, Vect3>(ownship.get_s(),ownship.get_v());
  } else if (instantaneous_bands(parameters)) {
    real vs = ownship.velocityXYZ().vs()+(dir?1:-1)*target_step*get_step(parameters);
    posvel = std::pair<Position, Velocity>(ownship.positionXYZ(),ownship.velocityXYZ().mkVs(vs));
  } else {
    posvel = ProjectedKinematics::vsAccel(ownship.positionXYZ(),ownship.velocityXYZ(),time,(dir?1:-1)*parameters.getVerticalAcceleration());
  }
  return std::pair<Vect3, Vect3>(ownship.pos_to_s(posvel.first),ownship.vel_to_v(posvel.first,posvel.second));
}

real DaidalusVsBands::max_delta_resolution(const DaidalusParameters& parameters) const {
  return parameters.getPersistencePreferredVerticalSpeedResolution();
}

}

