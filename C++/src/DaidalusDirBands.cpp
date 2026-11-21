/*
 * Copyright (c) 2015-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */

#include "Real.h"
#include "DaidalusDirBands.h"
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

DaidalusDirBands::DaidalusDirBands() : DaidalusRealBands(2*Pi), below_as_thr_(0) {}

DaidalusDirBands::DaidalusDirBands(const DaidalusDirBands& b) : DaidalusRealBands(b), below_as_thr_(0) {}

bool DaidalusDirBands::do_recovery(const DaidalusParameters& parameters) const {
  return parameters.isEnabledRecoveryHorizontalDirectionBands();
}

real DaidalusDirBands::get_step(const DaidalusParameters& parameters) const {
  return parameters.getHorizontalDirectionStep();
}

real DaidalusDirBands::get_min(const DaidalusParameters& parameters) const {
  return real(0.0);
}

real DaidalusDirBands::get_max(const DaidalusParameters& parameters) const {
  return get_mod();
}

real DaidalusDirBands::get_min_rel(const DaidalusParameters& parameters) const {
  return parameters.getLeftHorizontalDirection();
}

real DaidalusDirBands::get_max_rel(const DaidalusParameters& parameters) const {
  return parameters.getRightHorizontalDirection();
}

bool DaidalusDirBands::saturate_corrective_bands(const DaidalusParameters& parameters, const SpecialBandFlags& special_flags) const {
  return special_flags.dta_status > 0 && parameters.getDTALogic() < 0;
}

void DaidalusDirBands::set_special_configuration(const DaidalusParameters& parameters, const SpecialBandFlags& special_flags) {
  if (special_flags.below_as_thr) {
    below_as_thr_ = parameters.getHorizontalDirBandsBelowAirspeedThr();
  } else {
    below_as_thr_ = 0;
  }
}

bool DaidalusDirBands::instantaneous_bands(const DaidalusParameters& parameters) const {
  return (parameters.getTurnRate() == real(0.0) && parameters.getBankAngle() == real(0.0)) || below_as_thr_ > 0;
}

real DaidalusDirBands::own_val(const TrafficState& ownship) const {
  return ownship.velocityXYZ().compassAngle();
}

real DaidalusDirBands::time_step(const DaidalusParameters& parameters, const TrafficState& ownship) const {
	real gso = ownship_gs(parameters,ownship);
  real omega = parameters.getTurnRate() == real(0.0) ? Kinematics::turnRate(gso,parameters.getBankAngle()) : parameters.getTurnRate();
        if (omega == real(0.0)) {
			return real(0.0);
		} else {
      return get_step(parameters)/omega;
    }
}

real DaidalusDirBands::ownship_gs(const DaidalusParameters& parameters, const TrafficState& ownship) const {
  real gso = below_as_thr_ == 0 ? ownship.velocityXYZ().gs() :
    std::max(parameters.getAirspeedThreshold(),ownship.velocityXYZ().gs());
  return  std::max(parameters.getHorizontalSpeedStep(),gso);
}

Velocity DaidalusDirBands::ownship_vel(const DaidalusParameters& parameters, const TrafficState& ownship) const {
 	return ownship.velocityXYZ().mkGs(ownship_gs(parameters,ownship));
}

std::pair<Vect3, Vect3> DaidalusDirBands::trajectory(const DaidalusParameters& parameters, const TrafficState& ownship, real time, bool dir, int target_step, bool instantaneous) const {
  std::pair<Position,Velocity> posvel;
  Velocity ownship_velocityXYZ = ownship_vel(parameters,ownship);
  if (time == real(0.0) && target_step == real(0.0)) {
    return std::pair<Vect3, Vect3>(ownship.get_s(),ownship_velocityXYZ.vect3());
  } else if (instantaneous) {
    real trk = ownship_velocityXYZ.compassAngle()+(dir?1:-1)*target_step*get_step(parameters);
    posvel = std::pair<Position, Velocity>(ownship.positionXYZ(),ownship_velocityXYZ.mkTrk(trk));
  } else {
    real gso = ownship_velocityXYZ.gs();
    real bank = parameters.getTurnRate() == real(0.0) ? parameters.getBankAngle() : std::abs(Kinematics::bankAngle(gso,parameters.getTurnRate()));
    real R = Kinematics::turnRadius(gso,bank);
    posvel = ProjectedKinematics::turn(ownship.positionXYZ(),ownship_velocityXYZ,time,R,dir);
  }
  return std::pair<Vect3, Vect3>(ownship.pos_to_s(posvel.first),ownship.vel_to_v(posvel.first,posvel.second));
}

real DaidalusDirBands::max_delta_resolution(const DaidalusParameters& parameters) const {
  return parameters.getPersistencePreferredHorizontalDirectionResolution();
}

std::string DaidalusDirBands::rawString() const {
  std::string s(DaidalusRealBands::rawString());
	s += "below_as_thr_ = "+Fmi(below_as_thr_)+"\n";
	return s;
}

}

