/*
 * Copyright (c) 2015-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */

#include "Real.h"
#include "DaidalusAltBands.h"
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

DaidalusAltBands::DaidalusAltBands() {}

DaidalusAltBands::DaidalusAltBands(const DaidalusAltBands& b) : DaidalusRealBands(b) {}

bool DaidalusAltBands::do_recovery(const DaidalusParameters& parameters) const {
  return parameters.isEnabledRecoveryAltitudeBands();
}

real DaidalusAltBands::get_step(const DaidalusParameters& parameters) const {
  return parameters.getAltitudeStep();
}

real DaidalusAltBands::get_min(const DaidalusParameters& parameters) const {
  return parameters.getMinAltitude();
}

real DaidalusAltBands::get_max(const DaidalusParameters& parameters) const {
  return parameters.getMaxAltitude();
}

real DaidalusAltBands::get_min_rel(const DaidalusParameters& parameters) const {
  return parameters.getBelowRelativeAltitude();
}

real DaidalusAltBands::get_max_rel(const DaidalusParameters& parameters) const {
  return parameters.getAboveRelativeAltitude();
}

void DaidalusAltBands::set_special_configuration(const DaidalusParameters& parameters, const SpecialBandFlags& special_flags) {
  if (special_flags.dta_status > 0) {
    set_min_max_rel(0,-1);
  }
}

bool DaidalusAltBands::instantaneous_bands(const DaidalusParameters& parameters) const {
  return parameters.getVerticalRate() == 0 ||
      parameters.getVerticalAcceleration() == 0;
}

real DaidalusAltBands::own_val(const TrafficState& ownship) const {
  return ownship.positionXYZ().alt();
}

real DaidalusAltBands::time_step(const DaidalusParameters& parameters, const TrafficState& ownship) const {
  return 1;
}

namespace {

class LevelOutTrajectory : public DaidalusIntegerBands {
  real target_alt_;
  real end_time_;

public:
  const Tuple5<real,real,real,real,real> times;

  LevelOutTrajectory(const DaidalusParameters& parameters, const TrafficState& ownship, real target_alt, real time_step)
    : target_alt_(target_alt),
      times(Kinematics::vsLevelOutTimes(ownship.positionXYZ().alt(),ownship.velocityXYZ().vs(),
          parameters.getVerticalRate(),target_alt,parameters.getVerticalAcceleration(),-parameters.getVerticalAcceleration(),true)) {
    end_time_ = (times.first < 0 ? real(-1) : times.third)+time_step;
  }

  std::pair<Vect3,Vect3> trajectory(const DaidalusParameters&, const TrafficState& ownship, real time, bool, int, bool) const override {
    const Position& position = ownship.positionXYZ();
    const Velocity& velocity = ownship.velocityXYZ();
    if (time <= end_time_) {
      const auto vertical = Kinematics::vsLevelOutCalc(position.alt(),velocity.vs(),target_alt_,
          times.fourth,times.fifth,times.first,times.second,times.third,time);
      return std::make_pair(position.vect3().linear(velocity.vect3(),time).mkZ(vertical.first),
          velocity.mkVs(vertical.second).vect3());
    }
    return std::make_pair(position.linear(velocity,time).mkZ(target_alt_).vect3(),velocity.mkVs(0).vect3());
  }
};

}

std::pair<Vect3, Vect3> DaidalusAltBands::trajectory(const DaidalusParameters& parameters, const TrafficState& ownship, real time, bool dir, int target_step, bool instantaneous) const {
  real target_alt = get_min_val_()+target_step*get_step(parameters);
  if (instantaneous) {
    Position position = ownship.positionXYZ().mkZ(target_alt);
    return std::make_pair(ownship.pos_to_s(position),ownship.vel_to_v(position,ownship.velocityXYZ().mkVs(0)));
  }
  LevelOutTrajectory level_out(parameters,ownship,target_alt,time_step(parameters,ownship));
  return level_out.trajectory(parameters,ownship,time,dir,target_step,instantaneous);
}

// In PVS: alt_bands@conflict_free_traj_step
bool DaidalusAltBands::conflict_free_traj_step(const Detection3D& conflict_det, const Detection3D& recovery_det, real B, real T,
    const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic, int target_step, bool instantaneous) const {
  bool trajdir = true;
  if (instantaneous) {
    return no_CD_future_traj(conflict_det,recovery_det,B,T,trajdir,real(0.0),parameters,ownship,traffic,target_step,instantaneous);
  } else {
    real tstep = time_step(parameters,ownship);
    real target_alt = get_min_val_()+target_step*get_step(parameters);
    const LevelOutTrajectory level_out(parameters,ownship,target_alt,tstep);
    real tsqj1 = level_out.times.first;
    real tsqj2 = level_out.times.second;
    real tsqj3 = level_out.times.third+tstep;
    for (int i=0; i<=std::floor(tsqj1/tstep);++i) {
      real tsi = i*tstep;
      if ((B <= tsi && level_out.LOS_at(conflict_det,trajdir,tsi,parameters,ownship,traffic,target_step,instantaneous)) ||
          (recovery_det.isValid() && 0 <= tsi && tsi <= B &&
              level_out.LOS_at(recovery_det,trajdir,tsi,parameters,ownship,traffic,target_step,instantaneous))) {
        return false;
      }
    }
    real tsk1 = Util::max(tsqj1,real(0.0));
    real time_horizon_1 = parameters.isEnabledBandsAddTimeToManeuver() ? T : T+tsk1;
    if ((tsqj2 >= B &&
        level_out.CD_future_traj(conflict_det,B,std::min(time_horizon_1,tsqj2),trajdir,tsk1,parameters,ownship,traffic,target_step,instantaneous)) ||
        (recovery_det.isValid() && tsqj2 >= 0 &&
            level_out.CD_future_traj(recovery_det,0,Util::min(B,tsqj2),trajdir,tsk1,parameters,ownship,traffic,target_step,instantaneous))) {
      return false;
    }
    for (int i=(int)std::ceil(tsqj2/tstep); i<=std::floor(tsqj3/tstep);++i) {
      real tsi = i*tstep;
      if ((B <= tsi && level_out.LOS_at(conflict_det,trajdir,tsi,parameters,ownship,traffic,target_step,instantaneous)) ||
          (recovery_det.isValid() && 0 <= tsi && tsi <= B &&
              level_out.LOS_at(recovery_det,trajdir,tsi,parameters,ownship,traffic,target_step,instantaneous))) {
        return false;
      }
    }
    real tsk3 = Util::max(tsqj3,real(0.0));
    real time_horizon_3 = parameters.isEnabledBandsAddTimeToManeuver() ? T : T+tsk3;
    return level_out.no_CD_future_traj(conflict_det,recovery_det,B,time_horizon_3,trajdir,tsk3,parameters,ownship,traffic,target_step,instantaneous);
  }
}

// In PVS: alt_bands@alt_bands_generic
void DaidalusAltBands::alt_bands_generic(std::vector<Integerval>& l,
    const Detection3D& conflict_det, const Detection3D& recovery_det, real B, real T,
    int maxup, const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic, bool instantaneous) const {
  int d = -1; // Set to the first index with no conflict
  for (int k = 0; k <= maxup; ++k) {
    if (d >=0 && conflict_free_traj_step(conflict_det,recovery_det,B,T,parameters,ownship,traffic,k,instantaneous)) {
      continue;
    } else if (d >=0) {
      l.push_back(Integerval(d,k-1));
      d = -1;
    } else if (conflict_free_traj_step(conflict_det,recovery_det,B,T,parameters,ownship,traffic,k,instantaneous)) {
      d = k;
    }
  }
  if (d >= 0) {
    l.push_back(Integerval(d,maxup));
  }
}

void DaidalusAltBands::none_bands(IntervalSet& noneset, const Detection3D& conflict_det, const Detection3D& recovery_det,
    int epsh, int epsv, real B, real T, const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic) const {
  int maxup = (int)std::floor((get_max_val_()-get_min_val_())/get_step(parameters))+1;
  std::vector<Integerval> altint;
  alt_bands_generic(altint,conflict_det,recovery_det,B,T,maxup,parameters,ownship,traffic,instantaneous_bands(parameters));
  toIntervalSet(noneset,altint,get_step(parameters),get_min_val_());
}

bool DaidalusAltBands::any_red(const Detection3D& conflict_det, const Detection3D& recovery_det,
    int epsh, int epsv, real B, real T, const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic) const {
  return first_band_alt_generic(conflict_det,recovery_det,B,T,parameters,ownship,traffic,true,false,instantaneous_bands(parameters)) >= 0 ||
      first_band_alt_generic(conflict_det,recovery_det,B,T,parameters,ownship,traffic,false,false,instantaneous_bands(parameters)) >= 0;
}

bool DaidalusAltBands::all_red(const Detection3D& conflict_det, const Detection3D& recovery_det,
    int epsh, int epsv, real B, real T, const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic) const {
  return first_band_alt_generic(conflict_det,recovery_det,B,T,parameters,ownship,traffic,true,true,instantaneous_bands(parameters)) < 0 &&
      first_band_alt_generic(conflict_det,recovery_det,B,T,parameters,ownship,traffic,false,true,instantaneous_bands(parameters)) < 0;
}

int DaidalusAltBands::first_nat(int mini, int maxi, bool dir, const Detection3D& conflict_det, const Detection3D& recovery_det,
    real B, real T,  const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic, bool green, bool instantaneous) const {
  while (mini <= maxi) {
    if (dir && green == conflict_free_traj_step(conflict_det,recovery_det,B,T,parameters,ownship,traffic,mini,instantaneous)) {
      return mini;
    } else if (dir) {
      ++mini;
    } else {
      if (green == conflict_free_traj_step(conflict_det,recovery_det,B,T,parameters,ownship,traffic,maxi,instantaneous)) {
        return maxi;
      } else if (maxi == 0) {
        return -1;
      } else {
        --maxi;
      }
    }
  }
  return -1;
}

int DaidalusAltBands::first_band_alt_generic(const Detection3D& conflict_det, const Detection3D& recovery_det,
    real B, real T,
    const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic, bool dir, bool green, bool instantaneous) const {
  int upper = (int)(dir ? std::floor((get_max_val_()-get_min_val_())/get_step(parameters))+1 :
      std::floor((ownship.positionXYZ().alt()-get_min_val_())/get_step(parameters)));
  int lower = dir ? (int)(std::ceil(ownship.positionXYZ().alt()-get_min_val_())/get_step(parameters)) : 0;
  if (ownship.positionXYZ().alt() < get_min_val_() || ownship.positionXYZ().alt() > get_max_val_()) {
    return -1;
  } else {
    return first_nat(lower,upper,dir,conflict_det,recovery_det,B,T,parameters,ownship,traffic,green,instantaneous);
  }
}

real DaidalusAltBands::max_delta_resolution(const DaidalusParameters& parameters) const {
  return parameters.getPersistencePreferredAltitudeResolution();
}

}

