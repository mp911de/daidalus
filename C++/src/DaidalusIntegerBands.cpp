/*
 * Copyright (c) 2015-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */

#include "Real.h"
#include "DaidalusIntegerBands.h"
#include "CriteriaCore.h"
#include "TCASTable.h"
#include "Util.h"
#include <vector>
#include <string>

namespace larcfm {

TrafficState DaidalusIntegerBands::projected_ownship(const DaidalusParameters& parameters, const TrafficState& ownship, real tsk, bool trajdir, int target_step, bool instantaneous) const {
  std::pair<Vect3,Vect3> sovot = trajectory(parameters,ownship,tsk,trajdir,target_step,instantaneous);
  Vect3 sot = sovot.first;
  Vect3 vot = sovot.second;
  Vect3 sat = tsk == real(0.0) ? sot : vot.ScalAdd(-tsk,sot);
  TrafficState own = ownship;
  own.setPosition(Position::make(sat));
  own.setAirVelocity(Velocity::make(vot));
  return own;
}

/**
 * In PVS: int_bands@CD_future_traj
 */
bool DaidalusIntegerBands::CD_future_traj(const Detection3D& det, real B, real T, bool trajdir, real tsk,
    const DaidalusParameters& parameters,  const TrafficState& ownship, const TrafficState& traffic, int target_step, bool instantaneous) const {
  T = Util::min(parameters.getLookaheadTime(),T);
  if (tsk > T || B > T) return false;
  TrafficState own = projected_ownship(parameters,ownship,tsk,trajdir,target_step,instantaneous);
  return det.conflictWithTrafficState(own,traffic,Util::max(B,tsk),T);
}

/**
 * In PVS:
 * NOT CD_future_traj(CD,B,T,traj,kts,si,vi)) AND (NOT (useLOS2 AND CD_future_traj(CD2,traj,kts,si,vi)
 */
bool DaidalusIntegerBands::no_CD_future_traj(const Detection3D& conflict_det, const Detection3D& recovery_det, real B, real T,  bool trajdir, real tsk,
    const DaidalusParameters& parameters,  const TrafficState& ownship, const TrafficState& traffic, int target_step, bool instantaneous) const {
  T = Util::min(parameters.getLookaheadTime(),T);
  if (tsk > T || B > T) {
    return !(recovery_det.isValid() && CD_future_traj(recovery_det,0,B,trajdir,tsk,parameters,ownship,traffic,target_step,instantaneous));
  }
  TrafficState own = projected_ownship(parameters,ownship,tsk,trajdir,target_step,instantaneous);
  if (conflict_det.conflictWithTrafficState(own,traffic,Util::max(B,tsk),T)) return false;
  real recovery_end = Util::min(parameters.getLookaheadTime(),B);
  return !recovery_det.isValid() || tsk > recovery_end || 0 > recovery_end ||
      !recovery_det.conflictWithTrafficState(own,traffic,Util::max(real(0),tsk),recovery_end);
}

/**
 * In PVS:
 * LET kts = k*ts,
 * sot=traj(kts)`1,
 * vot=traj(kts)`2 IN
 * LOS_at(sot-kts*vot,vot,si,vi,kts)
 */
bool DaidalusIntegerBands::LOS_at(const Detection3D& det, bool trajdir, real tsk,
    const DaidalusParameters& parameters,  const TrafficState& ownship, const TrafficState& traffic, int target_step, bool instantaneous) const {
  if (tsk >= parameters.getLookaheadTime()) {
      return false;
  }
  std::pair<Vect3,Vect3> sovot = trajectory(parameters,ownship,tsk,trajdir,target_step,instantaneous);
  Vect3 sot = sovot.first;
  Vect3 vot = sovot.second;
  Vect3 sat = vot.ScalAdd(-tsk,sot);
  TrafficState own = ownship;
  own.setPosition(Position::make(sat));
  own.setAirVelocity(Velocity::make(vot));
  return det.violationAtWithTrafficState(own,traffic,tsk);
}

bool DaidalusIntegerBands::LOS_at_coast(const Detection3D& det, bool trajdir, real tsk, real tcoast,
      const DaidalusParameters& parameters, const TrafficState& ownship, 
      const TrafficState& traffic, int target_step, bool instantaneous) const {
  if (tsk >= parameters.getLookaheadTime()) {
		return false;
	}
	std::pair<Vect3,Vect3> sovot = trajectory(parameters,ownship,tsk,trajdir,target_step,instantaneous);
	Vect3 sot = sovot.first;
	Vect3 vot = sovot.second;
	Vect3 sat = vot.ScalAdd(-tsk,sot);
	TrafficState own = ownship;
	own.setPosition(Position::make(sat));
	own.setAirVelocity(Velocity::make(vot));
	return det.conflictWithTrafficState(own,traffic,tsk,tsk+tcoast); 
}

// In PVS: int_bands@first_los_step
int DaidalusIntegerBands::kinematic_first_los_step(const Detection3D& det, real tstep, bool trajdir,
    int min, int max,const DaidalusParameters& parameters,  const TrafficState& ownship, const TrafficState& traffic) const {
  for (int k=min; k<=max; ++k) {
    if (LOS_at_coast(det,trajdir,k*tstep,tstep,parameters,ownship,traffic,0,false)) {
      return k;
    }
  }
  return -1;
}

// In PVS: kinematic_bands@first_los_search_index
int DaidalusIntegerBands::kinematic_first_los_search_index(const Detection3D& conflict_det, const Detection3D& recovery_det, real tstep,
    real B,  bool trajdir, int max,
    const DaidalusParameters& parameters,  const TrafficState& ownship, const TrafficState& traffic) const {
  int FirstLosK = (int)std::ceil(B/tstep); // first k such that k*ts>=B
  int FirstLosN = Util::min((int)std::floor(parameters.getLookaheadTime()/tstep),max); // last k<=MaxN such that k*ts<=T
  int FirstLosK2 = 0;
  int FirstLosN2 = Util::min((int)std::floor(B/tstep),max);
  int FirstLosInit = recovery_det.isValid() ? kinematic_first_los_step(recovery_det,tstep,trajdir,FirstLosK2,FirstLosN2,parameters,ownship,traffic) : -1;
  int FirstLos = kinematic_first_los_step(conflict_det,tstep,trajdir,FirstLosK,FirstLosN,parameters,ownship,traffic);
  int LosInitIndex = FirstLosInit < 0 ? max+1 : FirstLosInit;
  int LosIndex = FirstLos < 0 ? max+1 : FirstLos;
  return Util::min(LosInitIndex,LosIndex);
}

// In PVS: kinematic_bands@bands_search_index
// epsh == epsv == 0, if traffic is not the repulsive aircraft
int DaidalusIntegerBands::kinematic_bands_search_index(const Detection3D& conflict_det, const Detection3D& recovery_det, real tstep,
    real B, bool trajdir, int max,const DaidalusParameters& parameters,  const TrafficState& ownship, const TrafficState& traffic,
    int epsh, int epsv) const {
  bool usehcrit = epsh != 0;
  bool usevcrit = epsv != 0;
  int FirstLos = kinematic_first_los_search_index(conflict_det,recovery_det,tstep,B,trajdir,max,parameters,ownship,traffic);
  int FirstNonHRep = !usehcrit || FirstLos == 0 ? FirstLos :
      kinematic_first_nonrepulsive_step(tstep,trajdir,FirstLos-1,parameters,ownship,traffic,epsh);
  int FirstProbHcrit = FirstNonHRep < 0 ? max+1 : FirstNonHRep;
  int FirstProbHL = Util::min(FirstLos,FirstProbHcrit);
  int FirstNonVRep = !usevcrit || FirstProbHL == 0 ? FirstProbHL :
      kinematic_first_nonvert_repul_step(tstep,trajdir,FirstProbHL-1,parameters,ownship,traffic,epsv);
  int FirstProbVcrit = FirstNonVRep < 0 ? max+1 : FirstNonVRep;
  return Util::min(FirstProbHL,FirstProbVcrit);
}

// In PVS: int_bands@traj_conflict_only_band, int_bands@nat_bands, and int_bands@nat_bands_rec

bool DaidalusIntegerBands::kinematic_traj_conflict_only_bands(std::vector<Integerval>& l,
    const Detection3D& conflict_det, const Detection3D& recovery_det, real tstep, real B, real T,
    bool trajdir, int max,const DaidalusParameters& parameters,  const TrafficState& ownship, const TrafficState& traffic, const BandTest* test) const {
  return collect_runs(l,max,trajdir,false,[&](int k) {
    real tsk = tstep*k;
    real horizon = parameters.isEnabledBandsAddTimeToManeuver() ? T : T+tsk;
    return no_CD_future_traj(conflict_det,recovery_det,B,horizon,trajdir,tsk,parameters,ownship,traffic,0,false);
  },test);
}

// In PVS: kinematic_bands@kinematic_bands
bool DaidalusIntegerBands::kinematic_bands(std::vector<Integerval>& l, const Detection3D& conflict_det, const Detection3D& recovery_det, real tstep,
    real B, real T,
    bool trajdir, int max,const DaidalusParameters& parameters,  const TrafficState& ownship, const TrafficState& traffic,
    int epsh, int epsv, const BandTest* test) const {
  l.clear();
  int bsi = kinematic_bands_search_index(conflict_det,recovery_det,tstep,B,trajdir,max,parameters,ownship,traffic,epsh,epsv);
  return bsi != 0 && kinematic_traj_conflict_only_bands(l,conflict_det,recovery_det,tstep,B,T,trajdir,bsi-1,parameters,ownship,traffic,test);
}

// In PVS: kinematic_bands_exist@first_green
int DaidalusIntegerBands::first_kinematic_green(const Detection3D& conflict_det, const Detection3D& recovery_det, real tstep,
    real B, real T,
    bool trajdir, int max,const DaidalusParameters& parameters,  const TrafficState& ownship, const TrafficState& traffic,
    int epsh, int epsv) const {
  bool usehcrit = epsh != 0;
  bool usevcrit = epsv != 0;
  for (int k=0; k <= max; ++k) {
    real tsk = tstep*k;
    real time_horizon = parameters.isEnabledBandsAddTimeToManeuver() ? T : T+tsk;
    if ((B <= tsk && LOS_at(conflict_det,trajdir,tsk,parameters,ownship,traffic,0,false)) ||
        (recovery_det.isValid() && 0 <= tsk && tsk <= B &&
            LOS_at(recovery_det,trajdir,tsk,parameters,ownship,traffic,0,false)) ||
            (usehcrit && !kinematic_repulsive_at(tstep,trajdir,k,parameters,ownship,traffic,epsh)) ||
            (usevcrit && !kinematic_vert_repul_at(tstep,trajdir,k,parameters,ownship,traffic,epsv))) {
      return -1;
    } else if (no_CD_future_traj(conflict_det,recovery_det,B,
      time_horizon,trajdir,tsk,parameters,ownship,traffic,0,false)) {
      return k;
    }
  }
  return -1;
}

// In PVS: int_bands@repulsive_at
bool DaidalusIntegerBands::kinematic_repulsive_at(real tstep, bool trajdir, int k,const DaidalusParameters& parameters,  const TrafficState& ownship, const TrafficState& traffic, int epsh) const {
  // k >= 0
  if (k==0) {
    return true;
  }
  auto previous = trajectory(parameters,ownship,(k-1)*tstep,trajdir,0,false);
  auto current = trajectory(parameters,ownship,k*tstep,trajdir,0,false);
  Vect2 vop = current.first.Sub(previous.first).Scal(1/tstep).vect2();
  Vect2 vi = traffic.get_v().vect2();
  Vect2 si = traffic.get_s().vect2();
  if (k == 1 && !CriteriaCore::horizontal_new_repulsive_criterion(previous.first.vect2().Sub(si),previous.second.vect2(),vi,vop,epsh)) return false;
  Vect3 next = trajectory(parameters,ownship,(k+1)*tstep,trajdir,0,false).first;
  Vect2 vok = next.Sub(current.first).Scal(1/tstep).vect2();
  Vect2 st = current.first.vect2().Sub(vi.ScalAdd(k*tstep,si));
  Vect2 vot = current.second.vect2();
  return CriteriaCore::horizontal_new_repulsive_criterion(st,vop,vi,vot,epsh) &&
      CriteriaCore::horizontal_new_repulsive_criterion(st,vot,vi,vok,epsh) &&
      CriteriaCore::horizontal_new_repulsive_criterion(st,vop,vi,vok,epsh);
}

// In PVS: int_bands@first_nonrepulsive_step
int DaidalusIntegerBands::kinematic_first_nonrepulsive_step(real tstep, bool trajdir, int max,const DaidalusParameters& parameters,  const TrafficState& ownship, const TrafficState& traffic, int epsh) const {
  for (int k=0; k <= max; ++k) {
    if (!kinematic_repulsive_at(tstep,trajdir,k,parameters,ownship,traffic,epsh)) {
      return k;
    }
  }
  return -1;
}

// In PVS: int_bands@vert_repul_at
bool DaidalusIntegerBands::kinematic_vert_repul_at(real tstep, bool trajdir, int k,const DaidalusParameters& parameters,  const TrafficState& ownship, const TrafficState& traffic, int epsv) const {
  // traffic is valid and k >= 0
  if (k==0) {
    return true;
  }
  auto previous = trajectory(parameters,ownship,(k-1)*tstep,trajdir,0,false);
  auto current = trajectory(parameters,ownship,k*tstep,trajdir,0,false);
  Vect3 vop = current.first.Sub(previous.first).Scal(1/tstep);
  Vect3 vi = traffic.get_v();
  Vect3 si = traffic.get_s();
  if (k == 1 && !CriteriaCore::vertical_new_repulsive_criterion(previous.first.Sub(si),previous.second,vi,vop,epsv)) return false;
  Vect3 next = trajectory(parameters,ownship,(k+1)*tstep,trajdir,0,false).first;
  Vect3 vok = next.Sub(current.first).Scal(1/tstep);
  Vect3 st = current.first.Sub(vi.ScalAdd(k*tstep,si));
  Vect3 vot = current.second;
  return CriteriaCore::vertical_new_repulsive_criterion(st,vop,vi,vot,epsv) &&
      CriteriaCore::vertical_new_repulsive_criterion(st,vot,vi,vok,epsv) &&
      CriteriaCore::vertical_new_repulsive_criterion(st,vop,vi,vok,epsv);
}

// In PVS: int_bands@first_nonvert_repul_step
int DaidalusIntegerBands::kinematic_first_nonvert_repul_step(real tstep, bool trajdir, int max,const DaidalusParameters& parameters,  const TrafficState& ownship, const TrafficState& traffic, int epsv) const {
  for (int k=0; k <= max; ++k) {
    if (!kinematic_vert_repul_at(tstep,trajdir,k,parameters,ownship,traffic,epsv)) {
      return k;
    }
  }
  return -1;
}

// In PVS: first_conflict_step(CD,B,T,traj,0,ts,si,vi,MaxN) >= 0
bool DaidalusIntegerBands::kinematic_any_conflict_step(const Detection3D& det, real tstep, real B, real T, bool trajdir, int max,
    const DaidalusParameters& parameters,  const TrafficState& ownship, const TrafficState& traffic) const {
  for (int k=0; k <= max; ++k) {
    real tsk = tstep*k;
    real time_horizon = parameters.isEnabledBandsAddTimeToManeuver() ? T : T+tsk;
    if (CD_future_traj(det,B,time_horizon,trajdir,tsk,parameters,ownship,traffic,0,false)) {
      return true;
    }
  }
  return false;
}

// In PVS: kinematic_bands_exist@red_band_exist
bool DaidalusIntegerBands::kinematic_red_band_exist(const Detection3D& conflict_det, const Detection3D& recovery_det, real tstep,
    real B, real T,
    bool trajdir, int max,const DaidalusParameters& parameters,  const TrafficState& ownship, const TrafficState& traffic,
    int epsh, int epsv) const {
  bool usehcrit = epsh != 0;
  bool usevcrit = epsv != 0;
  return (usehcrit && kinematic_first_nonrepulsive_step(tstep,trajdir,max,parameters,ownship,traffic,epsh) >= 0) ||
      (usevcrit && kinematic_first_nonvert_repul_step(tstep,trajdir,max,parameters,ownship,traffic,epsv) >= 0) ||
      kinematic_any_conflict_step(conflict_det,tstep,B,T,trajdir,max,parameters,ownship,traffic) ||
      (recovery_det.isValid() && kinematic_any_conflict_step(recovery_det,tstep,0,B,trajdir,max,parameters,ownship,traffic));
}

int DaidalusIntegerBands::first_instantaneous_green(const Detection3D& conflict_det, const Detection3D& recovery_det,
    real B, real T,
    bool trajdir, int max,const DaidalusParameters& parameters,  const TrafficState& ownship, const TrafficState& traffic,
    int epsh, int epsv) const {
  for (int k = 0; k <= max; ++k) {
    if (no_instantaneous_conflict(conflict_det,recovery_det,B,T,trajdir,parameters,ownship,traffic,epsh,epsv,k)) {
      return k;
    }
  }
  return -1;
}

//In PVS: inst_bands@conflict_free_track_step, inst_bands@conflict_free_gs_step, inst_bands@conflict_free_vs_step
bool DaidalusIntegerBands::no_instantaneous_conflict(const Detection3D& conflict_det, const Detection3D& recovery_det, real B, real T,
    bool trajdir,const DaidalusParameters& parameters,  const TrafficState& ownship, const TrafficState& traffic,
    int epsh, int epsv, int target_step) const {
  if (epsh == 0 && epsv == 0) {
    return no_CD_future_traj(conflict_det,recovery_det,B,T,trajdir,real(0.0),parameters,ownship,traffic,target_step,true);
  }
  bool usehcrit = epsh != 0;
  bool usevcrit = epsv != 0;
  std::pair<Vect3,Vect3> nsovo = trajectory(parameters,ownship,0,trajdir,target_step,true);
  Vect3 so = ownship.get_s();
  Vect3 vo = ownship.get_v();
  Vect3 si = traffic.get_s();
  Vect3 vi = traffic.get_v();
  Vect3 nvo = nsovo.second;
  Vect3 s = so.Sub(si);
  return
      (!usehcrit || CriteriaCore::horizontal_new_repulsive_criterion(s.vect2(),vo.vect2(),vi.vect2(),nvo.vect2(),epsh)) &&
      (!usevcrit || CriteriaCore::vertical_new_repulsive_criterion(s,vo,vi,nvo,epsv)) &&
      no_CD_future_traj(conflict_det,recovery_det,B,T,trajdir,real(0.0),parameters,ownship,traffic,target_step,true);
}

//In PVS: int_bands@nat_bands, int_bands@nat_bands_rec
bool DaidalusIntegerBands::instantaneous_bands(std::vector<Integerval>& l,
    const Detection3D& conflict_det, const Detection3D& recovery_det, real B, real T,
    bool trajdir, int max,const DaidalusParameters& parameters,  const TrafficState& ownship, const TrafficState& traffic,
    int epsh, int epsv, const BandTest* test) const {
  return collect_runs(l,max,trajdir,false,[&](int k) {
    return no_instantaneous_conflict(conflict_det,recovery_det,B,T,trajdir,parameters,ownship,traffic,epsh,epsv,k);
  },test);
}

bool DaidalusIntegerBands::instantaneous_red_band_exist(const Detection3D& conflict_det, const Detection3D& recovery_det,
    real B, real T,
    bool trajdir, int max,const DaidalusParameters& parameters,  const TrafficState& ownship, const TrafficState& traffic,
    int epsh, int epsv) const {
  for (int k = 0; k <= max; ++k) {
    if (!no_instantaneous_conflict(conflict_det,recovery_det,B,T,trajdir,parameters,ownship,traffic,epsh,epsv,k)) {
      return true;
    }
  }
  return false;
}

void DaidalusIntegerBands::append_intband(std::vector<Integerval>& l, std::vector<Integerval>& r) {
  // Append in place
  int last = l.size()-1;
  if (!l.empty() && !r.empty() && r[0].lb-l[last].ub <= 1) {
    l[last].ub = r[0].ub;
    r.erase(r.begin()); // side effect
  }
  l.insert(l.end(), r.begin(), r.end());
}

void DaidalusIntegerBands::neg(std::vector<Integerval>& l) {
  // Negate, flip, and reverse in place
  if (l.empty()) return;
  int mid = (l.size()-1)/2;
  bool odd = l.size() % 2 != 0;
  for (int i=0; i <= mid; ++i) {
    if (i == mid && odd) {
      int x = l[i].lb;
      l[i].lb = -l[i].ub;
      l[i].ub = -x;
    } else {
      int i2 = l.size()-i-1;
      int x = l[i].lb;
      l[i].lb = -l[i2].ub;
      l[i2].ub = -x;
      x = l[i].ub;
      l[i].ub = -l[i2].lb;
      l[i2].lb = -x;
    }
  }
}

bool DaidalusIntegerBands::all_kinematic_red(const Detection3D& conflict_det, const Detection3D& recovery_det, real tstep,
    real B, real T,
    int maxl, int maxr,const DaidalusParameters& parameters,  const TrafficState& ownship, const TrafficState& traffic,
    int epsh, int epsv, int dir) const {
  bool leftans = dir > 0 || first_kinematic_green(conflict_det,recovery_det,tstep,B,T,false,maxl,parameters,ownship,traffic,epsh,epsv) < 0;
  bool rightans = dir < 0 || first_kinematic_green(conflict_det,recovery_det,tstep,B,T,true,maxr,parameters,ownship,traffic,epsh,epsv) < 0;
  return leftans && rightans;
}

bool DaidalusIntegerBands::all_instantaneous_red(const Detection3D& conflict_det, const Detection3D& recovery_det,
    real B, real T,
    int maxl, int maxr,const DaidalusParameters& parameters,  const TrafficState& ownship, const TrafficState& traffic,
    int epsh, int epsv, int dir) const {
  bool leftans = dir > 0 || first_instantaneous_green(conflict_det,recovery_det,B,T,false,maxl,parameters,ownship,traffic,epsh,epsv) < 0;
  bool rightans = dir < 0 || first_instantaneous_green(conflict_det,recovery_det,B,T,true,maxr,parameters,ownship,traffic,epsh,epsv) < 0;
  return leftans && rightans;
}

bool DaidalusIntegerBands::any_kinematic_red(const Detection3D& conflict_det, const Detection3D& recovery_det, real tstep,
    real B, real T,
    int maxl, int maxr,const DaidalusParameters& parameters,  const TrafficState& ownship, const TrafficState& traffic,
    int epsh, int epsv, int dir) const {
  bool leftred = dir <= 0 && kinematic_red_band_exist(conflict_det,recovery_det,tstep,B,T,false,maxl,parameters,ownship,traffic,epsh,epsv);
  bool rightred = dir >= 0 && kinematic_red_band_exist(conflict_det,recovery_det,tstep,B,T,true,maxr,parameters,ownship,traffic,epsh,epsv);
  return leftred || rightred;
}

bool DaidalusIntegerBands::any_instantaneous_red(const Detection3D& conflict_det, const Detection3D& recovery_det,
    real B, real T,
    int maxl, int maxr,const DaidalusParameters& parameters,  const TrafficState& ownship, const TrafficState& traffic,
    int epsh, int epsv, int dir) const {
  bool leftred = dir <= 0 && instantaneous_red_band_exist(conflict_det,recovery_det,B,T,false,maxl,parameters,ownship,traffic,epsh,epsv);
  bool rightred = dir >= 0 && instantaneous_red_band_exist(conflict_det,recovery_det,B,T,true,maxr,parameters,ownship,traffic,epsh,epsv);
  return leftred || rightred;
}

// INTERFACE FUNCTIONS
void DaidalusIntegerBands::integer_bands_combine(std::vector<Integerval>& l, const Detection3D& conflict_det, const Detection3D& recovery_det, real tstep,
    real B, real T,
    int maxl, int maxr,const DaidalusParameters& parameters,  const TrafficState& ownship, const TrafficState& traffic,
    int epsh, int epsv) const {
  integer_bands_until(l,conflict_det,recovery_det,tstep,B,T,maxl,maxr,parameters,ownship,traffic,epsh,epsv,nullptr);
}

bool DaidalusIntegerBands::integer_bands_until(std::vector<Integerval>& l, const Detection3D& conflict_det, const Detection3D& recovery_det, real tstep,
    real B, real T,
    int maxl, int maxr,const DaidalusParameters& parameters,  const TrafficState& ownship, const TrafficState& traffic,
    int epsh, int epsv, const BandTest* test) const {
  auto scan = [&](std::vector<Integerval>& runs, bool dir, int max) {
    if (tstep == 0) {
      return instantaneous_bands(runs,conflict_det,recovery_det,B,T,dir,max,parameters,ownship,traffic,epsh,epsv,test);
    }
    return kinematic_bands(runs,conflict_det,recovery_det,tstep,B,T,dir,max,parameters,ownship,traffic,epsh,epsv,test);
  };
  std::vector<Integerval> right;
  if (scan(l,false,maxl) || scan(right,true,maxr)) return true;
  neg(l);
  append_intband(l,right);
  return false;
}

bool DaidalusIntegerBands::all_integer_red(const Detection3D& conflict_det, const Detection3D& recovery_det, real tstep,
    real B, real T,
    int maxl, int maxr,const DaidalusParameters& parameters,  const TrafficState& ownship, const TrafficState& traffic,
    int epsh, int epsv, int dir) const {
  return tstep == 0 ?
      all_instantaneous_red(conflict_det,recovery_det,
          B,T,maxl,maxr, parameters,ownship,traffic,
          epsh,epsv,dir)
          : all_kinematic_red(conflict_det,recovery_det,tstep,
              B,T,maxl,maxr, parameters,ownship,traffic,
              epsh,epsv,dir);
}

bool DaidalusIntegerBands::any_integer_red(const Detection3D& conflict_det, const Detection3D& recovery_det, real tstep,
    real B, real T,
    int maxl, int maxr,const DaidalusParameters& parameters,  const TrafficState& ownship, const TrafficState& traffic,
    int epsh, int epsv, int dir) const {
  return tstep == 0 ?
      any_instantaneous_red(conflict_det,recovery_det,
          B,T,maxl,maxr, parameters,ownship,traffic,
          epsh,epsv,dir)
          : any_kinematic_red(conflict_det,recovery_det,tstep,
              B,T,maxl,maxr, parameters,ownship,traffic,
              epsh,epsv,dir);
}

}
