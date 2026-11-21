/*
 * Copyright (c) 2015-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */
#ifndef DAIDALUSINTEGERBANDS_H_
#define DAIDALUSINTEGERBANDS_H_
#include "Real.h"

#include "Velocity.h"
#include "Vect2.h"
#include "Vect3.h"
#include "Position.h"
#include "Detection3D.h"
#include "Integerval.h"
#include "IntervalSet.h"
#include "CriteriaCore.h"
#include "DaidalusParameters.h"

#include <vector>
#include <string>
#include "TrafficState.h"

namespace larcfm {


class DaidalusIntegerBands {

public:
  class BandTest {
  public:
    virtual bool accepts(bool dir, int first, int last) const = 0;
  protected:
    ~BandTest() = default;
  };

protected:
  template <typename Green>
  static bool collect_runs(std::vector<Integerval>& runs, int max, bool dir, bool keep_last_single,
      const Green& green, const BandTest* test) {
    int first = -1;
    for (int k = 0; k <= max; ++k) {
      if (green(k)) {
        if (first < 0) first = k;
        else if (test && test->accepts(dir,k-1,k)) return true;
      } else if (first >= 0) {
        runs.emplace_back(first,k-1);
        first = -1;
      }
    }
    if (first >= 0 && (keep_last_single || first != max)) runs.emplace_back(first,max);
    return false;
  }

  // This class computes NONE bands

public:
  // trajdir == false is left/down
  // target_step is used by instantaneous_bands and altitude_bands
  virtual std::pair<Vect3,Vect3> trajectory(const DaidalusParameters& parameters, const TrafficState& ownship,
      real time, bool dir, int target_step, bool instantaneous) const = 0;

  virtual ~DaidalusIntegerBands() {}

  /*
   * In PVS: int_bands@CD_future_traj
   */
  bool CD_future_traj(const Detection3D& det, real B, real T, bool trajdir, real tsk,
      const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic, int target_step, bool instantaneous) const;

  /**
   * In PVS:
   * NOT CD_future_traj(CD,B,T,traj,kts,si,vi)) AND (NOT (useLOS2 AND CD_future_traj(CD2,B2,T2,traj,kts,si,vi)
   */
  bool no_CD_future_traj(const Detection3D& conflict_det, const Detection3D& recovery_det, real B, real T,  bool trajdir, real tsk,
      const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic, int target_step, bool instantaneous) const;

  /**
   * In PVS:
   * LET kts = k*ts,
   * sot=traj(kts)`1,
   * vot=traj(kts)`2 IN
   * LOS_at(sot-kts*vot,vot,si,vi,kts)
   */
  bool LOS_at(const Detection3D& det, bool trajdir, real tsk,
      const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic, int target_step, bool instantaneous) const;

private:
  TrafficState projected_ownship(const DaidalusParameters& parameters, const TrafficState& ownship, real tsk, bool trajdir, int target_step, bool instantaneous) const;

  bool LOS_at_coast(const Detection3D& det, bool trajdir, real tsk, real tcoast,
      const DaidalusParameters& parameters, const TrafficState& ownship, 
      const TrafficState& traffic, int target_step, bool instantaneous) const;	

    // In PVS: int_bands@first_los_step
  int kinematic_first_los_step(const Detection3D& det, real tstep, bool trajdir,
      int min, int max, const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic) const;

  // In PVS: kinematic_bands@first_los_search_index
  int kinematic_first_los_search_index(const Detection3D& conflict_det, const Detection3D& recovery_det, real tstep,
      real B, bool trajdir, int max,
      const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic) const;

  // In PVS: kinematic_bands@bands_search_index
  // epsh == epsv == 0, if traffic is not the repulsive aircraft
  int kinematic_bands_search_index(const Detection3D& conflict_det, const Detection3D& recovery_det, real tstep,
      real B, bool trajdir, int max, const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic,
      int epsh, int epsv) const;

  // In PVS: int_bands@traj_conflict_only_band, int_bands@nat_bands, and int_bands@nat_bands_rec
  bool kinematic_traj_conflict_only_bands(std::vector<Integerval>& l,
      const Detection3D& conflict_det, const Detection3D& recovery_det, real tstep, real B, real T,
      bool trajdir, int max, const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic, const BandTest* test) const;

  // In PVS: kinematic_bands@kinematic_bands
  bool kinematic_bands(std::vector<Integerval>& l, const Detection3D& conflict_det, const Detection3D& recovery_det, real tstep,
      real B, real T,
      bool trajdir, int max, const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic,
      int epsh, int epsv, const BandTest* test) const;

  // In PVS: kinematic_bands_exist@first_green
  int first_kinematic_green(const Detection3D& conflict_det, const Detection3D& recovery_det, real tstep,
      real B, real T,
      bool trajdir, int max, const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic,
      int epsh, int epsv) const;

  // In PVS: int_bands@repulsive_at
  bool kinematic_repulsive_at(real tstep, bool trajdir, int k, const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic, int epsh) const;

  // In PVS: int_bands@first_nonrepulsive_step
  int kinematic_first_nonrepulsive_step(real tstep, bool trajdir, int max, const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic, int epsh) const;

  // In PVS: int_bands@vert_repul_at
  bool kinematic_vert_repul_at(real tstep, bool trajdir, int k, const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic, int epsv) const;

  // In PVS: int_bands@first_nonvert_repul_step
  int kinematic_first_nonvert_repul_step(real tstep, bool trajdir, int max, const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic, int epsv) const;

  // In PVS: first_conflict_step(CD,B,T,traj,0,ts,si,vi,MaxN) >= 0
  bool kinematic_any_conflict_step(const Detection3D& det, real tstep, real B, real T, bool trajdir, int max,
      const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic) const;

  // In PVS: kinematic_bands_exist@red_band_exist
  bool kinematic_red_band_exist(const Detection3D& conflict_det, const Detection3D& recovery_det, real tstep,
      real B, real T,
      bool trajdir, int max, const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic,
      int epsh, int epsv) const;

  int first_instantaneous_green(const Detection3D& conflict_det, const Detection3D& recovery_det,
      real B, real T,
      bool trajdir, int max, const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic,
      int epsh, int epsv) const;

  //In PVS: inst_bands@conflict_free_track_step, inst_bands@conflict_free_gs_step, inst_bands@conflict_free_vs_step
  bool no_instantaneous_conflict(const Detection3D& conflict_det, const Detection3D& recovery_det, real B, real T,
      bool trajdir, const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic,
      int epsh, int epsv, int target_step) const;

  //In PVS: int_bands@nat_bands, int_bands@nat_bands_rec
  bool instantaneous_bands(std::vector<Integerval>& l,
      const Detection3D& conflict_det, const Detection3D& recovery_det, real B, real T,
      bool trajdir, int max, const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic,
      int epsh, int epsv, const BandTest* test) const;

  bool instantaneous_red_band_exist(const Detection3D& conflict_det, const Detection3D& recovery_det,
      real B, real T,
      bool trajdir, int max, const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic,
      int epsh, int epsv) const;

  bool all_kinematic_red(const Detection3D& conflict_det, const Detection3D& recovery_det, real tstep,
      real B, real T,
      int maxl, int maxr, const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic,
      int epsh, int epsv, int dir) const;

  bool all_instantaneous_red(const Detection3D& conflict_det, const Detection3D& recovery_det,
      real B, real T,
      int maxl, int maxr, const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic,
      int epsh, int epsv, int dir) const;

  bool any_kinematic_red(const Detection3D& conflict_det, const Detection3D& recovery_det, real tstep,
      real B, real T,
      int maxl, int maxr, const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic,
      int epsh, int epsv, int dir) const;

  bool any_instantaneous_red(const Detection3D& conflict_det, const Detection3D& recovery_det,
      real B, real T,
      int maxl, int maxr, const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic,
      int epsh, int epsv, int dir) const;

public:

  static void append_intband(std::vector<Integerval>& l, std::vector<Integerval>& r);

  static void neg(std::vector<Integerval>& l);

  // INTERFACE FUNCTIONS

  void integer_bands_combine(std::vector<Integerval>& l, const Detection3D& conflict_det, const Detection3D& recovery_det, real tstep,
      real B, real T,
      int maxl, int maxr, const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic,
      int epsh, int epsv) const;

  bool integer_bands_until(std::vector<Integerval>& l, const Detection3D& conflict_det, const Detection3D& recovery_det, real tstep,
      real B, real T,
      int maxl, int maxr, const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic,
      int epsh, int epsv, const BandTest* test) const;

  bool all_integer_red(const Detection3D& conflict_det, const Detection3D& recovery_det, real tstep,
      real B, real T,
      int maxl, int maxr, const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic,
      int epsh, int epsv, int dir) const;

  bool any_integer_red(const Detection3D& conflict_det, const Detection3D& recovery_det, real tstep,
      real B, real T,
      int maxl, int maxr, const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic,
      int epsh, int epsv, int dir) const;

};

}

#endif
