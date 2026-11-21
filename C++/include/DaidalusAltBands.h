/*
 * Copyright (c) 2015-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */
#ifndef DAIDALUSALTBANDS_H_
#define DAIDALUSALTBANDS_H_
#include "Real.h"

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

class DaidalusAltBands : public DaidalusRealBands {
  friend class DaidalusCore;
  DaidalusCore* level_out_core_ = nullptr;


public:
  DaidalusAltBands();

  DaidalusAltBands(const DaidalusAltBands& b);

  virtual void refresh(DaidalusCore& core);

  virtual bool do_recovery(const DaidalusParameters& parameters) const;

  virtual real get_step(const DaidalusParameters& parameters) const;

  virtual real get_min(const DaidalusParameters& parameters) const;

  virtual real get_max(const DaidalusParameters& parameters) const;

  virtual real get_min_rel(const DaidalusParameters& parameters) const;

  virtual real get_max_rel(const DaidalusParameters& parameters) const;

  virtual bool instantaneous_bands(const DaidalusParameters& parameters) const;
  
  virtual bool saturate_corrective_bands(const DaidalusParameters& parameters, const SpecialBandFlags& special_flags) const {
    return false;
  }

  virtual void set_special_configuration(const DaidalusParameters& parameters, const SpecialBandFlags& special_flags);

  virtual real own_val(const TrafficState& ownship) const;

  virtual real time_step(const DaidalusParameters& parameters, const TrafficState& ownship) const;

  virtual std::pair<Vect3, Vect3> trajectory(const DaidalusParameters& parameters, const TrafficState& ownship, real time, bool dir, int target_step, bool instantaneous) const;

  // In PVS: alt_bands@conflict_free_traj_step
private:
  bool conflict_free_traj_step(const Detection3D& conflict_det, const Detection3D& recovery_det, real B, real T,
      const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic, int target_step, bool instantaneous) const ;

  // In PVS: alt_bands@alt_bands_generic
  void alt_bands_generic(std::vector<Integerval>& l,
      const Detection3D& conflict_det, const Detection3D& recovery_det, real B, real T,
      int maxup, const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic, bool instantaneous) const;

  int first_nat(int mini, int maxi, bool dir, const Detection3D& conflict_det, const Detection3D& recovery_det,
      real B, real T,  const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic, bool green, bool instantaneous) const;

  int first_band_alt_generic(const Detection3D& conflict_det, const Detection3D& recovery_det, real B, real T,
      const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic, bool dir, bool green, bool instantaneous) const;

public:
  virtual void none_bands(IntervalSet& noneset, const Detection3D& conflict_det, const Detection3D& recovery_det,
      int epsh, int epsv, real B, real T, const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic) const;

  virtual bool any_red(const Detection3D& conflict_det, const Detection3D& recovery_det,
      int epsh, int epsv, real B, real T, const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic) const;

  virtual bool all_red(const Detection3D& conflict_det, const Detection3D& recovery_det,
      int epsh, int epsv, real B, real T, const DaidalusParameters& parameters, const TrafficState& ownship, const TrafficState& traffic) const;

  virtual real max_delta_resolution(const DaidalusParameters& parameters) const;

};

}

#endif
