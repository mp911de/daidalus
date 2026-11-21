/*
 * Copyright (c) 2015-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */
#ifndef DAIDALUSVSBANDS_H_
#define DAIDALUSVSBANDS_H_
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

class DaidalusVsBands : public DaidalusRealBands {

public:
  DaidalusVsBands();

  DaidalusVsBands(const DaidalusVsBands& b);

  virtual bool do_recovery(const DaidalusParameters& parameters) const;

  virtual real get_step(const DaidalusParameters& parameters) const;

  virtual real get_min(const DaidalusParameters& parameters) const;

  virtual real get_max(const DaidalusParameters& parameters) const;

  virtual real get_min_rel(const DaidalusParameters& parameters) const;

  virtual real get_max_rel(const DaidalusParameters& parameters) const;

  virtual bool saturate_corrective_bands(const DaidalusParameters& parameters, const SpecialBandFlags& special_flags) const {
    return false;
  }
 
  virtual void set_special_configuration(const DaidalusParameters& parameters, const SpecialBandFlags& special_flag);

  virtual bool instantaneous_bands(const DaidalusParameters& parameters) const;

  virtual real own_val(const TrafficState& ownship) const;

  virtual real time_step(const DaidalusParameters& parameters, const TrafficState& ownship) const;

  virtual std::pair<Vect3, Vect3> trajectory(const DaidalusParameters& parameters, const TrafficState& ownship, real time, bool dir, int target_step, bool instantaneous) const;

  virtual real max_delta_resolution(const DaidalusParameters& parameters) const;

};

}

#endif
