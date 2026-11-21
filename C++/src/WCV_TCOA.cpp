/*
 * Copyright (c) 2012-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */
#include "Real.h"
#include "WCV_TCOA.h"
#include "Vertical.h"
#include "Util.h"

namespace larcfm {

// Vertical Not Well Clear Violation
// ZTHR and TTHR are altitude and time thresholds
bool WCV_TCOA::vertical_WCV(real ZTHR, real TCOA, real sz, real vz) const {
  return std::abs(sz) <= ZTHR ||
      (vz != 0 && sz*vz <= 0 && Vertical::time_coalt(sz,vz) <= TCOA);
}

Interval WCV_TCOA::vertical_WCV_interval(real ZTHR, real TCOA, real B, real T, real sz, real vz) const {
  real time_in = B;
  real time_out = T;
  if (Util::almost_equals(vz,0) && std::abs(sz) <= ZTHR) // [CAM] Changed from == to almost_equals to mitigate numerical problems
    return Interval(time_in,time_out);
  if (Util::almost_equals(vz,0)) { // [CAM] Changed from == to almost_equals to mitigate numerical problems
    time_in = T;
    time_out = B;
    return Interval(time_in,time_out);
  }
  real act_H = Util::max(ZTHR,std::abs(vz)*TCOA);
  real tentry = Vertical::Theta_H(sz,vz,-1,act_H);
  real texit = Vertical::Theta_H(sz,vz,1,ZTHR);
  if (T < tentry || texit < B) {
    time_in = T;
    time_out = B;
    return Interval(time_in,time_out);
  }
  time_in = Util::max(B,tentry);
  time_out = Util::min(T,texit);
  return Interval(time_in,time_out);
}

/**
 * Returns a deep copy of this WCV_TCOA object, including any results that have been calculated.
 */
WCV_Vertical* WCV_TCOA::copy() const {
  WCV_TCOA* ret = new WCV_TCOA();
  return ret;
}


}
