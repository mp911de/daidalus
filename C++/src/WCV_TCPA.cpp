/*
 * Copyright (c) 2015-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */

#include "Real.h"
#include "WCV_TCPA.h"
#include "CDCylinder.h"
#include "WCV_TCOA.h"
#include "Vect3.h"
#include "Velocity.h"
#include "Horizontal.h"
#include "WCVTable.h"
#include "LossData.h"
#include "Util.h"
#include "format.h"
#include "string_util.h"

namespace larcfm {

/** Constructor that uses the default TCAS tables. */
WCV_TCPA::WCV_TCPA() : WCV_tvar(new WCV_TCOA()) {}

WCV_TCPA::WCV_TCPA(const WCV_TCPA& wcv) : WCV_tvar(wcv) {}

/**
 * @return one static WCV_TCPA
 */
const WCV_TCPA& WCV_TCPA::A_WCV_TCPA() {
  static WCV_TCPA dwc;
  return dwc;
}

real WCV_TCPA::horizontal_tvar(const Vect2& s, const Vect2& v) const {
  // Time variable is Modified Tau
  real TCPA = -1;
  real sdotv = s.dot(v);
  if (sdotv < 0)
    return (Util::sq(getDTHR())-s.sqv())/sdotv;
  return TCPA;
}

LossData WCV_TCPA::horizontal_WCV_interval(real T, const Vect2& s, const Vect2& v) const {
  real time_in = T;
  real time_out = 0;
  real sqs = s.sqv();
  real sqv = v.sqv();
  real sdotv = s.dot(v);
  real sqD = Util::sq(getDTHR());
  if (Util::almost_equals(sqv,0) && sqs <= sqD) { // [CAM] Changed from == to almost_equals to mitigate numerical problems
    time_in = 0;
    time_out = T;
    return LossData(time_in,time_out);
  }
  if (Util::almost_equals(sqv,0)) // [CAM] Changed from == to almost_equals to mitigate numerical problems
    return LossData(time_in,time_out);
  if (sqs <= sqD) {
    time_in = 0;
    time_out = Util::min(T,Horizontal::Theta_D(s,v,1,getDTHR()));
    return LossData(time_in,time_out);
  }
  if (sdotv > 0)
    return LossData(time_in,time_out);
  real tcpa = Horizontal::tcpa(s,v);
  if (v.ScalAdd(tcpa, s).norm() > getDTHR())
    return LossData(time_in,time_out);
  real Delta = Horizontal::Delta(s,v,getDTHR());
  if (Delta < 0 && tcpa - getTTHR() > T)
    return LossData(time_in,time_out);
  if (Delta < 0) {
    time_in = Util::max(real(0.0),tcpa-getTTHR());
    time_out = Util::min(T,tcpa);
    return LossData(time_in,time_out);
  }
  const auto times = Horizontal::Theta_D_times(s,v,getDTHR());
  real tmin = Util::min(times.first,tcpa-getTTHR());
  if (tmin > T)
    return LossData(time_in,time_out);
  time_in = Util::max(real(0.0),tmin);
  time_out = Util::min(T,times.second);
  return LossData(time_in,time_out);
}

Detection3D* WCV_TCPA::make() const {
  return new WCV_TCPA();
}

/**
 * Returns a deep copy of this WCV_TCPA object, including any results that have been calculated.
 */
Detection3D* WCV_TCPA::copy() const {
  return new WCV_TCPA(*this);
}

std::string WCV_TCPA::getSimpleClassName() const {
  return "WCV_TCPA";
}

bool WCV_TCPA::contains(const Detection3D& cd) const {
  if (larcfm::equals(getCanonicalClassName(), cd.getCanonicalClassName())) {
    return containsTable((WCV_tvar&)cd);
  }
  return false;
}

void WCV_TCPA::hazard_zone_far_end(std::vector<Position>& haz,
    const Position& po, const Velocity& v, const Vect3& pu, real T) const {
  Position npo = po.linear(v,getTTHR()+T);
  Velocity vu = Velocity::make(pu);
  haz.push_back(npo.linear(vu,-getDTHR()));
  real b = v.vect3().norm2D()*getTTHR();
  if (Util::almost_greater(getDTHR(),b)) {
      // Far end has the form of a cap
      real a = Util::sqrt_safe(Util::sq(getDTHR())-Util::sq(b));
      real alpha = Util::acos_safe(b/getDTHR());
      Vect3 vD = pu.ScalAdd(-a,v.vect3().Hat().Scal(b));
      CDCylinder::circular_arc(haz,po.linear(v,T),
              Velocity::make(vD),2*alpha,true);
  }
  haz.push_back(npo.linear(vu,getDTHR()));
}


}

