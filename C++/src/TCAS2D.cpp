/*
 * Copyright (c) 2012-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */

#include "Real.h"
#include "Vect2.h"
#include "TCAS2D.h"
#include "Util.h"
#include "Horizontal.h"

namespace larcfm {

// Compute modified tau
real TCAS2D::tau_mod(real DMOD, const Vect2& s, const Vect2& v) {
  real sdotv = s.dot(v);

  if (Util::almost_equals(sdotv,0)) // [CAM] Changed from == to almost_equals to mitigate numerical problems
    return 0;
  return (Util::sq(DMOD)-s.sqv())/sdotv;
}

bool TCAS2D::horizontal_RA(real DMOD, real Tau, const Vect2& s, const Vect2& v) {
  if (s.dot(v) >= 0) return s.norm() <= DMOD;
  else return s.norm() <= DMOD || tau_mod(DMOD,s,v) <= Tau;
}

bool TCAS2D::horizontal_RA_at(real DMOD, real Tau, const Vect2& s, const Vect2& v, real t) {
  const Vect2& sat = v.ScalAdd(t,s);
  return horizontal_RA(DMOD,Tau,sat,v);
}

real TCAS2D::nominal_tau(real B, real T, const Vect2& s, const Vect2& v, real rr) {
  if (v.isZero())
    return B;
  return Util::max(B,Util::min(T,-s.dot(v) / v.sqv()-rr/2));
}

real TCAS2D::time_of_min_tau(real DMOD, real B, real T, const Vect2& s, const Vect2& v) {
  if (v.ScalAdd(B,s).dot(v) >= 0)
    return B;
  real d = Horizontal::Delta(s,v,DMOD);
  real rr = 0;
  if (d < 0)
    rr = 2*std::sqrt(-d) / v.sqv();
  if (v.ScalAdd(T,s).dot(v) < 0)
    return T;
  return nominal_tau(B,T,s,v,rr);
}

bool TCAS2D::min_tau_undef(real DMOD, real B, real T, Vect2 s, Vect2 v) {
  return Horizontal::Delta(s,v,DMOD) >= 0 &&
      v.ScalAdd(B,s).dot(v) < 0 && v.ScalAdd(T,s).dot(v) >= 0;
}


bool TCAS2D::RA2D(real DMOD, real Tau, real B, real T, const Vect2& s, const Vect2& v) {
  if (min_tau_undef(DMOD,B,T,s,v))
    return true;
  real mt = time_of_min_tau(DMOD,B,T,s,v);
  return horizontal_RA_at(DMOD,Tau,s,v,mt);
}

void TCAS2D::RA2D_interval(real DMOD, real Tau, real B, real T, Vect2 s, Vect2 vo, Vect2 vi) {
  t_in = B;
  t_out = T;
  Vect2 v = vo.Sub(vi);
  real sqs = s.sqv();
  real sdotv = s.dot(v);
  real sqD = Util::sq(DMOD);
  if (vo.almostEquals(vi) && sqs <= sqD)
    return;
  real sqv = v.sqv();
  if (sqs <= sqD) {
    t_out = Util::root2b(sqv,sdotv,sqs-sqD,1);
    return;
  }
  real b = 2*sdotv+Tau*sqv;
  real c = sqs+Tau*sdotv-sqD;
  real discriminant = Util::discr(sqv,b,c);
  if (sdotv >= 0 || discriminant < 0) {
    t_in = T+1;
    t_out = 0;
    return;
  }
  if (Horizontal::Delta(s,v,DMOD) >= 0) {
    t_in = Util::root(sqv,b,c,-1);
    t_out = Util::root2b(sqv,sdotv,sqs-sqD,1);
  } else {
    const auto times = Util::roots(sqv,b,c);
    t_in = times.first;
    t_out = times.second;
  }
}

}
