/*
 * Copyright (c) 2012-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */

#ifndef TCAS2D_H_
#define TCAS2D_H_
#include "Real.h"

#include "Vect2.h"

namespace larcfm {

class TCAS2D {
  
public:

  real t_in;
  real t_out;

  // Compute modified tau
  static real tau_mod(real DMOD, const Vect2& s, const Vect2& v);
  
  static bool horizontal_RA(real DMOD, real Tau, const Vect2& s, const Vect2& v);
  
  static bool horizontal_RA_at(real DMOD, real Tau, const Vect2& s, const Vect2& v, real t);
  
  static real nominal_tau(real B, real T, const Vect2& s, const Vect2& v, real rr);
  
  static real time_of_min_tau(real DMOD, real B, real T, const Vect2& s, const Vect2& v);

  static bool min_tau_undef(real DMOD, real B, real T, Vect2 s, Vect2 v);

  static bool RA2D(real DMOD, real Tau, real B, real T, const Vect2& s, const Vect2& v);

  void RA2D_interval(real DMOD, real Tau, real B, real T, Vect2 s, Vect2 vo, Vect2 vi);

  
};
}
#endif
