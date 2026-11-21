/*
 * TangentLine.cpp
 * Release: ACCoRD++-2.b (08/22/10)
 *
 * Contact: Cesar Munoz (cesar.a.munoz@nasa.gov)
 * NASA LaRC
 * http://shemesh.larc.nasa.gov/people/cam/ACCoRD
 *
 * Computes a vector that is tangent to the protected zone.
 *
 * Copyright (c) 2011-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */

#include "Real.h"
#include "TangentLine.h"

namespace larcfm {

TangentLine::TangentLine() : Vect2(real(0.0),real(0.0)) {
  eps = 0;
  //x = y = 0.0;
}

Vect2 TangentLine::Q(const Vect2& s, const real D, const int eps) {
  real sq_s = s.sqv();
  real sq_D = sq(D);
  real delta = sq_s -sq_D;
  if (delta >= 0) { 
    real alpha = sq_D/sq_s;
    real beta  = D*sqrt_safe(delta)/sq_s;
    return Vect2(alpha*s.x+eps*beta*s.y,
		 alpha*s.y-eps*beta*s.x);   
    }
    return Vect2();
}

real InitTangentLineX (const Vect2& s, const real D, const int eps) {
  real x;
//  double y;
  real sq_s  = s.sqv();
  real sq_D  = sq(D);
  if (Util::almost_equals(sq_s,sq_D))  {
    x = eps*s.y;
//    y = -eps*s.x;
  } else {
    Vect2 q = TangentLine::Q(s,D,eps);
    x = q.x;
//    y = q.y;
    if (!q.isZero()) {
      x -= s.x;
//      y -= s.y;
    }
  }
  return x;
}

real InitTangentLineY (const Vect2& s, const real D, const int eps) {
//  double x;
  real y;
  real sq_s  = s.sqv();
  real sq_D  = sq(D);
  if (Util::almost_equals(sq_s,sq_D))  {
//    x = eps*s.y;
    y = -eps*s.x;
  } else {
    Vect2 q = TangentLine::Q(s,D,eps);
//    x = q.x;
    y = q.y;
    if (!q.isZero()) {
//      x -= s.x;
      y -= s.y;
    }
  }
  return y;
}

TangentLine::TangentLine(const Vect2& s, const real D, const int eps) : Vect2(InitTangentLineX(s, D, eps), InitTangentLineY(s, D, eps) ) {
  this->eps = eps;
}

bool TangentLine::horizontal_criterion(const Vect2& v) const {
  return eps*det(v) >= 0;    
}

int TangentLine::get_eps() const {
    return eps;
}

}
