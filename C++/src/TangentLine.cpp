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

static Vect2 tangent_point(const Vect2& s, real alpha, real beta, int eps) {
  return Vect2(alpha*s.x+eps*beta*s.y,alpha*s.y-eps*beta*s.x);
}

Vect2 TangentLine::Q(const Vect2& s, const real D, const int eps) {
  real sq_s = s.sqv();
  real sq_D = sq(D);
  real delta = sq_s-sq_D;
  if (!(delta >= 0)) return Vect2();
  return tangent_point(s,sq_D/sq_s,D*sqrt_safe(delta)/sq_s,eps);
}

static Vect2 tangent_direction(const Vect2& s, const Vect2& q) {
  return q.isZero() ? q : q.Sub(s);
}

std::pair<Vect2,Vect2> TangentLine::directions(const Vect2& s, real D) {
  real sq_s = s.sqv();
  real sq_D = sq(D);
  if (Util::almost_equals(sq_s,sq_D)) {
    return std::make_pair(Vect2(-s.y,s.x),Vect2(s.y,-s.x));
  }
  real delta = sq_s-sq_D;
  if (!(delta >= 0)) return std::make_pair(Vect2(),Vect2());
  real alpha = sq_D/sq_s;
  real beta = D*sqrt_safe(delta)/sq_s;
  return std::make_pair(tangent_direction(s,tangent_point(s,alpha,beta,-1)),
      tangent_direction(s,tangent_point(s,alpha,beta,1)));
}

TangentLine::TangentLine(const Vect2& s, const real D, const int eps)
  : Vect2(Util::almost_equals(s.sqv(),sq(D)) ? Vect2(eps*s.y,-eps*s.x) :
      tangent_direction(s,Q(s,D,eps))), eps(eps) {}

bool TangentLine::horizontal_criterion(const Vect2& v) const {
  return eps*det(v) >= 0;    
}

int TangentLine::get_eps() const {
    return eps;
}

}
