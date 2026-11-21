/*
 * Vertical.cpp
 * Release: ACCoRD++-2.b (08/22/10)
 *
 * Contact: Cesar Munoz (cesar.a.munoz@nasa.gov)
 * NASA LaRC
 * http://shemesh.larc.nasa.gov/people/cam/ACCoRD
 *
 * Copyright (c) 2011-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */

#include "Real.h"
#include "Vertical.h"
#include "Horizontal.h"
#include <limits>
#include "format.h"

namespace larcfm {

Vertical::Vertical() : z(NaN), udef(true) {}

Vertical::Vertical(const real vz) : z(vz), udef(false) {}

bool Vertical::operator == (const Vertical& v) const {  // strict equality
  return ! udef && ! v.udef && z==v.z;
}

bool Vertical::operator != (const Vertical& v) const {  // strict disequality
  return !(*this == v);
}

const Vertical& Vertical::NoVerticalSolution() {
  const static Vertical v;
  return v;
}

Vertical Vertical::add_this(const real vz) {
  if (!udef) {
    z += vz;
  }
  return (*this);
}

bool Vertical::undef() const {
  return udef;
}

bool Vertical::almost_vertical_los(const real sz, const real H) {
  real absz = std::abs(sz);
  return !Util::almost_equals(absz,H) && absz < H;
}

// Computes times when sz,vz intersects rectangle of half height H
// eps = -1 : larcfm::Entry
// eps =  1 : larcfm::Exit
real Vertical::Theta_H(const real sz, const real vz, const int eps, const real H) {
  if (vz == 0) return NaN;
  return (eps*sign(vz)*H-sz)/vz;
}

real Vertical::time_coalt(real sz, real vz) {
  if (sz == 0) return 0;
  if (vz == 0) return NaN;
  return -sz/vz;
}

/* Vertical miss distance within lookahead time */
real Vertical::vmd(real sz, real vz, real T) {
  if (sz*vz < 0) {
    // aircraft are vertically converging
    if (time_coalt(sz,vz) <= T) {
      return 0;
    } else {
      return std::abs(sz+T*vz);
    }
  }
  return std::abs(sz);
}

Vertical vs_at(const real sz, const real t,
    const int eps, const real H) {
  if (t == 0)
    return Vertical::NoVerticalSolution();
  return Vertical((eps*H-sz)/t);
}

Vertical vs_only(const Vect3& s, const Vect2& v, const real t,
    const int eps, const real D, const real H) {

  const Vect2& s2 = s.vect2();

  if (eps*s.z() < H && s2.sqv() > sq(D) && Horizontal::Delta(s2,v,D) > 0)
    return vs_at(s.z(),Horizontal::Theta_D(s2,v,larcfm::Entry,D),eps,H);
  else if (eps*s.z() >= H && t > 0)
    return vs_at(s.z(),t,eps,H);
  return Vertical::NoVerticalSolution();
}

/* Solve the following equation on vz:
 *   sz+t*vz = eps*H,
 *
 * where t = Theta_D(s,v,epsp).
 * eps determines the bottom, i.e.,-1, or top, i.e., 1, circle.
 */
Vertical Vertical::vs_circle(const Vect3& s, const Vect3& vo, const Vect3& vi,
    const int eps, const real D, const real H) {

  const Vect2& s2 = s.vect2();
  const Vect2& vo2 = vo.vect2();
  const Vect2& vi2 = vi.vect2();
  Vect2 v2 = vo2-vi2;

  if (vo2.almostEquals(vi2) && eps == sign(s.z()))
    return Vertical(vi.z());
  else if (Horizontal::Delta(s2,v2,D) > 0)
    return vs_only(s,v2,Horizontal::Theta_D(s2,v2,larcfm::Exit,D),eps,D,H).add_this(vi.z());
  return NoVerticalSolution();
}

Vertical Vertical::vs_circle_at(const real sz, const real viz,
    const real t, const int eps, const int dir, const real H) {
  if (t > 0 && dir*eps*sz <= dir*H)
    return vs_at(sz,t,eps,H).add_this(viz);
  return NoVerticalSolution();
}

Vertical Vertical::vs_los_recovery(const Vect3& s, const Vect3& vo, const Vect3& vi,
    const real H, const real t, int epsv) {
  if (t <= 0)
    return NoVerticalSolution();
  Vect3 v = vo-vi;
  real nvz = (epsv*H - s.z())/t;
  if (s.z()*v.z() >= 0 && std::abs(v.z()) >= std::abs(nvz))
    return Vertical(vo.z());
  else
    return Vertical(nvz+vi.z());
}

std::string Vertical::toString() {
  return " udef = ["+bool2str(udef)+" z = "+FmPrecision(z)+"]";
}
}
