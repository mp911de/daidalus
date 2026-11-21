/*
 * Vertical.h
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

#ifndef VERTICAL_H_
#define VERTICAL_H_
#include "Real.h"

#include "Consts.h"
#include "Vect3.h"

namespace larcfm {

/** Vertical solution */
class Vertical {
public:

  /** z solution */
  real  z;

  /** Construct a vertical solution */
  Vertical();

  /** Construct a vertical solution */
  explicit Vertical(const real vz);

  /** Is this solution undefined? */
  bool undef() const;

  /** method */
  Vertical add_this(const real vz);

  /** Exact equality */
  bool operator == (const Vertical& v) const;
  /** Is any component not exactly equal */
  bool operator != (const Vertical& v) const;

  /** "Solution" indicating no solution */
  static const Vertical& NoVerticalSolution();

  /** */
  static bool almost_vertical_los(const real sz,const real H);

  /** */
  static real Theta_H(const real sz, const real v, const int eps, const real H);

  static real time_coalt(real sz, real vz);

  /* Vertical miss distance within lookahead time */
  static real vmd(real sz, real vz, real T);

  /** Solve the following equation on vz:
   *   sz+t*vz = eps*H,
   *
   * where t = Theta_D(s,v,eps).
   * eps determines the bottom, i.e.,-1, or top, i.e., 1, circle.
   */
  static Vertical vs_circle(const Vect3& s, const Vect3& vo, const Vect3& vi,
      const int eps, const real D, const real H);

  /** */
  static Vertical vs_circle_at(const real sz, const real viz,
      const real t, const int eps, const int dir, const real H);

  /** */
  static Vertical vs_los_recovery(const Vect3& s, const Vect3& vo, const Vect3& vi,
      const real H, const real t, int epsv);

  std::string toString();

private:
  bool udef;

};


}

#endif /* VERTICAL_H_ */
