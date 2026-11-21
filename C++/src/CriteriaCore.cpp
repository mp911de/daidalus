/**
 * CriteriaCore.cpp
 *
 * State-based Implicit Criteria
 *
 * Copyright (c) 2011-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */


#include "Real.h"
#include "CriteriaCore.h"
#include "format.h"
#include "Util.h"
#include "Kinematics.h"
#include "ACCoRDConfig.h"
#include <string>
#include <float.h>

namespace larcfm {

int CriteriaCore::horizontalCoordination(const Vect2& s, const Vect2& v) {
  return sign(v.det(s));
}

int CriteriaCore::horizontalCoordination(const Vect3& s, const Vect3& v) {
  return horizontalCoordination(s.vect2(), v.vect2());
}

/**
 * An algebraic way to determine a*sqrt(b) = e.
 * @return true, if a*sqrt(b) = e
 */
bool eq(real a,real b,real e) {
  return a*e >= 0 && sq(a)*b == sq(e);
}

/**
 * An algebraic way to determine a*sqrt(b) > e.
 * @return true, if a*sqrt(b) > e
 */
bool gt(real a,real b,real e) {
  return a >= 0 ? e < 0 || sq(a)*b > sq(e):
      e < 0 && sq(a)*b < sq(e);
}

int CriteriaCore::verticalCoordination(const Vect3& s, const Vect3& vo, const Vect3& vi, real D, real H, const std::string& ownship, const std::string& traffic) {
  if (CD3D::LoS(s,D,H)) {
    return verticalCoordinationLoS(s,vo,vi,ownship,traffic);
  } else {
    Vect3 v = vo.Sub(vi);
    return verticalCoordinationConflict(s,v,D,ownship,traffic);
  }
}

int CriteriaCore::verticalCoordinationConflict(const Vect3& s, const Vect3& v, real D, const std::string& ownship, const std::string& traffic) {
  Vect2 s2 = s.vect2();
  Vect2 v2 = v.vect2();
  real a = v2.sqv();
  real b = s2*v2;
  real c = s2.sqv()-sq(D);
  real d = sq(b)-a*c;
  if (Util::almost_equals(v.z(),0) || v2.isZero() ||
      d < 0 || eq(v.z(),Util::sq(b)-a*c,s.z()*a-v.z()*b)) {
    return breakSymmetry(s, ownship, traffic);
  }
  if (gt(v.z(),Util::sq(b)-a*c,s.z()*a-v.z()*b)) {
    return -1;
  }
  return 1;
}

int CriteriaCore::verticalCoordinationLoS(const Vect3& s, const Vect3& vo, const Vect3& vi, const std::string& ownship, const std::string& traffic) {
  int epsv;
  epsv = losr_vs_dir(s,vo,vi,ACCoRDConfig::NMAC_D, ACCoRDConfig::NMAC_H,ownship,traffic);
  return epsv;
}

real R(const Vect2& sp, real D) {
  return std::sqrt(sp.sqv()-sq(D))/D;
}

bool CriteriaCore::horizontal_criterion(const Vect2& sp, const  Vect2& v, real D, int epsh) {
  return sp.dot(v) >= R(sp,D)*epsh*sp.det(v);
}

bool  closed_region_3D(const Vect3& s,const Vect3& p, int eps, int dir, const Vect3& v, real D, real H) {
  Vect2 s2 = s.vect2();
  Vect2 v2 = v.vect2();
  Vect2 p2 = p.vect2();
  real t = (sq(D)-s2.dot(p2))/(v2.dot(p2));
  return v2.dot(p2) != 0 && sign(v2.dot(p2)) == dir && t >= 0 &&  eps*(s.z() + t*v.z()) >= H &&
      ((std::abs(s.z()) >= H && dir == eps*sign(s.z())) || (std::abs(s.z()) < H && dir == -1));
}

bool  CriteriaCore::vertical_criterion(int eps, const Vect3& s, const Vect3& v, const Vect3& nv, real D, real H) {
  int Entry = -1;
  int dir = Entry;
  if (std::abs(s.z()) >= H) dir = eps*sign(s.z());
  Vect2 s2 = s.vect2();
  Vect2 v2 = v.vect2();
  Vect3 p = Vect3(s2.Add(v2.Scal(Horizontal::Theta_D(s2,v2,dir,D))), eps*H);
  return ((v.x() == 0 && v.y() == 0) && eps*nv.z() >= 0 && eps*s.z() >= H)
      || (Horizontal::Delta(s2,v2,D) > 0 && Horizontal::Theta_D(s2,v2,dir,D) > 0
          && closed_region_3D(s,p,eps,dir,nv,D,H));
}

bool horizontal_los(const Vect2& s, real D) {
  return s.sqv() < sq(D);
}

bool vertical_los(real sz, real H) {
  return std::abs(sz) < H;
}

// sp = relative position
//  v = relative velocity
bool CriteriaCore::horizontal_criterion_0(const Vect2& sp, int eps, const Vect2& v, real D) {
  Vect2 vv = v;  // copy?
  if (Util::almost_equals(v.norm(),real(0.0),PRECISION13)) vv = Vect2::ZERO();        // criterion_3D's last term passes in nv - v which is exactly 0 for vertical solutions
  return horizontal_criterion(sp,vv,D,eps);
}

// from PVS
bool CriteriaCore::criterion_3D(const Vect3& sp, const Velocity&  v, int epsH, int epsV, const Velocity&  nv, real D, real H) {
  return (Horizontal::horizontal_sep(sp.vect2(),D) && horizontal_criterion_0(sp.vect2(),epsH,nv.vect2(),D)) ||
      (vertical_criterion(epsV,sp,v.vect3(),nv.vect3(),D,H) &&
          (horizontal_los(sp.vect2(),D) ||
              horizontal_criterion_0(sp.vect2(),epsH,nv.vect2().Sub(v.vect2()),D)));
}

// ***************************************** Repulsive Criteria ******************************************

// caD and caH are the diameter and height of a collision avoidance zone, e.g., 350ft x 100 ft
Vect3 CriteriaCore::vertical_decision_vect(const Vect3& s, const Vect3& vo, const Vect3& vi, real caD, real caH) {
  Vect2 s2 = s.vect2();
  Vect2 vo2 = vo.vect2();
  Vect2 vi2 = vi.vect2();
  Vect3 v = vo.Sub(vi);
  if ((!s.isZero() && CD3D::cd3d(s,vo,vi,caD,caH)) || Util::almost_equals(vo.z(),vi.z())) {
    return s;
  } else if (vo2.almostEquals(vi2) || s.isZero())  {
    return v;
  } else if (s2.dot(v.vect2()) <= 0) {
    return s.AddScal(CD2D::tcpa(s2,vo2,vi2),v);
  } else {
    return s;
  }
}

// Compute an absolute repulsive vertical direction
int CriteriaCore::losr_vs_dir(const Vect3& s, const Vect3& vo, const Vect3& vi,
    real caD, real caH, const std::string& ownship, const std::string& traffic) {
  int rtn = breakSymmetry(vertical_decision_vect(s,vo,vi,caD,caH),ownship,traffic);
  return rtn;
}

bool CriteriaCore::horizontal_los_criterion(const Vect2& s, const Vect2& vo, const Vect2& vi,const Vect2& nvo, int eps) {
  Vect2 v = vo.Sub(vi);
  Vect2 nv = nvo.Sub(vi);
  bool rtn = !s.isZero() && !nv.isZero() &&
      eps*s.det(v) <= 0 && eps*s.det(nv) <= 0
      && ((s.dot(v) < 0 &&  eps*nv.det(v) < 0)
          || (s.dot(v) >= 0 && (!v.isZero() || s.dot(nv) >= 0) &&
              (v.isZero() || s.dot(nv) > s.dot(v))));
  return rtn;
}

bool CriteriaCore::horizontalRepulsiveCriteria(const Vect3& s3, const Velocity& vo3, const Velocity& vi3, const Velocity& nvo3, int eps) {
  return horizontal_los_criterion(s3.vect2(),vo3.vect2(), vi3.vect2(),nvo3.vect2(),eps);
}

bool CriteriaCore::vs_bound_crit(const Vect3& v, const Vect3& nv, int eps) {
  Vect2 v2 = v.vect2();
  bool rtn;
  if (eps*v.z() > 0)
    rtn = eps*nv.z() > eps*v.z() && -eps*v.z()*nv.vect2().dot(v2) + eps*nv.z()*v2.sqv() >= 0;
  else
    rtn = eps*nv.z() >= 0;
  return rtn;
}

real CriteriaCore::min_rel_vert_speed(real vz, int eps, real minrelvs) {
  real rtn;
  if (eps*vz <= 0)
    rtn = minrelvs;
  else
    rtn = Util::max(minrelvs,std::abs(vz));
  return rtn;
}

bool CriteriaCore:: vertical_los_criterion(const Vect3& s, const Vect3& v, const Vect3& nv, int eps, real H, real minrelvs) {
  return std::abs(s.z()) < H &&
      vs_bound_crit(v,nv,eps) &&
      eps*nv.z() >= min_rel_vert_speed(v.z(),eps,minrelvs);
}

// caD and caH are the dimater and height of a collision avoidance zone, e.g., 350ft x 100 ft
// minrelvs is a minium relative vertical speed
bool CriteriaCore::verticalRepulsiveCriterion(const Vect3& s, const Vect3& vo, const Vect3& vi, const Vect3& nvo,
    real H, real minrelvs, int epsv) {
  bool rtn =  vertical_los_criterion(s,vo.Sub(vi),nvo.Sub(vi),epsv,H,minrelvs);
  return rtn;
}

/** Perform a symmetry calculation */
int CriteriaCore::breakSymmetry(const Vect3& s, const std::string& ownship, const std::string& traffic) {
 if (Util::almost_equals(s.z(),0)) {
    return Util::less_or_equal(ownship, traffic) ? 1 : -1;
  } else if (s.z() > 0) {
    return 1;
  } else {
    return -1;
  }
}

static bool trkChanged(const Velocity& vo, const Velocity& nvo) {
  return std::abs(vo.trk() - nvo.trk()) > Units::from("deg",real(0.001));
}

static bool gsChanged(const Velocity& vo, const Velocity& nvo) {
  return std::abs(vo.gs() - nvo.gs()) > Units::from("kn",real(0.001));
}

static bool vsChanged(const Velocity& vo, const Velocity& nvo) {
  return std::abs(vo.vs() - nvo.vs()) > Units::from("fpm",real(0.001));
}

bool CriteriaCore::criteria(const Vect3& s, const Velocity&  vo, const Velocity&  vi, const Velocity& nvo,
    real minRelVs, real D, real H, int epsh, int epsv) {
  if (horizontal_los(s.vect2(),D) && vertical_los(s.z(),H)) {
    bool horizChange = trkChanged(vo,nvo) || gsChanged(vo,nvo);
    bool vertChange = vsChanged(vo,nvo);
    bool vlc, hlc;
    vlc = verticalRepulsiveCriterion(s,vo.vect3(),vi.vect3(),nvo.vect3(),H, minRelVs, epsv);
    hlc = horizontalRepulsiveCriteria(s,  vo, vi, nvo, epsh);
    if (horizChange && vertChange) return hlc && vlc;
    else if (horizChange) return hlc;
    else if (vertChange) return vlc;
    else return (hlc || vlc);
  } else {
    Velocity v = vo.Sub(vi.vect3());
    Velocity nv = nvo.Sub(vi.vect3());
    return criterion_3D(s,v,epsh,epsv,nv,D,H);
  }
}

// -------------------- repulsive iterative search directions --------------------------------

// Compute a new track only vector that is one step to the dir of vo
// dir = 1 is right
Vect2 CriteriaCore::incr_trk_vect(const Vect2& vo, real step, int dir) {
  Vect2 rtn = vo.Scal(std::cos(dir*step)).Add(vo.PerpR().Scal(std::sin(dir*step)));
  return rtn;
}

int CriteriaCore::losr_trk_iter_dir(const Vect2& s, const Vect2& vo, const Vect2& vi, real step, int eps) {
  if (horizontal_los_criterion(s,vo,vi,incr_trk_vect(vo,step,1),eps))
    return 1;
  if (horizontal_los_criterion(s,vo,vi,incr_trk_vect(vo,step,-1),eps))
    return -1;
  return 0;
}

int CriteriaCore::trkSearchDirection(const Vect3& s, const Vect3& vo, const Vect3& vi, int eps) {
  return losr_trk_iter_dir(s.vect2(),vo.vect2(),vi.vect2(),Units::from("deg",1), eps);
}

Vect2 CriteriaCore::incr_gs_vect(const Vect2& vo, real step, int dir) {
  real normvo = vo.norm();
  return vo.Scal((normvo+dir*step)/normvo);
}

int CriteriaCore::losr_gs_iter_dir(const Vect2& s, const Vect2& vo, const Vect2& vi, real mings, real maxgs, real step, int eps) {
  if (vo.norm() + step <= maxgs && horizontal_los_criterion(s,vo,vi,incr_gs_vect(vo,step,1),eps))
    return 1;
  if (vo.norm() - step >= mings && horizontal_los_criterion(s,vo,vi,incr_gs_vect(vo,step,-1),eps))
    return -1;
  return 0;
}

int CriteriaCore::gsSearchDirection(const Vect3& s, const Vect3& vo, const Vect3& vi, int eps) {
  real mings = 0; // Units::from("kn",150);
  real maxgs = std::numeric_limits<real>::max(); // Units::from("kn",700);
  return losr_gs_iter_dir(s.vect2(),vo.vect2(),vi.vect2(),mings, maxgs, Units::from("kn",1), eps);
}

int CriteriaCore::vsSearchDirection(int epsv) {
  return epsv;
}

int CriteriaCore::dataVsRateEpsilon(int epsv, real vsRate){
  int trafSrchDir = vsSearchDirection(epsv);
  int absDir = -1;
  if (vsRate >= 0) absDir = 1;
  if (absDir == trafSrchDir) return epsv;
  else return -epsv;
}

int CriteriaCore::dataTurnEpsilon(const Vect3& s, const Velocity& vo, const Velocity& vi, int epsh, real trackRate){
  int trafSrchDir = trkSearchDirection(s.Neg(),vi.vect3(),vo.vect3(), epsh);
  int absDir = -1;
  if (trackRate >= 0) absDir = 1;
  if (absDir == trafSrchDir) return epsh;
  else return -epsh;
}

bool CriteriaCore::horizontal_new_repulsive_criterion(const Vect2& s,  const Vect2& vo, const Vect2& vi, const Vect2& nvo, int eps) {
  Vect2 v = vo.Sub(vi);
  Vect2 nv = nvo.Sub(vi);
  bool rtn = !s.isZero() && !nv.isZero() &&
      eps*s.det(v) <= 0 && eps*s.det(nv) < 0
      && ((s.dot(v) < 0 &&  eps*nv.det(v) < 0)
          || (s.dot(v) >= 0
              && (!v.isZero() || s.dot(nv) >= 0)
              && (v.isZero() || s.dot(nv) > s.dot(v))
              && eps*nv.det(v) <= 0));
  return rtn;
}

bool CriteriaCore::horizontal_new_repulsive_criterion(const Vect3& s, const Vect3& vo, const Vect3& vi, const Vect3& nvo, int eps) {
  Vect2 s2 = s.vect2();
  Vect2 v = vo.Sub(vi).vect2();
  Vect2 nv = nvo.Sub(vi).vect2();
  bool rtn = !s.isZero() && !nv.isZero() &&
      eps*s2.det(v) <= 0 && eps*s2.det(nv) < 0
      && ((s2.dot(v) < 0 &&  eps*nv.det(v) < 0)
          || (s2.dot(v) >= 0
              && (!v.isZero() || s2.dot(nv) >= 0)
              && (v.isZero() || s2.dot(nv) > s2.dot(v))
              && eps*nv.det(v) <= 0));
  return rtn;
}

bool CriteriaCore::vertical_new_repulsive_criterion(const Vect3& s, const Vect3& vo, const Vect3& vi, const Vect3& nvo, int eps) {
  Vect3 v = vo.Sub(vi);
  Vect3 nv = nvo.Sub(vi);
  Vect2 v2 = v.vect2();
  return eps*nv.z() > eps*v.z() && -eps*v.z()*nv.vect2().dot(v2) + eps*nv.z()*v2.sqv() >= 0;
}

}
