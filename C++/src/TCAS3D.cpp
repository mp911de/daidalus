/*
 * Copyright (c) 2012-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */
#include "Real.h"
#include "TCAS3D.h"
#include "TCAS2D.h"
#include "CD2D.h"
#include "CD3D.h"
#include "Vect3.h"
#include "TCASTable.h"
#include "Util.h"
#include "Vertical.h"
#include "ConflictData.h"
#include "Triple.h"
#include "CDCylinder.h"
#include "WCV_TAUMOD.h"
#include <cfloat>
#include <limits>
#include <cmath>

namespace larcfm {

TCAS3D::TCAS3D() : table_(TCASTable::make_TCASII_Table(true)), id("") {}

TCAS3D::TCAS3D(const TCASTable& tab) : table_(tab), id("") {}

/**
 * @return one static TCAS3D
 */
const TCAS3D& TCAS3D::A_TCAS3D() {
  static TCAS3D dwc;
  return dwc;
}

/**
 * @return one static TCASII_RA
 */
const TCAS3D& TCAS3D::TCASII_RA() {
  static TCAS3D ra = make_TCASII_RA();
  return ra;
}

/**
 * @return one static TCASII_TA
 */
const TCAS3D& TCAS3D::TCASII_TA() {
  static TCAS3D ta = make_TCASII_TA();
  return ta;
}

/** Make TCAS3D object with empty Table **/
TCAS3D TCAS3D::make_Empty() {
  TCAS3D tcas3d;
  tcas3d.table_.clear();
  return tcas3d;
}

/** Make TCAS3D object with an RA Table **/
TCAS3D TCAS3D::make_TCASII_RA() {
  TCAS3D tcas3d;
  return tcas3d;
}

/** Make TCAS3D objec with a TA Table **/
TCAS3D TCAS3D::make_TCASII_TA() {
  TCAS3D tcas3d;
  tcas3d.table_.setDefaultTCASIIThresholds(false);
  return tcas3d;
}

TCASTable& TCAS3D::getTCASTable() {
  return table_;
}

/**
 * Set table to TCASII Thresholds (RA Table when ra is true, TA Table when ra is false)
 */
void TCAS3D::setDefaultTCASIIThresholds(bool ra) {
  table_.setDefaultTCASIIThresholds(ra);
}

LossData TCAS3D::conflictInterval(const Vect3& so, const Vect3& vo, const Vect3& si, const Vect3& vi, real B, real T) const {
  return RA3D_times(so,vo,si,vi,B,T,nullptr);
}

LossData TCAS3D::conflictIntervalWithTrafficState(const TrafficState& ownship, const TrafficState& intruder, real B, real T) const {
  return conflictInterval(ownship.get_s(),ownship.get_v(),intruder.get_s(),intruder.get_v(),B,T);
}

ConflictData TCAS3D::conflictDetection(const Vect3& so, const Vect3& vo, const Vect3& si, const Vect3& vi, real B, real T) const {
  return RA3D(so,vo,si,vi,B,T);
}

// pointer to new instance of this object
TCAS3D* TCAS3D::make() const {
  return new TCAS3D();
}

// pointer to deep copy of current object
TCAS3D* TCAS3D::copy() const {
  return new TCAS3D(*this);
}

bool TCAS3D::vertical_RA(real sz, real vz, real ZTHR, real TCOA) {
  if (std::abs(sz) <= ZTHR) return true;
  if (Util::almost_equals(vz,0)) return false; // [CAM] Changed from == to almost_equals to mitigate numerical problems
  real tcoa = Vertical::time_coalt(sz,vz);
  return 0 <= tcoa && tcoa <= TCOA;
}

bool TCAS3D::cd2d_TCAS_after(real HMD, Vect2 s, Vect2 vo, Vect2 vi, real t) {
  Vect2 v = vo.Sub(vi);
  return
      (vo.almostEquals(vi) && s.sqv() <= Util::sq(HMD)) ||
      (v.sqv() > 0 && Horizontal::Delta(s,v,HMD) >= 0 &&
          Horizontal::Theta_D(s,v,1,HMD) >= t);
}

bool TCAS3D::cd2d_TCAS(real HMD, Vect2 s, Vect2 vo, Vect2 vi) {
  return cd2d_TCAS_after(HMD,s,vo,vi,0);
}

// if true, then ownship has a TCAS resolution advisory at current time
bool TCAS3D::TCASII_RA(const Vect3& so, const Vect3& vo, const Vect3& si, const Vect3& vi) const {

  Vect2 so2 = so.vect2();
  Vect2 si2 = si.vect2();
  Vect2 s2 = so2.Sub(si2);
  Vect2 vo2 = vo.vect2();
  Vect2 vi2 = vi.vect2();
  Vect2 v2 = vo2.Sub(vi2);
  int sl = table_.getSensitivityLevel(so.z());
  bool usehmdf = table_.getHMDFilter();
  real TAU  = table_.getTAU(sl);
  real TCOA = table_.getTCOA(sl);
  real DMOD = table_.getDMOD(sl);
  real HMD  = table_.getHMD(sl);
  real ZTHR = table_.getZTHR(sl);

  return (!usehmdf || cd2d_TCAS(HMD,s2,vo2,vi2)) &&
      TCAS2D::horizontal_RA(DMOD,TAU,s2,v2) &&
      vertical_RA(so.z()-si.z(),vo.z()-vi.z(),ZTHR,TCOA);
}

// if true, within lookahead time interval [B,T], the ownship has a TCAS resolution advisory (effectively conflict detection)
// B must be non-negative and T > B

ConflictData TCAS3D::RA3D(const Vect3& so, const Vect3& vo, const Vect3& si, const Vect3& vi, real B, real T) const {
  real tmin;
  LossData interval = RA3D_times(so,vo,si,vi,B,T,&tmin);
  Vect3 s = so.Sub(si);
  Vect3 v = vo.Sub(vi);
  int sl = table_.getMaxSensitivityLevel();
  real dmin = s.linear(v,tmin).cyl_norm(table_.getDMOD(sl),table_.getZTHR(sl));
  return ConflictData(interval,tmin,dmin,s,v);
}

LossData TCAS3D::RA3D_times(const Vect3& so, const Vect3& vo, const Vect3& si, const Vect3& vi, real B, real T, real* critical_time) const {
  Vect2 so2 = so.vect2();
  Vect2 vo2 = vo.vect2();
  Vect2 si2 = si.vect2();
  Vect2 vi2 = vi.vect2();

  real tin = INFINITY;
  real tout = -INFINITY;
  real tmin = INFINITY;
  int sl_first = table_.getSensitivityLevel(so.z()+B*vo.z());
  int sl_last = table_.getSensitivityLevel(so.z()+T*vo.z());
  if (sl_first == sl_last || Util::almost_equals(vo.z(),real(0.0))) {
    Triple<real,real,real> ra3dint = RA3D_interval(sl_first,so2,so.z(),vo2,vo.z(),si2,si.z(),vi2,vi.z(),B,T,critical_time != nullptr);
    tin = ra3dint.first;
    tout = ra3dint.second;
    tmin = ra3dint.third;
  } else {
    int sl = sl_first;
    for (real t_B = B; t_B < T; sl = sl_first < sl_last ? sl+1 : sl-1) {
      if (table_.isValidSensitivityLevel(sl)) {
        real level = sl_first < sl_last ? table_.getLevelAltitudeUpperBound(sl) :table_.getLevelAltitudeLowerBound(sl);
        real t_level = !ISFINITE(level) ? INFINITY :(level-so.z())/vo.z();
        Triple<real,real,real> ra3dint = RA3D_interval(sl,so2,so.z(),vo2,vo.z(),si2,si.z(),vi2,vi.z(),t_B,Util::min(t_level,T),critical_time != nullptr);
        if (Util::almost_less(ra3dint.first,ra3dint.second)) {
          tin = Util::min(tin,ra3dint.first);
          tout = Util::max(tout,ra3dint.second);
        }
        tmin = Util::min(tmin,ra3dint.third);
        t_B = t_level;
        if (sl == sl_last) {
          break;
        }
      }
    }
  }
  if (critical_time) *critical_time = tmin;
  return LossData(tin,tout);
}

// Assumes 0 <= B < T
Triple<real,real,real> TCAS3D::RA3D_interval(int sl, const Vect2& so2, real soz, const Vect2& vo2, real voz,
    const Vect2& si2, real siz, const Vect2& vi2, real viz, real B, real T) const {
  return RA3D_interval(sl,so2,soz,vo2,voz,si2,siz,vi2,viz,B,T,true);
}

Triple<real,real,real> TCAS3D::RA3D_interval(int sl, const Vect2& so2, real soz, const Vect2& vo2, real voz,
    const Vect2& si2, real siz, const Vect2& vi2, real viz, real B, real T, bool severity) const {
  Vect2 s2 = so2.Sub(si2);
  Vect2 v2 = vo2.Sub(vi2);
  real sz = soz-siz;
  real vz = voz-viz;
  bool usehmdf = table_.getHMDFilter();
  real TAU  = table_.getTAU(sl);
  real TCOA = table_.getTCOA(sl);
  real DMOD = table_.getDMOD(sl);
  real HMD  = table_.getHMD(sl);
  real ZTHR = table_.getZTHR(sl);

  auto result = [&](real in, real out, real begin, real end) {
    real critical = severity ? TCAS2D::time_of_min_tau(DMOD,begin,end,s2,v2) : PINFINITY;
    return Triple<real,real,real>(in,out,critical);
  };
  bool level = Util::almost_equals(voz,viz);
  if ((usehmdf && !cd2d_TCAS_after(HMD,s2,vo2,vi2,B)) || (level && std::abs(sz) > ZTHR)) {
    return result(T,B,B,T);
  }
  real entry = B;
  real exit = T;
  if (!level) {
    real height = Util::max(ZTHR,std::abs(vz)*TCOA);
    entry = Vertical::Theta_H(sz,vz,-1,height);
    exit = Vertical::Theta_H(sz,vz,1,ZTHR);
  }
  if (exit < B || T < entry) return result(T,B,B,T);

  real begin = Util::max(B,entry);
  real end = Util::min(T,exit);
  TCAS2D ra;
  ra.RA2D_interval(DMOD,TAU,begin,end,s2,vo2,vi2);
  if (ra.t_in > ra.t_out || ra.t_out < begin || ra.t_in > end) return result(T,B,B,T);
  real in = Util::max(begin,Util::min(end,ra.t_in));
  real out = Util::max(begin,Util::min(end,ra.t_out));
  if (usehmdf && HMD < DMOD) {
    Vect2 at_entry = v2.ScalAdd(entry,s2);
    if (at_entry.dot(v2) >= 0 && !(at_entry.sqv() <= Util::sq(HMD))) return result(T,B,B,T);
    real exit_theta = T;
    if (v2.sqv() > 0) exit_theta = Util::max(B,Util::min(Horizontal::Theta_D(s2,v2,1,HMD),T));
    out = Util::min(out,exit_theta);
    if (!(in <= out)) return result(in,out,B,T);
  }
  return result(in,out,in,out);
}


/**
 * Returns TAU threshold for sensitivity level sl in seconds
 */
real TCAS3D::getTAU(int sl) const {
  return table_.getTAU(sl);
}

/**
 * Returns TCOA threshold for sensitivity level sl in seconds
 */
real TCAS3D::getTCOA(int sl) const {
  return table_.getTCOA(sl);
}

/**
 * Returns DMOD for sensitivity level sl in internal units.
 */
real TCAS3D::getDMOD(int sl) const {
  return table_.getDMOD(sl);
}

/**
 * Returns DMOD for sensitivity level sl in u units.
 */
real TCAS3D::getDMOD(int sl, const std::string& u) const {
  return table_.getDMOD(sl,u);
}

/**
 * Returns Z threshold for sensitivity level sl in internal units.
 */
real TCAS3D::getZTHR(int sl) const {
  return table_.getZTHR(sl);
}

/**
 * Returns Z threshold for sensitivity level sl in u units.
 */
real TCAS3D::getZTHR(int sl,const std::string& u) const {
  return table_.getZTHR(sl,u);
}

/**
 * Returns HMD for sensitivity level sl in internal units.
 */
real TCAS3D::getHMD(int sl) const {
  return table_.getHMD(sl);
}

/**
 * Returns HMD for sensitivity level sl in u units.
 */
real TCAS3D::getHMD(int sl, const std::string& u) const {
  return table_.getHMD(sl,u);
}

/** Modify the value of Tau Threshold for a given sensitivity level (2-8)
 * Parameter val is given in seconds
 */
void TCAS3D::setTAU(int sl, real val) {
  table_.setTAU(sl,val);
}

/** Modify the value of TCOA Threshold for a given sensitivity level (2-8)
 * Parameter val is given in seconds
 */
void TCAS3D::setTCOA(int sl, real val) {
  table_.setTCOA(sl,val);
}

/** Modify the value of DMOD for a given sensitivity level (2-8)
 * Parameter val is given in internal units
 */
void TCAS3D::setDMOD(int sl, real val) {
  table_.setDMOD(sl, val);
}

/** Modify the value of DMOD for a given sensitivity level (2-8)
 * Parameter val is given in u units
 */
void TCAS3D::setDMOD(int sl, real val, const std::string& u) {
  table_.setDMOD(sl,val,u);
}

/** Modify the value of ZTHR for a given sensitivity level (2-8)
 * Parameter val is given in internal units
 */
void TCAS3D::setZTHR(int sl, real val) {
  table_.setZTHR(sl,val);
}

/** Modify the value of ZTHR for a given sensitivity level (2-8)
 * Parameter val is given in u units
 */
void TCAS3D::setZTHR(int sl, real val, const std::string& u) {
  table_.setZTHR(sl,val,u);
}

/**
 * Modify the value of HMD for a given sensitivity level (2-8)
 * Parameter val is given in internal units
 */
void TCAS3D::setHMD(int sl, real val) {
  table_.setHMD(sl,val);
}

/**
 * Modify the value of HMD for a given sensitivity level (2-8)
 * Parameter val is given in u units
 */
void TCAS3D::setHMD(int sl, real val, const std::string& u) {
  table_.setHMD(sl,val,u);
}

void TCAS3D::setHMDFilter(bool flag) {
  table_.setHMDFilter(flag);
}

bool TCAS3D::getHMDFilter() {
  return table_.getHMDFilter();
}

std::string TCAS3D::toString() const {
  return (id == "" ? "" : id+" = ")+getSimpleClassName()+": {"+table_.toString()+"}";
}

std::string TCAS3D::toPVS() const {
  return getSimpleClassName()+"("+table_.toPVS()+")";
}

ParameterData TCAS3D::getParameters() const {
  ParameterData p;
  updateParameterData(p);
  return p;
}

void TCAS3D::updateParameterData(ParameterData& p) const {
  table_.updateParameterData(p);
  p.set("id",id);
}

void TCAS3D::setParameters(const ParameterData& p) {
  table_.setParameters(p);
  if (p.contains("id")) {
    id = p.getString("id");
  }
}

std::string TCAS3D::getSimpleClassName() const {
  return "TCAS3D";
}


std::string TCAS3D::getIdentifier() const {
  return id;
}

void TCAS3D::setIdentifier(const std::string& s) {
  id = s;
}

bool TCAS3D::equals(const Detection3D& cd) const {
  if (!larcfm::equals(getCanonicalClassName(), cd.getCanonicalClassName())) return false;
  if (!larcfm::equals(id, cd.getIdentifier())) return false;
  if (!table_.equals(((TCAS3D&)cd).table_)) return false;
  return true;
}

bool TCAS3D::contains(const Detection3D& cd) const {
  if (larcfm::equals(getCanonicalClassName(), cd.getCanonicalClassName())) {
    const TCAS3D& td = (TCAS3D&)cd;
    return table_.contains(td.table_);
  }
  return false;
}

void TCAS3D::horizontalHazardZone(std::vector<Position>& haz, const TrafficState& ownship, const TrafficState& intruder,
    real T) const {
  int sl = table_.getSensitivityLevel(ownship.altitude());
  bool usehmdf = table_.getHMDFilter();
  real TAUMOD  = table_.getTAU(sl);
  real DMOD = Util::max(table_.getDMOD(sl),table_.getHMD(sl));
  haz.clear();
  Position po = ownship.getPosition();
  Velocity v = ownship.getAirVelocity().Sub(intruder.getAirVelocity().vect3());
  if (Util::almost_equals(TAUMOD+T,0) || Util::almost_equals(v.vect3().norm2D(),0)) {
    CDCylinder::circular_arc(haz,po,Velocity::mkVxyz(DMOD,0,0),2*Pi,false);
  } else {
    Vect3 sD = Horizontal::unit_perpL(v.vect3()).Scal(DMOD);
    Velocity vD = Velocity::make(sD);
    CDCylinder::circular_arc(haz,po,vD,Pi,usehmdf);
    Position TAU_center = WCV_TAUMOD::TAU_center(po,v,TAUMOD,T);
    Vect3 vC = v.vect3().Scal(real(0.5)*TAUMOD);     // TAUMOD Center (relative)
    if (usehmdf) {
      Vect3 vDC = vC.Sub(vD.vect3()); // Far end point opposite to -vD (vC-relative);
      Vect3 nvDC = vC.Add(vD.vect3()); // Far end point opposite to vD (vC-relative);
      real sqa = vDC.sqv2D();
      real alpha = Util::atan2_safe(vDC.det2D(nvDC)/sqa,vDC.dot2D(nvDC)/sqa);
      Velocity velDC = Velocity::make(vDC);
      CDCylinder::circular_arc(haz,TAU_center,velDC,alpha,true);
    } else {
      Vect3 nsCD=sD.Neg().Sub(vC);
      Vect3 sCD=sD.Sub(vC);
      real sqa = sCD.sqv2D();
      Velocity nvCD = Velocity::make(nsCD);
      if (Util::almost_equals(T,0)) { // Two circles: DMOD and TAUMO. They intersect at +/-vD
        real alpha = Util::atan2_safe(nsCD.det2D(sCD)/sqa,nsCD.dot2D(sCD)/sqa);
        CDCylinder::circular_arc(haz,TAU_center,nvCD,alpha,false);
      } else { // Two circles: DMOD and TAUMOD. They intersect at +/- vD.
        Vect3 sT = Horizontal::unit_perpL(v.vect3()).Scal(std::sqrt(sqa));
        Velocity vT = Velocity::make(sT);
        Vect3 nsT = sT.Neg();
        Velocity nvT = Velocity::make(nsT);
        real alpha = Util::atan2_safe(nsCD.det2D(nsT)/sqa,nsCD.dot2D(nsT)/sqa);
        Position TAU_center_0 = WCV_TAUMOD::TAU_center(po,v,TAUMOD,0);
        CDCylinder::circular_arc(haz,TAU_center_0,nvCD,alpha,true);
        CDCylinder::circular_arc(haz,TAU_center,nvT,Pi,true);
        CDCylinder::circular_arc(haz,TAU_center_0,vT,alpha,false);                }
    }
  }
}


}
