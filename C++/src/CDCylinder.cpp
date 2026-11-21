/*
 * CD3D.h 
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
#include "CDCylinder.h"
#include "LossData.h"
#include "LossData.h"
#include "ConflictData.h"
#include "DaidalusParameters.h"
#include <cmath>

#include "CD3D.h"
#include "format.h"

namespace larcfm {

/**
 * Instantiates a new CD3D object.
 */
CDCylinder::CDCylinder(const std::string& s) : 
            D_(Units::from("nmi", real(5.0))),
            H_(Units::from("ft", real(1000.0))),
            id(s) {
  units_["D"] = "nmi";
  units_["H"] = "ft";
}

CDCylinder::CDCylinder(const CDCylinder& cdc) : 
            D_(cdc.D_), 
            H_(cdc.H_), 
            units_(cdc.units_), 
            id(cdc.id) {}

CDCylinder::CDCylinder(real d, real h) :
            D_(std::abs(d)), 
            H_(std::abs(h)),
            id("") {
  units_["D"] = "m";
  units_["H"] = "m";
}

CDCylinder::CDCylinder(real d, const std::string& dunit, real h, const std::string& hunit) :
            D_(Units::from(dunit,std::abs(d))),
            H_(Units::from(hunit,std::abs(h))),
            id("") {
  units_["D"] = dunit;
  units_["H"] = hunit;
}

CDCylinder CDCylinder::make(real distance, const std::string& dUnits, real height, const std::string& hUnits) {
  return CDCylinder(distance,dUnits,height,hUnits);
}

CDCylinder CDCylinder::mk(real distance, real height) {
  return CDCylinder(distance,height);
}

/**
 * @return one static CDCylinder
 */
const CDCylinder& CDCylinder::A_CDCylinder() {
  const static CDCylinder cd3d;
  return cd3d;
}

/**
 * @return CDCylinder thresholds, i.e., D=5nmi, H=1000ft.
 */
const CDCylinder& CDCylinder::CD3DCylinder() {
  return A_CDCylinder();
}

std::string CDCylinder::getUnits(const std::string& key) const {
  std::map<std::string,std::string>::const_iterator got = units_.find(key);
  if (got == units_.end()) {
    return "unspecified";
  }
  return DaidalusParameters::get_overriden_unit(got->second);
}

real CDCylinder::getHorizontalSeparation() const {
  return D_;
}

void CDCylinder::setHorizontalSeparation(real d) {
  D_ = std::abs(d);
}

real CDCylinder::getVerticalSeparation() const {
  return H_;
}

void CDCylinder::setVerticalSeparation(real h) {
  H_ = std::abs(h);
}

real CDCylinder::getHorizontalSeparation(const std::string& u) const {
  return Units::to(u, D_);
}

void CDCylinder::setHorizontalSeparation(real d, const std::string& u) {
  setHorizontalSeparation(Units::from(u,d));
  units_["D"]=u;
}

real CDCylinder::getVerticalSeparation(const std::string& u) const {
  return Units::to(u, H_);
}

void CDCylinder::setVerticalSeparation(real h, const std::string& u) {
  setVerticalSeparation(Units::from(u,h));
  units_["H"]=u;
}

LossData CDCylinder::detection(const Vect3& s, const Vect3& vo, const Vect3& vi, const real D, const real H, const real B, const real T) const {
  return CD3D::detection(s,vo,vi,D,H,B,T);
}

LossData CDCylinder::detection(const Vect3& s, const Vect3& vo, const Vect3& vi, const real D, const real H, const real T) const {
  return detection(s,vo,vi,D,H,real(0.0),T);
}

LossData CDCylinder::detection(const Vect3& s, const Vect3& vo, const Vect3& vi, const real D, const real H) const {
  return detection(s,vo,vi,D,H,real(0.0),PINFINITY);
}

ConflictData CDCylinder::conflict_detection(const Vect3& so, const Vect3& vo, const Vect3& si, const Vect3& vi, real D, real H, real B, real T) {
  Vect3 s = so.Sub(si);
  Vect3 v = vo.Sub(vi);
  real t_tca = CD3D::tccpa(s, vo, vi, D, H, B, T);
  real dist_tca = s.linear(v,t_tca).cyl_norm(D, H);
  LossData ld = CD3D::detection(s,vo,vi,D,H,B,T);
  return ConflictData(ld,t_tca,dist_tca,s,v);
}

ConflictData CDCylinder::conflictDetection(const Vect3& so, const Vect3& vo, const Vect3& si, const Vect3& vi, real B, real T) const {
  return conflict_detection(so,vo,si,vi,D_,H_,B,T);
}

real CDCylinder::time_of_closest_approach(const Vect3& so, const Vect3& vo, const Vect3& si, const Vect3& vi, real D, real H, real B, real T) {
  return CD3D::tccpa(so.Sub(si),vo,vi,D,H,B,T);
}

real CDCylinder::timeOfClosestApproach(const Vect3& so, const Vect3& vo, const Vect3& si, const Vect3& vi, real B, real T) const {
  return time_of_closest_approach(so,vo,si,vi,D_,H_,B,T);
}

CDCylinder* CDCylinder::make() const {
  return new CDCylinder();
}

CDCylinder* CDCylinder::copy() const {
  CDCylinder* cd = new CDCylinder(*this);
  return cd;
}

ParameterData CDCylinder::getParameters() const {
  ParameterData p;
  updateParameterData(p);
  return p;
}

void CDCylinder::updateParameterData(ParameterData& p) const {
  p.setInternal("D",D_,getUnits("D"));
  p.setInternal("H",H_,getUnits("H"));
  p.set("id",id);
}

void CDCylinder::setParameters(const ParameterData& p) {
  if (p.contains("D")) {
    setHorizontalSeparation(p.getValue("D"));
    units_["D"]=p.getUnit("D");
  }
  if (p.contains("H")) {
    setVerticalSeparation(p.getValue("H"));
    units_["H"]=p.getUnit("H");
  }
  if (p.contains("id")) {
    id = p.getString("id");
  }
}

std::string CDCylinder::getSimpleClassName() const {
  return "CDCylinder";
}

std::string CDCylinder::toString() const {
  return (id == "" ? "" : id+" : ")+getSimpleClassName()+" = {D = "+
      Units::str(getUnits("D"),D_)+", H = "+Units::str(getUnits("H"),H_)+"}";
}

std::string CDCylinder::toPVS() const {
  return getSimpleClassName()+"((# D:= "+FmPrecision(D_)+", H:= "+FmPrecision(H_)+" #))";
}

std::string CDCylinder::getIdentifier() const {
  return id;
}

void CDCylinder::setIdentifier(const std::string& s) {
  id = s;
}

bool CDCylinder::equals(const Detection3D& det) const {
  if (!larcfm::equals(getCanonicalClassName(), det.getCanonicalClassName())) return false;
  if (!larcfm::equals(id,det.getIdentifier())) return false;
  const CDCylinder& cd = (CDCylinder&)det;
  if (D_ != cd.D_) return false;
  if (H_ != cd.H_) return false;
  return true;
}

bool CDCylinder::contains(const Detection3D& det) const {
  if (larcfm::equals(getCanonicalClassName(), det.getCanonicalClassName())) {
    const CDCylinder& cd = (CDCylinder&)det;
    return D_ >= cd.D_ && H_ >= cd.H_;
  }
  return false;
}

/* Return a list of point representing a counter-clockwise circular arc centered at pc,
 * whose first point is pc+v(0), and the last one is pc+v(alpha), where alpha is
 * in [0,2*pi].
 */
void CDCylinder::circular_arc(std::vector<Position>& arc, const Position& pc, const Velocity& v,
    real alpha, bool include_last) {
  alpha = Util::almost_equals(alpha,2*Pi,DaidalusParameters::ALMOST_) ? alpha :  Util::to_2pi(alpha);
  real step = Pi/180;
  arc.push_back(pc.linear(v,1));
  real current_trk = v.trk();
  for (real a = step; Util::almost_less(a,alpha,DaidalusParameters::ALMOST_); a += step) {
    arc.push_back(pc.linear(v.mkTrk(current_trk-a),1));
  }
  if (include_last) {
    arc.push_back(pc.linear(v.mkTrk(current_trk-alpha),1));
  }
}

void CDCylinder::horizontalHazardZone(std::vector<Position>& haz,
    const TrafficState& ownship, const TrafficState& intruder, real T) const {
  haz.clear();
  const Position& po = ownship.getPosition();
  Velocity v = ownship.getAirVelocity().Sub(intruder.getAirVelocity().vect3());
  if (Util::almost_equals(T,0) || Util::almost_equals(v.vect3().norm2D(),0)) {
    circular_arc(haz,po,Velocity::mkVxyz(D_,0,0),2*Pi,false);
  } else {
    Vect3 sD = Horizontal::unit_perpL(v.vect3()).Scal(D_);
    Velocity vD = Velocity::make(sD);
    Velocity vnD = Velocity::make(sD.Neg());
    circular_arc(haz,po,vD,Pi,true);
    circular_arc(haz,po.linear(v,T),vnD,Pi,true);
  }
}

}

