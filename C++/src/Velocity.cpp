/*
 * Velocity.cpp
 *
 * Contact: Cesar Munoz (cesar.a.munoz@nasa.gov)
 * NASA LaRC
 * http://research.nianet.org/fm-at-nia/ACCoRD
 *
 * NOTES: 
 * Track is True North/clockwise
 *
 * Copyright (c) 2011-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */

#include "Real.h"
#include "Velocity.h"
#include "Units.h"
#include "Util.h" 
#include "format.h"
#include "Constants.h"
#include "string_util.h"
#include <limits>
#include <vector>

namespace larcfm {

Velocity::Velocity(const real trk, const real gs,
		const real vx, const real vy, const real vz) : trk_(trk),
		gs_(gs),v_(vx,vy,vz) {}

Velocity::Velocity(const real vx, const real vy, const real vz) :
	trk_(Util::atan2_safe(vx,vy)),
	gs_(Util::sqrt_safe(vx*vx+vy*vy)),
	v_(vx,vy,vz) {}

Velocity::Velocity() : trk_(real(0.0)),gs_(real(0.0)),v_(Vect3::ZERO()) {}

bool Velocity::isZero() const {
	return v_.isZero();
}

bool Velocity::isInvalid() const {
	return v_.isInvalid();
}

Velocity Velocity::make(const Vect3& v) {
	return Velocity(v.x(),v.y(),v.z());
}

Velocity Velocity::make(const Vect2& v) {
	return Velocity(v.x,v.y,real(0.0));
}

Velocity Velocity::mkVxyz(const real vx, const real vy, const real vz) {
	return Velocity(vx,vy,vz);
}

Velocity Velocity::makeVxyz(const real vx, const real vy, const std::string& uvxy,
		const real vz, const std::string& uvz) {
	return Velocity(Units::from(uvxy,vx),Units::from(uvxy,vy),Units::from(uvz,vz));
}

Velocity Velocity::mkTrkGsVs(const real trk, const real gs, const real vs){
	return Velocity(trk,gs,trkgs2vx(trk,gs),trkgs2vy(trk,gs),vs);
}

Velocity Velocity::makeTrkGsVs(const real trk, const std::string& utrk,
		const real gs, const std::string& ugs,
		const real vs, const std::string& uvs) {
	return mkTrkGsVs(Units::from(utrk,trk), Units::from(ugs,gs),Units::from(uvs,vs));
}

Velocity Velocity::mkVel(const Vect3& p1,const Vect3& p2, real speed) {
	return make(p2.Sub(p1).Hat().Scal(speed));
}

real Velocity::track(const Vect3& p1,const Vect3& p2) {
	return Util::atan2_safe(p2.x()-p1.x(),p2.y()-p1.y());
}

Velocity Velocity::Neg() const {
	return Velocity(to_pi(trk_+real(M_PI)),gs_,-v_.x(),-v_.y(),-v_.z());
}

Velocity Velocity::Add(const Vect3& v) const {
    if (Util::almost_equals(v_.x(),-v.x()) && Util::almost_equals(v_.y(),-v.y())) {
        // Set to zero but maintain the original track
		return Velocity(trk_,real(0.0),real(0.0),real(0.0),v_.z()+v.z());
    }
    return mkVxyz(v_.x()+v.x(),v_.y()+v.y(),v_.z()+v.z());
}

Velocity Velocity::Sub(const Vect3& v) const {
    if (Util::almost_equals(v_.x(),v.x()) && Util::almost_equals(v_.y(),v.y())) {
        // Set to zero but maintain the original track
		return Velocity(trk_,real(0.0),real(0.0),real(0.0),v_.z()-v.z());
    }
    return mkVxyz(v_.x()-v.x(),v_.y()-v.y(),v_.z()-v.z());
}

/**
 * Make a unit 2D vector from the velocity vector. 
 * @return the unit 2D vector
 */
Vect2 Velocity::Hat2D() const {
	return Vect2(std::sin(trk_),std::cos(trk_));
}

Velocity Velocity::genVel(const Vect3& p1, const Vect3& p2, real dt) {
	return make(p2.Sub(p1).Scal(1/dt));
}

Velocity Velocity::mkAddTrk(real atrk) const {
	real s = std::sin(atrk);
	real c = std::cos(atrk);
    return Velocity(to_pi(trk_+atrk),gs_,v_.x()*c+v_.y()*s,-v_.x()*s+v_.y()*c,v_.z());
}

real Velocity::trkgs2vx(real trk, real gs) {
	return gs * std::sin(trk);
}

real Velocity::trkgs2vy(real trk, real gs) {
	return gs * std::cos(trk);
}

Vect2 Velocity::trkgs2v(real trk, real gs) {
	return Vect2(trkgs2vx(trk,gs), trkgs2vy(trk,gs));
}

real Velocity::angle() const {
	return to_pi(real(M_PI)/real(2.0)-trk_);
}

real Velocity::angle(const std::string& uangle) const {
	return Units::to(uangle,angle());
}

real Velocity::track(const std::string& utrk) const {
	return Units::to(utrk,trk_);
}

/**
 * Compass angle in radians in the range [<code>0</code>, <code>2*Math.PI</code>).
 * Convention is clockwise with respect to north.
 *
 * @return the compass angle [rad]
 */
real Velocity::compassAngle() const {
	return to_2pi(trk_);
}

/**
 * Compass angle in explicit units in corresponding range [<code>0</code>, <code>2*Math.PI</code>).
 * Convention is clockwise with respect to north.
 *
 *  @param u the explicit units of compass angle
 *
 *  @return the compass angle [u]
 */
real Velocity::compassAngle(const std::string& u) const {
	return Units::to(u,compassAngle());
}

real Velocity::groundSpeed(const std::string& ugs) const {
	return Units::to(ugs,gs_);
}

real Velocity::vs() const {
	return v_.z();
}

real Velocity::verticalSpeed(const std::string& uvs) const {
	return Units::to(uvs,v_.z());
}

bool Velocity::compare(const Velocity& v, real maxTrk, real maxGs, real maxVs) const {
	if (Util::turnDelta(v.trk(),trk()) > maxTrk) return false;
	if (std::abs(v.gs() - gs()) > maxGs) return false;
	if (std::abs(v.vs() - vs()) > maxVs) return false;
	return true;
}

bool Velocity::compare(const Velocity& v, real horizDelta, real vertDelta) const {
	return std::abs(v_.z()-v.z()) <= vertDelta && vect2().Sub(v.vect2()).norm() <= horizDelta;
}

const Velocity& Velocity::ZERO() {
	const static Velocity v;
	return v;
}

const Velocity& Velocity::INVALID() {
	const static Velocity v(NaN,NaN,NaN,NaN,NaN);
	return v;
}

/**
 * New velocity from existing velocity, changing only the track
 * @param trk track angle [rad]
 * @return new velocity
 */
Velocity Velocity::mkTrk(real trk) const {
	return mkTrkGsVs(trk, gs_,v_.z());
}

/**
 * New velocity from existing velocity, changing only the track
 * @param trk track angle [u]
 * @param u units
 * @return new velocity
 */
Velocity Velocity::mkTrk(real trk, const std::string& u) const {
	return mkTrk(Units::from(u,trk));
}

/**
 * New velocity from existing velocity, changing only the ground speed
 * @param gs [m/s]
 * @return
 */
Velocity Velocity::mkGs(real ags) const {
	if (ags < 0) {
		return INVALID();
	}
        if (gs_ > real(0.0)) {
        real scal = ags/gs_;
    	return Velocity(trk_,ags,v_.x()*scal,v_.y()*scal,v_.z());
    }
    return mkTrkGsVs(trk_,ags,v_.z());  
}

/**
 * New velocity from existing velocity, changing only the ground speed
 * @param ags [u]
 * @param u unit
 * @return
 */
Velocity Velocity::mkGs(real ags, const std::string& u) const {
	return mkGs(Units::from(u,ags));
}

/**
}
 * New velocity from existing velocity, changing only the vertical speed
 * @param vs [m/s]
 * @return
 */
Velocity Velocity::mkVs(real vs) const {
	return Velocity(trk_,gs_,v_.x(), v_.y(),vs);
}

/**
 * New velocity from existing velocity, changing only the vertical speed
 * @param vs [u]
 * @param u units
 * @return
 */
Velocity Velocity::mkVs(real vs, const std::string& u) const {
	return mkVs(Units::from(u,vs));
}

Velocity Velocity::zeroSmallVs(real threshold) const {
	if (std::abs(v_.z()) < threshold) {
		return mkVs(real(0.0));
	}
	return Velocity(trk_,gs_,v_.x(),v_.y(),v_.z());
}

Velocity Velocity::parseXYZ(const std::string& str) {
	return Velocity::make(Vect3::parse(str));
}

std::string Velocity::toString() const {
	return toString(Constants::get_output_precision());
}

std::string Velocity::toString(int prec) const {
	 return"("+Units::str("deg",compassAngle(),prec)+", "+Units::str("knot",gs_,prec)+", "+Units::str("fpm",v_.z(),prec)+")";
}

std::string Velocity::toStringUnits() const {
	return toStringUnits("deg","knot","fpm");
}

std::string Velocity::toStringUnits(const std::string& trkUnits, const std::string& gsUnits, const std::string& vsUnits) const {
	return  "("+Units::str(trkUnits,compassAngle())+", "+ Units::str(gsUnits,gs_)+", "+ Units::str(vsUnits,v_.z())+")";
}

std::string Velocity::toStringXYZ() const {
	return toStringXYZ(Constants::get_output_precision());
}

std::string Velocity::toStringXYZ(int prec) const {
	return "("+FmPrecision(Units::to("knot", v_.x()),prec)+", "+FmPrecision(Units::to("knot", v_.y()),prec)+", "+FmPrecision(Units::to("fpm", v_.z()),prec)+")";
}

std::vector<std::string> Velocity::toStringList() const {
	std::vector<std::string> ret;
	if (isInvalid()) {
		ret.push_back("-");
		ret.push_back("-");
		ret.push_back("-");
	} else {
		ret.push_back(Fm12(Units::to("deg", compassAngle())));
		ret.push_back(Fm12(Units::to("knot", gs_)));
		ret.push_back(Fm12(Units::to("fpm", v_.z())));
	}
	return ret;
}

std::vector<std::string> Velocity::toStringList(int precision) const {
	std::vector<std::string> ret;
	if (isInvalid()) {
		ret.push_back("-");
		ret.push_back("-");
		ret.push_back("-");
	} else {
		ret.push_back(FmPrecision(Units::to("deg", compassAngle()),precision));
		ret.push_back(FmPrecision(Units::to("knot", gs_),precision));
		ret.push_back(FmPrecision(Units::to("fpm", v_.z()),precision));
	}
	return ret;
}

std::vector<std::string> Velocity::toStringXYZList() const {
	std::vector<std::string> ret;
	if (isInvalid()) {
		ret.push_back("-");
		ret.push_back("-");
		ret.push_back("-");
	} else {
		ret.push_back(Fm12(Units::to("knot", v_.x())));
		ret.push_back(Fm12(Units::to("knot", v_.y())));
		ret.push_back(Fm12(Units::to("fpm", v_.z())));
	}
	return ret;
}

std::vector<std::string> Velocity::toStringXYZList(int precision) const {
	std::vector<std::string> ret;
	if (isInvalid()) {
		ret.push_back("-");
		ret.push_back("-");
		ret.push_back("-");
	} else {
		ret.push_back(FmPrecision(Units::to("knot", v_.x()),precision));
		ret.push_back(FmPrecision(Units::to("knot", v_.y()),precision));
		ret.push_back(FmPrecision(Units::to("fpm", v_.z()),precision));
	}
	return ret;
}

std::string Velocity::toStringNP() const {
	return toStringNP(Constants::get_output_precision());
}

std::string Velocity::toStringNP(int precision) const {
	return FmPrecision(Units::to("deg", compassAngle()), precision)+", "+FmPrecision(Units::to("knot", gs_),precision)+", "+FmPrecision(Units::to("fpm", v_.z()),precision);
}

std::string Velocity::toStringNP(const std::string& utrk, const std::string& ugs, const std::string& uvs, int precision) const {
    return FmPrecision(Units::to(utrk, compassAngle()), precision)+", "+FmPrecision(Units::to(ugs,gs_),precision)+", "+FmPrecision(Units::to(uvs, v_.z()), precision);
}

std::string Velocity::toStringNP(const std::string& utrk, const std::string& ugs, const std::string& uvs) const {
    return toStringNP(utrk,ugs,uvs,Constants::get_output_precision());
}

}
