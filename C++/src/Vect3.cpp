/*
 * Vect3.cpp
 *
 * Contact: Cesar Munoz (cesar.a.munoz@nasa.gov)
 * NASA LaRC
 *
 * 3-D vectors.
 *
 * Copyright (c) 2011-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */

#include "Real.h"
#include "Vect3.h"
#include "Units.h"
#include <limits>
#include "string_util.h"
#include "format.h"
#include "Constants.h"
#include "Util.h" // NaN def

namespace larcfm {

Vect3::Vect3(real xx, real yy, real zz) : v2_(xx,yy), z_(zz) {}

Vect3::Vect3(const Vect2&v, real vz) : v2_(v),z_(vz) {}

Vect3 Vect3::makeXYZ(real x, const std::string& ux, real y, const std::string& uy, real z, const std::string& uz) {
	return Vect3(Units::from(ux,x),Units::from(uy,y),Units::from(uz,z));
}

Vect3 Vect3::mkXYZ(real x, real y, real z) {
	return Vect3(x,y,z);
}

Vect3 Vect3::mkX(real nx) {
	return Vect3(nx,y(),z());
}
Vect3 Vect3::mkY(real ny) {
	return Vect3(x(),ny,z());
}
Vect3 Vect3::mkZ(real nz) {
	return Vect3(x(),y(),nz);
}

bool Vect3::isZero() const {
	return x() == real(0.0) && y() == real(0.0) && z() == real(0.0);
}

bool Vect3::almostEquals(const Vect3& v) const {
	return Util::almost_equals(x(),v.x()) && Util::almost_equals(y(),v.y()) && Util::almost_equals(z(),v.z());
}

bool Vect3::almostEquals(const Vect3& v, INT64FM maxUlps) const {
	return Util::almost_equals(x(),v.x(),maxUlps) && Util::almost_equals(y(),v.y(),maxUlps) && Util::almost_equals(z(),v.z(),maxUlps);
}

bool Vect3::almostEquals2D(const Vect3& v, real horizEps) const {
	return (vect2().Sub(v.vect2())).norm() < horizEps;
}

bool Vect3::within_epsilon(const Vect3& v2, real epsilon) const {
	if (std::abs(x() - v2.x()) > epsilon) return false;
	if (std::abs(y() - v2.y()) > epsilon) return false;
	if (std::abs(z() - v2.z()) > epsilon) return false;
	return true;
}


Vect3 Vect3::operator + (const Vect3& v) const {
	return this->Add(v);
}

Vect3 Vect3::operator - (const Vect3& v) const {
	return this->Sub(v);
}

Vect3 Vect3::operator - () const {
	return Vect3(-x(),-y(),-z());
}

Vect3 Vect3::operator * (const real k) const {
	return this->Scal(k);
}

real Vect3::operator * (const Vect3& v) const { // Dot product
	return dot(v.x(),v.y(),v.z());
}

bool Vect3::operator == (const Vect3& v) const {  // strict equality
	return x() == v.x() && y() == v.y() && z() == v.z();
}

bool Vect3::operator != (const Vect3& v) const {  // strict disequality
	return x() !=v.x() || y() !=v.y() || z() != v.z();
}

Vect3 Vect3::Hat() const {
	real n = norm();
	if ( n == real(0.0) ) { // this is only checking the divide by zero case, so an exact comparison is correct.
		return ZERO();
	}
	return Vect3(x()/n, y()/n, z()/n);
}

Vect3 Vect3::cross(const Vect3& v) const {
	return Vect3(y()*v.z() - z()*v.y(), z()*v.x() - x()*v.z(), x()*v.y() - y()*v.x());
}

bool Vect3:: parallel(const Vect3& v) const {
	return cross(v).almostEquals(ZERO());
}

Vect3 Vect3::Add(const Vect3& v) const{
	return Vect3(x()+v.x(), y()+v.y(), z()+v.z());
}

Vect3 Vect3::Sub(const Vect3& v) const {
	return Vect3(x()-v.x(),y()-v.y(),z()-v.z());
}

Vect3 Vect3::Neg() const {
	return Vect3(-x(),-y(),-z());
}

Vect3 Vect3::Scal(real k) const {
	return Vect3(k*x(),k*y(),k*z());
}

Vect3 Vect3::ScalAdd(const real k, const Vect3& v) const {
	return Vect3(k*x()+v.x(), k*y()+v.y(), k*z()+v.z());
}

Vect3 Vect3::AddScal(real k, const Vect3& v) const {
	return Vect3(x()+k*v.x(), y()+k*v.y(), z()+k*v.z());
}

Vect3 Vect3::PerpR() const {
	return Vect3(y(),-x(), 0);
}

Vect3 Vect3::PerpL() const {
	return Vect3(-y(),x(), 0);
}

Vect3 Vect3::linear(const Vect3& v, real t) const {
	return Vect3(x()+v.x()*t, y()+v.y()*t, z()+v.z()*t);
}

Vect3 Vect3::linearByDist2D(real track, real d) const {
	return Vect3(x()+d*std::sin(track), y()+d*std::cos(track), z());
}


real Vect3::dot(const real xx, const real yy, const real zz) const {
	return x()*xx + y()*yy + z()*zz;
}

real Vect3::dot(const Vect3& v) const {
	return dot(v.x(), v.y(), v.z());
}

real Vect3::sqv() const {
	return dot(x(),y(),z());
}

real Vect3::norm() const {
	return std::sqrt(sqv());
}

real Vect3::cyl_norm(const real d, const real h) const {
	return Util::max(vect2().sqv()/sq(d),sq(z()/h));
}

real Vect3::distanceH(const Vect3& w) const {
	return (vect2()-w.vect2()).norm();
}
real Vect3::distanceV(const Vect3& w) const {
	return z() - w.z();
}

std::string Vect3::toString() const {
	return toString(Constants::get_output_precision());
}

std::string Vect3::toString(int precision) const {
	return formatXYZ(precision,"(",", ",")");
}

std::string Vect3::toStringNP(const std::string& xunit, const std::string& yunit, const std::string& zunit) const {
	return toStringNP(xunit,yunit,zunit,Constants::get_output_precision());
}

std::string Vect3::toStringNP(const std::string& xunit, const std::string& yunit, const std::string& zunit, int prec) const {
	return  FmPrecision(Units::to(xunit,x()),prec)+", "
           +FmPrecision(Units::to(yunit,y()),prec)+", " 
		   +FmPrecision(Units::to(zunit,z()),prec);
}

std::string Vect3::formatXYZ(const std::string& pre, const std::string& mid, const std::string& post) const {
	return formatXYZ(Constants::get_output_precision(),pre,mid,post);
}

std::string Vect3::formatXYZ(int prec, const std::string& pre, const std::string& mid, const std::string& post) const {
	return pre+FmPrecision(x(),prec)+mid+FmPrecision(y(),prec)+mid+FmPrecision(z(),prec)+post;
}

std::string Vect3::toPVS() const {
	return toPVS(Constants::get_output_precision());
}

std::string Vect3::toPVS(int precision) const {
	return "(# x:= "+FmPrecision(x(),precision)+", y:= "+FmPrecision(y(),precision)+", z:= "+FmPrecision(z(),precision)+" #)";
}

/** 3-D time of closest point of approach
 * if time is negative or velocities are parallel returns 0
 * @param so position of one
 * @param vo velocity of one
 * @param si position of two
 * @param vi velocity of two
 * @return time of closest point of approach
 */
real Vect3::tcpa(const Vect3& so, const Vect3& vo, const Vect3& si, const Vect3& vi) {
	real t;
	Vect3 s = so.Sub(si);
	Vect3 v = vo.Sub(vi);
	real nv = v.sqv();
	if (nv > 0)
		t = Util::max(real(0.0),-s.dot(v)/nv);
	else
		t = 0;
	return t;
}

/**
 * Returns true if the current vector has an "invalid" value
 */
bool Vect3::isInvalid() const {
	return ISNAN(x()) || ISNAN(y()) || ISNAN(z()); 
}

const Vect3& Vect3::ZERO() {
	const static Vect3 v(real(0.0),real(0.0),real(0.0));
	return v;
}

const Vect3& Vect3::INVALID() {
	const static Vect3 v(NaN, NaN, NaN);
	return v;
}

Vect3 Vect3::parse(const std::string& s) {
	std::vector<std::string> fields = split(s, Constants::wsPatternParens);
	while (fields.size() > 0 && equals(fields[0], "")) {
		fields.erase(fields.begin());
	}
	if (fields.size() == 3) {
		return Vect3(Util::parse_double(fields[0]),Util::parse_double(fields[1]),Util::parse_double(fields[2]));
	} else if (fields.size() == 6) {
		return Vect3(Units::from(Units::clean(fields[1]),Util::parse_double(fields[0])),
				Units::from(Units::clean(fields[3]),Util::parse_double(fields[2])),
				Units::from(Units::clean(fields[5]),Util::parse_double(fields[4])));
	}
	return Vect3::INVALID();
}

/*
 * Two dimensional calculations on Vect3s.  z components will be ignored or set to zero.
 */

real Vect3::det2D(const Vect3& v) const {
	return x()*v.y() - y()*v.x();
}

real Vect3::dot2D(const Vect3& v) const {
	return x()*v.x() + y()*v.y();
}

real Vect3::sqv2D() const {
	return x()*x()+y()*y();
}

real Vect3::norm2D() const {
	return std::sqrt(sqv2D());
}

Vect3 Vect3::Hat2D() const {
	real n = norm2D();
	if (n == real(0.0)) { // this is only checking the divide by zero case, so an exact comparison is correct.
		return ZERO();
	}
	return Vect3(x()/n, y()/n, real(0.0));
}

}
