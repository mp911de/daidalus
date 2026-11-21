/*
 * Vect4.cpp
 *
 * 4-D vectors.
 *
 * Copyright (c) 2011-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */

#include "Real.h"
#include "Vect4.h"

namespace larcfm {

  Vect4::Vect4(const real xx, const real yy, const real zz, const real tt) : x(xx), y(yy), z(zz), t(tt) {}

  Vect4::Vect4(const Vect3& v, real vt) : x(v.x()), y(v.y()), z(v.z()), t(vt) {}

  bool Vect4::isZero() const {
    return x == real(0.0) && y == real(0.0) && z == real(0.0) && t == real(0.0);
  }

  bool Vect4::almostEquals(const Vect4& v) const {
    return Util::almost_equals(x,v.x) && Util::almost_equals(y,v.y) && Util::almost_equals(z,v.z) && Util::almost_equals(t,v.t);
  }

  bool Vect4::almostEquals(const Vect4& v, INT64FM maxUlps) const {
    return Util::almost_equals(x,v.x,maxUlps) && Util::almost_equals(y,v.y,maxUlps) && Util::almost_equals(z,v.z,maxUlps) && Util::almost_equals(t,v.t,maxUlps);
  }
  
  Vect4 Vect4::operator + (const Vect4& v) const {
	  return this->Add(v);
  }

  Vect4 Vect4::operator - (const Vect4& v) const {
	  return this->Sub(v);
  }

  Vect4 Vect4::operator - () const {
    return Vect4(-x,-y,-z,-t);
  }

  Vect4 Vect4::operator * (const real k) const {
	  return this->Scal(k);
  }

  real Vect4::operator * (const Vect4& v) const { // Dot product
    return dot(v.x,v.y,v.z,v.t);
  }

  bool Vect4::operator == (const Vect4& v) const {  // strict equality
    return x==v.x && y==v.y && z==v.z && t==v.t;
  }

  bool Vect4::operator != (const Vect4& v) const {  // strict disequality
    return x!=v.x || y!=v.y || z!=v.z || t!=v.t;
  }

  Vect2 Vect4::vect2() const {
    return Vect2(x,y);
  }


  Vect3 Vect4::vect3() const {
    return Vect3(x,y,z);
  }

  Vect4 Vect4::Add(const Vect4& v) const {
	  return Vect4(x+v.x,y+v.y,z+v.z,t+v.t);
  }

  Vect4 Vect4::Sub(const Vect4& v) const {
	  return Vect4(x-v.x,y-v.y,z-v.z,t-v.t);
  }

  Vect4 Vect4::Neg() const {
      return Vect4(-x,-y,-z,-t);
  }

  Vect4 Vect4::Scal(real k) const {
	  return Vect4(k*x,k*y,k*z,k*t);
  }

  Vect4 Vect4::ScalAdd(const real k, const Vect4& v) const {
	  return Vect4(k*x+v.x,k*y+v.y,k*z+v.z,k*t+v.t);
  }

  real Vect4::dot(const real x, const real y, const real z, const real t) const {
    return this->x*x + this->y*y + this->z*z + this->t*t;
  }

  real Vect4::sqv() const {
    return dot(x,y,z,t);
  }

  real Vect4::norm() const {
    return sqrt_safe(sqv());
  }
    

  real Vect4::lat() const {
    return y;
  }
  
  
  real Vect4::lon() const {
    return x;
  }
  
  
  std::string Vect4::toString() const {
    std::stringstream temp;
    temp << "(" << x << ", " << y << ", " << z << "," << t << ")";
    return temp.str();
}

  const Vect4& Vect4::ZERO() {
		const static Vect4 v(real(0.0),real(0.0),real(0.0),real(0.0));
		return v;
  }


}
