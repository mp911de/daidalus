/*
 * SimpleProjection.cpp
 *
 * Contact: Jeff Maddalon (j.m.maddalon@nasa.gov)
 * NASA LaRC
 * 
 * Copyright (c) 2011-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */

#include "Real.h"
#include "SimpleProjection.h"
#include "LatLonAlt.h"
#include "GreatCircle.h"
#include "Util.h"
#include "Constants.h"
#include "format.h"
#include <stdexcept>
#include <limits>


namespace larcfm {
    using std::string;
    using std::cout;
    using std::endl;
    using std::runtime_error;
  
    static const real tranLat = Units::from("deg", real(85.0));
    
    
    SimpleProjection::SimpleProjection() :
      projLat(real(0.0)),
      projLon(real(0.0)),
      projAlt(real(0.0)),
      projNorth(false) {}
    
    SimpleProjection::SimpleProjection(const LatLonAlt& lla) :
        projLat(lla.lat()),
        projLon(lla.lon()),
        projAlt(lla.alt()),
        projNorth(projLat >= real(0.0)) {}
    
    SimpleProjection::SimpleProjection(real lat, real lon, real alt) :
        projLat(lat),
        projLon(lon),
        projAlt(alt),
        projNorth(projLat >= real(0.0)) {}
    
    SimpleProjection SimpleProjection::makeNew(const LatLonAlt& lla) const {
      return SimpleProjection(lla);
    }

    SimpleProjection SimpleProjection::makeNew(real lat, real lon, real alt) const {
      return SimpleProjection(lat, lon, alt);
    }

//    SimpleProjection SimpleProjection::makeNewCopy() const {
//      return SimpleProjection(projLat, projLon, projAlt);
//    }
  
    real SimpleProjection::conflictRange(real lat, real accuracy) const {
      real degs = Units::to("deg",lat);
      if (accuracy < Units::from("nm", real(0.1))) { //0.01 -- 35 nm - 3 nm
        if (degs < 30) {
          return Units::from("nm",10);
        } else {
          return Units::from("nm",5);
        }
      } else if (accuracy < Units::from("nm", real(0.5))) {  // 0.1 nm -- 185 nm - 4 nm
        if (degs < 20) {
          return Units::from("nm", 40);
        } else if (degs < 50) {
          return Units::from("nm", 25);
        } else if (degs < 70) {
          return Units::from("nm", 15);
        } else if (degs < 80) {
          return Units::from("nm", 10);
        } else {
          return Units::from("nm", 5);
        }    		
      } else { // 0.5 nm -- 330 nm - 13 nm
        if (degs < 20) {
          return Units::from("nm", 95);
        } else if (degs < 50) {
          return Units::from("nm", 50);
        } else if (degs < 70) {
          return Units::from("nm", 20);
        } else if (degs < 80) {
          return Units::from("nm", 10);
        } else {
          return Units::from("nm", 5);
        }
      }
    }
    
    real SimpleProjection::maxRange() const{
      return std::numeric_limits<real>::max();
    }
    
    LatLonAlt SimpleProjection::getProjectionPoint() const {
      return LatLonAlt::mk(projLat, projLon, projAlt);
    }
    
    bool SimpleProjection::isPolar() const {
      return (std::abs(projLat) > tranLat);
    }
    
    Vect2 SimpleProjection::polar_xy(const LatLonAlt& lla, bool north) {
       	int sgn = 1;
    	if (! north) sgn = -1;
        real a = sgn*Units::from("deg",real(90.0)) - lla.lat();
        real r = std::abs(GreatCircle::distance_from_angle(a,real(0.0)));
    	//f.pln("^^ polarXY: a = "+Units::to("deg",a)+" r = "+Units::to("nm",r));
        return Vect2(r*std::sin(lla.lon()),r*std::cos(lla.lon()));
    }
    
    LatLonAlt SimpleProjection::polar_inverse(const Vect2& v, real alt, bool north) {
        real lon = to_pi(v.compassAngle());
        real d = v.norm();
        real a = GreatCircle::angle_from_distance(d,real(0.0));
    	//f.pln("^^ polarLL: a = "+Units::to("deg",a)+" d = "+Units::to("nm",d));
    	int sgn = 1;
    	if (! north) sgn = -1;
        real lat = sgn*(Units::from("deg",real(90.0)) - a);
    	//f.pln("^^ polarLL: lat = "+Units::to("deg",lat)+" lon = "+Units::to("deg",lon));
    	return LatLonAlt::mk(lat,lon,alt);
    }

    Vect2 SimpleProjection::project2(const LatLonAlt& lla) const {
      if (std::abs(projLat) >  tranLat)
    	   return polar_xy(lla, projNorth);
    	else
           return SimpleProjection::projectXY(projLat,projLon,lla.lat(),lla.lon());
           
    }

    Vect3 SimpleProjection::project(const LatLonAlt& lla) const {
       return Vect3(project2(lla),lla.alt() - projAlt);
    }      
    
    Vect3 SimpleProjection::project(const Position& sip) const {
    	Vect3 si;
    	if (sip.isLatLon()) {
    		si = project(sip.lla());
    	} else {
    		si = sip.vect3();
    	}
    	return si;
    }

    Vect3 SimpleProjection::projectPoint(const Position& sip) const {
    	return project(sip);
    }

    LatLonAlt SimpleProjection::inverse(const Vect2& xy, real alt) const {
      if (std::abs(projLat) >  tranLat)
    		return polar_inverse(xy, alt + projAlt, projNorth);
    	else {
           real a = xy.compassAngle();
           real d = xy.norm();
    	   LatLonAlt origin = LatLonAlt::mk(projLat,projLon,alt + projAlt);
           return GreatCircle::linear_initial(origin,a,d);
    	}
    }

    LatLonAlt SimpleProjection::inverse(const Vect3& xyz) const {  
    	return inverse(xyz.vect2(), xyz.z());
    }
    
  Velocity SimpleProjection::projectVelocity(const LatLonAlt& lla, const Velocity& v) const {
    real timeStep = real(10.0);
    LatLonAlt ll2 = GreatCircle::linear_initial(lla,v,timeStep);
    Vect3 se = project(lla);
    Vect3 s2 = project(ll2);
    Vect3 vn = s2.Sub(se).Scal(1/timeStep);
    return Velocity::make(vn);
  }
  
  Velocity SimpleProjection::projectVelocity(const Position& ss, const Velocity& v) const {
    Velocity vn;
    if (ss.isLatLon()) {
      vn = projectVelocity(ss.lla(),v);
    } else {
      vn = v;
    }
    return vn;
  }

  // The user needs to keep track of whether to translate back (i.e. whether original was LatLon())
  Velocity SimpleProjection::inverseVelocity(const Vect3& s, const Velocity& v, bool toLatLon) const {
    if (toLatLon) {
      real timeStep = real(10.0);
      Vect3 s2 = s.linear(v.vect3(),timeStep);
      LatLonAlt lla1 = inverse(s);
      LatLonAlt lla2 = inverse(s2);
      Velocity nv = GreatCircle::velocity_initial(lla1,lla2,timeStep);
      return nv;
    } else {
      return v;
    }
  }

  std::pair<Vect3,Velocity> SimpleProjection::project(const Position& p, const Velocity& v) const {
   	return std::pair<Vect3,Velocity>(project(p),projectVelocity(p,v));
   }

  std::pair<Vect3,Velocity> SimpleProjection::project(const LatLonAlt& p, const Velocity& v) const {
   	return std::pair<Vect3,Velocity>(project(p),projectVelocity(p,v));
   }

  std::pair<Position,Velocity> SimpleProjection::inverse(const Vect3& p, const Velocity& v, bool toLatLon) const {
    if (toLatLon) {
      return std::pair<Position,Velocity>(Position(inverse(p)),inverseVelocity(p,v,true));
    } else {
      return std::pair<Position,Velocity>(Position(p),v);
    }
   }


  Vect2 SimpleProjection::projectXY(real lat0, real lon0, real lat1, real lon1) {
    real ad = GreatCircle::angular_distance(lat0, lon0, lat1, lon1);
//std::cout << "SimpleProjection::projectXY ad=" << Fm16(ad) << std::endl;;
    if (Constants::almost_equals_radian(ad)) {
      // If the two points are about 1 meter apart, then count them as the
      // same.
      return Vect2(real(0.0), real(0.0));
    }

    LatLonAlt p0 = LatLonAlt::mk(lat0, lon0, real(0.0));
    LatLonAlt p1 = LatLonAlt::mk(lat1, lon1, real(0.0));

    real a;
    if (Constants::almost_equals_radian(ad / real(2.0))) {
      // The two points are about 2 meters apart, so skip the midpoint
      a = GreatCircle::initial_course(p0, p1);
//std::cout << "SimpleProjection::projectXY a(1)=" << Fm16(a) << std::endl;;
    } else {
      LatLonAlt mid = GreatCircle::interpolate(p0, p1, real(0.5));
      a = GreatCircle::initial_course(mid, p1);
//std::cout << "SimpleProjection::projectXY mid=" << mid.toString15() << std::endl;;
//std::cout << "SimpleProjection::projectXY p1 =" << p1.toString15() << std::endl;;
//std::cout << "SimpleProjection::projectXY a(2)=" << Fm16(a) << std::endl;;
    }
    real d = GreatCircle::distance_from_angle(ad,real(0.0));
//std::cout << "SimpleProjection::projectXY d=" << d << std::endl;;
//std::cout << "SimpleProjection::projectXY sin(a)=" << sin(a) << std::endl;;
//std::cout << "SimpleProjection::projectXY cos(a)=" << cos(a) << std::endl;;
    return Vect2(d * std::sin(a), d * std::cos(a));
  }

  Vect3 SimpleProjection::projectXYZ(const LatLonAlt& p0, const LatLonAlt& p1) {
    Vect2 r = projectXY(p0.lat(), p0.lon(), p1.lat(), p1.lon());
    return Vect3(r.x, r.y, p1.alt() - p0.alt());
  }


  //**************************************
  // deprecated functions:

  Vect2 polar_xy(const LatLonAlt& lla, bool north) {
  	return SimpleProjection::polar_xy(lla, north);
  }

  LatLonAlt polar_inverse(const Vect2& v, real alt, bool north) {
  	return SimpleProjection::polar_inverse(v, alt, north);
  }




}

