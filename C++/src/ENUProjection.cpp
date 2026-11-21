/*
 * ENUProjection.cpp
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
#include "ENUProjection.h"
#include "EuclideanProjection.h"
#include "LatLonAlt.h"
#include "GreatCircle.h"
#include "Util.h"
#include "Units.h"
//#include "UnitSymbols.h"
#include "Constants.h"
#include <stdexcept>

// For latitudes up to 86 degrees this should provide accuracies within 0.01 NM for legs up to 55 NM in length (ownship as reference point)
// "   "         "  "  89 degrees "    "      "       "          "      "    "  "   "    "  "  25 NM "  "       "       "  "         "
// "   "            >  89 degrees "    "      "       "          "      "    "  "   "    "  "   5 NM "  "       "       "  "         "
// For latitudes up to 86 degrees this should provide accuracies within 0.1 NM for legs up to 120 NM in length (ownship as reference point)
// "   "         "  "  89 degrees "    "      "       "          "      "    "  "   "    "  "  29 NM "  "       "       "  "         "
// "   "            >  89 degrees "    "      "       "          "      "    "  "   "    "  "   7 NM "  "       "       "  "         "
// This should provide accuracies within 0.5 NM for legs up to 190 NM in length (ownship as reference point)


namespace larcfm {
    using std::string;
    using std::cout;
    using std::endl;
    using std::runtime_error;
  
  // Static functions
  
  static Vect3 spherical2xyz(real lat, real lon) {
    real r = GreatCircle::spherical_earth_radius;
    // convert latitude to 0-PI
    real theta = Pi/2 - lat;
    real phi = Pi - lon;
    real rsintheta = r*std::sin(theta);
    real x = rsintheta*std::cos(phi);
    real y = rsintheta*std::sin(phi);
    real z = r*std::cos(theta);
    return Vect3(x,y,z);
  }
  
  static LatLonAlt xyz2spherical(const Vect3& v, real alt) {
    real r = GreatCircle::spherical_earth_radius;
    real theta = acos_safe(v.z()/r);
    real phi = atan2_safe(v.y(), v.x());
    real lat = Pi/2 - theta;
    real lon = to_pi(Pi - phi);
    return LatLonAlt::mk(lat, lon, alt);
  }
  
  static Vect3 vect3_orthog_toy(const Vect3& v) {
    if (!Util::within_epsilon(v.x(), Constants::GPS_LIMIT_HORIZONTAL) || !Util::within_epsilon(v.y(), Constants::GPS_LIMIT_HORIZONTAL)) {
      return Vect3(v.y(), -v.x(), 0);
    } else {
      return Vect3(1,0,0);
    }
  }
  
  static Vect3 equator_map(const Vect3& ref, const Vect3& p) {
    Vect3 xmult = ref.Hat();
    Vect3 ymult = vect3_orthog_toy(ref).Hat();
    Vect3 zmult = ref.cross(vect3_orthog_toy(ref)).Hat();
    return Vect3(xmult.dot(p), ymult.dot(p), zmult.dot(p));
  }
  
  static Vect3 equator_map_inv(const Vect3& ref, const Vect3& p) {
    Vect3 xmult = ref.Hat();
    Vect3 ymult = vect3_orthog_toy(ref).Hat();
    Vect3 zmult = ref.cross(vect3_orthog_toy(ref)).Hat();
    Vect3 xmultInv = Vect3(xmult.x(), ymult.x(), zmult.x());
    Vect3 ymultInv = Vect3(xmult.y(), ymult.y(), zmult.y());
    Vect3 zmultInv = Vect3(xmult.z(), ymult.z(), zmult.z());
    return  Vect3(xmultInv.dot(p), ymultInv.dot(p), zmultInv.dot(p));
  }
  
  static Vect2 sphere_to_plane(const Vect3& ref, const Vect3& p) {
    Vect3 v = equator_map(ref,p);
    return Vect2(v.y(), -v.z());
  }
  
  static Vect3 plane_to_sphere(const Vect2& v) {
    real r = GreatCircle::spherical_earth_radius;
    real x = std::sqrt(r*r - v.x*v.x - v.y*v.y);
    return Vect3(x, v.x, -v.y);
  }
  
  
    // ENUProjection
    
    ENUProjection::ENUProjection() : 
    	ref(Vect3()),
    	llaRef(LatLonAlt::ZERO()) {
      projAlt = 0;
      //ref = Vect3();
      //llaRef = LatLonAlt::ZERO();
    }
    
    ENUProjection::ENUProjection(const LatLonAlt& lla) :
    	ref(spherical2xyz(lla.lat(),lla.lon())),
    	llaRef(lla) {
        projAlt = lla.alt();
        //ref = spherical2xyz(lla.lat(),lla.lon());
        //llaRef = lla;
    }
    
    ENUProjection::ENUProjection(real lat, real lon, real alt) :
    	ref(spherical2xyz(lat,lon)),
    	llaRef(LatLonAlt::mk(lat, lon, alt)) {
        projAlt = alt;
        //ref = spherical2xyz(lat,lon);
        //llaRef = LatLonAlt::mk(lat, lon, alt);
    }
    
    ENUProjection ENUProjection::makeNew(const LatLonAlt& lla) const {
      return ENUProjection(lla);
    }

    ENUProjection ENUProjection::makeNew(real lat, real lon, real alt) const {
      return ENUProjection(lat, lon, alt);
    }

    real ENUProjection::conflictRange(real lat, real accuracy) const {
//      if (accuracy < Units::from("NM", 0.01)) { //~0.001 nm accuracy
//        return Units::from("NM", 18);
//      } else if (accuracy < Units::from("NM", 0.1)) {	//0.01 nm accuracy
//        return Units::from("NM", 50);
//      } else if (accuracy < Units::from("NM", 0.5)) { //0.1 nm accuracy
//        return Units::from("NM", 110);
//      } else {  //0.5 nm accuracy
//        return Units::from("NM", 205);
//      }
      return Units::from("NM", std::floor(real(243.0)*std::pow(Units::to("NM",std::ceil(accuracy)),real(1.0)/real(3.0))));
   }
    
    real ENUProjection::maxRange() const{
      return Units::from("NM", 3400);
    }
    
	LatLonAlt ENUProjection::getProjectionPoint() const {
	  return llaRef;
	}

    Vect2 ENUProjection::project2(const LatLonAlt& lla) const {
      return sphere_to_plane(ref, spherical2xyz(lla.lat(),lla.lon()));
    }

    Vect3 ENUProjection::project(const LatLonAlt& lla) const {
       return Vect3(project2(lla),lla.alt() - projAlt);
    }      
    
    Vect3 ENUProjection::project(const Position& sip) const {
    	Vect3 si;
    	if (sip.isLatLon()) {
    		si = project(sip.lla());
    	} else {
    		si = sip.vect3();
    	}
    	return si;
    }

    Vect3 ENUProjection::projectPoint(const Position& sip) const {
    	return project(sip);
    }

    LatLonAlt ENUProjection::inverse(const Vect2& xy, real alt) const {
      return xyz2spherical(equator_map_inv(ref, plane_to_sphere(xy)), alt + projAlt);
    }

    LatLonAlt ENUProjection::inverse(const Vect3& xyz) const {  
    	return inverse(xyz.vect2(), xyz.z());
    }
    
  Velocity ENUProjection::projectVelocity(const LatLonAlt& lla, const Velocity& v) const {
    real timeStep = real(10.0);
    LatLonAlt ll2 = GreatCircle::linear_initial(lla,v,timeStep);
    Vect3 se = project(lla);
    Vect3 s2 = project(ll2);
    Vect3 vn = s2.Sub(se).Scal(1/timeStep);
    return Velocity::make(vn);
  }
  
  Velocity ENUProjection::projectVelocity(const Position& ss, const Velocity& v) const {
	  if (ss.isLatLon()) {
        return projectVelocity(ss.lla(),v);
	  } else {
		  return v;
	  }
  }

  // The user needs to keep track of whether to translate back (i.e. whether original was LatLon())
  Velocity ENUProjection::inverseVelocity(const Vect3& s, const Velocity& v, bool toLatLon) const {
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

  std::pair<Vect3,Velocity> ENUProjection::project(const Position& p, const Velocity& v) const {
   	return std::pair<Vect3,Velocity>(project(p),projectVelocity(p,v));
   }

  std::pair<Vect3,Velocity> ENUProjection::project(const LatLonAlt& p, const Velocity& v) const {
   	return std::pair<Vect3,Velocity>(project(p),projectVelocity(p,v));
   }

  std::pair<Position,Velocity> ENUProjection::inverse(const Vect3& p, const Velocity& v, bool toLatLon) const {
    if (toLatLon) {
      return std::pair<Position,Velocity>(Position(inverse(p)),inverseVelocity(p,v,true));
    } else {
      return std::pair<Position,Velocity>(Position(p),v);
    }
   }


}

