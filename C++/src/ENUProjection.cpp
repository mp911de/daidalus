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
  
  static Vect3 plane_to_sphere(const Vect2& v) {
    real r = GreatCircle::spherical_earth_radius;
    real x = std::sqrt(r*r - v.x*v.x - v.y*v.y);
    return Vect3(x, v.x, -v.y);
  }
  
  
    // ENUProjection
    
    ENUProjection::ENUProjection() :
      projAlt(0), radial{0,0,0}, east{1,0,0}, south{0,0,0}, llaRef(LatLonAlt::ZERO()) {}

    ENUProjection::ENUProjection(const LatLonAlt& lla) :
      projAlt(lla.alt()), llaRef(lla) {
      const Vect3 ref = spherical2xyz(lla.lat(),lla.lon());
      const Vect3 orthogonal = vect3_orthog_toy(ref);
      auto axis = [](const Vect3& v) {
        const Vect3 unit = v.Hat();
        return Axis{unit.x(),unit.y(),unit.z()};
      };
      radial = axis(ref);
      east = axis(orthogonal);
      south = axis(ref.cross(orthogonal));
    }

    ENUProjection::ENUProjection(real lat, real lon, real alt) :
      ENUProjection(LatLonAlt::mk(lat,lon,alt)) {}

    std::string ENUProjection::toString() const {
      const Vect3 ref = radial.x == 0 && radial.y == 0 && radial.z == 0 ?
          Vect3() : spherical2xyz(llaRef.lat(),llaRef.lon());
      return "ENUProjection "+to_string(projAlt)+" "+ref.toString();
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
      const Vect3 p = spherical2xyz(lla.lat(),lla.lon());
      return Vect2(p.dot(east.x,east.y,east.z), -p.dot(south.x,south.y,south.z));
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
      const Vect3 p = plane_to_sphere(xy);
      return xyz2spherical(Vect3(p.dot(radial.x,east.x,south.x),
                                 p.dot(radial.y,east.y,south.y),
                                 p.dot(radial.z,east.z,south.z)),alt + projAlt);
    }

    LatLonAlt ENUProjection::inverse(const Vect3& xyz) const {  
    	return inverse(xyz.vect2(), xyz.z());
    }
    
  Velocity ENUProjection::projectVelocity(const LatLonAlt& lla, const Velocity& v) const {
    return project(lla,v).second;
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
    return toLatLon ? inverse(s,v,true).second : v;
  }

  std::pair<Vect3,Velocity> ENUProjection::project(const Position& p, const Velocity& v) const {
    return p.isLatLon() ? project(p.lla(),v) : std::make_pair(p.vect3(),v);
   }

  std::pair<Vect3,Velocity> ENUProjection::project(const LatLonAlt& p, const Velocity& v) const {
    const real timeStep = 10;
    const Vect3 start = project(p);
    const Vect3 end = project(GreatCircle::linear_initial(p,v,timeStep));
    return std::make_pair(start,Velocity::make(end.Sub(start).Scal(1/timeStep)));
   }

  std::pair<Position,Velocity> ENUProjection::inverse(const Vect3& p, const Velocity& v, bool toLatLon) const {
    if (toLatLon) {
      const real timeStep = 10;
      const LatLonAlt start = inverse(p);
      const LatLonAlt end = inverse(p.linear(v.vect3(),timeStep));
      return std::make_pair(Position(start),GreatCircle::velocity_initial(start,end,timeStep));
    } else {
      return std::pair<Position,Velocity>(Position(p),v);
    }
   }


}

