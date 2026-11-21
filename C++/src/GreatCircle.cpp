/* 
 * GreatCircle
 *
 * Contact: Jeff Maddalon
 * Organization: NASA/Langley Research Center
 *
 * Copyright (c) 2011-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */

#include "Real.h"
#include "Util.h"
#include "GreatCircle.h"
#include "Triple.h"

#include <cmath>
namespace larcfm {



static const real EPS = real(1.0e-15); // small number, about machine

const real GreatCircle::spherical_earth_radius = real(6366707.0194937070000000000); // Units::from(_m, 1.0 / angle_from_distance(1.0));

const real GreatCircle::minDt = real(1E-5);

real GreatCircle::angle_from_distance(real distance) {
	return Units::to("nm", distance) * Pi / (real(180.0) * real(60.0));
}

real GreatCircle::angle_from_distance(real distance, real h) {
	return angle_from_distance(distance * spherical_earth_radius
			/ (spherical_earth_radius + h));
}

real GreatCircle::distance_from_angle(real angle, real h) {
	return (spherical_earth_radius+h) * angle;
}

real GreatCircle::angular_distance(real lat1, real lon1, real lat2, real lon2) {
	//	  if (Util::almost_equals(lat1,lat2) && Util::almost_equals(lon1,lon2)) return 0.0;
	real rtn= asin_safe(sqrt_safe(sq(std::sin((lat1 - lat2) / 2))
			+ std::cos(lat1) * std::cos(lat2)
			* sq(std::sin((lon1 - lon2) / 2)))) * real(2.0);
	return rtn;
}

real GreatCircle::angular_distance(const LatLonAlt& p1, const LatLonAlt& p2) {
	return angular_distance(p1.lat(), p1.lon(), p2.lat(), p2.lon());
}

real GreatCircle::distance(real lat1, real lon1, real lat2, real lon2) {
	return distance_from_angle(angular_distance(lat1, lon1, lat2, lon2), real(0.0));
}

real GreatCircle::distance(const LatLonAlt& p1, const LatLonAlt& p2) {
	//std::cout << p1.toString() << " , " << p2.toString() << std::endl;
	//return distance_from_angle(angular_distance(p1, p2),(p1.alt() + p2.alt())/2.0);
	return distance_from_angle(angular_distance(p1, p2), real(0.0));
}


bool GreatCircle::almost_equals(real lat1, real lon1, real lat2, real lon2) {
	return Constants::almost_equals_radian(angular_distance(lat1, lon1, lat2, lon2));
}


bool GreatCircle::almost_equals(real lat1, real lon1, real lat2,
		real lon2, real epsilon) {
	return Util::within_epsilon(distance(lat1, lon1, lat2, lon2), epsilon);
}

  bool GreatCircle::almostEquals(const LatLonAlt& b, const LatLonAlt& a) {
	return almost_equals(b.lat(), b.lon(), a.lat(), a.lon())
	&& Constants::almost_equals_alt(b.alt(), a.alt());
}

bool GreatCircle::almostEquals(const LatLonAlt& b, const LatLonAlt& a, real horizEps, real vertEps) {
	return almost_equals(b.lat(), b.lon(), a.lat(), a.lon(), horizEps)
	&& Util::within_epsilon(b.alt(), a.alt(), vertEps);
}

bool GreatCircle::almostEquals2D(const LatLonAlt& b, const LatLonAlt& a, real horizEps) {
  return almost_equals(b.lat(), b.lon(), a.lat(), a.lon(), horizEps);
}

static real initial_course_impl(const LatLonAlt& p1, const LatLonAlt& p2) {
	real lat1 = p1.lat();
	real lon1 = p1.lon();
	real lat2 = p2.lat();
	real lon2 = p2.lon();

	if (std::cos(lat1) < EPS) { // EPS a small number, about machine precision
		if (lat1 > 0) {
			return Pi; // starting from North pole, all directions are south
		} else {
			return real(2.0) * Pi; // starting from South pole, all directions are
			// north. JMM: why not 0?
		}
	}

	real tc1 = to_2pi(atan2_safe(std::sin(lon2-lon1)*std::cos(lat2), std::cos(lat1)*std::sin(lat2)-std::sin(lat1)*std::cos(lat2)*std::cos(lon2-lon1)));
	return tc1;
}



real GreatCircle::initial_course(real lat1, real lon1, real lat2, real lon2) {
	LatLonAlt p1 = LatLonAlt::mk(lat1, lon1, real(0.0));
	LatLonAlt p2 = LatLonAlt::mk(lat2, lon2, real(0.0));
	//double d = angular_distance(lat1, lon1, lat2, lon2);
	return initial_course_impl(p1, p2); //, d);
}

real GreatCircle::initial_course(const LatLonAlt& p1, const LatLonAlt& p2) {
	return initial_course(p1.lat(),p1.lon(),
			p2.lat(),p2.lon());
}


real GreatCircle::final_course(const LatLonAlt& p1, const LatLonAlt& p2) {
	return initial_course(p2, p1)+real(M_PI);
}


real GreatCircle::max_latitude_gc(real lat1, real lon1, real lat2, real lon2) {
	return max_latitude_gc_course(lat1, initial_course(lat1,lon1,lat2,lon2));
}


real GreatCircle::max_latitude_gc_course(real lat1, real trk) {
	return Util::acos_safe(std::abs(std::sin(trk)*std::cos(lat1)));
}

real GreatCircle::max_latitude_gc(const LatLonAlt& p1, const LatLonAlt& p2) {
	return max_latitude(p1.lat(), p1.lon(), p2.lat(), p2.lon());
}

real GreatCircle::min_latitude_gc(real lat1, real lon1, real lat2, real lon2) {
	return min_latitude_gc_course(lat1, initial_course(lat1,lon1,lat2,lon2));
}

real GreatCircle::min_latitude_gc_course(real lat1, real trk) {
	return -Util::acos_safe(std::abs(std::sin(trk)*std::cos(lat1)));
}

real GreatCircle::min_latitude_gc(const LatLonAlt& p1, const LatLonAlt& p2) {
	return min_latitude(p1.lat(), p1.lon(), p2.lat(), p2.lon());
}

// Given a great circle defined by Point 1 and 2, find the longitude of where it
// crosses the latitude defined by lat3.
//
// from Aviation Formulary
// longitude sign is reversed from the formulary!
real GreatCircle::lonCross(real lat1, real lon1, real lat2, real lon2, real lat3) {
	//f.pln("lonCross 1 "+i);
	//double lat1 = ac.point(i).lat();
	//double lon1 = ac.point(i).lon();
	//double lat2 = ac.point(i+1).lat();
	//double lon2 = ac.point(i+1).lon();
	//double tc = ac.initialVelocity(i).compassAngle();
	real tc = initial_course(lat1,lon1,lat2,lon2);
	bool NW = (tc > real(M_PI)/2 && tc <= real(M_PI)) || tc >= 3*real(M_PI)/2;
	real l12;

	if (NW) l12 = lon1-lon2;
	else l12 = lon2-lon1;

	real A = std::sin(lat1)*std::cos(lat2)*std::cos(lat3)*std::sin(l12);
	real B = std::sin(lat1)*std::cos(lat2)*std::cos(lat3)*std::cos(l12) - std::cos(lat1)*std::sin(lat2)*std::cos(lat3);
	real lon;

	if (NW) lon = lon1 + Util::atan2_safe(B,A) + real(M_PI);
	else lon = lon1 - Util::atan2_safe(B,A) - real(M_PI);

	if (lon >= 2*real(M_PI)) lon = lon-2*real(M_PI);

	if (NW) {
		lon = lon-real(M_PI);
	} else {
		lon = real(M_PI)+lon;
	}

	if (lon < -real(M_PI)) lon = 2*real(M_PI)+lon;
	if (lon > real(M_PI)) lon = -2*real(M_PI)+lon;

	//f.pln("lonCross 2 "+i);
	return lon;
}


real GreatCircle::max_latitude(real lat1, real lon1, real lat2, real lon2) {
	if (lat1 <= 0 && lat2 <= 0) return Util::max(lat1, lat2);  // the segment will, at most, curve north
	real dist = distance(lat1,lon1,lat2,lon2);
	if (Util::sign(lat1) != Util::sign(lat2) && dist < spherical_earth_radius*Pi/2) return Util::max(lat1, lat2);
//	double maxLat = max_latitude_gc(lat1,lon1,lat2,lon2);
	real maxLat;
	// this block is new --- GEH
	if (dist < GreatCircle::spherical_earth_radius*Pi) {
		// special case: we are in the northern hemisphere and we are heading south, then the first point is the max
		if (lat1 >= lat2) {
			real trk = initial_course(lat1, lon1, lat2, lon2);
			if (lat1 > 0 && (trk >= real(0.5)*Pi && trk <= real(1.5)*Pi)) {
				return lat1;
			} else {
				maxLat = max_latitude_gc_course(lat1,trk);
			}
		} else {
			real trk = initial_course(lat2, lon2, lat1, lon1);
			if (lat2 > 0 && (trk >= real(0.5)*Pi && trk <= real(1.5)*Pi)) {
				return lat2;
			} else {
				maxLat = max_latitude_gc_course(lat2,trk);
			}
		}
		// END BLOCK
	} else {
		maxLat = max_latitude_gc(lat1,lon1,lat2,lon2);
	}
	real maxLon = lonCross(lat1,lon1,lat2,lon2,maxLat);
	if (Util::max(distance(lat1,lon1,maxLat,maxLon),distance(lat2,lon2,maxLat,maxLon)) < distance(lat1,lon1,lat2,lon2)) {
		return maxLat;
	}
	return Util::max(lat1, lat2);
}

real GreatCircle::max_latitude(const LatLonAlt& p1, const LatLonAlt& p2) {
	return max_latitude(p1.lat(), p1.lon(), p2.lat(), p2.lon());
}

real GreatCircle::min_latitude(real lat1, real lon1, real lat2, real lon2) {
	if (lat1 >= 0 && lat2 >= 0) return Util::min(lat1, lat2);  // the segment will, at most, curve north
	real dist = distance(lat1,lon1,lat2,lon2);
	if (Util::sign(lat1) != Util::sign(lat2) && dist < spherical_earth_radius*Pi/2) return Util::min(lat1, lat2);
//	double minLat = min_latitude_gc(lat1,lon1,lat2,lon2);
	real minLat;
	// this block is new --- GEH
	if (dist < GreatCircle::spherical_earth_radius*Pi) {
		// special case: we are in the southern hemisphere and we are heading north, then the first point is the max
		if (lat1 <= lat2) {
			real trk = initial_course(lat1, lon1, lat2, lon2);
			if (lat1 < 0 && (trk <= real(0.5)*Pi || trk >= real(1.5)*Pi)) {
				return lat1;
			} else {
				minLat = min_latitude_gc_course(lat1,trk);
			}
		} else {
			real trk = initial_course(lat2, lon2, lat1, lon1);
			if (lat2 < 0 && (trk <= real(0.5)*Pi || trk >= real(1.5)*Pi)) {
				return lat2;
			} else {
				minLat = min_latitude_gc_course(lat2,trk);
			}
		}
		// END BLOCK
	} else {
		minLat = min_latitude_gc(lat1,lon1,lat2,lon2);
	}
	real minLon = lonCross(lat1,lon1,lat2,lon2,minLat);
	if (Util::max(distance(lat1,lon1,minLat,minLon),distance(lat2,lon2,minLat,minLon)) < distance(lat1,lon1,lat2,lon2)) {
		return minLat;
	}
	return Util::min(lat1, lat2);
}

real GreatCircle::min_latitude(const LatLonAlt& p1, const LatLonAlt& p2) {
	return min_latitude(p1.lat(), p1.lon(), p2.lat(), p2.lon());
}


// parameter d is the angular distance between lat/long #1 and #2
static LatLonAlt interpolate_impl(const LatLonAlt& p1, const LatLonAlt& p2, real d, real f, real alt) {
	if (Constants::almost_equals_radian(d) ) {
		return p1.mkAlt(alt);
		// if the two points are almost the same, then consider the two
		// points the same and arbitrarily return one of them (in this case
		// p1)
		// with the altitude that was provided
	}

	real lat1 = p1.lat();
	real lon1 = p1.lon();
	real lat2 = p2.lat();
	real lon2 = p2.lon();

	real a = std::sin((1 - f) * d) / std::sin(d);
	real b = std::sin(f * d) / std::sin(d);
	real x = a * std::cos(lat1) * std::cos(lon1) + b * std::cos(lat2) * std::cos(lon2);
	real y = a * std::cos(lat1) * std::sin(lon1) + b * std::cos(lat2) * std::sin(lon2);
	real z = a * std::sin(lat1) + b * std::sin(lat2);
	return LatLonAlt::mk(atan2_safe(z, std::sqrt(x * x + y * y)), // lat
			atan2_safe(y, x), // longitude
			alt); // long
}

real GreatCircle::representative_course(real lat1, real lon1, real lat2, real lon2) {
	LatLonAlt p1 = LatLonAlt::mk(lat1, lon1, real(0.0));
	LatLonAlt p2 = LatLonAlt::mk(lat2, lon2, real(0.0));
	real d = angular_distance(lat1, lon1, lat2, lon2);
	LatLonAlt midPt = interpolate_impl(p1, p2, d, real(0.5), real(0.0));
	return initial_course_impl(midPt, p2); //, d / 2.0);
}


real GreatCircle::representative_course(const LatLonAlt& p1, const LatLonAlt& p2) {
	return representative_course(p1.lat(),p1.lon(),p2.lat(),p2.lon());
}

LatLonAlt GreatCircle::interpolate(const LatLonAlt& p1, const LatLonAlt& p2, real f) {
	real d = angular_distance(p1, p2);
	return interpolate_impl(p1, p2, d, f, (p2.alt() - p1.alt())*f + p1.alt());
}

LatLonAlt GreatCircle::interpolateEst(const LatLonAlt& p1, const LatLonAlt& p2, real f) {
	return LatLonAlt::mk((p2.lat() - p1.lat()) * f + p1.lat(),
			(p2.lon() - p1.lon()) * f + p1.lon(),
			(p2.alt() - p1.alt()) * f + p1.alt());
}

static LatLonAlt linear_initial_impl(const LatLonAlt& s, real track, real d, real vertical) {
	real lat = asin_safe(std::sin(s.lat())*std::cos(d)+std::cos(s.lat())*std::sin(d)*std::cos(track));
	real dlon = atan2_safe(std::sin(track)*std::sin(d)*std::cos(s.lat()),std::cos(d)-std::sin(s.lat())*std::sin(lat));
	// slightly different from aviation formulary because I use
	// "east positive" convention
	real lon = to_pi(s.lon() + dlon);

	return LatLonAlt::mk(lat, lon, s.alt()+vertical);
}

static LatLonAlt linear_rhumb_impl(const LatLonAlt& s, real track, real d, real vertical) {
	// -- Based on the calculation in the "Rhumb line" section of the
	//    Aviation Formulary v1.44
	// -- Weird things happen to rhumb lines that go to the poles, therefore
	//    force any polar latitudes to be "near" the pole

	static const real eps = real(1e-15);
	real s_lat = Util::max(Util::min(s.lat(), Pi/2-eps), -Pi/2+eps);
	real lat = s_lat + d * std::cos(track);
	lat = Util::max(Util::min(lat, Pi/2-eps), -Pi/2+eps);

	real q;
	if ( Constants::almost_equals_radian(lat, s_lat) ) {
		// (std::abs(lat - lat1) < EPS) {
		q = std::cos(s_lat);
	} else {
		real dphi = std::log(std::tan(lat / 2 + Pi / 4)
				/ std::tan(s_lat / 2 + Pi / 4));
		q = (lat - s_lat) / dphi;
	}
	real dlon = -d * std::sin(track) / q;

	// slightly different from aviation formulary because I use
	// "east positive" convention
	real lon = to_pi(s.lon() - dlon);
	return LatLonAlt::mk(lat, lon, s.alt()+vertical);
}

LatLonAlt GreatCircle::linear_gcgs(const LatLonAlt& p1, const LatLonAlt& p2, const Velocity& v, real t) {
	real d = GreatCircle::angular_distance(p1, p2);
	if ( Constants::almost_equals_radian(d) ) {
		// If the two points are about 1 meter apart, then count them as the
		// same.
		return p1;
	}
	real f = angle_from_distance(v.gs() * t) / d;
	return interpolate_impl(p1, p2, d, f, p1.alt() + v.z()*t);
}

LatLonAlt GreatCircle::linear_gc(const LatLonAlt& p1, const LatLonAlt& p2, real d) {
	//return GreatCircle.linear_initial(p1, initial_course(p1,p2), d);
	real dist = angular_distance(p1,p2);
	real f = angle_from_distance(d)/dist;
	return interpolate_impl(p1, p2, dist, f, (p2.alt() - p1.alt())*f + p1.alt());
}

LatLonAlt GreatCircle::linear_rhumb(const LatLonAlt& s, const Velocity& v, real t) {
	return linear_rhumb_impl(s, v.trk(), GreatCircle::angle_from_distance(v.gs() * t), v.z()*t);
}

LatLonAlt GreatCircle::linear_rhumb(const LatLonAlt& s, real track, real dist) {
	return linear_rhumb_impl(s, track, GreatCircle::angle_from_distance(dist), real(0.0));
}


/**
 * Solve the spherical triangle when one has a side (in angular distance), another side, and an angle between sides.
 * The angle is <b>not</b> between the sides.  The sides are labeled a, b, and c.  The angles are labelled A, B, and
 * C.  Side a is opposite angle A, and so forth.<p>
 *
 * Given these constraints, in some cases two solutions are possible.  To
 * get one solution set the parameter firstSolution to true, to get the other set firstSolution to false.  A firstSolution == true
 * will return a smaller angle, B, than firstSolution == false.
 *
 * @param b one side (in angular distance)
 * @param a another side (in angular distance)
 * @param A the angle opposite the side a
 * @param firstSolution select which solution to use
 * @return a Triple of angles B and C, and the side c.
 */
Triple<real,real,real> GreatCircle::side_side_angle(real b, real a, real A, bool firstSolution) {
	// This function follows the convention of "Spherical Trigonometry" by Todhunter, Macmillan, 1886
	//   Note, angles are labelled counter-clockwise a, b, c

	// Law of sines
	real B = Util::asin_safe(std::sin(b)*std::sin(A)/std::sin(a));  // asin returns [-pi/2,pi/2]
	if ( ! firstSolution) {
		B = real(M_PI) - B;
	}

	// one of Napier's analogies
	real c = 2 * Util::atan2_safe(std::sin(real(0.5)*(a+b))*std::cos(real(0.5)*(A+B)),std::cos(real(0.5)*(a+b))*std::cos(real(0.5)*(A-B)));

	// Law of cosines
	real C = Util::acos_safe(-std::cos(A)*std::cos(B)+std::sin(A)*std::sin(B)*std::cos(c));

	if ( gauss_check(a,b,c,A,B,C)) {
		return Triple<real,real,real>(Util::to_pi(B),C,Util::to_2pi(c));
	} else {
		return Triple<real,real,real>(real(0.0),real(0.0),real(0.0));
	}
}

/**
 * Solve the spherical triangle when one has a side (in angular distance), and two angles.
 * The side is <b>not</b> between the angles.  The sides are labeled a, b, and c.  The angles are labelled A, B, and
 * C.  Side a is opposite angle A, and so forth.<p>
 *
 * Given these constraints, in some cases two solutions are possible.  To
 * get one solution set the parameter firstSolution to true, to get the other set firstSolution to false.  A firstSolution == true
 * will return a smaller side, b, than firstSolution == false.
 *
 * @param a one side (in angular distance)
 * @param A the angle opposite the side a
 * @param B another angle
 * @param firstSolution select which solution to use
 * @return a Triple of side b, angle C, and the side c.
 */
Triple<real,real,real> GreatCircle::side_angle_angle(real a, real A, real B, bool firstSolution) {
	// This function follows the convention of "Spherical Trigonometry" by Todhunter, Macmillan, 1886
	//   Note, angles are labelled counter-clockwise a, b, c

	// Law of sines
	real b = Util::asin_safe(std::sin(a)*std::sin(B)/std::sin(A));  // asin returns [-pi/2,pi/2]
	if ( ! firstSolution) {
		b = real(M_PI) - b;
	}

	// one of Napier's analogies
	real c = 2 * Util::atan2_safe(std::sin(real(0.5)*(a+b))*std::cos(real(0.5)*(A+B)),std::cos(real(0.5)*(a+b))*std::cos(real(0.5)*(A-B)));

	// Law of cosines
	real C = Util::acos_safe(-std::cos(A)*std::cos(B)+std::sin(A)*std::sin(B)*std::cos(c));

	if ( gauss_check(a,b,c,A,B,C)) {
		return Triple<real,real,real>(Util::to_2pi(b),Util::to_2pi(C),Util::to_2pi(c));
	} else {
		return Triple<real,real,real>(real(0.0),real(0.0),real(0.0));
	}
}

/**
 * This implements the spherical cosine rule to complete a triangle on the unit sphere
 * @param a side a (angular distance)
 * @param C angle between sides a and b
 * @param b side b (angular distance)
 * @return triple of A,B,c (angle opposite a, angle opposite b, side opposite C)
 */
Triple<real,real,real> GreatCircle::side_angle_side(real a, real C, real b) {
	real c = Util::acos_safe(std::cos(a)*std::cos(b)+std::sin(a)*std::sin(b)*std::cos(C));
	real cRatio = std::sin(C)/std::sin(c);
	real A = Util::asin_safe(std::sin(a)*cRatio);
	real B = Util::asin_safe(std::sin(b)*cRatio);
	return Triple<real,real,real>(A,B,c);
}

/**
 * This implements the supplemental (polar triangle) spherical cosine rule to complete a triangle on the unit sphere
 * @param A angle A
 * @param c side between A and B (angular distance
 * @param B angle B
 * @return triple of a,b,C (side opposite A, side opposite B, angle opposite c)
 */
Triple<real,real,real> GreatCircle::angle_side_angle(real A, real c, real B) {
	real C = Util::acos_safe(-std::cos(A)*std::cos(B)+std::sin(A)*std::sin(B)*std::cos(c));
	real cRatio = std::sin(c)/std::sin(C);
	real a = Util::asin_safe(std::sin(A)*cRatio);
	real b = Util::asin_safe(std::sin(B)*cRatio);
	return Triple<real,real,real>(a,b,C);
}

bool GreatCircle::gauss_check(real a, real b, real c, real A, real B, real C) {
	// This function follows the convention of "Spherical Trigonometry" by Todhunter, Macmillan, 1886
	//   Note, angles are labelled counter-clockwise a, b, c
	A = Util::to_pi(A);
	B = Util::to_pi(B);
	C = Util::to_pi(C);
	a = Util::to_2pi(a);
	b = Util::to_2pi(b);
	c = Util::to_2pi(c);
	if (A==real(0.0) || A==real(M_PI) || B==real(0.0) || B==real(M_PI) || C==real(0.0) || C==real(M_PI)) return false;
	if (a==real(0.0) || b==real(0.0) || c==real(0.0)) return false;
	//		f.pln("gauss "+std::cos(0.5*(A+B))*std::cos(0.5*c)+" "+std::cos(0.5*(a+b))*std::sin(0.5*C));
	return Util::almost_equals(std::cos(real(0.5)*(A+B))*std::cos(real(0.5)*c),std::cos(real(0.5)*(a+b))*std::sin(real(0.5)*C),PRECISION13);
}



LatLonAlt GreatCircle::linear_initial(const LatLonAlt& s, const Velocity& v, real t) {
	return linear_initial_impl(s, v.trk(), GreatCircle::angle_from_distance(v.gs() * t), v.z()*t);
}

LatLonAlt GreatCircle::linear_initial(const LatLonAlt& s, real track, real dist) {
	return linear_initial_impl(s, track, GreatCircle::angle_from_distance(dist), real(0.0));
}

real GreatCircle::cross_track_distance(const LatLonAlt& p1, const LatLonAlt& p2, const LatLonAlt& offCircle) {
	real dist_p1oc = angular_distance(p1,offCircle);
	real trk_p1oc = initial_course_impl(p1,offCircle); //,dist_p1oc);
	real trk_p1p2 = initial_course(p1,p2);
	// This is a direct application of the "spherical law of sines"
	return distance_from_angle(Util::asin_safe(std::sin(dist_p1oc)*std::sin(trk_p1oc-trk_p1p2)), (p1.alt()+p2.alt()+offCircle.alt())/real(3.0));
}

bool GreatCircle::collinear(const LatLonAlt& p1, const LatLonAlt& p2, const LatLonAlt& p3) {
	real epsilon = real(1E-7);
	return Util::within_epsilon(cross_track_distance(p1,p2,p3),epsilon);
}


/**
 * This returns the point on the great circle running through p1 and p2 that is closest to point x.
 * The altitude of the output is the same as x.<p>
 * If p1 and p2 are the same point, then every great circle runs through them, thus x is on one of these great circles.  In this case, x will be returned.
 * This assumes any 2 points will be within 90 degrees of each other (angular distance).
 * @param p1 the starting point of the great circle
 * @param p2 another point on the great circle
 * @param x point to determine closest segment point to.
 * @return the LatLonAlt point on the segment that is closest (horizontally) to x
 */
LatLonAlt GreatCircle::closest_point_circle(const LatLonAlt& p1, const LatLonAlt& p2, const LatLonAlt& x) {
	// almost same point or antipode:
	if ((Util::almost_equals(p1.lat(),p2.lat()) && Util::almost_equals(p1.lon(),p2.lon())) ||
			(Util::almost_equals(p1.lat(),-p2.lat()) && Util::almost_equals(p1.lon(),Util::to_pi(p2.lon()+Pi)))) return x;
	Vect3 a = spherical2xyz(p1.lat(), p1.lon());
	Vect3 b = spherical2xyz(p2.lat(), p2.lon());
	if (a.almostEquals(b) ||a.almostEquals(b.Neg())) return x;
	Vect3 c = a.cross(b);
	Vect3 p = spherical2xyz(x.lat(), x.lon());
	Vect3 g = p.Sub(c.Scal(p.dot(c)/c.sqv()));
	real v = spherical_earth_radius/g.norm();
	return xyz2spherical(g.Scal(v)).mkAlt(x.alt()); // return to x's altitude
}

/**
 * This returns the point on the great circle segment running through p1 and p2 that is closest to point x.
 * This will return either p1 or p2 if the actual closest point is outside the segment.
 * @param p1 the starting point of the great circle
 * @param p2 another point on the great circle
 * @param x point to determine closest segment point to.
 * @return the LatLonAlt point on the segment that is closest (horizontally) to x
 */
LatLonAlt GreatCircle::closest_point_segment(const LatLonAlt& p1, const LatLonAlt& p2, const LatLonAlt& x) {
	LatLonAlt c = closest_point_circle(p1,p2,x);
	real d12 = distance(p1,p2);
	real d1c = distance(p1,c);
	real d2c = distance(p2,c);
	if (d1c < d12 && d2c < d12) {
		return c;
	}
	if (d1c < d2c) {
		return p1;
	} else {
		return p2;
	}
}


//  LatLonAlt GreatCircle::closest_point_circle(const LatLonAlt& p1, const LatLonAlt& p2, const LatLonAlt& x) {
//		double a = angular_distance(x,p2);
//		double b = angular_distance(p1,p2);
//		double c = angular_distance(p1,x);
//		double A = angle_between(p2,p1,x);
//		double B = angle_between(p1,x,p2);
//		double C = angle_between(x,p2,p1);
//		return closest_point_circle(p1,p2,x,a,b,c,A,B,C);
//  }

LatLonAlt GreatCircle::closest_point_circle(const LatLonAlt& p1, const LatLonAlt& p2, const LatLonAlt& x, real a, real b, real c, real A, real B, real C) {
	//       x (B)
	//      / \.
	// (A) p1--p2 (C)
	//		double a = angular_distance(x,p2);
	//		double b = angular_distance(p1,p2);
	//		double c = angular_distance(p1,x);
	//		double A = angle_between(p2,p1,x);
	//		double B = angle_between(p1,x,p2);
	//		double C = angle_between(x,p2,p1);

	//		if (Util::almost_equals(A, 0.0) || Util::almost_equals(C, 0.0) || Util::almost_equals(A, M_PI) || Util::almost_equals(C, M_PI)) {			// collinear
	if (Util::within_epsilon(A, real(0.000001)) || Util::within_epsilon(C, real(0.000001)) || Util::within_epsilon(real(M_PI)-A, real(0.000001)) || Util::within_epsilon(real(M_PI)-C, real(0.000001))) {
		return x;
	}
	if (Util::almost_equals(b,real(0.0))) {
		return x;   // if p1==p2, every great circle runs through them, thus x is on one of these great circles
	}
	if (A+B+C < real(M_PI) || A+B+C >= real(M_PI)*3) {
		//			fpln("GreatCircle.closestPoint ERROR: not a triangle p1="+p1+"p2="+p2+"x="+x+" A+B+C="+(A+B+C)+" = "+Units.to("deg", A+B+C)+" deg");
		// if the triangle is relatively small, it is probably collinear
		if (a < real(M_PI)/2 && b < real(M_PI)/2 && c < real(M_PI)/2) {
			return x;
		}
		return LatLonAlt::INVALID();
	}
	if (almostEquals(p1,x) || Util::almost_equals(A, real(M_PI)/2)) {
		return p1.mkAlt(x.alt());
	}
	if (almostEquals(p2,x) || Util::almost_equals(C, real(M_PI)/2)) {
		return p2.mkAlt(x.alt());
	}

	// general case p1 @ A, x @ B, p2 @ C
	//		double d1 = 0;
	//f.pln("GreatCircle.closest_point_circle DEG A="+Units.to("deg", A)+" B="+Units.to("deg", B)+" C="+Units.to("deg", C)+" / a="+Units.to("deg", a)+" b="+Units.to("deg", b)+" c="+Units.to("deg", c));
	//f.pln("GreatCircle.closest_point_circle RAD A="+A+" B="+B+" C="+C+" / a="+a+" b="+b+" c="+c);
	//f.pln("GreatCircle.closest_point_circle A/a="+(std::sin(A)/std::sin(a))+" B/b="+(std::sin(B)/std::sin(b))+" C/c="+(std::sin(C)/std::sin(c)));
	if (A <= real(M_PI)/2 && C <= real(M_PI)/2) {
		//   B       C1
		//  / \     / |.
		// A---C   B1-A1
		if (A < C) {
			real a1 = c;
			real A1 = real(M_PI)/2;
			real B1 = A;
			real c1 = side_angle_angle(a1,A1,B1, true).third;
			real ff = (c1/b);
			//f.pln("GreatCircle.closest_point_circle a1) ff="+ff);
			return interpolate(p1,p2,ff);
		} else {
			//   B     C1
			//  / \    | \.
			// A---C   A1-B1
			real a1 = a;
			real A1 = real(M_PI)/2;
			real B1 = C;
			real c1 = side_angle_angle(a1,A1,B1, true).third;
			real ff = (c1/b);
			//f.pln("GreatCircle.closest_point_circle a2) ff="+ff);
			return interpolate(p2,p1,ff);
		}

		//			d1 = side_angle_angle(a,M_PI/2,C,true).third;
		//f.pln("GreatCircle.closest_point_circle #1 d1="+d1+" "+Units.to("deg", d1)+Units.degreeStr);
		//			double ff = 1-(d1/b);
		//f.pln("GreatCircle.closest_point_circle p1="+p1+" p2="+p2+" ff="+ff);
		//			return interpolate(p1, p2, ff);
	} else if (A <= real(M_PI)/2 && C > real(M_PI)/2) {
		//    -- B    C1
		//  /   /    / |.
		// A---C    B1-A1
		real a1 = a;
		real A1 = real(M_PI)/2;
		real B1 = real(M_PI)-C;
		real c1 = side_angle_angle(a1,A1,B1, true).third;
		real ff = 1+(c1/b);
		//f.pln("GreatCircle.closest_point_circle b) ff="+ff);
		return interpolate(p1,p2,ff);

		//			d1 = side_angle_angle(a,M_PI/2,M_PI-C,true).third;
		//f.pln("GreatCircle.closest_point_circle #2 d="+d1+" "+Units.to("deg", d1)+Units.degreeStr);
		//			return linear_initial(p1,initial_course(p1,p2),distance_from_angle(b+d1,0)).mkAlt(x.alt());
	} else if (A > real(M_PI)/2 && C <= real(M_PI)/2) {
		// B--		  C1
		//  \   \	  | \.
		//   A---C	  A1-B1
		real a1 = c;
		real A1 = real(M_PI)/2;
		real B1 = real(M_PI)-A;
		real c1 = side_angle_angle(a1,A1,B1, true).third;
		real ff = 1+(c1/b);
		//f.pln("GreatCircle.closest_point_circle c) ff="+ff);
		return interpolate(p2,p1,ff);
		//			d1 = side_angle_angle(a,M_PI/2,M_PI-A,true).third;
		//f.pln("GreatCircle.closest_point_circle #3 d="+d1+" "+Units.to("deg", d1)+Units.degreeStr);
		//			return linear_initial(p2,initial_course(p2,p1),distance_from_angle(b+d1,0)).mkAlt(x.alt());
	}
	//f.pln("GreatCircle.closest_point_circle INVALID: weird triangle");
	return LatLonAlt::INVALID(); // weird triangle
}

/**
 * Given two great circles defined by a1,a2 and b1,b2, return the intersection poin that is closest a1.  Use LatLonAlt.antipode() to get the other value.
 * This assumes that the arc distance between a1,a2 < 90 and b1,b2 < 90
 * This returns an INVALID value if both segments are collinear
 * EXPERIMENTAL
 */
LatLonAlt GreatCircle::intersection(const LatLonAlt& a1, const LatLonAlt& a2, const LatLonAlt& b1, const LatLonAlt& b2) {
	Vect3 va = spherical2xyz(a1.lat(), a1.lon()).cross(spherical2xyz(a2.lat(), a2.lon()));
	Vect3 vb = spherical2xyz(b1.lat(), b1.lon()).cross(spherical2xyz(b2.lat(), b2.lon()));
	real r = GreatCircle::spherical_earth_radius;
	Vect3 vavb = va.cross(vb);
	if (vavb.almostEquals(Vect3::ZERO())) {
		return LatLonAlt::INVALID();
	}
	Vect3 v1 = vavb.Scal(r / vavb.norm());
	Vect3 v2 = vavb.Scal(-r / vavb.norm());
	LatLonAlt x1 = xyz2spherical(v1).mkAlt(a1.alt());
	LatLonAlt x2 = xyz2spherical(v2).mkAlt(a1.alt());
	if (distance(a1,x1) < distance(a1,x2)) {
		return x1;
	} else {
		return x2;
	}
}

LatLonAlt GreatCircle::intersectSegments(const LatLonAlt& so, const LatLonAlt& so2, const LatLonAlt& si, const LatLonAlt& si2) {
	LatLonAlt interSec = GreatCircle::intersection(so,so2, si, si2);
	if (interSec.isInvalid()) return LatLonAlt::INVALID();
	real fco = GreatCircle::final_course(so,so2);
	real fco_int = GreatCircle::final_course(so,interSec);
	real turnDelta_o = Util::turnDelta(fco,fco_int);
	bool before_o = turnDelta_o > Pi/2;
	if (before_o) return LatLonAlt::INVALID();                      // BEFORE [so,so2]
	real gso = GreatCircle::distance(so,so2);
	real oFrac = GreatCircle::distance(so,interSec)/gso;
	if (oFrac > 1) return LatLonAlt::INVALID();                     // AFTER [so,so2]
	real fci = GreatCircle::final_course(si,si2);
	real fci_int = GreatCircle::final_course(si,interSec);
	real turnDelta_i = Util::turnDelta(fci,fci_int);
	bool before_i = turnDelta_i > Pi/2;                 // BEFORE [si,si2]
	if (before_i) return LatLonAlt::INVALID();
	real gsi = GreatCircle::distance(si,si2);
	real iFrac = GreatCircle::distance(si,interSec)/gsi;
	if (iFrac > 1) return LatLonAlt::INVALID();                    // AFTER [si,si2]
    return interSec;
}

LatLonAlt GreatCircle::intersectionSegment(real T, const LatLonAlt& so, const Velocity& vo, const LatLonAlt& si, const LatLonAlt& si2) {
	LatLonAlt so2 = linear_initial(so, vo,T);
	return intersectSegments(so, so2, si, si2);
}


/**
 * Given two great circles defined by so, so2 and si, si2 return the intersection point that is closest to so.
 * (Note. because on a sphere there are two intersection points)
 *  Calculate altitude of intersection using linear extrapolation from line (so,so2)
 *
 * @param so     first point of line o
 * @param so2    second point of line o
 * @param dto    the delta time between point so and point so2.
 * @param si     first point of line i
 * @param si2    second point of line i
 * @return a pair: intersection point and the delta time from point "so" to the intersection, can be negative if intersect
 *                 point is in the past. If intersection point is invalid then the returned delta time is -1
 */
std::pair<LatLonAlt,real>  GreatCircle::intersectionExtrapAlt(const LatLonAlt& so, const LatLonAlt& so2, real dto, const LatLonAlt& si, const LatLonAlt& si2) {
	LatLonAlt lgc = GreatCircle::intersection(so, so2, si, si2);
	if (lgc.isInvalid()) return std::pair<LatLonAlt,real>(lgc,-real(1.0));
	real gso = distance(so,so2)/dto;
	real intTm = distance(so,lgc)/gso;  // relative to so
	bool behind = GreatCircle::behind(lgc, so, GreatCircle::velocity_average(so, so2, real(1.0)));
	if (behind) intTm = -intTm;
	// compute a better altitude
	real vs = (so2.alt() - so.alt())/dto;
	real nAlt = so.alt() + vs*(intTm);
	LatLonAlt pgc = LatLonAlt::mk(lgc.lat(),lgc.lon(),nAlt);
	return std::pair<LatLonAlt,real>(pgc,intTm);
}

std::pair<LatLonAlt,real>  GreatCircle::intersectionAvgAlt(const LatLonAlt& so, const LatLonAlt& so2, real dto, const LatLonAlt& si, const LatLonAlt& si2) {
	LatLonAlt interSec = GreatCircle::intersection(so, so2, si, si2);
	//f.pln(" %%% GreatCircle.intersection: lgc = "+lgc.toString(15));
	if (interSec.isInvalid()) return std::pair<LatLonAlt,real>(interSec,-real(1.0));
	real gso = distance(so,so2)/dto;
	real intTm = distance(so,interSec)/gso;  // relative to so
	//f.pln(" ## gso = "+Units.str("kn", gso)+" distance(so,lgc) = "+Units.str("NM",distance(so,lgc)));
	bool behind = GreatCircle::behind(interSec, so, GreatCircle::velocity_average(so, so2, real(1.0))); //TODO: initial?
	//		f.pln("behind="+behind+" interSec="+interSec+" so="+so+" vo="+GreatCircle.velocity_average(so, so2, 1.0));
	if (behind) intTm = -intTm;
	// compute a better altitude using average of near points
	real do1 = distance(so,interSec);
	real do2 = distance(so2,interSec);
	real alt_o = so.alt();
	if (do2 < do1) alt_o = so2.alt();
	real di1 = distance(si,interSec);
	real di2 = distance(si2,interSec);
	real alt_i = si.alt();
	if (di2 < di1) alt_i = si2.alt();
	real nAlt = (alt_o + alt_i)/real(2.0);
	//    	f.pln(" $$ LatLonAlt.intersection: so.alt() = "+Units.str("ft",so.alt())+" so2.alt() = "+Units.str("ft",so2.alt())+
	//    			" si.alt() = "+Units.str("ft",si.alt())+" si2.alt() = "+Units.str("ft",si2.alt())+
	//    			" nAlt() = "+Units.str("ft",nAlt));
	//f.pln(" $$ LatLonAlt.intersection: intTm = "+intTm+" vs = "+Units.str("fpm",vs)+" nAlt = "+Units.str("ft",nAlt)+" "+behind);
	LatLonAlt pgc = LatLonAlt::mk(interSec.lat(),interSec.lon(),nAlt);
	return std::pair<LatLonAlt,real>(pgc,intTm);

}

std::pair<LatLonAlt,real> GreatCircle::intersection(const LatLonAlt& so, const Velocity& vo, const LatLonAlt& si, const Velocity& vi) {
	LatLonAlt so2 = linear_initial(so, vo, 1000);
	LatLonAlt si2 = linear_initial(si, vi, 1000);
	LatLonAlt i = intersection(so, so2, si, si2);
	if (i.isInvalid()) return std::pair<LatLonAlt,real>(LatLonAlt::INVALID(),-real(1.0)); // collinear or (nearly) same position or cross in the past
	real dt = distance(so,i)/vo.gs();
	if (behind(i, so, vo)) dt = -dt;
	return std::pair<LatLonAlt,real>(i,dt);
}

real GreatCircle::angleBetween(const LatLonAlt& a1, const LatLonAlt& a2, const LatLonAlt& b1, const LatLonAlt& b2) {
	Vect3 va = spherical2xyz(a1.lat(), a1.lon()).cross(spherical2xyz(a2.lat(), a2.lon())).Hat(); // normal 1
	Vect3 vb = spherical2xyz(b1.lat(), b1.lon()).cross(spherical2xyz(b2.lat(), b2.lon())).Hat(); // normal 2
	real cosx = va.dot(vb);
	return Util::acos_safe(cosx);
}


real GreatCircle::angle_betweenOLD(const LatLonAlt& a, const LatLonAlt& b, const LatLonAlt& c) {
	real a1 = angular_distance(c,b);
	real b1 = angular_distance(a,c);
	real c1 = angular_distance(b,a);
	real d = std::sin(c1)*std::sin(a1);
	if (d == real(0.0)) {
		return real(M_PI);
	}
	return Util::acos_safe( (std::cos(b1)-std::cos(c1)*std::cos(a1)) / d );
}

real GreatCircle::angle_between(const LatLonAlt& a, const LatLonAlt& b, const LatLonAlt& c) {
//	double ang1 = initial_course(b,a);
//	double ang2 = initial_course(b,c);
//	double angleBetw = Util::turnDelta(ang1,ang2);
//	return angleBetw;
	return angleBetween(b,a,b,c);
}

real GreatCircle::angle_between(const LatLonAlt& a, const LatLonAlt& b, const LatLonAlt& c, int dir) {
	real trk1 = GreatCircle::final_course(a, b);
	real trk2 = GreatCircle::initial_course(b,c);
	int calcdir = Util::turnDir(trk1, trk2);
		real theta = angleBetween(b,a,b,c);
		if (dir==calcdir) return theta;
		return 2*Pi - theta; // note that this is not quite the same as theta+PI
}


/**
 * Return true if x is "behind" ll, considering its current direction of travel, v.
 * "Behind" here refers to the hemisphere aft of ll.
 * That is, x is within the region behind the perpendicular line to v through ll.
 * @param ll aircraft position
 * @param v aircraft velocity
 * @param x intruder positino
 * @return
 */
bool GreatCircle::behind(const LatLonAlt& x, const LatLonAlt& ll, const Velocity& v) {
	Velocity v2 = velocity_initial(ll, x, 100);
	return Util::turnDelta(v.trk(), v2.trk()) > real(M_PI)/real(2.0);
}

/**
 * Returns values describing if the ownship state will pass in front of or behind the intruder (from a horizontal perspective)
 * @param so ownship position
 * @param vo ownship velocity
 * @param si intruder position
 * @param vi intruder velocity
 * @return 1 if ownship will pass in front (or collide, from a horizontal sense), -1 if ownship will pass behind, 0 if collinear or parallel or closest intersection is behind you
 */
int GreatCircle::passingDirection(const LatLonAlt& so, const Velocity& vo, const LatLonAlt& si, const Velocity& vi) {
	std::pair<LatLonAlt,real> p = intersection(so,vo,si,vi);
	if (p.second < 0) return 0;
	LatLonAlt si3 = linear_initial(si,vi,p.second); // intruder position at time of intersection
	if (behind(p.first, si3, vi)) return -1;
	return 1;
}

//    int dirForBehind(LatLonAlt so, Velocity vo, LatLonAlt si, Velocity vi) {
//    	LatLonAlt so2 = linear_initial(so, vo, 1000);
//    	LatLonAlt si2 = linear_initial(si, vi, 1000);
//    	LatLonAlt i = intersection(so, so2, si, si2);
//    	if (i.isInvalid() || behind(so,vo,i) || behind(si,vi,i)) return 0; // collinear or (nearly) same position or cross in the past
//    	double tso = distance(so,i)/vo.gs();
//    	if (behind(so,vo,i)) tso = -tso;
//    	LatLonAlt siXP = linear_initial(si,vi,tso);
//    	int ahead = (behind(siXP,vi,i) ? -1 : 1);
//    	int onRight = Util::sign(cross_track_distance(so,i,siXP));
//fpln("ahead="+ahead+"  onRight="+onRight+" siXP="+siXP.toStringNP(8)+" i="+i.toStringNP(8));
//    	return ahead*onRight;
//    }


int GreatCircle::dirForBehind(const LatLonAlt& so, const Velocity& vo, const LatLonAlt& si, const Velocity& vi) {
	LatLonAlt so2 = linear_initial(so, vo, 1000);
	LatLonAlt si2 = linear_initial(si, vi, 1000);
	LatLonAlt i = intersection(so, so2, si, si2);
	if (i.isInvalid() || behind(i,so,vo) || behind(i,si,vi)) return 0; // collinear or (nearly) same position or cross in the past
	int onRight = Util::sign(cross_track_distance(si,si2,so));
	return -onRight;
}



Velocity GreatCircle::velocity_initial(const LatLonAlt& p1, const LatLonAlt& p2, real t) {
	// p1 is the source position, p2 is another point to form a great circle
	// positive time is moving from p1 toward p2
	// negative time is moving from p1 away from p2
	if (std::abs(t) < minDt || Util::almost_equals(std::abs(t) + minDt, minDt,
			PRECISION7)) {
		// time is negative or very small (less than 1 ms)
		return Velocity::ZERO();
	}
	real d = angular_distance(p1, p2);
	if (Constants::almost_equals_radian(d)) {
		if (Constants::almost_equals_alt(p1.alt(), p2.alt())) {
			// If the two points are about 1 meter apart, then count them as
			// the same.
			return Velocity::ZERO();
		} else {
			return Velocity::ZERO().mkVs((p2.alt() - p1.alt()) / t);
		}
	}
	real gs = GreatCircle::distance_from_angle(d, real(0.0)) / t;
	real crs = initial_course_impl(p1, p2); //, d);
	return Velocity::mkTrkGsVs(crs, gs, (p2.alt() - p1.alt()) / t);
}

Velocity GreatCircle::velocity_average(const LatLonAlt& p1, const LatLonAlt& p2, real t) {
	// p1 is the source position, p2 is another point on that circle
	// positive time is moving from p1 toward p2
	// negative time is moving from p1 away from p2
	if (t >= real(0.0)) {
		return GreatCircle::velocity_initial(GreatCircle::interpolate(p1, p2, real(0.5)), p2, t / real(2.0));
	} else {
		return GreatCircle::velocity_average(p1, GreatCircle::interpolate(p1, p2, -real(1.0)), -t);
	}
}

Velocity GreatCircle::velocity_average_speed(const LatLonAlt& s1, const LatLonAlt& s2, real speed) {
	real dist = GreatCircle::distance(s1, s2);
	real dt = dist/speed;
	return GreatCircle::velocity_average(s1, s2, dt);
}


Velocity GreatCircle::velocity_final(const LatLonAlt& p1, const LatLonAlt& p2, real t) {
	// p1 is the source position, p2 is another point on that circle
	// positive time is moving from p1 toward p2
	// negative time is moving from p1 away from p2 (final velocity is the
	// velocity at that time)
	if (t >= real(0.0)) {
		return GreatCircle::velocity_initial(p2, p1, -t);
	} else {
		return GreatCircle::velocity_initial(GreatCircle::interpolate(p1, p2, -real(1.0)), p1, t);
	}
}



Vect3 GreatCircle::spherical2xyz(real lat, real lon) {
	real r = GreatCircle::spherical_earth_radius;
	// convert latitude to 0-PI
	real theta = real(M_PI)/2 - lat;
	real phi = lon; //M_PI - lon;
	real x = r*std::sin(theta)*std::cos(phi);
	real y = r*std::sin(theta)*std::sin(phi);
	real z = r*std::cos(theta);
	return Vect3(x,y,z);
}

Vect3 GreatCircle::spherical2xyz(const LatLonAlt& lla) {
	return spherical2xyz(lla.lat(), lla.lon());
}

LatLonAlt GreatCircle::xyz2spherical(const Vect3& v) {
	real r = GreatCircle::spherical_earth_radius;
	real theta = Util::acos_safe(v.z()/r);
	real phi = Util::atan2_safe(v.y(), v.x());
	real lat = real(M_PI)/2 - theta;
	real lon = Util::to_pi(phi); //M_PI - phi);
	return LatLonAlt::mk(lat, lon, 0);
}




real GreatCircle::chord_distance(real lat1, real lon1, real lat2, real lon2) {
	Vect3 v1 = spherical2xyz(lat1,lon1);
	Vect3 v2 = spherical2xyz(lat2,lon2);
	return v1.Sub(v2).norm();
}


real GreatCircle::chord_distance(real surface_dist) {
	real theta = angle_from_distance(surface_dist,real(0.0));
	return real(2.0)*std::sin(theta/real(2.0))*GreatCircle::spherical_earth_radius;
}

real GreatCircle::chord_distance(const LatLonAlt& lla1, const LatLonAlt& lla2) {
	return chord_distance(lla1.lat(), lla1.lon(), lla2.lat(), lla2.lon());
}



real GreatCircle::surface_distance(real chord_distance) {
	real theta = real(2.0)*Util::asin_safe(chord_distance*real(0.5) / GreatCircle::spherical_earth_radius);
	return distance_from_angle(theta,real(0.0));
}

real  GreatCircle::to_chordal_radius(real surface_radius) {
	return chord_distance(surface_radius*real(2.0))/real(2.0);
}

real  GreatCircle::to_surface_radius(real chord_radius) {
	return surface_distance(chord_radius*real(2.0))/real(2.0);
}



LatLonAlt GreatCircle::small_circle_rotation(const LatLonAlt& so, const LatLonAlt& center, real angle) {
	if (Util::almost_equals(angle, 0)) return so;
	real R = angular_distance(so, center);
	Triple<real,real,real>ABc = side_angle_side(R, angle, R);
	real A = ABc.first;
	real c = distance_from_angle(ABc.third, real(0.0));
	real crs = initial_course(so, center);
	if (crs > real(M_PI)) crs = crs-2*real(M_PI);
	real trk = Util::to_2pi(crs - A);
	LatLonAlt ret = linear_initial(so, trk, c);
	return ret;
}

/**
 * Accurately calculate the linear distance of an arc on a small circle (turn) on the sphere.
 * @param radius along-surface radius of small circle
 * @param arcAngle angular (radian) length of the arc.  This is the angle between two great circles that intersect at the small circle's center.
 * @return linear distance of the small circle arc
 * Note: A 100 km radius turn over 60 degrees produces about 4.3 m error.
 */
real GreatCircle::small_circle_arc_length(real radius, real arcAngle) {
	// use the chord of the diameter to determine the radius in the ECEF Euclidean frame
	real r2 = chord_distance(radius*2)/2;
	// because this is a circle in a Euclidean frame, use the normal calculations
	return arcAngle*r2;
}

/**
 * Accurately calculate the angular distance of an arc on a small circle (turn) on the sphere.
 * @param radius along-surface radius of small circle
 * @param arcLength linear (m) length of the arc.  This is the along-line length of the arc.
 * @return Angular distance of the arc around the small circle (from 0 o 2pi)
 * Note: A 100 km radius turn over 100 km of turn produces about 0.0024 degrees of error.
 */
real GreatCircle::small_circle_arc_angle(real radius, real arcLength) {
	// use the chord of the diameter to determine the radius in the ECEF Euclidean frame
	real r2 = chord_distance(radius*2)/2;
	if (r2 == real(0.0)) return real(0.0);
	// because this is a circle in a Euclidean frame, use the normal calculations
	return arcLength/r2;
}

}
