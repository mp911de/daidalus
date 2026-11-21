/* KinematicsLatLon
 *
 * Authors:  Ricky Butler              NASA Langley Research Center
 *           George Hagen              NASA Langley Research Center
 *           Jeff Maddalon             NASA Langley Research Center
  *
 * Copyright (c) 2011-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */

#include "Real.h"
#include "Util.h"
#include "Units.h"
#include "Vect2.h"
#include "Vect3.h"
#include "Triple.h"
#include "Velocity.h"
#include "Kinematics.h"
#include "StateVector.h"
#include "Constants.h"
#include "VectFuns.h"
#include "KinematicsLatLon.h"
#include "GreatCircle.h"
#include "Quad.h"
#include "Tuple5.h"
#include "format.h"
#include <cmath>
#include <float.h>


namespace larcfm {


   bool KinematicsLatLon::chordalSemantics = false;

std::pair<LatLonAlt,Velocity> KinematicsLatLon::linear(const LatLonAlt& so, const Velocity& vo, real t) {
	LatLonAlt sn = GreatCircle::linear_initial(so, vo, t);
	return std::pair<LatLonAlt,Velocity>(sn,vo);
}

std::pair<LatLonAlt,Velocity> KinematicsLatLon::linear(std::pair<LatLonAlt,Velocity> sv0, real t) {
	LatLonAlt s0 = sv0.first;
	Velocity v0 = sv0.second;
	//fpln(" $$$$$$$ in linear at time: "+t);
	return std::pair<LatLonAlt,Velocity>(linear(s0,v0,t).first,v0);
}

  real KinematicsLatLon::turnRadiusByRate(real speed, real omega) {
    real R = Kinematics::turnRadiusByRate(speed, omega);
    return GreatCircle::surface_distance(R*2)/real(2.0);
  }
  
std::pair<LatLonAlt,Velocity> KinematicsLatLon::turnOmega(const LatLonAlt& so, const Velocity& vo, real t, real omega) {
	real currentTrk = vo.trk();
	real perpTrk;
	if (omega > 0) perpTrk = currentTrk+real(M_PI)/2;
	else perpTrk = currentTrk-real(M_PI)/2;
	real radius = Kinematics::turnRadiusByRate(vo.gs(), omega);
	LatLonAlt center = GreatCircle::linear_initial(so, perpTrk, radius);
	real alpha = omega*t;
	real vFinalTrk = GreatCircle::initial_course(center,so);
	real nTrk = vFinalTrk + alpha;
	LatLonAlt sn = GreatCircle::linear_initial(center, nTrk, radius);
	sn = sn.mkAlt(so.alt() + vo.z()*t);
	//double finalTrk = currentTrk+theta;
	real final_course = GreatCircle::final_course(center,sn);   // TODO: THIS IS PROBABLY BETTER
	real finalTrk = final_course + Util::sign(omega)*real(M_PI)/2;
	Velocity vn = vo.mkTrk(finalTrk);
	return std::pair<LatLonAlt,Velocity>(sn,vn);
}

std::pair<LatLonAlt,Velocity> KinematicsLatLon::turnOmega(std::pair<LatLonAlt,Velocity> sv0, real t, real omega) {
	if (Util::almost_equals(omega,0))
		return linear(sv0,t);
	LatLonAlt s0 = sv0.first;
	Velocity v0 = sv0.second;
	return turnOmega(s0,v0,t,omega);
}

//  std::pair<LatLonAlt,Velocity> KinematicsLatLon::turnRadius(const LatLonAlt& so, const Velocity& vo, double t, double signedRadius) {
//    double currentTrk = vo.trk();
//    double dir = Util::sign(signedRadius);
//    double perpTrk = currentTrk+dir*M_PI/2;
//    double radius = std::abs(signedRadius);
//    LatLonAlt center = GreatCircle::linear_initial(so, perpTrk, radius);
//    double pathDist = vo.gs()*t;
//    double theta = pathDist/radius;
//    //TODO: theta = pathDist/radius, assumes radius is a chord radius, when it is actually a great circle radius.
//    //      for small distances the difference is not that big, but still...
//    //Note: The other problem is that this assumes a constant speed and constant ground speed through the turn.
//    //      this may or may not be true.
//    double vFinalTrk = GreatCircle::initial_course(center,so);
//    double nTrk = vFinalTrk + dir*theta;
//    LatLonAlt sn = GreatCircle::linear_initial(center, nTrk, radius);
//    sn = sn.mkAlt(so.alt() + vo.z*t);
//    double final_course = GreatCircle::final_course(center,sn);
//    double finalTrk = final_course + dir*M_PI/2;
//    Velocity vn = vo.mkTrk(finalTrk);
//    return std::pair<LatLonAlt,Velocity>(sn,vn);
//  }
//

  
std::pair<LatLonAlt,Velocity> KinematicsLatLon::turnByDist2D(const LatLonAlt& so, const LatLonAlt& center, int dir, real d, real gsAtd) {
	LatLonAlt sn = turnByDist2D(so, center, dir, d);
	real final_course = GreatCircle::final_course(center,sn);
	real finalTrk = final_course + dir*real(M_PI)/2;
    Velocity vn = Velocity::mkTrkGsVs(finalTrk,gsAtd,real(0.0));
	return std::pair<LatLonAlt,Velocity>(sn,vn);
}

LatLonAlt KinematicsLatLon::turnByDist2D(const LatLonAlt& so, const LatLonAlt& center, int dir, real d) {
    real R = GreatCircle::distance(so, center);
    real theta;
    if (chordalSemantics) {
        real chordalRadius = GreatCircle::to_chordal_radius(R);
    	theta = dir*d/chordalRadius;
    } else {
    	theta = dir*d/R;
    }
	real vFinalTrk = GreatCircle::initial_course(center,so);
	real nTrk = vFinalTrk + theta;
	LatLonAlt sn = GreatCircle::linear_initial(center, nTrk, R);
	sn = sn.mkAlt(real(0.0));
	return sn;
}

LatLonAlt KinematicsLatLon::turnByAngle2D(const LatLonAlt& so, const LatLonAlt& center, real alpha) {
    real R = GreatCircle::distance(so, center);
	real trkFromCenter = GreatCircle::initial_course(center,so);
	real nTrk = trkFromCenter + alpha;
	LatLonAlt sn = GreatCircle::linear_initial(center, nTrk, R);
	sn = sn.mkAlt(real(0.0));
	return sn;
}

std::pair<LatLonAlt,Velocity> KinematicsLatLon::turn(LatLonAlt s0, Velocity v0, real t, real R, bool turnRight) {
	if (Util::almost_equals(R,0))
		return std::pair<LatLonAlt,Velocity>(s0,v0);
	int dir = -1;
	if (turnRight) dir = 1;
	real omega = dir*v0.gs()/R;
	return turnOmega(s0,v0,t,omega);
}

std::pair<LatLonAlt,Velocity> KinematicsLatLon::turn(std::pair<LatLonAlt,Velocity> sv0, real t, real R, bool turnRight) {
	return turn(sv0.first,sv0.second,t,R,turnRight);
}


std::pair<LatLonAlt,Velocity> KinematicsLatLon::turnUntilTimeOmega(std::pair<LatLonAlt,Velocity> svo, real t, real turnTime, real omega) {
	std::pair<LatLonAlt,Velocity> tPair;
	if (t <= turnTime) {
		tPair = turnOmega(svo, t, omega);
	} else {
		tPair = turnOmega(svo, turnTime, omega);
		tPair = linear(tPair,t-turnTime);
	}
	return tPair;
}

std::pair<LatLonAlt,Velocity> KinematicsLatLon::turnUntilTimeOmega(const LatLonAlt& so, const Velocity& vo, real t, real turnTime, real omega) {
	return turnUntilTimeOmega(std::pair<LatLonAlt,Velocity>(so,vo), t, turnTime, omega);
}


std::pair<LatLonAlt,Velocity> KinematicsLatLon::turnUntil(std::pair<LatLonAlt,Velocity> svo, real t, real goalTrack, real maxBank) {
	real omega = Kinematics::turnRateGoal(svo.second, goalTrack, maxBank);
	real turnTime = Kinematics::turnTime(svo.second, goalTrack, maxBank);
	return turnUntilTimeOmega(svo,t,turnTime,omega);
}


std::pair<LatLonAlt,Velocity> KinematicsLatLon::turnUntil(const LatLonAlt& so, const Velocity& vo, real t, real goalTrack, real maxBank) {
	return turnUntil(std::pair<LatLonAlt,Velocity>(so,vo), t, goalTrack, maxBank);
}

std::pair<LatLonAlt,Velocity> KinematicsLatLon::turnUntilTimeRadius(const std::pair<LatLonAlt,Velocity>& svo, real t, real turnTime, real R, bool turnRight) {
	std::pair<LatLonAlt,Velocity> tPair;
	if (t <= turnTime) {
		tPair = turn(svo, t, R, turnRight);
	} else {
		tPair = turn(svo, turnTime, R, turnRight);
		tPair = linear(tPair,t-turnTime);
	}
	return tPair;
}


LatLonAlt KinematicsLatLon::center(const LatLonAlt& s0, real trk, real radius, int dir) {
	  return GreatCircle::linear_initial(s0, trk+dir*real(M_PI)/2, radius);
}


LatLonAlt KinematicsLatLon::center(const LatLonAlt& s0, const Velocity& v0, real omega) {
	  real v = v0.gs();
    real R = v/omega;
    return center(s0, v0.trk(), R, Util::sign(omega));
}

real KinematicsLatLon::closestTimeOnTurn(const LatLonAlt& s0, const Velocity& v0, real omega, const LatLonAlt& x, real endTime) {
	LatLonAlt cent = center(s0,v0,omega);
	if (GreatCircle::almostEquals(x.mkAlt(0), cent.mkAlt(0))) return -real(1.0);
	real ang1 = GreatCircle::initial_course(cent,s0);
	real ang2 = GreatCircle::initial_course(cent,x);
	real delta = Util::turnDelta(ang1, ang2, Util::sign(omega));
	real t = std::abs(delta/omega);
	if (endTime > 0 && (t < 0 || t > endTime)) {
		real maxTime = 2*real(M_PI)/std::abs(omega);
		if (t > (maxTime + endTime) / 2) {
			return real(0.0);
		} else {
			return endTime;
		}
	}
	return t;

}


real KinematicsLatLon::closestDistOnTurn(const LatLonAlt& s0, const Velocity& v0, real R, int dir, const LatLonAlt& x, real maxDist) {
	LatLonAlt cent = center(s0, v0.trk(), R, dir);
	if (GreatCircle::almostEquals(x.mkAlt(0),cent.mkAlt(0))) return -real(1.0);
	real ang1 = GreatCircle::initial_course(cent,s0);
	real ang2 = GreatCircle::initial_course(cent,x);
	real delta = Util::turnDelta(ang1, ang2, dir);
	real t = GreatCircle::small_circle_arc_length(R, delta);
	if (maxDist > 0 && (t < 0 || t > maxDist)) {
		real maxD = 2*real(M_PI)*R;
		if (t > (maxD + maxDist) / 2) {
			return real(0.0);
		} else {
			return maxDist;
		}
	}
	return t;
}


std::pair<LatLonAlt,Velocity> KinematicsLatLon::gsAccel(const LatLonAlt& so, const Velocity& vo,  real t, real a) {
	real dist = vo.gs()*t + real(0.5)*a*t*t;
	real currentTrk = vo.trk();
	LatLonAlt sn = GreatCircle::linear_initial(so, currentTrk, dist);
	sn = sn.mkAlt(so.alt() + vo.z()*t);
	real vnGs = vo.gs() + a*t;
	Velocity vn = vo.mkGs(vnGs);
	//fpln(" $$$$$ gsAccel: sn = "+sn+" vn = "+vn);
	return std::pair<LatLonAlt,Velocity>(sn,vn);
}


std::pair<LatLonAlt,Velocity> KinematicsLatLon::gsAccelUntil(const LatLonAlt& so, const Velocity& vo, real t, real goalGS, real gsAccel_d) {
	if (gsAccel_d < 0 ) {
		fpln("Kinematics::gsAccelUntil: user supplied negative gsAccel!!");
		gsAccel_d = -gsAccel_d;                              // make sure user supplies positive value
	}
	real accelTime = Kinematics::gsAccelTime(vo,goalGS,gsAccel_d);
	int sgn = 1;
	if (goalGS < vo.gs()) sgn = -1;
	real a = sgn*gsAccel_d;
	std::pair<LatLonAlt, Velocity> nsv = gsAccel(so, vo, accelTime, a);
	if (t<=accelTime) return gsAccel(so, vo, t, a);
	else  return gsAccel(nsv.first, nsv.second, t-accelTime, 0);
}

std::pair<LatLonAlt,Velocity> KinematicsLatLon::gsAccelUntil(std::pair<LatLonAlt,Velocity> sv0, real t, real goalGs, real gsAccel) {
	return gsAccelUntil(sv0.first, sv0.second, t, goalGs, gsAccel);
}

// ****************************** Vertical Speed KINEMATIC CALCULATIONS *******************************

std::pair<LatLonAlt,Velocity> KinematicsLatLon::vsAccel(const LatLonAlt& so, const Velocity& vo,  real t, real a) {
	real dist = vo.gs()*t;
	real currentTrk = vo.trk();
	LatLonAlt sn = GreatCircle::linear_initial(so, currentTrk, dist);
	real nsz = so.alt() + vo.z()*t + real(0.5)*a*t*t;
	sn = sn.mkAlt(nsz);
	Velocity  vn = vo.mkVs(vo.z() + a*t);
	return std::pair<LatLonAlt,Velocity>(sn,vn);
}

std::pair<LatLonAlt,Velocity>  KinematicsLatLon::vsAccel(std::pair<LatLonAlt,Velocity> svo,  real t, real a) {
	return vsAccel(svo.first, svo.second, t, a);
}

std::pair<LatLonAlt,Velocity> KinematicsLatLon::vsAccelUntil(const LatLonAlt& so, const Velocity& vo, real t, real goalVs, real vsAccel_d) {
	if (vsAccel_d < 0 ) {
		fpln("Kinematics::vsAccelUntil: user supplied negative vsAccel!!");
		vsAccel_d = -vsAccel_d;                              // make sure user supplies positive value
	}
	real accelTime = Kinematics::vsAccelTime(vo,goalVs, vsAccel_d);
	int sgn = 1;
	if (goalVs < vo.vs()) sgn = -1;
	//LatLonAlt ns = LatLonAlt.ZERO;
	if (t <= accelTime)
		return vsAccel(so, vo, t, sgn*vsAccel_d);
	else {
		LatLonAlt posEnd = vsAccel(so,vo,accelTime,sgn*vsAccel_d).first;
		Velocity nvo = Velocity::mkVxyz(vo.x(),vo.y(), goalVs);
		return linear(posEnd,nvo,t-accelTime);
	}
}

std::pair<LatLonAlt,Velocity> KinematicsLatLon::vsAccelUntil(std::pair<LatLonAlt,Velocity> sv0, real t, real goalVs, real vsAccel) {
	return vsAccelUntil(sv0.first, sv0.second,t,goalVs, vsAccel);
}

std::pair<LatLonAlt, Velocity> KinematicsLatLon::vsLevelOutCalculation(std::pair<LatLonAlt,Velocity> sv0,
		real targetAlt, real a1, real a2, real t1, real t2, real t3,  real t) {
	LatLonAlt s0 = sv0.first;
	Velocity v0 = sv0.second;
	real soz = s0.alt();
	real voz = v0.z();
	std::pair<real, real> vsL = Kinematics::vsLevelOutCalc(soz,voz, targetAlt, a1, a2, t1, t2, t3, t);
	real nz = vsL.first;
	real nvs = vsL.second;
	Velocity nv = v0.mkVs(nvs);
	LatLonAlt ns = linear(s0,v0,t).first.mkAlt(nz);
	return std::pair<LatLonAlt, Velocity>(ns,nv);
}

std::pair<LatLonAlt, Velocity> KinematicsLatLon::vsLevelOut(std::pair<LatLonAlt, Velocity> sv0, real t, real climbRate,
		real targetAlt, real accelUp, real accelDown, bool allowClimbRateChange) {
	Tuple5<real,real,real,real,real> LevelParams = Kinematics::vsLevelOutTimes(sv0.first.alt(), sv0.second.vs(), climbRate, targetAlt, accelUp, accelDown, allowClimbRateChange);
	return vsLevelOutCalculation(sv0, targetAlt, LevelParams.fourth, LevelParams.fifth, LevelParams.first, LevelParams.second, LevelParams.third, t);
}

std::pair<LatLonAlt, Velocity> KinematicsLatLon::vsLevelOut(std::pair<LatLonAlt, Velocity> sv0, real t, real climbRate,
		real targetAlt, real a, bool allowClimbRateChange) {
	return vsLevelOut(sv0, t, climbRate, targetAlt, a, -a, allowClimbRateChange);
}

std::pair<LatLonAlt, Velocity> KinematicsLatLon::vsLevelOut(std::pair<LatLonAlt, Velocity> sv0, real t, real climbRate,
		real targetAlt, real a) {
	return vsLevelOut(sv0, t, climbRate, targetAlt, a, -a, true);
}

}
