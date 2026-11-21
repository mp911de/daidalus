/* Kinematics
 *
 * Authors:  Ricky Butler              NASA Langley Research Center
 *           George Hagen              NASA Langley Research Center
 *           Jeff Maddalon             NASA Langley Research Center
 *           Cesar Munoz               NASA Langley Research Center
 *           Anthony Narkawicz         NASA Langley Research Center
 *           Aaron Dutle               NASA Langley Research Center
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
#include "Quad.h"
#include "Tuple5.h"
#include "Debug.h"
#include "format.h"
#include <cmath>
#include <float.h>


namespace larcfm {


real Kinematics::turnRadius(real speed, real bank, real g) {
	real abank = std::abs(bank);
	if ( g <= real(0.0) ) {
		g = Units::gn; // I could flag this as an error, but I will assume sea-level gravity.
	}
	abank = Util::min(real(M_PI)/real(2.0),abank);
	if ( abank == real(0.0) ) {
		return MAXDOUBLE;
	}
	return speed*speed/(g*std::tan(abank));
}

real Kinematics::turnRadius(real speed, real bank) {
	return Kinematics::turnRadius(speed, bank, Units::gn);
}


real Kinematics::turnRadiusByRate(real speed, real omega) {
	  if (Util::almost_equals(omega,real(0.0))) return MAXDOUBLE;
	  return std::abs(speed/omega);
 }


real Kinematics::speedOfTurn(real R, real bank) {
	real abank = std::abs(bank);
	abank = Util::min(real(M_PI)/real(2.0),abank);
	R = Util::max(real(0.0), R);
	return std::sqrt(Units::gn*std::tan(abank)*R);
}

real Kinematics::turnRate(real speed, real bankAngle) {
	if (Util::almost_equals(speed,0)) return real(0.0);
	bankAngle = Util::min(real(M_PI)/real(2.0),Util::max(bankAngle, -real(M_PI)/real(2.0)));
	return Units::gn*std::tan(bankAngle)/speed;
}

real Kinematics::turnRateByRadius(real speed, real R) {
	if (Util::almost_equals(R,real(0.0))) return real(0.0);
	return speed/R;
}


real Kinematics::bankAngleFromRadius(real speed, real R) {
	if (R <= real(0.0)) return real(0.0);
	return Util::atan2_safe(speed*speed,(R*Units::gn));
}

real Kinematics::bankAngleFromRadius(real speed, real R, bool turnRight) {
	return Util::sign(turnRight)*bankAngleFromRadius(speed, R);
}

real Kinematics::bankAngleGoal(real track, real goalTrack, real signedBank) {
	return Util::turnDir(track, goalTrack)*signedBank;
}

real Kinematics::bankAngle(real speed, real turnRate) {
	real tanb = turnRate*speed/Units::gn;
	//if (std::abs(tanb) >= M_PI/2) return -1;
	real b = std::atan(tanb);
	return b;
}

real Kinematics::turnRateGoal(const Velocity& vo, real goalTrack, real signedBank) {
	real bankAngle = bankAngleGoal(vo.trk(), goalTrack, signedBank);
	return turnRate(vo.gs(), bankAngle);
}


bool Kinematics::turnDone(real currentTrack, real targetTrack, bool turnRight) {
	if (Util::turnDelta(currentTrack,targetTrack) < real(0.0001)) return true;
	if (turnRight) return !Util::clockwise(currentTrack,targetTrack);
	else return Util::clockwise(currentTrack,targetTrack);
}

real Kinematics::turnTime(const Velocity& v0, real goalTrack, real signedBank, int turnDir) {
	real deltaTrk = Util::turnDelta(v0.trk(), goalTrack, turnDir);
	return turnTime(v0.gs(), deltaTrk, signedBank);
}

real Kinematics::turnTime(const Velocity& v0, real goalTrack, real signedBank, bool turnRight) {
	real deltaTrk = Util::turnDelta(v0.trk(), goalTrack, turnRight);
	return turnTime(v0.gs(), deltaTrk, signedBank);
}

real Kinematics::turnTime(const Velocity& v0, real goalTrack, real signedBank) {
	real deltaAng = Util::signedTurnDelta(v0.trk(), goalTrack);
	return turnTime(v0.gs(), deltaAng, signedBank);
}


real Kinematics::turnTime(real groundSpeed, real deltaTrack, real bankAngle) {
	real omega = Kinematics::turnRate(groundSpeed, bankAngle);
	return turnTime(deltaTrack, omega);
}

real Kinematics::turnTime(real deltaTrack, real omega) {
	if (omega == 0) return MAXDOUBLE;
	return std::abs(deltaTrack/omega);
}


std::pair<Vect3,Velocity> Kinematics::linear(const std::pair<Vect3,Velocity>& sv0, real t) {
	return linear(sv0.first, sv0.second, t);
}

std::pair<Vect3,Velocity> Kinematics::linear(const Vect3& so, const Velocity& vo, real t) {
	return std::pair<Vect3,Velocity>(so.linear(vo.vect3(),t),vo);
}

/**
 * Position/Velocity
 * @param s0          starting position
 * @param center
 * @param d           distance into turn
 * @param gsAt_d
 * @return Position/Velocity after t time
 */
std::pair<Vect3,Velocity> Kinematics::turnByDist2D(const Vect3& so, const Vect3& center, int dir, real d, real gsAt_d) {
    real R = so.distanceH(center);
    if (R==real(0.0)) return std::pair<Vect3,Velocity>(so,Velocity::INVALID());
	real alpha = dir*d/R;
	real trkFromCenter = Velocity::track(center, so);
	real nTrk = trkFromCenter + alpha;
	Vect3 sn = center.linearByDist2D(nTrk, R);
	sn = sn.mkZ(real(0.0));
	real finalTrk = nTrk + dir*real(M_PI)/2;
    Velocity vn = Velocity::mkTrkGsVs(finalTrk,gsAt_d,real(0.0));
	return std::pair<Vect3,Velocity>(sn,vn);
}

Vect3 Kinematics::turnByDist2D(const Vect3& so, const Vect3& center, int dir, real d) {
    real R = so.distanceH(center);
    if (R==real(0.0)) return so;
	real alpha = dir*d/R;
	real trkFromCenter = Velocity::track(center, so);
	real nTrk = trkFromCenter + alpha;
	Vect3 sn = center.linearByDist2D(nTrk, R);
	return sn;
}


Vect3 Kinematics::turnByAngle2D(const Vect3& so, const Vect3& center, real alpha) {
    real R = so.distanceH(center);
	real trkFromCenter = Velocity::track(center, so);
	real nTrk = trkFromCenter + alpha;
	Vect3 sn = center.linearByDist2D(nTrk, R);
	return sn;
}


std::pair<Vect3,Velocity> Kinematics::turnOmega(const Vect3& s0, const Velocity& v0, real t, real omega) {
	if (Util::almost_equals(omega,0)) {
		return linear(s0, v0, t);
	}
	// New implementation avoids calculating track and groundspeed,
	// reduces trig functions to one sine and one cosine.
	Velocity nv = v0.mkAddTrk(omega*t);
	real xT = s0.x() + (v0.y()-nv.y())/omega;
	real yT = s0.y() + (-v0.x()+nv.x())/omega;
	real zT = s0.z() + v0.z()*t;
	Vect3 ns = Vect3(xT,yT,zT);
	return std::pair<Vect3,Velocity>(ns,nv);
}

std::pair<Vect3,Velocity> Kinematics::turn(const Vect3& s0, const Velocity& v0, real t, real R,  int dir) {
	// if (Util::almost_equals(R,0)) {
	// 	return std::pair<Vect3,Velocity>(s0,v0);
	// }
	// double omega = dir*v0.gs()/R;
	real omega = dir*turnRateByRadius(v0.gs(),R);
	return turnOmega(s0,v0,t,omega);
}

std::pair<Vect3,Velocity> Kinematics::turn(const Vect3& so, const Velocity& vo, real t, real R,  bool turnRight) {
	return turn(so, vo, t, R, Util::sign(turnRight));
}


// std::pair<Vect3,Velocity> Kinematics::turn(const Vect3& s0, const Velocity& v0, double t, double bank) {
// 	// if (Util::almost_equals(bank,0)) {
// 	// 	return std::pair<Vect3,Velocity>(s0.linear(v0,t),v0);
// 	// } else {
// 	// 	double R = turnRadius(v0.gs(),bank);
// 	// 	bool turnRight = (bank >= 0);
// 	// 	return turn(s0,v0,t,R,turnRight);
// 	// }
// 	double omega = turnRate(v0.gs(),bank);
// 	return turnOmega(s0,v0,t,omega);
// }

std::pair<Vect3,Velocity>  Kinematics::turnUntilTimeRadius(const Vect3& so, const Velocity& vo, real t, real turnTime, real R, int dir) {
	  std::pair<Vect3,Velocity> tPair;
	  if (t <= turnTime) {
		  tPair = Kinematics::turn(so, vo, t, R, dir);
	  } else {
		  tPair = Kinematics::turn(so, vo, turnTime, R, dir);
		  tPair = linear(tPair,t-turnTime);
	  }
	  return tPair;
}


std::pair<Vect3,Velocity> Kinematics::turnUntilTimeOmega(const Vect3& so, const Velocity& vo, real t, real turnTime, real omega) {
	std::pair<Vect3,Velocity> tPair;
	if (t <= turnTime) {
		tPair = Kinematics::turnOmega(so, vo, t, omega);
	} else {
		tPair = Kinematics::turnOmega(so, vo, turnTime, omega);
		tPair = linear(tPair,t-turnTime);
	}
	return tPair;
}


std::pair<Vect3,Velocity> Kinematics::turnUntilTrack(const Vect3& so, const Velocity& vo, real t, real goalTrack, real signedBank) {
	real omega = turnRateGoal(vo, goalTrack, signedBank);
	real turnTime = Kinematics::turnTime(vo, goalTrack, signedBank);

	return turnUntilTimeOmega(so,vo,t,turnTime,omega);
}



Vect3 Kinematics::turnTrack(const Vect3& so, const Velocity& vo, real goalTrack, real signedBank) {
	real omega = turnRateGoal(vo, goalTrack, signedBank);
	real turnTime = Kinematics::turnTime(vo, goalTrack, signedBank);
	return turnOmega(so, vo, turnTime, omega).first;
}


real Kinematics::turnTimeClosest(const Vect3& s0, const Velocity& v0, real omega, const Vect3& x, real endTime) {
	Vect2 center2 = centerOfTurn(s0,v0,omega);
	if (x.vect2().almostEquals(center2)) return -real(1.0);
	real trk1 = s0.vect2().Sub(center2).trk();
	real trk2 = x.vect2().Sub(center2).trk();
	real delta = Util::turnDelta(trk1, trk2, Util::sign(omega));
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

real Kinematics::turnDistClosest(const Vect3& s0, const Velocity& v0, real R, int dir, const Vect3& x, real maxDist) {
	Vect2 center = centerOfTurn(s0.vect2(), v0.vect2(), R, dir);
	if (x.vect2().almostEquals(center)) return -real(1.0);
	real trk1 = s0.vect2().Sub(center).trk();
	real trk2 = x.vect2().Sub(center).trk();
	real delta = Util::turnDelta(trk1, trk2, dir);
	real d = delta*R;
	if (maxDist > 0 && (d < 0 || d > maxDist)) {
		real maxD = 2*real(M_PI)*R;
		if (d > (maxD + maxDist) / 2) {
			return real(0.0);
		} else {
			return maxDist;
		}
	}
	return d;
	// double omega = turnRateByRadius(v0.gs(),R*dir);
	// if (omega*R == 0.0) return 0.0;
	// double t = turnTimeClosest(s0, v0, omega, x, maxDist/(omega*R));
	// return omega*t*R;
}


/** Q function from ACCoRD.TangentLine. Returns the point on the circle (with 0,0 origin)  that is tangent with the outside point
 * @param s vertex point
 * @param D radius of the circle
 * @param eps direction of turn (+1 is turn LEFT, -1 is turn RIGHT in he absolute frame)
 * @return tangent point on the circle, or an INVALID vector if the vertex is within the circle
 */
Vect2 Q(const Vect2& s, real D, int eps) {
	real sq_s = s.sqv();
	real sq_D = Util::sq(D);
	real delta = sq_s -sq_D;
	if (delta >= 0) {
		real alpha = sq_D/sq_s;
		real beta  = D*Util::sqrt_safe(delta)/sq_s;
		return Vect2(alpha*s.x+eps*beta*s.y,
				alpha*s.y-eps*beta*s.x);
	}
	return Vect2::INVALID();
}


/**
 * true iff the trajectory following "from" to "to" is a left turn
 * "from" is your current vector of travel, "to" is your goal vector of travel
 */
int isRightTurn(const Vect2& from, const Vect2& to) {
	real detv = to.det(from);
	if (detv < 0) return -1;
	return 1;
}


std::pair<Vect2,Vect2> Kinematics::directTo(const Vect2& bot, const Vect2& v0, const Vect2& goal, real R) {
	//f.pln("@@@ p0 = "+f.sStr(p0)+"  v0 = "+f.vStr(v0)+" np = "+f.sStr(np)+"  R = "+Units.str("nm",R));
	Vect2 newV = goal.Sub(bot);					// vector to goal
	int eps = -isRightTurn(v0,newV);	 
	Vect2 vperp;
	if (eps > 0)    // Turn left
		vperp = v0.PerpL().Hat();    // unit velocity vector (perpendicular to initial velocity)
	else
		vperp = v0.PerpR().Hat();    // unit velocity vector (perpendicular to initial velocity)
	Vect2 center = bot.Add(vperp.Scal(R));		// center of turn
	// Shift coordinate system so that center is located at (0,0) Use ACCoRD tangent point Q calculation
	Vect2 s = goal.Sub(center);					// from center to goal
	Vect2 rop = Q(s,R,eps);						// tangent in relative frame (wrt center of circle)
	Vect2 EOT = rop.Add(center);				// return from relative (translate tangent point back to absolute frame)
	return std::pair<Vect2,Vect2>(EOT,center);
}

Quad<Vect3,Velocity,real,int> Kinematics::directToPoint(const Vect3& so, const Velocity& vo, const Vect3& wp, real R){
	Vect2 EOT = directTo(so.vect2(),vo.vect2(),wp.vect2(),R).first;
	if (EOT.isInvalid()) return Quad<Vect3,Velocity,real,int>(Vect3::INVALID(), Velocity::INVALID(), -real(1.0), 0);
	real finalTrack = wp.vect2().Sub(EOT).trk();
	// this should not be based on final track direction, but rather on the actual turn taken.
	//		double turnDelta = Util.signedTurnDelta(vo.trk(),finalTrack);
	real turnDir = isRightTurn(vo.vect2(),wp.Sub(so).vect2());
	real turnDelta = Util::turnDelta(vo.trk(), finalTrack, turnDir > 0);	// angle change in that direction
	real omega = turnDir*vo.gs()/R;
	real turnTime = std::abs(turnDelta/omega);
	std::pair<Vect3,Velocity> p2 = turnOmega(so,vo,turnTime,omega);
	return Quad<Vect3,Velocity,real,int>(p2.first, p2.second, turnTime, (int)turnDir);
}

Triple<Vect3,real,real> Kinematics::genDirectToVertex(const Vect3& so, const Velocity& vo, const Vect3& wp, real bankAngle, real timeBeforeTurn){
	//std::pair<Vect2,Vect2> eot = directTo(Vect2 bot, Vect2 v0, Vect2 goal, double R) {
	Vect3 soPlus = so.Add(vo.vect3().Scal(timeBeforeTurn));
	real R = Kinematics::turnRadius(vo.gs(), bankAngle);
	//public Triple<Vect3,Velocity,double> directToPoint(Vect3 soPlus, Velocity vo, Vect3 wp, double R) {
	Quad<Vect3,Velocity,real,int> dtp = directToPoint(soPlus,vo,wp,R);
	if (dtp.third < 0) {
		return Triple<Vect3,real,real>(Vect3::INVALID(), -real(1.0), -real(1.0));
	}
	Vect3 si = dtp.first;
	Velocity vi = dtp.second;
	std::pair<Vect3,real> ipPair = VectFuns::intersection(soPlus,vo.vect3(),si,vi.vect3());
	if (ipPair.second != ipPair.second) { // NaN
		return Triple<Vect3,real,real>(Vect3::INVALID(), -real(1.0), -real(1.0));
	}
	Vect3 ip = ipPair.first;
	return Triple<Vect3,real,real>(ip,ipPair.second+timeBeforeTurn,dtp.third+timeBeforeTurn);
}



Vect2 Kinematics::centerOfTurn(const Vect2& so, const Vect2& vo, real R, int dir) {
	Vect2 vperp;
	if (dir > 0)    // turn Right
		vperp = vo.PerpR().Hat();    // unit velocity vector
	else
		vperp = vo.PerpL().Hat();    // unit velocity vector
	Vect2 center = so.Add(vperp.Scal(R));
	//fpln("%% center = "+f.sStr(center));
	return center;
}

Vect2 Kinematics::centerOfTurn(const Vect3& so, const Velocity& vo, real omega) {
	real v = vo.gs();
	real theta = vo.trk();
	real R = v/omega;
	return Vect2(so.x() + R*std::cos(theta),so.y() - R*std::sin(theta));
}



bool Kinematics::testLoSTrk(const Vect3& so, const Velocity& vo, const Velocity& nvo, const Vect3& si, const Velocity& vi,
		real bankAngOwn, real stopTime, real D, real H) {
	real step = real(1.0);
	bool rtn = false;
	for (real t = 0; t < stopTime; t = t + step) {
		Vect3 soAtTm = turnUntilTrack(so, vo, t, nvo.trk(), bankAngOwn).first;
		Vect3 siAtTm = si.linear(vi.vect3(),t);
		real distH = soAtTm.Sub(siAtTm).vect2().norm();
		real distV = std::abs(soAtTm.Sub(siAtTm).z());
		if (distH < D && distV< H) {
			rtn =true;
			break;
		}
	}
	return rtn;
}

// ****************************** Ground Speed KINEMATIC CALCULATIONS *******************************


Vect3 Kinematics::gsAccelPos(const Vect3& so3, const Velocity& vo3, real tm, real a) {
	Vect2 so = so3.vect2();
	Vect2 sK = so.Add(vo3.Hat2D().Scal(vo3.gs()*tm+real(0.5)*a*tm*tm));
	real nz = so3.z() + vo3.z()*tm;
	return Vect3(sK,nz);
}

std::pair<Vect3,Velocity> Kinematics::gsAccel(const Vect3& so, const Velocity& vo,  real t, real a) {
	real nvoGs = vo.gs() + a*t;
	Velocity nvo = vo.mkGs(nvoGs);
	return std::pair<Vect3,Velocity>(gsAccelPos(so,vo,t,a),nvo);
}


real Kinematics::accelTime(real gs0,real goalGs, real gsAccel) {
	real deltaGs = gs0 - goalGs;
	if (deltaGs == real(0.0) || gsAccel == real(0.0)) return 0;
	real rtn = std::abs(deltaGs/gsAccel);
	return rtn;
}

real Kinematics::gsAccelTime(const Velocity& vo,real goalGs, real gsAccel) {
    return accelTime(vo.gs(),goalGs,gsAccel);
}


Triple<Vect3,Velocity,real> Kinematics::gsAccelGoal(const Vect3& so, const Velocity& vo, real goalGs, real gsAccel) {
	real accelTime = gsAccelTime(vo, goalGs, gsAccel);
	gsAccel = std::abs(gsAccel);
	int sgn = 1;
	if (goalGs < vo.gs()) sgn = -1;
	Vect3 nso = gsAccelPos(so, vo, accelTime, sgn*gsAccel);
	Velocity nvo = vo.mkGs(goalGs);
	return Triple<Vect3,Velocity,real>(nso,nvo,accelTime);
}


std::pair<Vect3,Velocity> Kinematics::gsAccelUntil(const Vect3& so, const Velocity& vo, real t, real goalGs, real gsAccel_d) {
	real accelTime = gsAccelTime(vo,goalGs,gsAccel_d);
	gsAccel_d = std::abs(gsAccel_d);
	int sgn = 1;
	if (goalGs < vo.gs()) sgn = -1;
	if (t <= accelTime) {
		return gsAccel(so, vo, t, sgn*gsAccel_d);
	} else {
		std::pair<Vect3,Velocity> nsv = gsAccel( so, vo, accelTime, sgn*gsAccel_d);
		return linear(nsv.first, nsv.second, t-accelTime);
	}
}



// std::pair<Vect3,Velocity> Kinematics::gsAccelUntil(const std::pair<Vect3,Velocity>& sv0, double t, double goalGs, double gsAccel) {
// 	return gsAccelUntil(sv0.first, sv0.second, t, goalGs, gsAccel);
// }

std::pair<real,real> Kinematics::accelToDist(real gsIn, real dist, real gsAccel) {
	// if (gsIn < 0 || dist < 0 || (gsAccel < 0 && dist < -0.5*gsIn*gsIn/gsAccel)) {
	// 	return std::pair<double,double>(0.0,-1.0);
	// }

	// double A = 0.5*gsAccel;
	// double B = gsIn;
	// double C = -dist;

	// double ta = (-B+std::sqrt(B*B-4*A*C))/(2*A); // try first root
	// double tb = (-B-std::sqrt(B*B-4*A*C))/(2*A);
	// double t = -1;
	// if (ta >= 0) {
	// 	t = ta;
	// } else if (tb >= 0) {
	// 	t = tb;
	// } else {
	// 	return std::pair<double,double>(0.0, -1.0);
	// }
	// return std::pair<double,double>(gsIn+gsAccel*t, t);

	if (gsIn < 0 || dist < 0 || (gsAccel < 0 && dist < -real(0.5)*gsIn*gsIn/gsAccel)) {
		return std::pair<real,real>(real(0.0),-real(1.0));
	}

        real ta = Util::root(real(0.5) * gsAccel, gsIn, -dist, 1);

	if (ta >= 0) {
		return std::pair<real,real>(gsIn+gsAccel*ta, ta);
	} else {
                real tb = Util::root(real(0.5) * gsAccel, gsIn, -dist, -1);
			
		if (tb >= 0) {
			return std::pair<real,real>(gsIn+gsAccel*tb, tb);
		} else {
			return std::pair<real,real>(real(0.0), -real(1.0)); // current ground speed turns negative before distance is reached
		}
	} 
}


// double Kinematics::distanceToGsAccelTime(double gs, double gsAccel, double dist) {
// 	double t1 = Util::root(0.5 * gsAccel, gs, -dist, 1);
// 	double t2 = Util::root(0.5 * gsAccel, gs, -dist, -1);
// 	//f.pln("$$$ Kinematics.distanceToGsAccelTime t1="+t1+" t2="+t2);		
// 	double dt = ISNAN(t1) || t1 < 0 ? t2 : (ISNAN(t2) || t2 < 0 ? t1 : Util::min(t1, t2));
// 	return dt;
// }


real Kinematics::accel(real gs1, real gs2, real a) {
	// if (gs1 == gs2 || a == 0.0) return 0.0;
	// int sign = Util::sign(gs2 - gs1);
	// double t = (gs2 - gs1) / (sign * a);
	// return gs1 * t + sign * 0.5 * a * t * t;

	if (a == real(0.0)) return real(0.0);
	real accelTime = std::abs((gs1 - gs2) / a);
	real neededDist = accelTime * (gs1 + gs2) / real(2.0);
	return neededDist;
}

std::pair<real, real> Kinematics::accelUntil(real gs0, real gsTarget, real gsAccel, real dt) {
	real ds = -1;
	if (Util::almost_equals(gsAccel, 0)) {
		return std::pair<real, real>(gs0 * dt, gs0);
	}
	real deltaGs = gsTarget - gs0;
	real t0 = std::abs(deltaGs / gsAccel);  // time to reach gsTarget
	real a = Util::sign(deltaGs) * std::abs(gsAccel);	// sign of acceleration
	real gsFinal = gsTarget;
	if (dt < t0) {
		ds = gs0 * dt + real(0.5) * a * dt * dt;
		gsFinal = gs0 + a * dt;
	}
	else ds = gs0 * t0 + real(0.5) * a * t0 * t0 + (dt - t0) * gsTarget;
	return std::pair<real, real>(ds, gsFinal);
}


/**
 * The time required to cover distance "dist" if initial speed is "gs" and acceleration is "a_gs"
 *
 * @param gs       initial ground speed
 * @param a_gs     signed ground speed acceleration
 * @param dist     non-negative distance
 * @return         time required to cover distance, -1 if speed will reach zero before achieving the distance
 *
 */
real Kinematics::timeToDistance(real gs, real a_gs, real dist) {
	real t1 = Util::rootNegC(real(0.5) * a_gs, gs, -dist);
	//f.pln(" $$ timeToDistance: RETURN t1 = "+t1+" gs = "+gs+", a_gs = "+a_gs+", dist = "+dist+", discr = "+discr);
	return t1;
}



std::pair<real,real> Kinematics::accelSpeedToRTA(real gsIn, real dist, real rta, real gsAccel) {

	real avgGs = dist/rta;
	int sign = 1;
	if (avgGs < gsIn) {
		sign = -1;
	}
	real a = gsAccel*sign;
	real A = real(0.5)*a;
	real B = -a*rta;
	real C = dist - gsIn*rta;
	real z = B*B-4*A*C;
	if (z < real(0.0)) {
		return std::pair<real,real>(-1,-1);
	}
	real ta = (-B+std::sqrt(z))/(2*A); // try first root
	real tb = (-B-std::sqrt(z))/(2*A);
	real t = -1;
	if (ta < rta && ta > 0) {
		t = ta;
	} else if (tb < rta && tb > 0) {
		t = tb;
	}
	return std::pair<real,real>(gsIn + a*t, t);
}




bool Kinematics::testLoSGs(const Vect3& so, const Velocity& vo, const Velocity& nvo, const Vect3& si, const Velocity& vi,
		real gsAccelOwn, real stopTime, real D, real H) {
	//		fpln(" $$$$ testLoSGs: vo = "+vo.toString()+" vi = "+vi.toString()+"  nvo = "+nvo.toString()+" stopTime="+Fm3(stopTime));
	real step = real(1.0);
	bool rtn = false;
	for (real t = 0; t < stopTime; t = t + step) {
		Vect3 soAtTm = gsAccelUntil(so, vo, t, nvo.gs(), gsAccelOwn).first;
		Vect3 siAtTm = si.linear(vi.vect3(),t);
		real distH = soAtTm.Sub(siAtTm).vect2().norm();
		real distV = std::abs(soAtTm.Sub(siAtTm).z());
		if (distH < D && distV < H) rtn =true;
	}
	return rtn;
}



// ****************************** Vertical Speed KINEMATIC CALCULATIONS *******************************

/**
 * Return the elevation angle (alternatively the negative glide-slope angle) for a climb (descent)
 * @return elevation angle [radians]
 */
real Kinematics::elevationAngle(const Velocity& v) {
	return Util::atan2_safe(v.vs(), v.gs());
}


Vect3 Kinematics::vsAccelPos(const Vect3& so3, const Velocity& vo3, real t, real a) {
	return Vect3(so3.x() + t*vo3.x(), so3.y() + t*vo3.y(), so3.z() + vo3.z()*t + real(0.5)*a*t*t);
}

std::pair<Vect3,Velocity> Kinematics::vsAccel(const Vect3& so3, const Velocity& vo3,  real t, real a) {
	real nvoVs = vo3.vs() + a*t;
	Velocity nvo = vo3.mkVs(nvoVs);
	return std::pair<Vect3,Velocity>(vsAccelPos(so3,vo3,t,a),nvo);
}

// std::pair<Vect3,Velocity> Kinematics::vsAccel(const std::pair<Vect3,Velocity>& sv0,  double t, double a) {
// 	  return vsAccel(sv0.first, sv0.second,t,a);
// }


real Kinematics::vsAccelTime(const Velocity& vo,real goalVs, real vsAccel) {
	return accelTime(vo.vs(),goalVs, vsAccel);;
}

// double Kinematics::vsAccelTime(double vs, double goalVs, double vsAccel) {
// 	double deltaVs = vs - goalVs;
// 	double rtn = std::abs(deltaVs/vsAccel);
// 	//f.pln("#### vsAccelTime: vs() = "+Units.str("fpm",vs)+" deltaVs = "+Units.str("fpm",deltaVs)+" rtn = "+rtn);
// 	return rtn;
// }


Triple<Vect3,Velocity,real> Kinematics::vsAccelGoal(const Vect3& so, const Velocity& vo, real goalVs, real vsAccel) {
	int sgn = 1;
	if (goalVs < vo.vs()) sgn = -1;
	real accelTime = vsAccelTime(vo, goalVs, vsAccel);
	Vect3 nso = vsAccelPos(so, vo, accelTime,  sgn*vsAccel);
	Velocity nvo = Velocity::mkVxyz(vo.x(),vo.y(),goalVs);
	return Triple<Vect3,Velocity,real>(nso,nvo,accelTime);
}


std::pair<Vect3,Velocity> Kinematics::vsAccelUntil(const Vect3& so, const Velocity& vo, real t, real goalVs, real vsAccel_d) {
	vsAccel_d = std::abs(vsAccel_d);
	// if (vsAccel_d < 0 ) {
	// 	std::cout << "Kinematics::vsAccelUntil: user supplied negative vsAccel!!" << std::endl;
	// 	vsAccel_d = -vsAccel_d;                              // make sure user supplies positive value
	// }
	real accelTime = vsAccelTime(vo,goalVs, vsAccel_d);
	int sgn = 1;
	if (goalVs < vo.vs()) sgn = -1;
	//Vect3 ns = Vect3::ZERO();
	if (t <= accelTime) {
		return vsAccel(so,vo,t,sgn*vsAccel_d);
	} else {
		Vect3 posEnd = vsAccelPos(so,vo,accelTime,sgn*vsAccel_d);
		Velocity nvo = Velocity::mkVxyz(vo.x(),vo.y(), goalVs);
		return linear(posEnd,nvo,t-accelTime);
	}
}


//  std::pair<Vect3,Velocity> Kinematics::vsAccelUntil(const std::pair<Vect3,Velocity>& sv0, double t, double goalVs, double vsAccel) {
//  	return vsAccelUntil(sv0.first, sv0.second, t, goalVs, vsAccel);
//  }



bool Kinematics::testLoSVs(const Vect3& so, const Velocity& vo, const Velocity& nvo, const Vect3& si, const Velocity& vi,
		real vsAccelOwn, real stopTime, real D, real H) {
	//fpln(" $$$$ testLoSTrk: vo = "+vo+" vi = "+vi+"  nvo = "+nvo+" nvi = "+nvi);
	real step = real(1.0);
	bool rtn = false;
	for (real t = 0; t < stopTime; t = t + step) {
		Vect3 soAtTm = vsAccelUntil(so, vo, t, nvo.vs(), vsAccelOwn).first;
		Vect3 siAtTm = si.linear(vi.vect3(),t);
		real distH = soAtTm.Sub(siAtTm).vect2().norm();
		real distV = std::abs(soAtTm.Sub(siAtTm).z());
		if (distH < D && distV < H) rtn =true;
	}
	return rtn;
}

real Kinematics::timeNeededForFLC(real deltaZ, real vsFLC, real vsAccel, bool kinematic) {
	real rtn;
	if (kinematic) {
		rtn = std::abs(deltaZ / vsFLC) + std::abs(vsFLC / vsAccel);
	} else {
		rtn = std::abs(deltaZ / vsFLC);
	}
	return rtn;
}

//static double antiDer1(double voz, double t, double a1) {   // alpha
//	//fpln(" $$$$ alpha: t ="+t+" a1 = "+a1+" voz = "+Units.str("fpm",voz)+"   return:"+(voz*t + 0.5*a1*t*t));
//	return voz*t + 0.5*a1*t*t;
//}
//
//static double antiDer2(double climbRate, double t) { // beta
//	return  climbRate*t;  // voz*t+ a1*T1*t;
//}
//
//static double antiDer3(double climbRate, double t, double a2) {  // gamma
//	return climbRate*t + 0.5*a2*t*t;
//}
//
////// Computes time for the vsLevelOut method,
//// Note: if T2 < T1, there is no constant vertical speed phase, If T1 < 0, target altitude is not achieveable
////
//// allowClimbRateChange   if true, the climbRate can be reduced, otherwise T1 might be set to -1, which
////                                 indicates failure
//// return <T1 = end of first accel ,T2 = start of constant vertical speed phase, T3 = start of deceleration, climbRate'>
////
//Quad<double,double,double,double> Kinematics::vsLevelOutTimesRWB(const std::pair<Vect3,Velocity>& sv0, double climbRate, double targetAlt,
//		double a1, double a2, bool allowClimbRateIncrease, bool allowClimbRateDecrease) {
//	Vect3    s0 = sv0.first;
//	Velocity v0 = sv0.second;
//	double soz = s0.z;
//	double voz = v0.z;
//	int altDir = -1;
//	if (targetAlt >= soz) altDir = 1;
//	climbRate = altDir*std::abs(climbRate);
//	int vsDir1 = -1;
//	if (climbRate >= voz) vsDir1 = 1;
//	int vsDir2 = -altDir;
//	a1 = vsDir1*std::abs(a1);
//	a2 = vsDir2*std::abs(a2);
//	//fpln(" vsLevelOutTimes: altDir = "+altDir+" vsDir1 = "+vsDir1);
//	//fpln(" #### vsLevelOutTimes: a = "+a+" climbRate = "+Units.str("fpm",climbRate)+" targetAlt = "+Units.str("ft",targetAlt));
//	double T1 = (climbRate - voz)/a1;
//	double S = (targetAlt-soz);
//	double S1 = antiDer1(voz,T1,a1);
//	double t3 = -climbRate/a2;
//	double S3 = antiDer3(climbRate,t3,a2);
//	double S2 = S - S1 - S3;
//	//fpln(" #### S = "+Units.str("ft",S)+" S1 = "+Units.str("ft",S1)+" S3 = "+Units.str("ft",S3)+" S2 = "+Units.str("ft",S2));
//	double T2 = S2/climbRate + T1;
//	double T3 = -climbRate/a2 + T2;
//	double cc = 0;
//	if (T2 < T1) {
//		//fpln("  vsLevelOutTimes: Case 2: no constant vertical speed phase! T1 = "+T1+" T2 ="+T2);
//		if (allowClimbRateDecrease) {
//			//double aa = a;
//			//double bb = 2*voz;
//			//cc = voz*voz/(2*a) - S;
//			double aa = 0.5*a1*(1 - a1/a2);
//			double bb = voz - (a1/a2)*voz;
//			cc = -voz*voz/(2*a2) - S;
//			double root1 = Util::root(aa,bb,cc,1);
//			double root2 = Util::root(aa,bb,cc,-1);
//			//fpln(" root1 = "+root1+" root2 = "+root2);
//			if (root1 >= 0) T1 = root1;
//			else if (root2 >= 0) T1 = root2;
//			else {
//				fpln(" vsLevelOut: Both roots are negative!  root1 = "+Fm2(root1)+" root2 = "+Fm2(root2));
//				T1 = -1;
//				return Quad<double,double,double,double>(-1.0,-1.0,-1.0,climbRate);
//			}
//			T2 = T1;
//			climbRate = voz + a1*T1;
//			T3 = -climbRate/a2 + T1;
//		} else {
//			T1 = -1;  // FAILURE
//		}
//	}
//	// Deal with special case where current vertical speed already exceeds climbRate (CHANGE climbRate)
//	if (vsDir1 != altDir && allowClimbRateIncrease) {
//		climbRate = voz;  // increase climb rate
//		//fpln(" vsLevelOutTimes: recompute climbRate = "+Units.str("fpm",climbRate) );
//		Quad<double,double,double,double> qTemp = vsLevelOutTimesRWB(sv0, climbRate, targetAlt, a1, a2, false, false);
//		T1 = qTemp.first;
//		T2 = qTemp.second;
//		T3 = qTemp.third;
//	}
//	//fpln(" T1 = "+T1+"  T2 = "+T2+" T3 = "+T3+" t3 = "+t3+" climbRate = "+Units.str("fpm",climbRate));
//	return Quad<double,double,double,double>(T1,T2,T3,climbRate);
//}
//
//Quad<double,double,double,double> Kinematics::vsLevelOutTimesRWB(const std::pair<Vect3,Velocity>& sv0, double climbRate, double targetAlt,
//		double a, bool allowClimbRateChange) {
//	return vsLevelOutTimesRWB(sv0, climbRate, targetAlt, a, -a, allowClimbRateChange, allowClimbRateChange);
//}
//
//Quad<double,double,double,double> Kinematics::vsLevelOutTimesRWB(const std::pair<Vect3,Velocity>& sv0, double climbRate, double targetAlt,
//		double a) {
//	return vsLevelOutTimesRWB(sv0, climbRate, targetAlt, a, -a, false, true);
//}
//
//
//double Kinematics::vsLevelOutTimeRWB(const std::pair<Vect3,Velocity>& sv0, double climbRate, double targetAlt, double a, bool allowClimbRateChange) {
//	Quad<double,double,double,double> qV = vsLevelOutTimesRWB(sv0,climbRate,targetAlt,a, allowClimbRateChange);
//	if (qV.first < 0) return -1;
//	else return qV.third;
//}
//
///** returns Pair that contains position and velocity at time t due to level out maneuver
// *
// * @param sv0        current position and velocity vectors
// * @param t          time point of interest
// * @param T1         end of first accel
// * @param T2         end of constant vertical speed phase
// * @param T3         end of deceleration (Level after this)
// * @param climbRate  climb rate
// * @param targetAlt  target altitude
// * @param a1         first acceleration
// * @param a2         second acceleration
// * @return
// */
//std::pair<Vect3, Velocity> Kinematics::vsLevelOutCalculationRWB(const std::pair<Vect3,Velocity>& sv0, double t,  double T1, double T2,
//		double T3, double climbRate, double targetAlt, double a1, double a2) {
//	Vect3    s0 = sv0.first;
//	Velocity v0 = sv0.second;
//	//fpln(f.Fm1(t)+"  #### vsLevelOutCalculation: s0 = "+f.sStr(s0)+" v0 = "+v0+" a = "+a+" climbRate = "+Units.str("fpm",climbRate)+" targetAlt = "+Units.str("ft",targetAlt));
//	//fpln(f.Fm1(t)+"  #### vsLevelOutCalculation: T1 ="+T1+" T2 = "+T2+" climbRate = "+Units.str("fpm",climbRate)+" targetAlt = "+Units.str("ft",targetAlt));
//	double soz = s0.z;
//	double voz = v0.z;
//	int altDir = -1;
//	if (targetAlt >= soz) altDir = 1;
//	climbRate = altDir*std::abs(climbRate);
//	int vsDir1 = -1;
//	if (climbRate >= voz) vsDir1 = 1;
//	int vsDir2 = -altDir;
//	a1 = vsDir1*std::abs(a1);
//	a2 = vsDir2*std::abs(a2);
//	Velocity nv = Velocity::ZEROV;
//	Vect3    ns = Vect3::ZERO();
//	if (t <= T1) {
//		nv = v0.mkVs(voz + a1*t);
//		ns = s0.linear(v0,t).mkZ(soz + antiDer1(voz,t,a1));
//		//fpln(t+" T<=T1: vsDir1 = "+vsDir1+" soz = "+Units.str("ft",soz)+" antiDer1(voz,t,a1) = "+Units.str("ft",antiDer1(voz,t,a1)));
//	} else if (t <= T2) {
//		nv = v0.mkVs(climbRate);
//		//fpln(t+" T<=T2: soz = "+Units.str("ft",soz)+" antiDer1(voz,T1,a1) = "+Units.str("ft",antiDer1(voz,T1,a1))+" antiDer2(climbRate,t-T1) = "+Units.str("ft",antiDer2(climbRate,t-T1)));
//		ns = s0.linear(v0,t).mkZ(soz + antiDer1(voz,T1,a1)+antiDer2(climbRate,t-T1));
//	} else if (t <= T3) {
//		nv = v0.mkVs(climbRate + a2*(t-T2));
//		//fpln("t<=T3: soz = "+Units.str("ft",soz)+" alpha(voz,T1,a) = "+Units.str("ft",antiDer1(voz,T1,a))
//		//		+" beta(climbRate,T2-T1,T1,a) = "+Units.str("ft",antiDer2(climbRate,t-T1,a))+" gamma(voz,t-T2,T1,T2,a) = "+Units.str("ft",antiDer3(voz,t-T2,a)));
//		ns = s0.linear(v0,t).mkZ(soz + antiDer1(voz,T1,a1)+ antiDer2(climbRate,T2-T1) +antiDer3(climbRate,t-T2,a2));
//	} else {
//		nv = v0.mkVs(0);
//		ns = s0.linear(v0,t).mkZ(soz + antiDer1(voz,T1,a1)+ antiDer2(climbRate,T2-T1) + antiDer3(climbRate,T3-T2,a2));
//	}
//	//fpln(f.Fm1(t)+"  #### vsLevelOutCalculation: vsDir = "+vsDir+" T2 = "+T2+"  ns = "+f.sStr(ns));
//	return std::pair<Vect3, Velocity>(ns,nv);
//}
//
//
//std::pair<Vect3,Velocity> Kinematics::vsLevelOutRWB(const std::pair<Vect3,Velocity>& sv0, double t, double climbRate, double targetAlt, double a, bool allowClimbRateChange) {
//	Quad<double,double,double,double> qV =  vsLevelOutTimesRWB(sv0,climbRate,targetAlt,a, allowClimbRateChange);
//	double T1 = qV.first;
//	double T2 = qV.second;
//	double T3 = qV.third;
//	if (T1 < 0) {
//		return std::pair<Vect3,Velocity>(Vect3::INVALID,Velocity::INVALIDV);
//	}
//	return vsLevelOutCalculationRWB(sv0, t, T1, T2, T3, qV.fourth, targetAlt, a,-a);
//}
//
//StateVector Kinematics::vsLevelOutFinalRWB(const std::pair<Vect3,Velocity>& sv0, double climbRate, double targetAlt, double a, bool allowClimbRateChange) {
//	Quad<double,double,double,double> qV =  vsLevelOutTimesRWB(sv0,climbRate,targetAlt,a, allowClimbRateChange);
//	double T1 = qV.first;
//	double T3 = qV.third;
//	if (T1 < 0) {         //  overshoot case
//		return StateVector(Vect3::INVALID,Velocity::INVALIDV,-1.0);
//	}
//	return StateVector(vsLevelOutCalculationRWB(sv0, qV.third, T1, qV.second, qV.third, qV.fourth,targetAlt,a,-a), T3);
//}
//
//
///** This version prioritizes being able to reach an altitude at all and then achieving the specified climb rate.
// *  The returned Statevector contains position,velocity, and time to reach target altitude.
// * @param sv0
// * @param climbRate
// * @param targetAlt
// * @param a
// * @return
// */
//StateVector Kinematics::vsLevelOutFinalRWB(const std::pair<Vect3,Velocity>& sv0, double climbRate, double targetAlt, double a) {
//	Quad<double,double,double,double> qV =  vsLevelOutTimesRWB(sv0,climbRate,targetAlt,a);
//	double T1 = qV.first;
//	double T3 = qV.third;
//	if (T1 < 0) {         //  overshoot case
//		return StateVector(Vect3::INVALID,Velocity::INVALIDV,-1.0);
//	}
//	return StateVector(vsLevelOutCalculationRWB(sv0, qV.third, T1, qV.second, qV.third, qV.fourth,targetAlt,a,-a), T3);
//}
//
//double Kinematics::vsLevelOutClimbRateRWB(const std::pair<Vect3, Velocity>& svo, double climbRate, double targetAlt,
//		double a1, double a2, bool allowClimbRateChange) {
//	Quad<double,double,double,double> ntp = vsLevelOutTimesRWB(svo, climbRate, targetAlt, a1, a2, false, false);
//	//fpln(" $$$ vsLevelOutTimes: "+ntp.first+" "+ ntp.second+" "+ ntp.third+" "+ntp.fourth+" "+ntp.fifth);
//	return  vsLevelOutCalculationRWB(svo, ntp.first, ntp.first, ntp.second, ntp.third, ntp.fourth, targetAlt, a1, a2).second.z;
//}


real V1(real voz, real a1, real t) {   // alpha
	//fpln(" $$$$ alpha: t ="+t+" a = "+a+" voz = "+Units::str("fpm",voz)+"   return:"+(voz*t + 0.5*a*t*t));
	return voz + a1*t;
}

real S1(real voz, real a1, real t) {   // alpha
	//fpln(" $$$$ alpha: t ="+t+" a = "+a+" voz = "+Units::str("fpm",voz)+"   return:"+(voz*t + 0.5*a*t*t));
	return voz*t + real(0.5)*a1*t*t;
}

real T3(real voz, real a1) {   // alpha
	//fpln(" $$$$ alpha: t ="+t+" a = "+a+" voz = "+Units::str("fpm",voz)+"   return:"+(voz*t + 0.5*a*t*t));
	return -voz/a1;
}

real S3(real voz, real a1) {   // alpha
	//fpln(" $$$$ alpha: t ="+t+" a = "+a+" voz = "+Units::str("fpm",voz)+"   return:"+(voz*t + 0.5*a*t*t));
	return S1(voz, a1, T3(voz, a1));
}


StateVector Kinematics::vsLevelOutFinal(const std::pair<Vect3,Velocity>& sv0, real climbRate, real targetAlt, real a, bool allowClimbRateChange) {
	Tuple5<real, real,real,real,real> qV =  vsLevelOutTimes(sv0,climbRate,targetAlt,a, allowClimbRateChange);
	real T1 = qV.first;
	real T3 = qV.third;
	if (T1 < 0) {         //  overshoot case
		//fpln(" $$$$$$ vsLevelOutFinal: T1 < 0,      targetAlt = "+Units::str("ft",targetAlt)+" currentAlt = "+Units::str("ft",sv0.first.z()));
		return StateVector(Vect3::INVALID(),Velocity::INVALID(),-real(1.0));
	}
	return StateVector(vsLevelOutCalculation(sv0, targetAlt, qV.fourth, qV.fifth, T1, qV.second, T3, T3),T3);
}

StateVector Kinematics::vsLevelOutFinal(const std::pair<Vect3,Velocity>& sv0, real climbRate, real targetAlt, real a) {
	Tuple5<real,real,real,real,real> qV = vsLevelOutTimesBase(sv0.first.z(),sv0.second.z(),climbRate,targetAlt,a,-a,true);
	real T1 = qV.first;
	real T3 = qV.third;
	if (T1 < 0) {         //  overshoot case
		//f.pln(" $$$$$$ vsLevelOutFinal: T1 < 0,      targetAlt = "+Units.str("ft",targetAlt)+" currentAlt = "+Units.str("ft",sv0.first.z()));
		return StateVector(Vect3::INVALID(),Velocity::INVALID(),-real(1.0));
	}
	return StateVector(vsLevelOutCalculation(sv0, targetAlt, qV.fourth, qV.fifth, T1, qV.second, T3, T3),T3);
}


// bool Kinematics::overShoot(const std::pair<Vect3, Velocity>& svo, double climbRate, double targetAlt, double accelup,
// 		double acceldown, bool allowClimbRateChange){
// 	double a2 = vsLevelOutTimes(svo, climbRate, targetAlt, accelup, acceldown, allowClimbRateChange).fifth;
// 	int sgnv = -1;
// 	if (svo.second.z>=0) sgnv =1;
// 	int altDir = -1;
// 	if (targetAlt-svo.first.z>=0) altDir = 1;
// 	if (sgnv==altDir && std::abs(targetAlt-svo.first.z)< std::abs(S3(svo.second.z, a2))) return true;
// 	else return false;
// }


///**
// * Returns true if time t is within the constant velocity segment of the climb
// * All values are in internal units
// */
//bool Kinematics::inConstantClimb(double sz, double vz, double t, double climbRate, double targetAlt, double a) {
//	std::pair<Vect3,Velocity> sv0 = std::pair<Vect3,Velocity>(Vect3::mkXYZ(0,0,sz),Velocity::mkVxyz(0, 0, vz));
//	Quad<double,double,double,double> qV =  vsLevelOutTimesRWB(sv0,climbRate,targetAlt,a,true);
//	return t > qV.first && t < qV.second;
//}
//


// /** Helper function for vsLevelOutTimesAD.
//  *  Note: This could be integrated into the function vsLevelOutTimesAD as a recursive call if desired.
//  *
//  * @param s0z          initial vertical position
//  * @param v0z		   initial vertical velocity
//  * @param climbRate    desired vertical speed for the climb/descent (positive), sign calculated in code
//  * @param targetAlt    target altitude
//  * @param accelup      maximum positive acceleration
//  * @param acceldown    maximum negative acceleration
//  * @param allowClimbRateChange	if true, if the current velocity is of greater magnitude than the specified climb rate,
//  * 										then continue at the current velocity (prioritize achieving the desired altitude).
//  * 										If false, first achieve the goal climb rate (prioritize achieving the indicated vs)
//  *
//  *
//  * @return <T1 = end of first accel ,T2 = end of constant vertical speed phase, T3 = end of deceleration, a1 = acceleration for phase 1, a2 =acceleration for phase 2>
//  */
Tuple5<real,real,real,real,real> Kinematics::vsLevelOutTimesBase(real s0z, real v0z, real climbRate, real targetAlt,
		real accelup, real acceldown, bool allowClimbRateChange) {

	int altDir = -1;
	if (targetAlt >= s0z) altDir = 1;
	climbRate = altDir*std::abs(climbRate);
	if (allowClimbRateChange) climbRate = altDir*(Util::max(std::abs(climbRate), std::abs(v0z)));
	real S = targetAlt-s0z;
	real a1 = acceldown;
	if (climbRate>=v0z) a1 = accelup;
	real a2 = accelup;
	if (targetAlt>=s0z) a2 = acceldown;
	real T1 = (climbRate - v0z)/a1;

	if (std::abs(S)>= std::abs(S1(v0z, a1, T1)+S3(V1(v0z, a1, T1), a2))) {
		real T2 = (S - S1(v0z, a1, T1)-S3(V1(v0z, a1, T1), a2))/climbRate;
		//fpln("times1 case1");
		return Tuple5<real,real,real,real,real>(T1, T1+T2, T1+T2+T3(climbRate, a2), a1, a2);
	} else {
		real aa = real(0.5)*a1*(1 - a1/a2);
		real bb = v0z*(1- (a1/a2));
		real cc = -v0z*v0z/(2*a2) - S;
		real root1 = Util::root(aa,bb,cc,1);
		real root2 = Util::root(aa,bb,cc,-1);
		if (root1<0)  T1 = root2;
		else if (root2<0) T1 = root1;
		else
			T1= Util::min(root1, root2);
		//fpln("times1 case2");
		return Tuple5<real, real,real,real,real>(T1, T1, T1+T3(V1(v0z, a1, T1), a2),a1,a2);
	}
}

Tuple5<real,real,real,real,real> Kinematics::vsLevelOutTimes(real s0z, real v0z, real climbRate, real targetAlt,
		real accelup, real acceldown, bool allowClimbRateChange) {

	int sgnv = -1;
	if (v0z >= 0) sgnv = 1;
	int altDir = -1;
	if (targetAlt >= s0z) altDir = 1;
	real S = targetAlt-s0z;
	real a1 = acceldown;
	if (targetAlt>=s0z) a1 = accelup;
	real a2 = accelup;
	if (targetAlt>=s0z) a2 = acceldown;


	if (sgnv==altDir || Util::almost_equals(v0z, real(0.0))) {
		if (std::abs(S)>=std::abs(S3(v0z, a2))) {
			//fpln(" ##times Case1.1");
			return vsLevelOutTimesBase(s0z, v0z, climbRate, targetAlt, accelup, acceldown, allowClimbRateChange);
		} else {
			Tuple5<real,real,real, real, real> ot = vsLevelOutTimesBase(s0z+S3(v0z, a2), real(0.0), climbRate, targetAlt, accelup, acceldown, allowClimbRateChange);
			//fpln("times Case1.2");
			return Tuple5<real, real,real,real,real>(-v0z/a2+ot.first, -v0z/a2+ot.second, -v0z/a2+ot.third , ot.fourth, ot.fifth);
		}
	} else {
		Tuple5<real,real,real, real, real> ot = vsLevelOutTimesBase(s0z+ S3(v0z, a1), real(0.0), climbRate, targetAlt, accelup, acceldown, allowClimbRateChange);
		//fpln("times Case 2");
		return Tuple5<real,real,real,real,real>(-v0z/a1+ot.first, -v0z/a1+ot.second, -v0z/a1+ot.third , ot.fourth, ot.fifth);
	}
}


Tuple5<real,real,real,real,real> Kinematics::vsLevelOutTimes(const std::pair<Vect3, Velocity>& svo, real climbRate, real targetAlt,
		real accelup, real acceldown, bool allowClimbRateChange) {
	real s0z = svo.first.z();
	real v0z = svo.second.z();
	return vsLevelOutTimes(s0z, v0z, climbRate, targetAlt, accelup, acceldown, allowClimbRateChange);
}

Tuple5<real,real,real,real,real> Kinematics::vsLevelOutTimes(const std::pair<Vect3, Velocity>& svo, real climbRate, real targetAlt,
		real a, bool allowClimbRateChange) {
	real s0z = svo.first.z();
	real v0z = svo.second.z();
	return vsLevelOutTimes(s0z, v0z, climbRate, targetAlt, a, -a, allowClimbRateChange);
}

// Tuple5<double,double,double,double,double> Kinematics::vsLevelOutTimes(const std::pair<Vect3, Velocity>& svo, double climbRate, double targetAlt, double a) {
// 	double s0z = svo.first.z;
// 	double v0z = svo.second.z;
// 	return vsLevelOutTimes(s0z, v0z, climbRate, targetAlt, a, -a, true);
// }


real Kinematics::vsLevelOutClimbRate(const std::pair<Vect3, Velocity>& svo, real climbRate, real targetAlt,
		real accelup, real acceldown, bool allowClimbRateChange) {
	Tuple5<real,real,real,real,real> ntp = vsLevelOutTimes(svo, climbRate, targetAlt, accelup, acceldown, allowClimbRateChange);
	//fpln(" $$$ vsLevelOutTimes: "+ntp.first+" "+ ntp.second+" "+ ntp.third+" "+ntp.fourth+" "+ntp.fifth);
	return vsLevelOutCalculation(svo, targetAlt, ntp.fourth, ntp.fifth, ntp.first, ntp.second, ntp.third, ntp.first).second.z();

}


real Kinematics::vsLevelOutTime(const std::pair<Vect3,Velocity>& sv0, real climbRate, real targetAlt, real a, bool allowClimbRateChange) {
	Tuple5<real,real,real,real,real> qV = vsLevelOutTimes(sv0,climbRate,targetAlt,a, -a, allowClimbRateChange);
	if (qV.first < 0) return -1;
	else return qV.third;
}


real Kinematics::vsLevelOutTime(const std::pair<Vect3,Velocity>& sv0, real climbRate, real targetAlt, real a) {
	return vsLevelOutTime(sv0, climbRate, targetAlt, a, true);
}

std::pair<real, real> Kinematics::vsLevelOutCalc(real soz, real voz, real targetAlt, real a1, real a2, real t1, real t2, real t3,  real t) {
	real nz = 0;
	real nvs = 0;
	if (t <= t1) {
		nvs = (voz + a1*t);
		nz = (soz + S1(voz,a1, t));
		//	fpln("#Phase1, time "+t+" Vel "+nv.z+" Pos "+ns.z);
	} else if (t <= t2) {
		nvs = (voz+a1*t1);
		nz = (soz + S1(voz,a1, t1)+ V1(voz, a1, t1)*(t-t1));
		//	fpln("#Phase2, time "+t+" Vel "+nv.z+" Pos "+ns.z);
	} else if (t <= t3) {
		nvs = (voz+a1*t1+a2*(t-t2));
		nz = (soz + S1(voz,a1, t1)+ V1(voz, a1, t1)*(t2-t1) + S1(V1(voz, a1, t1),a2, t-t2));
		//    fpln("#Phase3, time "+t+" Vel "+nv.z+" Pos "+ns.z);
	} else {
		nvs = 0;
		nz = targetAlt;
		//	fpln("#Phase4, time "+t+" Vel "+nv.z+" Pos "+ns.z);
	}
	return std::pair<real, real>(nz,nvs);
}



/** returns Pair that contains position and velocity at time t due to level out maneuver based on vsLevelOutTimesAD
 *
 * @param sv0        			current position and velocity vectors
 * @param t          			time point of interest
 * @param climbRate  			climb rate
 * @param targetAlt  			target altitude
 * @param a1         			first acceleration
 * @param a2         			second acceleration
 * @param allowClimbRateChange allows climbRate to change to initial velocity if it can help.
 * @return
 */
std::pair<Vect3, Velocity> Kinematics::vsLevelOutCalculation(const std::pair<Vect3,Velocity>& sv0,
		real targetAlt, real a1, real a2, real t1, real t2, real t3,  real t) {
	Vect3 s0 = sv0.first;
	Velocity v0 = sv0.second;
	real soz = s0.z();
	real voz = v0.z();
	std::pair<real, real> vsL = vsLevelOutCalc(soz,voz, targetAlt, a1, a2, t1, t2, t3, t);
	real nz = vsL.first;
	real nvs = vsL.second;
	Velocity nv = v0.mkVs(nvs);
	Vect3 ns = s0.linear(v0.vect3(),t).mkZ(nz);
	return std::pair<Vect3, Velocity>(ns,nv);
}

std::pair<Vect3, Velocity> Kinematics::vsLevelOut(const std::pair<Vect3, Velocity>& sv0, real t, real climbRate,
		real targetAlt, real accelUp, real accelDown, bool allowClimbRateChange) {
	Tuple5<real,real,real,real,real> LevelParams = vsLevelOutTimes(sv0, climbRate, targetAlt, accelUp, accelDown, allowClimbRateChange);
	return vsLevelOutCalculation(sv0, targetAlt, LevelParams.fourth, LevelParams.fifth, LevelParams.first, LevelParams.second, LevelParams.third, t);
}

std::pair<Vect3, Velocity> Kinematics::vsLevelOut(const std::pair<Vect3, Velocity>& sv0, real t, real climbRate,
		real targetAlt, real a, bool allowClimbRateChange) {
	return vsLevelOut(sv0, t, climbRate, targetAlt, a, -a, allowClimbRateChange);
}

//std::pair<Vect3, Velocity> Kinematics::vsLevelOut(const std::pair<Vect3, Velocity>& sv0, double t, double climbRate,
//		double targetAlt, double a) {
//	return vsLevelOut(sv0, t, climbRate, targetAlt, a, -a, true);
//}

real Kinematics::trackFrom(const Vect3& p1, const Vect3& p2) {
	return p2.Sub(p1).vect2().trk();
}



Vect4 Kinematics::minDistBetweenTrk(const Vect3& so, const Velocity& vo, const Velocity& nvo, const Vect3& si, const Velocity& vi, const Velocity& nvi,
		real bankAngOwn, real stopTime) {
	real minDist =  MAXDOUBLE;
	real minDistH =  MAXDOUBLE;
	real minDistV =  MAXDOUBLE;
	real minT = -1;
	real step = real(1.0);
	for (real t = 0; t < stopTime; t = t + step) {
		//Vect3 soAtTm = turnUntilPosition(so, vo, nvo.track(), bankAngOwn, t, turnRightOwn);
		//Vect3 siAtTm = turnUntilPosition(si, vi, nvi.track(), bankAngTraf, t, turnRightTraf);
		std::pair<Vect3,Velocity> psv = Kinematics::turnUntilTrack(so, vo, t, nvo.trk(), bankAngOwn);
		Vect3 soAtTm = psv.first;
		Vect3 vown = psv.second.vect3();
		std::pair<Vect3,Velocity> psvi = Kinematics::turnUntilTrack(si, vi, t, nvi.trk(), bankAngOwn);
		Vect3 siAtTm = psvi.first;
		Vect3 vtraf = psvi.second.vect3();
		//fpln(" $$$$ minDistBetweenTrk: soAtTm = "+f.sStr(soAtTm)+" siAtTm = "+f.sStr(siAtTm));
		Vect3 s = soAtTm.Sub(siAtTm);
		real dist = s.norm();
		//fpln(" $$$$ minDistBetweenTrk: t = "+t+"  dist = "+Units.str("nm",dist));
		if (dist < minDist) {               // compute distances at TCA in 3D
			real distH = s.vect2().norm();
			real distV = std::abs(s.z());
			minDist = dist;
			minDistH = distH;
			minDistV = distV;
			minT = t;
		}
        bool divg = s.dot(vown.Sub(vtraf)) > 0;
      	if (divg) break;

	}
	return Vect4(minDistH,minDist,minDistV,minT);
}


/** Minimum distance between two aircraft when BOTH aircraft gs accelerate, compute trajectories up to time stopTime
 *
 * @param so    initial position of ownship
 * @param vo    initial velocity of ownship
 * @param nvo   the target velocity of ownship (i.e. after turn maneuver complete)
 * @param si    initial position of traffic
 * @param vi    initial velocity of traffic
 * @param nvi           target velocity of traffic (i.e. after acceleration maneuver complete)
 * @param gsAccelOwn    ground speed acceleration of the ownship
 * @param gsAccelTraf   ground speed acceleration of the traffic
 * @param stopTime         the duration of the turns
 * @return                 minimum distance data packed in a Vect4
 */
Vect4 Kinematics::minDistBetweenGs(const Vect3& so, const Velocity& vo, const Velocity& nvo, const Vect3& si, const Velocity& vi,  const Velocity& nvi,
		real gsAccelOwn, real gsAccelTraf, real stopTime) {
	real minDist = MAXDOUBLE;
	real minDistH = MAXDOUBLE;
	real minDistV = MAXDOUBLE;
	//fpln(" $$$$ minDistBetween: vo = "+vo+" vi = "+vi+"  nvo = "+nvo+" nvi = "+nvi);
	real step = real(1.0);
	real minT = -1;
	//fpln(" $$$$$$$$$$$$$$$$$$$$ minDistBetweenTrk: step = "+step);
	for (real t = 0; t < stopTime; t = t + step) {
		Vect3 soAtTm = Kinematics::gsAccelUntil(so, vo, t, nvo.gs(), gsAccelOwn).first;
		Vect3 siAtTm = Kinematics::gsAccelUntil(si, vi, t, nvi.gs(), gsAccelTraf).first;
		Vect3 s = soAtTm.Sub(siAtTm);
		//fpln(" $$$$ minDistBetweenTrk: soAtTm = "+f.sStr(soAtTm)+" siAtTm = "+f.sStr(siAtTm));
		real dist = s.norm();
		//fpln(" $$$$ minDistBetweenTrk: dist = "+Units.str("nm",dist));
		if (dist < minDist) {               // compute distances at TCA in 3D
			real distH = s.vect2().norm();
			real distV = std::abs(s.z());
			minDist = dist;
			minDistH = distH;
			minDistV = distV;
			minT = t;
		}
	}
	return Vect4(minDistH,minDist,minDistV,minT);
}

Vect4 Kinematics::minDistBetweenGs(const Vect3& so, const Velocity& vo, const Velocity& nvo, const Vect3& si, const Velocity& vi,
		real gsAccelOwn, real stopTime) {
	real minDist = MAXDOUBLE;
	real minDistH = MAXDOUBLE;
	real minDistV = MAXDOUBLE;
	//fpln(" $$$$ minDistBetween: vo = "+vo+" vi = "+vi+"  nvo = "+nvo+" nvi = "+nvi);
	real step = real(1.0);
	real minT = -1;
	//fpln(" $$$$$$$$$$$$$$$$$$$$ minDistBetweenTrk: step = "+step);
	for (real t = 0; t < stopTime; t = t + step) {
		Vect3 soAtTm = Kinematics::gsAccelUntil(so, vo, t, nvo.gs(), gsAccelOwn).first;
		Vect3 siAtTm = si.linear(vi.vect3(),t);
		Vect3 s = soAtTm.Sub(siAtTm);
		//fpln(" $$$$ minDistBetweenTrk: soAtTm = "+f.sStr(soAtTm)+" siAtTm = "+f.sStr(siAtTm));
		real dist = s.norm();
		//fpln(" $$$$ minDistBetweenTrk: dist = "+Units.str("nm",dist));
		if (dist < minDist) {               // compute distances at TCA in 3D
			real distH = s.vect2().norm();
			real distV = std::abs(s.z());
			minDist = dist;
			minDistH = distH;
			minDistV = distV;
			minT = t;
		}
	}
	return Vect4(minDistH,minDist,minDistV,minT);
}


Vect4 Kinematics::minDistBetweenVs(const Vect3& so, const Velocity& vo, const Velocity& nvo, const Vect3& si, const Velocity& vi,  const Velocity& nvi,
		real vsAccelOwn, real vsAccelTraf, real stopTime) {
	real minDist = MAXDOUBLE;
	real minDistH = MAXDOUBLE;
	real minDistV = MAXDOUBLE;
	real minT = -1;
	//fpln(" $$$$ minDistBetweenVs: vo = "+vo+" vi = "+vi+"  nvo = "+nvo+" nvi = "+nvi);
	real step = real(1.0);
	for (real t = 0; t < stopTime; t = t + step) {
		Vect3 soAtTm = Kinematics::vsAccelUntil(so, vo, t, nvo.vs(), vsAccelOwn).first;
		Vect3 siAtTm = Kinematics::vsAccelUntil(si, vi, t, nvi.vs(), vsAccelTraf).first;
		//fpln(" $$$$ minDistBetweenVs: soAtTm = "+f.sStr(soAtTm)+" siAtTm = "+f.sStr(siAtTm));
		Vect3 s = soAtTm.Sub(siAtTm);
		real dist = s.norm();
		//fpln(" $$$$ minDistBetweenVs: dist = "+Units.str("nm",dist));
		if (dist < minDist) {               // compute distances at TCA in 3D
			real distH = s.vect2().norm();
			real distV = std::abs(s.z());
			minDist = dist;
			minDistH = distH;
			minDistV = distV;
			minT = t;
		}
	}
	return Vect4(minDistH,minDist,minDistV,minT);
}
Vect4 Kinematics::minDistBetweenVs(const Vect3& so, const Velocity& vo, const Velocity& nvo, const Vect3& si, const Velocity& vi,
		real vsAccelOwn, real stopTime) {
	real minDist = MAXDOUBLE;
	real minDistH = MAXDOUBLE;
	real minDistV = MAXDOUBLE;
	//fpln(" $$$$ minDistBetween: vo = "+vo+" vi = "+vi+"  nvo = "+nvo+" nvi = "+nvi);
	real step = real(1.0);
	real minT = -1;
	for (real t = 0; t < stopTime; t = t + step) {
		Vect3 soAtTm = Kinematics::vsAccelUntil(so, vo, t, nvo.vs(), vsAccelOwn).first;
		Vect3 siAtTm = si.linear(vi.vect3(),t);
		//fpln(" $$$$ minDistBetweenVs: soAtTm = "+f.sStr(soAtTm)+" siAtTm = "+f.sStr(siAtTm));
		Vect3 s = soAtTm.Sub(siAtTm);
		real dist = s.norm();
		//fpln(" $$$$ minDistBetweenVs: distV = "+Units.str("ft",dist));
		if (dist < minDist) {               // compute distances at TCA in 3D
			real distH = s.vect2().norm();
			real distV = std::abs(s.z());
			minDist = dist;
			minDistH = distH;
			minDistV = distV;
			//fpln(" $$$$ minDistBetweenVs: minDistV = "+Units.str("ft",minDistV));
			minT = t;
		}
	}
	return Vect4(minDistH,minDist,minDistV,minT);
}

real Kinematics::gsTimeConstantAccelFromDist(real gs1, real a, real dist) {
	real t1 = Util::root(real(0.5)*a,  gs1,  -dist, 1);
	real t2 = Util::root(real(0.5)*a,  gs1,  -dist, -1);
	return t1 < 0 ? t2 : t1;
}

}
