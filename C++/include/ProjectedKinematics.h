/*
 * Kinematics.h
 * 
 * Copyright (c) 2011-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */

#ifndef PROJECTEDKINEMATICS_H_
#define PROJECTEDKINEMATICS_H_
#include "Real.h"

#include "Position.h"
#include "Velocity.h"
#include "Kinematics.h"
#include "EuclideanProjection.h"
#include "Triple.h"
#include "Quad.h"


namespace larcfm {

/**
 * This class contains versions of the Kinematics functions that have been lifted to deal with Position objects instead of Vect3 objects.
 *
 */
class ProjectedKinematics {
public:
	/**
	 * Calculate the angle of a constant-radius turn from two points and the radius
	 * 
	 * @param s1
	 * @param s2
	 * @param R
	 * @return the turn angle
	 */
	static real turnAngle(const Position& s1, const Position& s2, real R);

  /**
   * Horizontal distance covered in a turn
   * 
   * @param s1
   * @param s2
   * @param R
   * @return the turn distance
   */
	static real turnDistance(const Position& s1, const Position& s2, real R);

	/**
	 * Given two points on a turn and the velocity (direction) at the first point, determine the direction for the shortest turn going through the second point,
	 * returning true if that relative direction is to the right
	 *
	 * @param s1
	 * @param v1
	 * @param s2
	 * @return true if clockwise turn
	 */
	static bool clockwise(const Position& s1, const Velocity& v1, const Position& s2);

	//static double closestTimeOnTurn(const Position& turnstart, const Velocity& v1, double omega, const Position& x, double endTime);

	/**
	 *  Return distance along turn where x is closest to.
	 *  This returns -1 if x is at the turn's center
	 *
	 * @param turnstart starting position of turn
	 * @param v1        velocity
	 * @param R         turn radius
	 * @param dir       +1 = turn right, -1  turn left
	 * @param x         a position
	 * @param endDist   end distance of turn
	 * @return          distance along turn
	 */
	static real closestDistOnTurn(const Position& turnstart, const Velocity& v1, real R, int dir, const Position& x, real endDist);

	static std::pair<Position,Velocity> turn(const Position& so, const Velocity& vo, real t, real R, bool turnRight);
	static std::pair<Position,Velocity> gsAccel(const Position& so, const Velocity& vo, real t, real accel);

	static std::pair<Position,Velocity> gsAccelUntil(const Position& so, const Velocity& vo, real t, real goalGs, real a);
	static std::pair<Position,Velocity> vsAccel(const Position& so, const Velocity& vo, real t, real accel);

	static std::pair<Position,Velocity> vsAccelUntil(const Position& so, const Velocity& vo, real t, real goalVs,  real accel);

	/**
	 *  Position and velocity after t time units turning in direction "dir" with radius R.  This is a wrapper around turnPosition and
	 *  turnVelocity for Position objects,and uses the projection defined in the static Projection class.
	 * @param so  starting position
	 * @param vo  initial velocity
	 * @param R   turn radius
	 * @param t   time of turn [secs]
	 * @param turnRight true iff only turn direction is to the right
	 * @return Position and Velocity after t time
	 */
	static std::pair<Position,Velocity> turnOmega(const Position& so, const Velocity& vo, real t, real omega);

	static std::pair<Position,Velocity> turnOmega(const Position& so, const Velocity& vo, real t, real omega, const EuclideanProjection& proj);


	/**
	 *  Position and velocity after t time units turning in direction "dir" with radius R.  This is a wrapper around turnUntil
	 *  and uses the projection defined in the static Projection class.
	 * @param so  starting position
	 * @param vo  initial velocity
	 * @param goalTrack the target track angle
	 * @param bankAngle the aircraft's bank angle
	 * @param t   time of turn [secs]
	 * @param turnRight true iff only turn direction is to the right
	 * @return Position and Velocity after t time
	 */
	static std::pair<Position,Velocity> turnUntil(const Position& so, const Velocity& vo, real t, real goalTrack, real bankAngle);

	static std::pair<Position,Velocity> turnUntil(const std::pair<Position,Velocity>& sv, real t, real goalTrack, real bankAngle);

	static std::pair<Position,Velocity> linear(const std::pair<Position,Velocity>& p, real t);

	static std::pair<Position,Velocity> linear(const Position& so, const Velocity& p, real t);

	static  std::pair<Position,real> intersection(const Position& so, const Velocity& vo, const Position& si, const Velocity& vi);


	static real timeOfintersection(const Position& so, const Velocity& vo, const Position& si, const Velocity& vi);


	static Quad<Position,Velocity,real,int> directToPoint(const Position& so, const Velocity& vo, const Position& wp, real R);


	static Triple<Position,real,real> genDirectToVertex(const Position& sop, const Velocity& vo, const Position& wpp, real bankAngle, real timeBeforeTurn);

	static Triple<Position,Velocity,real> vsLevelOutFinal(const Position& so, const Velocity& vo, real climbRate, real targetAlt, real a);

	static real vsLevelOutTime(const Position& so, const Velocity& vo, real climbRate, real targetAlt, real a);

	static std::pair<Position,Velocity> vsLevelOut(const Position& so, const Velocity& vo, real t, real climbRate, real targetAlt, real a);


};

}

#endif //PROJECTEDKINEMATICS_H_
