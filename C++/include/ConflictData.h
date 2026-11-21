/*
 * Copyright (c) 2015-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */
/*
 * LossData.h
 *
 *  Created on: Nov 21, 2014
 *      Author: rbutler
 */

#ifndef CONFLICTDATA_H_
#define CONFLICTDATA_H_
#include "Real.h"

#include "LossData.h"
#include "Vect3.h"
#include "Horizontal.h"
#include "Vertical.h"
#include "Interval.h"

#include <string>

namespace larcfm {

class ConflictData : public LossData {

public:

  real  time_crit;   // relative time to critical point
  real  dist_crit;   // distance or severity at critical point (0 is most critical, +inf is least severe)

private:
  Vect3 s_; // Relative position
  Vect3 v_; // Relative velocity

  static real squircle(real x, real y);

public:

	enum SLoWCClassification {
    SLoWC_Catastrophic, SLoWC_Hazardous, SLoWC_Major, SLoWC_Minor, SLoWC_Minimal, SLoWC_None
	};

  std::string toString() const;

  ConflictData(real t_in, real t_out, real t_crit, real d_crit, const Vect3& s, const Vect3& v);
  ConflictData(const LossData& ld, real t_crit, real d_crit, const Vect3& s, const Vect3& v);
  ConflictData();

  static const ConflictData& EMPTY();

  /**
   * Returns internal vector representation of relative aircraft position.
   */
  Vect3 const & get_s() const;

  /**
   * Returns internal vector representation of relative aircraft velocity.
   */
  Vect3 const & get_v() const;

  /**
   * Returns HMD, in internal units, within lookahead time t, in seconds, assuming straight line trajectory.
   */
  real HMD(real T) const;

  /**
   * Returns HMD, in specified units, within lookahead time t, in seconds, assuming straight line trajectory.
   */
  real HMD(const std::string& u, real T) const;

  /**
   * Returns VMD, in internal units, within lookahead time t, in seconds, assuming straight line trajectory.
   */
  real VMD(real T) const;

  /**
   * Returns VMD, in specified units, within lookahead time t, in seconds, assuming straight line trajectory.
   */
  real VMD(const std::string& u, real T) const;

  /**
   * Horizontal separation
   * @return Horizontal separation in internal units at current time
   */
  real horizontalSeparation() const;

  /**
   * Horizontal separation
   * @param u units
   * @return Horizontal separation in specified units u at current time
   */
  real horizontalSeparation(const std::string& u) const;

  /**
   * Horizontal separation
   * @param time in seconds
   * @return Horizontal separation in internal units at given time
   */
  real horizontalSeparationAtTime(real time) const;

  /**
   * Horizontal separation at given time
   * @param time time in seconds
   * @param u units
   * @return Horizontal separation in specified units at given time
   */
  real horizontalSeparationAtTime(const std::string& u, real time) const;

  /**
   * Vertical separation
   * @return Vertical separation in internal units at current time
   */
  real verticalSeparation() const;

  /**
   * Vertical separation
   * @param u units
   * @return Vertical separation in specified units at current time
   */
  real verticalSeparation(const std::string& u) const;

  /**
   * Vertical separation at given time
   * @param time time in seconds
   * @return Vertical separation in internal units at given time
   */
  real verticalSeparationAtTime(real time) const;

  /**
   * Vertical separation at given time
   * @param time time in seconds
   * @param u units
   * @return Vertical separation in specified units at given time
   */
  real verticalSeparationAtTime(const std::string& u, real time) const;

  /**
   * Time to horizontal closest point of approach in seconds.
   * When aircraft are diverging, tcpa is defined as 0.
   */
  real tcpa2D() const;

  /**
   * Time to 3D closest point of approach in seconds.
   * When aircraft are diverging, tcpa is defined as 0
   */
  real tcpa3D() const;

  /**
   * Time to co-altitude.
   * @return time to co-altitude in seconds. Returns NaN is v_.z is zero.
   */
  real tcoa() const;

  /**
   * Horizontal relative speed (ownship wrt intruder)
   * @return Horizontal relative speed in internal units at current time
   */
  real horizontalRelativeSpeed() const;

  /**
   * Horizontal relative speed (ownship wrt intruder)
   * @param u units
   * @return Horizontal relative speed in specified units u at current time
   */
  real horizontalRelativeSpeed(const std::string& u) const;

  /**
   * Horizontal closure rate
   * @return Horizontal closure rate in internal units at current time
   * If positive, aircraft are horizontally converging.
   */
  real horizontalClosureRate() const;

  /**
   * Horizontal closure rate
   * @param u units
   * @return Horizontal closure rate in specified units u at current time
   * If positive, aircraft are horizontally converging.
   */
  real horizontalClosureRate(const std::string& u) const;

  /**
   * Vertical relative speed (ownship wrt intruder)
   * @return Vertical relative speed in internal units at current time
   */
  real verticalRelativeSpeed() const;

  /**
   * Vertical relative speed (ownship wrt intruder)
   * @param u units
   * @return Vertical relative speed in specified units at current time
   */
  real verticalRelativeSpeed(const std::string& u) const;

  /**
   * Vertical closure rate
   * @return Vertical closure rate in internal units at current time
   * If positive, aircraft are vertically converging.
   */
  real verticalClosureRate() const;

  /**
   * Vertical closure rate
   * @param u units
   * @return Vertical closure rate in specified units at current time
   * If positive, aircraft are vertically converging.
   */
  real verticalClosureRate(const std::string& u) const;

  /**
   * return SLoWC classification for SC-228 SLoWC value
   **/
  static SLoWCClassification classifySLoWC(real slowc);

  static std::string SLoWCClassificationToString(SLoWCClassification slowc_class);

  /**
	 * Severity of Loss of Well Clear calculation for one time instant, 0.0 (minimal) - 100.0 (catastrophic).
	 * For an encounter, it should take max of every instance (SC-228)
	 **/
	real SLoWC(real DMOD, real TAUMOD, real ZTHR) const;

  /**
   * @return A time in seconds that can be used to compare severity of conflicts for arbitrary
   * well-clear volumes. This time is not necessarily TCPA. ** Don't use it as TCPA. **
   */
  real getCriticalTimeOfConflict() const;

  /**
   * @return A non-negative scalar that can be used to compare severity of conflicts for arbitrary
   * well-clear volumes. This scalar is a distance in the mathematical way. It is 0 when aircraft are
   * at the same poistion, but it isn't a distance in the physical way. In particular, this distance
   * is unitless. ** Don't use as CPA **
   */
  real getDistanceAtCriticalTime() const;

};

} /* namespace larcfm */

#endif /* SRC_LOSSDATA_H_ */
