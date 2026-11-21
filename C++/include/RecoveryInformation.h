/*
 * Copyright (c) 2015-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */
#ifndef RECOVERYINFORMATION_H_
#define RECOVERYINFORMATION_H_
#include "Real.h"

#include <string>

namespace larcfm {

class RecoveryInformation {

  private:
  real time_;  // Time to recovery in seconds. NaN: recovery no computed, Negative Infinity: recovery not available
  int nfactor_;  // Number of times the recovery volume was reduced
   real horizontal_distance_; // Minimum horizontal distance in internal units, i.e., [m]
   real vertical_distance_; // Minimum vertical distance in internal units, i.e., [m]

  public:
   RecoveryInformation(real t, int nfactor, real hs, real vs);

  /**
   * @return Number of time the recovery volume was reduced
   */
  int nFactor() const;


  /**
   * @return Time to recovery in seconds
   */
  real timeToRecovery() const;

  /**
   * @return Time to recovery in given units
   */
  real timeToRecovery(const std::string& u) const;

  /**
   * @return Recovery horizontal distance in internal units, i.e., [m]
   */
  real recoveryHorizontalDistance() const;

  /**
   * @return Recovery horizontal distance in given units
   */
   real recoveryHorizontalDistance(const std::string& u) const;

  /**
   * @return Recovery vertical distance in internal units, i.e., [m]
   */
   real recoveryVerticalDistance() const;

  /**
   * @return Recovery vertical distance in given units
   */
   real recoveryVerticalDistance(const std::string& u) const;

  /**
   * @return True if recovery bands are computed.
   */
   bool recoveryBandsComputed() const;

  /**
   * @return True if recovery are computed, but they are saturated.
   */
   bool recoveryBandsSaturated() const;

   std::string toString() const;

   std::string toStringUnits(const std::string& hunits,const std::string& vunits) const;

   std::string toPVS() const;
};
}

#endif
