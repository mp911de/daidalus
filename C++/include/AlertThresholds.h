/*
 * Copyright (c) 2015-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */
#ifndef ALERTTHRESHOLDS_H_
#define ALERTTHRESHOLDS_H_
#include "Real.h"

#include "Detection3D.h"
#include "CDCylinder.h"
#include "Vect3.h"
#include "Velocity.h"
#include "Interval.h"
#include "BandsRegion.h"
#include "ParameterData.h"
#include "ParameterTable.h"

#include <map>
#include <memory>

namespace larcfm {

class AlertThresholds : public ParameterTable {

private:
  std::unique_ptr<Detection3D> detector_; // State-based detector
  real  alerting_time_; // Alerting_time
  // If alerting_time > 0, alert is based on detection
  // If alerting_time = 0, alert is based on violation
  real early_alerting_time_; // Early alerting time (for maneuver guidance). If zero, same as alerting_time
  BandsRegion::Region region_;  // Guidance region for this alert
 	// Spread alerting -- EXPERIMENTAL FEATURE
  real spread_hdir_;  // Alert when direction band within spread (non-negative value)
  real spread_hs_;   // Alert when horizontal speed band within spread (non-negative value)
  real spread_vs_;   // Alert when vertical speed band within speed (non-negative value)
  real spread_alt_;  // Alert when altitude  band within spread (non-negative value)
  std::map<std::string,std::string> units_;

  void copyFrom(const AlertThresholds& athr);

public:

  static const AlertThresholds& INVALID();

  AlertThresholds();

  /** 
   * Beware of this constructor. From this point on, AlertThresholds owns
   * this Detection3D pointer.
   */
  explicit AlertThresholds(Detection3D* det);

  /**
   * Creates an alert threholds object. Parameter det is a detector,
   * alerting_time is a non-negative alerting time (possibly positive infinity),
   * early_alerting_time is a early alerting time >= at (for maneuver guidance),
   * region is the type of guidance
   */
  AlertThresholds(const Detection3D& det,
      real alerting_time, real early_alerting_time,
      BandsRegion::Region region);

  AlertThresholds(const AlertThresholds& athr);

  ~AlertThresholds();

  // needed because of pointer
  AlertThresholds& operator=(const AlertThresholds& athr);

  bool isValid() const;

  /**
   * Get detector
   */
  const Detection3D& getCoreDetection() const;

  /**
   * Set detector
   */
  void setCoreDetection(const Detection3D& det);

  /**
   * Return alerting time in seconds.
   */
  real getAlertingTime() const;

  /**
   * Return alerting time in specified units.
   */
  real getAlertingTime(const std::string& u) const;

  /**
   * Set alerting time in seconds. Alerting time is non-negative number.
   */
  void setAlertingTime(real val);

  /**
   * Set alerting time in specified units. Alerting time is non-negative number.
   */
  void setAlertingTime(real val, const std::string& u);

  /**
   * Return early alerting time in seconds.
   */
  real getEarlyAlertingTime() const;

  /**
   * Return early alerting time in specified units.
   */
  real getEarlyAlertingTime(const std::string& u) const;

  /**
   * Set early alerting time in seconds. Early alerting time is a positive number >= alerting time
   */
  void setEarlyAlertingTime(real val);

  /**
   * Set early alerting time in specified units. Early alerting time is a positive number >= alerting time
   */
  void setEarlyAlertingTime(real val, const std::string& u);

  /**
   * Return guidance region.
   */
  BandsRegion::Region getRegion() const;

  /**
   * Set guidance region.
   */
  void setRegion(BandsRegion::Region region);

  /**
   * Get track spread in internal units [rad] -- EXPERIMENTAL FEATURE
   * Spread is relative to ownship's track
   */
  real getHorizontalDirectionSpread() const;

  /**
   * Get track spread in given units [u] -- EXPERIMENTAL FEATURE
   * Spread is relative to ownship's track
   */
  real getHorizontalDirectionSpread(const std::string& u) const;

  /**
   * Set track spread in internal units -- EXPERIMENTAL FEATURE 
   * Spread is relative to ownship's track and is expected to be in [0,pi].
   */
  void setHorizontalDirectionSpread(real spread);

  /**
   * Set track spread in given units -- EXPERIMENTAL FEATURE 
   * Spread is relative to ownship's track and is expected to be in [0,pi] [u].
   */
  void setHorizontalDirectionSpread(real spread, const std::string& u);

  /**
   * Get ground speed spread in internal units [m/s] -- EXPERIMENTAL FEATURE
   * Spread is relative to ownship's ground speed
   */
  real getHorizontalSpeedSpread() const;

  /**
   * Get ground speed spread in given units -- EXPERIMENTAL FEATURE 
   * Spread is relative to ownship's ground speed
   */
  real getHorizontalSpeedSpread(const std::string& u) const;

  /**
   * Set ground speed spread in internal units [m/s] -- EXPERIMENTAL FEATURE
   * Spread is relative to ownship's ground speed and is expected to be non-negative
   */
  void setHorizontalSpeedSpread(real spread);

  /**
   * Set ground speed spread in given units -- EXPERIMENTAL FEATURE
   * Spread is relative to ownship's ground speed and is expected to be non-negative
   */
  void setHorizontalSpeedSpread(real spread, const std::string& u);

  /**
   * Get vertical speed spread in internal units [m/s] -- EXPERIMENTAL FEATURE
   * Spread is relative to ownship's vertical speed
   */
  real getVerticalSpeedSpread() const;

  /**
   * Get vertical speed spread in given units -- EXPERIMENTAL FEATURE 
   * Spread is relative to ownship's vertical speed
   */
  real getVerticalSpeedSpread(const std::string& u) const;

  /**
   * Set vertical speed spread in internal units [m/s] -- EXPERIMENTAL FEATURE
   * Spread is relative to ownship's vertical speed and is expected to be non-negative
   */
  void setVerticalSpeedSpread(real spread);

  /**
   * Set vertical speed spread in given units -- EXPERIMENTAL FEATURE
   * Spread is relative to ownship's vertical speed and is expected to be non-negative
   */
  void setVerticalSpeedSpread(real spread, const std::string& u);

  /**
   * Get altitude spread in internal units [m] -- EXPERIMENTAL FEATURE 
   * Spread is relative to ownship's altitude
   */
  real getAltitudeSpread() const;

  /**
   * Get altitude spread in given units -- EXPERIMENTAL FEATURE 
   * Spread is relative to ownship's altitude
   */
  real getAltitudeSpread(const std::string& u) const;

  /**
   * Set altitude spread in internal units [m] -- EXPERIMENTAL FEATURE
   * Spread is relative to ownship's altitude and is expected to be non-negative
   */
  void setAltitudeSpread(real spread);

  /**
   * Set altitude spread in given units -- EXPERIMENTAL FEATURE
   * Spread is relative to ownship's altitude and is expected to be non-negative
   */
  void setAltitudeSpread(real spread, const std::string& u);

  std::string toString() const;

  std::string toPVS() const;

  ParameterData getParameters() const;

  void updateParameterData(ParameterData& p) const;

  void setParameters(const ParameterData& p);

  std::string getUnits(const std::string& key) const;

};
}
#endif
