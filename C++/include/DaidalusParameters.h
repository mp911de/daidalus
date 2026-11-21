/*
 * Copyright (c) 2015-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */
#ifndef DAIDALUSPARAMETERS_H_
#define DAIDALUSPARAMETERS_H_
#include "DaidalusFeatures.h"
#include "Real.h"

#include "ParameterData.h"
#include "ParameterAcceptor.h"
#include "ErrorReporter.h"
#include "Alerter.h"
#include "Function.h"

#include <map>
#include <vector>

namespace larcfm {

class DaidalusParameters : public ErrorReporter, public ParameterAcceptor {

private:

	/**
	 * This variable overrides the unit returned by the method getUnitsOf so that it returns 
	 * internal units (SI) instead of the specified unit in the DaidalusParameter object.
	 * Index 0: angle or angular rate
	 * Index 1: distance
	 * Index 2: speed
	 * Index 3: acceleration
	 * Index 4: time 
	 * Setting the value of an index to false disables the overriding.
	 */
	static bool SI_[];

  std::vector<Alerter> alerters_;

  ErrorLog error;

  // Bands
  real lookahead_time_; // [s] Lookahead time
  real distance_filter_; // Distance filter (respect to the ownship)
	real altitude_filter_; // Altitude filter (respect to the ownship)
  real left_hdir_;      // Left horizontal direction [0 - pi]
  real right_hdir_;     // Right horizontal direction [0 - pi]
  real min_hs_;         // Minimum horizontal speed
  real max_hs_;         // Maximum horizontal speed
  real min_vs_;         // Minimum vertical speed
  real max_vs_;         // Maximum vertical speed
  real min_alt_;        // Minimum altitude
  real max_alt_;        // Maximum altitude

  // Relative bands
  // The following values specify above and below values for the computation of bands
  // relative to the current value. In general, if above >= 0 and below >= 0 the bands range is
  // [val-below,val+above], where val is current value. The following are special cases,
  // if below < 0 && above >= 0 then [min,val+above]
  // if below >=0 && above < 0  then [val-below,max]
  // if below == 0 && above == 0] then [min,max]
  real below_relative_hs_;
  real above_relative_hs_;
  real below_relative_vs_;
  real above_relative_vs_;
  real below_relative_alt_;
  real above_relative_alt_;

  // Kinematic bands
  real step_hdir_; // Direction step
  real step_hs_;  // Horizontal speed step
  real step_vs_;  // Vertical speed step
  real step_alt_; // Altitude step
  real horizontal_accel_; // Horizontal acceleration
  real vertical_accel_; // Vertical acceleration
  real turn_rate_; // Turn rate
  real bank_angle_; // Bank angles (only used when turn_rate is 0)
  real vertical_rate_; // Vertical rate

  // Recovery bands
  real min_horizontal_recovery_; // Horizontal distance protected during recovery. TCAS RA DMOD is used this value is 0
  real min_vertical_recovery_; // Vertical distance protected during recovery. TCAS RA ZTHR is used when this value is 0
  bool recovery_hdir_;
  bool recovery_hs_;
  bool recovery_vs_;
  bool recovery_alt_;

  // Collision avoidance bands 
  bool   ca_bands_;  // When true, compute recovery bands until NMAC
  real ca_factor_; // Reduction factor when computing CA bands. It;s a value in [0,1]
  real horizontal_nmac_; // Horizontal Near Mid-Air Collision
  real vertical_nmac_; // Vertical Near Mid-Air Collision

  // Hysteresis and persistence
  real recovery_stability_time_;    // Recovery stability time
  real hysteresis_time_;            // Resolution hysteresis time
  real persistence_time_;           // Persistence time
  bool   bands_persistence_;          // Bands persistence enabled/disabled
  real persistence_preferred_hdir_; // Persistence for preferred horizontal direction resolution
  real persistence_preferred_hs_;   // Persistence for preferred horizontal speed resolution
  real persistence_preferred_vs_;   // Persistence for preferred Max delta for vertical speed resolution
  real persistence_preferred_alt_;  // Persistence for preferred altitude resolution resolution
  int    alerting_m_; // Alerting M parameter of M of N
  int    alerting_n_; // Alerting N parameter of M of N

  // Implicit coordination -- EXPERIMENTAL FEATURE
  bool conflict_crit_; /* Use criteria for conflict bands */
  bool recovery_crit_; /* Use criteria for recovery bands */

  // Sensor Uncertainty Mitigation
  real h_pos_z_score_; // Number of horizontal position standard deviations
  real h_vel_z_score_min_; // Minimum number of horizontal velocity standard deviations
  real h_vel_z_score_max_; // Maximum number of horizontal velocity standard deviations
  real h_vel_z_distance_; // Distance at which h_vel_z_score scales from min to max as range decreases
  real v_pos_z_score_; // Number of vertical position standard deviations
  real v_vel_z_score_; // Number of vertical velocity standard deviations

  // Contours
  real contour_thr_; // Horizontal threshold, specified as an angle to the left/right of current aircraft direction,
  // for computing horizontal contours. A value of 0 means only conflict contours. A value of pi means all contours.

  // DAA Terminal Area (DTA)

  /**
   * DTA Logic:
   * 0: Disabled
   * 1: Enabled special DTA maneuver guidance. Horizontal recovery is enabled,
   * but vertical recovery blocks down resolutions when alert is higher than corrective.
   * -1: Enabled special DTA maneuver guidance. Horizontal recovery is disabled,
   * vertical recovery blocks down resolutions when alert is higher than corrective.
   * NOTE:
   * When DTA logic is enabled, DAIDALUS automatically switches to DTA alerter and to
   * special maneuver guidance, when aircraft enters DTA volume (depending on ownship- vs
   * intruder-centric logic).
   */
  int dta_logic_;
  real dta_latitude_;
  real dta_longitude_;
  Position dta_position_;
  real dta_radius_;
  real dta_height_;
  int dta_alerter_;

  // Airspeed threshold for calculation of horizontal direction bands
	real airspeed_thr_;
	/** 
	 * 0: Special horizontal dir bands logic when airspeed is below airspeed_thr is disabled
	 * 1: Instantaneous horizontal direction bands computed assuming airspeed_thr
	 * -1: Kinematic horizontal direction bands computed assumming airspeed_thr
	 */
	int hdir_bands_below_as_thr_;

  // Alerting logic
  // When true, alerting logic is ownship-centric. Otherwise, it is intruder-centric.
  bool ownship_centric_alerting_; 

	// By default the parameter bands_add_time_to_maneuver is set to false. When set to true, 
	// the time to maneuver is considerd against alerting time when computing peripheral bands. 
  // The latter was the behavior in v.2.0.1 and v1.
	bool bands_add_time_to_maneuver_;

  // Corrective region for all alerters.
  // (IMPORTANT: all alerters should declare an alert level with a corrective region!)
  BandsRegion::Region corrective_region_; 

  std::map<std::string,std::string> units_;

  static bool initialized;
  static void init();

  /** This method is needed because WCV_TAUMOD_SUM doesn't require the
   *  user to initialize SUM parameters, which may be specified globally.
   */
  void set_alerter_with_SUM_parameters(const Alerter& alerter) ;

  /** This method is needed because WCV_TAUMOD_SUM doesn't require the
   *  user to initialize SUM parameters, which may be specified globally.
   */
  void updateDTAPosition();
  void set_alerters_with_SUM_parameters();

  std::string units_of_daa_key(const std::string& key) const;

  std::string cannonical_key(const std::string& keystr) const;

public:

  /**
   * DAIDALUS version
   */
  static const std::string VERSION;
  static const INT64FM ALMOST_;

  /**
	 * Experimental features are turned off. They are not maintained and will
	 * eventually get removed from the code base.
	 */
	static const bool EXPERIMENTAL_FEATURES;

  static void overrideUnitsWithSI(bool angle, bool distance, bool speed, bool acceleration, bool time);

  static void overrideUnitsWithSI();

  static std::string get_overriden_unit(const std::string& units);

  /**
   * The following method set default output precision and enable/disable trailing zeros.
   * It doesn't affect computations
   */
  static void setDefaultOutputPrecision(int precision, bool trailing_zeros = false);

  DaidalusParameters();

  ~DaidalusParameters();

  /* Setting for WC Definitions RTCA DO-365 */

  /*
   * Set DaidalusParameters object such that
   * - Configure two alerters (Phase I and Phase II) as defined as in RTCA DO-365A
   * - Maneuver guidance logic assumes kinematic maneuvers
   * - Turn rate is set to 3 deg/s, when type is true, and to  1.5 deg/s
   *   when type is false.
   * - Configure Sensor Uncertainty Migitation (SUM) when sum is true
   * - Bands don't saturate until NMAC
   */
  void set_DO_365A(bool type=true, bool sum=true);

  /*
   * Set DaidalusParameters object such that
   * - Configure two alerters (Phase I, Phase II, and Non-Cooperative) as defined as in RTCA DO-365B
   * - Maneuver guidance logic assumes kinematic maneuvers
   * - Turn rate is set to 3 deg/s, when type is true, and to  1.5 deg/s
   *   when type is false.
   * - Configure Sensor Uncertainty Migitation (SUM) when sum is true
   * - Bands don't saturate until NMAC
   */
  void set_DO_365B(bool type=true, bool sum=true);

  /*
   * Set DaidalusParameters object such that
   * - Alerting thresholds are buffered
   * - Maneuver guidance logic assumes kinematic maneuvers
   * - Turn rate is set to 3 deg/s, when type is true, and to  1.5 deg/s
   *   when type is false.
   * - Bands don't saturate until NMAC
   */
  void set_Buffered_WC_DO_365(bool type);

  /* Set DaidalusParameters object such that alerting logic and maneuver guidance corresponds to
   * ACCoRD's CD3D, i.e.,
   * - Separation is given by a cylinder of of diameter 5nm and height 1000ft
   * - Lookahead time and alerting time is 180s
   * - Only 1 alert level
   * - Instantaneous maneuvers */
  void set_CD3D();

  /* Set DaidalusParameters object such that alerting logic and maneuver guidance corresponds to
   * ideal TCASII */
  void set_TCASII();

  typedef std::map< std::string,std::vector<std::string> > aliasestype;
  static const aliasestype& getAliases(); // List of key aliases

  int numberOfAlerters() const;

  const Alerter& getAlerterAt(int i) const;

  int getAlerterIndex(const std::string& id) const;

  void clearAlerters();

  /**
   * Add alerter (if id of alerter already exists, replaces alerter with new one).
   * Return index of added alerter
   */
  int addAlerter(const Alerter& alerter);

  /** 
   * @return lookahead time in seconds.
   */
  real getLookaheadTime() const;

  /** 
   * @return lookahead time in specified units [u].
   */
  real getLookaheadTime(const std::string& u) const;
  
  /** 
   * @return distance filter in meters.
   */
	real getDistanceFilter() const;

  /** 
   * @return distance filter in specified units [u].
   */
	real getDistanceFilter(const std::string& u) const;

  /** 
   * @return altitude filter in meters.
   */
	real getAltitudeFilter() const;

  /** 
   * @return altitude filter in specified units [u].
   */
	real getAltitudeFilter(const std::string& u) const;

  real getLeftHorizontalDirection() const;

  real getLeftHorizontalDirection(const std::string& u) const;

  real getRightHorizontalDirection() const;

  real getRightHorizontalDirection(const std::string& u) const;

  /** 
   * @return airspeed threshold in internal units [m/s].
   */
  real getAirspeedThreshold() const;

  /** 
   * @return airspeed threshold in specified units [u].    
   */
  real getAirspeedThreshold(const std::string& u) const;

  real getMinHorizontalSpeed() const;

  real getMinHorizontalSpeed(const std::string& u) const;

  real getMaxHorizontalSpeed() const;

  real getMaxHorizontalSpeed(const std::string& u) const;

  real getMinVerticalSpeed() const;

  real getMinVerticalSpeed(const std::string& u) const;

  real getMaxVerticalSpeed() const;

  real getMaxVerticalSpeed(const std::string& u) const;

  real getMinAltitude() const;

  real getMinAltitude(const std::string& u) const;

  real getMaxAltitude() const;

  real getMaxAltitude(const std::string& u) const;

  real getBelowRelativeHorizontalSpeed() const;

  real getBelowRelativeHorizontalSpeed(const std::string& u) const;

  real getAboveRelativeHorizontalSpeed() const;

  real getAboveRelativeHorizontalSpeed(const std::string& u) const;

  real getBelowRelativeVerticalSpeed() const;

  real getBelowRelativeVerticalSpeed(const std::string& u) const;

  real getAboveRelativeVerticalSpeed() const;

  real getAboveRelativeVerticalSpeed(const std::string& u) const;

  real getBelowRelativeAltitude() const;

  real getBelowRelativeAltitude(const std::string& u) const;

  real getAboveRelativeAltitude() const;

  real getAboveRelativeAltitude(const std::string& u) const;

  real getHorizontalDirectionStep() const;

  real getHorizontalDirectionStep(const std::string& u) const;

  real getHorizontalSpeedStep() const;

  real getHorizontalSpeedStep(const std::string& u) const;

  real getVerticalSpeedStep() const;

  real getVerticalSpeedStep(const std::string& u) const;

  real getAltitudeStep() const;

  real getAltitudeStep(const std::string& u) const;

  real getHorizontalAcceleration() const;

  real getHorizontalAcceleration(const std::string& u) const;

  real getVerticalAcceleration() const;

  real getVerticalAcceleration(const std::string& u) const;

  real getTurnRate() const;

  real getTurnRate(const std::string& u) const;

  real getBankAngle() const;

  real getBankAngle(const std::string& u) const;

  real getVerticalRate() const;

  real getVerticalRate(const std::string& u) const;

  real getHorizontalNMAC() const;

  real getHorizontalNMAC(const std::string& u) const;

  real getVerticalNMAC() const;

  real getVerticalNMAC(const std::string& u) const;

  real getRecoveryStabilityTime() const;

  real getRecoveryStabilityTime(const std::string& u) const;

  /**
   * @return hysteresis time in seconds.
   */
  real getHysteresisTime() const;

  /**
   * @return hysteresis time in specified units [u]
   */
  real getHysteresisTime(const std::string& u) const;

  real getPersistenceTime() const;

  real getPersistenceTime(const std::string& u) const;

  bool   isEnabledBandsPersistence() const;

  void   setBandsPersistence(bool flag);

  void   enableBandsPersistence();

  void   disableBandsPersistence();

  real getPersistencePreferredHorizontalDirectionResolution() const;

  real getPersistencePreferredHorizontalDirectionResolution(const std::string& u) const;

  real getPersistencePreferredHorizontalSpeedResolution() const;

  real getPersistencePreferredHorizontalSpeedResolution(const std::string& u) const;

  real getPersistencePreferredVerticalSpeedResolution() const;

  real getPersistencePreferredVerticalSpeedResolution(const std::string& u) const;

  real getPersistencePreferredAltitudeResolution() const;

  real getPersistencePreferredAltitudeResolution(const std::string& u) const;

  /**
   * @return Alerting parameter m of "M of N" strategy
   */
  int getAlertingParameterM() const;

  /**
   * @return Alerting parameter m of "M of N" strategy
   */
  int getAlertingParameterN() const;

  real getMinHorizontalRecovery() const;

  real getMinHorizontalRecovery(const std::string& u) const;

  real getMinVerticalRecovery() const;

  real getMinVerticalRecovery(const std::string& u) const;

  bool setLookaheadTime(real val);

  bool setLookaheadTime(real val, const std::string& u);

  /** 
   * Set distance filter to value in meters.
   */
	bool setDistanceFilter(real val);

  /** 
   * Set distance filter to value in specified units [u].
   */
  bool setDistanceFilter(real val, const std::string& u);
  
  /** 
   * Set altitude filter to value in meters.
   */
	bool setAltitudeFilter(real val);

  /** 
   * Set altitude filter to value in specified units [u].
   */
  bool setAltitudeFilter(real val, const std::string& u);

  bool setLeftHorizontalDirection(real val);

  bool setLeftHorizontalDirection(real val, const std::string& u);

  bool setRightHorizontalDirection(real val);

  bool setRightHorizontalDirection(real val, const std::string& u);

	/** 
   * Set airspeed threshold to value in internal units [m/s].
   */
	bool setAirspeedThreshold(real val);

  /** 
   * Set airspeed threshold to value in specified units [u].
   */
	bool setAirspeedThreshold(real val, const std::string& u);
    
  /**
   * Set minimum horizontal speed to value in internal units [m/s].
   * Minimum horizontal speed must be non-negative.
  */  
  bool setMinHorizontalSpeed(real val);

  /**
   * Set minimum horizontal speed to value in specified units.
   * Minimum horizontal speed must be non-negative.
  */
  bool setMinHorizontalSpeed(real val, const std::string& u);

  bool setMaxHorizontalSpeed(real val);

  bool setMaxHorizontalSpeed(real val, const std::string& u);

  bool setMinVerticalSpeed(real val);

  bool setMinVerticalSpeed(real val, const std::string& u);

  bool setMaxVerticalSpeed(real val);

  bool setMaxVerticalSpeed(real val, const std::string& u);

  bool setMinAltitude(real val);

  bool setMinAltitude(real val, const std::string& u);

  bool setMaxAltitude(real val);

  bool setMaxAltitude(real val, const std::string& u);

  bool setBelowRelativeHorizontalSpeed(real val);

  bool setBelowRelativeHorizontalSpeed(real val, const std::string& u);

  bool setAboveRelativeHorizontalSpeed(real val);

  bool setAboveRelativeHorizontalSpeed(real val, const std::string& u);

  bool setBelowRelativeVerticalSpeed(real val);

  bool setBelowRelativeVerticalSpeed(real val, const std::string& u);

  bool setAboveRelativeVerticalSpeed(real val);

  bool setAboveRelativeVerticalSpeed(real val, const std::string& u);

  bool setBelowRelativeAltitude(real val);

  bool setBelowRelativeAltitude(real val, const std::string& u);

  bool setAboveRelativeAltitude(real val);

  bool setAboveRelativeAltitude(real val, const std::string& u);

  void setBelowToMinRelativeHorizontalSpeed();

  void setAboveToMaxRelativeHorizontalSpeed();

  void setBelowToMinRelativeVerticalSpeed();

  void setAboveToMaxRelativeVerticalSpeed();

  void setBelowToMinRelativeAltitude();

  void setAboveToMaxRelativeAltitude();

  void disableRelativeHorizontalSpeedBands();

  void disableRelativeVerticalSpeedBands();

  void disableRelativeAltitude();

  bool setHorizontalDirectionStep(real val);

  bool setHorizontalDirectionStep(real val, const std::string& u);

  bool setHorizontalSpeedStep(real val);

  bool setHorizontalSpeedStep(real val, const std::string& u);

  bool setVerticalSpeedStep(real val);

  bool setVerticalSpeedStep(real val, const std::string& u);

  bool setAltitudeStep(real val);

  bool setAltitudeStep(real val, const std::string& u);

  bool setHorizontalAcceleration(real val);

  bool setHorizontalAcceleration(real val, const std::string& u);

  bool setVerticalAcceleration(real val);

  bool setVerticalAcceleration(real val, const std::string& u);

private:
  bool set_turn_rate(real val);

  bool set_bank_angle(real val);

  bool set_alerting_parameterM(int m);

  bool set_alerting_parameterN(int n);

public:
  bool setTurnRate(real val);

  bool setTurnRate(real val, const std::string& u);

  bool setBankAngle(real val);

  bool setBankAngle(real val, const std::string& u);

  bool setVerticalRate(real val);

  bool setVerticalRate(real val, const std::string& u);

  bool setHorizontalNMAC(real val);

  bool setHorizontalNMAC(real val, const std::string& u);

  bool setVerticalNMAC(real val);

  bool setVerticalNMAC(real val, const std::string& u);

  bool setRecoveryStabilityTime(real val);

  bool setRecoveryStabilityTime(real val, const std::string& u);

  bool setHysteresisTime(real val);

  bool setHysteresisTime(real val, const std::string& u);

  bool setPersistenceTime(real val);

  bool setPersistenceTime(real val, const std::string& u);

  bool setPersistencePreferredHorizontalDirectionResolution(real val);

  bool setPersistencePreferredHorizontalDirectionResolution(real val, const std::string& u);

  bool setPersistencePreferredHorizontalSpeedResolution(real val);

  bool setPersistencePreferredHorizontalSpeedResolution(real val, const std::string& u);

  bool setPersistencePreferredVerticalSpeedResolution(real val);

  bool setPersistencePreferredVerticalSpeedResolution(real val, const std::string& u);

  bool setPersistencePreferredAltitudeResolution(real val);

  bool setPersistencePreferredAltitudeResolution(real val, const std::string& u);

  /**
   * Set alerting parameters of M of N strategy
   */
  bool setAlertingMofN(int m, int n);

  bool setMinHorizontalRecovery(real val);

  bool setMinHorizontalRecovery(real val, const std::string& u);

  bool setMinVerticalRecovery(real val);

  bool setMinVerticalRecovery(real val, const std::string& u);

  /**
   * @return true if repulsive criteria is enabled for conflict bands -- EXPERIMENTAL FEATURE
   */
  bool isEnabledConflictCriteria() const;

  /**
   * Enable/disable repulsive criteria for conflict bands -- EXPERIMENTAL FEATURE
   */
  void setConflictCriteria(bool flag);

  /**
   * Enable repulsive criteria for conflict bands -- EXPERIMENTAL FEATURE
   */
  void enableConflictCriteria();

  /**
   * Disable repulsive criteria for conflict bands -- EXPERIMENTAL FEATURE
   */
  void disableConflictCriteria();

  /**
   * @return true if repulsive criteria is enabled for recovery bands -- EXPERIMENTAL FEATURE
   */
  bool isEnabledRecoveryCriteria() const;

  /**
   * Enable/disable repulsive criteria for recovery bands -- EXPERIMENTAL FEATURE
   */
  void setRecoveryCriteria(bool flag);

  /**
   * Enable repulsive criteria for recovery bands -- EXPERIMENTAL FEATURE
   */
  void enableRecoveryCriteria();
  /**
   * Disable repulsive criteria for recovery bands -- EXPERIMENTAL FEATURE
   */
  void disableRecoveryCriteria();

  /**
   * Enable/disable repulsive criteria for conflict and recovery bands -- EXPERIMENTAL FEATURE
   */
  void setRepulsiveCriteria(bool flag);

  /**
   * Enable repulsive criteria for conflict and recovery bands -- EXPERIMENTAL FEATURE
   */
  void enableRepulsiveCriteria();

  /**
   * Disable repulsive criteria for conflict and recovery bands -- EXPERIMENTAL FEATURE
   */
  void disableRepulsiveCriteria();

  /**
   * @return recovery bands flag for track bands.
   */
  bool isEnabledRecoveryHorizontalDirectionBands() const;

  /**
   * @return recovery bands flag for ground speed bands.
   */
  bool isEnabledRecoveryHorizontalSpeedBands() const;

  /**
   * @return true if recovery bands for vertical speed bands is enabled.
   */
  bool isEnabledRecoveryVerticalSpeedBands() const;

  /**
   * @return true if recovery bands for altitude bands is enabled.
   */
  bool isEnabledRecoveryAltitudeBands() const;

  /**
   * Enable/disable recovery bands for track, ground speed, vertical speed, and altitude.
   */
  void setRecoveryBands(bool flag);

  /**
   * Enable all recovery bands for track, ground speed, vertical speed, and altitude.
   */ 
  void enableRecoveryBands();

  /** 
   * Disable all recovery bands for track, ground speed, vertical speed, and altitude.
   */
  void disableRecoveryBands();

  /**
   * Sets recovery bands flag for track bands to value.
   */ 
  void setRecoveryHorizontalDirectionBands(bool flag);

  /** 
   * Sets recovery bands flag for ground speed bands to value.
   */ 
  void setRecoveryHorizontalSpeedBands(bool flag);

  /** 
   * Sets recovery bands flag for vertical speed bands to value.
   */ 
  void setRecoveryVerticalSpeedBands(bool flag);

  /**
   * Sets recovery bands flag for altitude bands to value.
   */
  void setRecoveryAltitudeBands(bool flag);

  /**
   * @return true if collision avoidance bands are enabled.
   */
  bool isEnabledCollisionAvoidanceBands() const;

  /**
   * Enable/disable collision avoidance bands.
   */
  void setCollisionAvoidanceBands(bool flag);

  /**
   * Enable collision avoidance bands.
   */
  void enableCollisionAvoidanceBands();

  /**
   * Disable collision avoidance bands.
   */
  void disableCollisionAvoidanceBands();

  /**
   * @return get factor for computing collision avoidance bands. Factor value is in (0,1]
   */
  real getCollisionAvoidanceBandsFactor() const;

  /**
   * @return set factor for computing collision avoidance bands. Factor value is in (0,1]
   */
  bool setCollisionAvoidanceBandsFactor(real val);

  real getHorizontalPositionZScore() const;

  bool setHorizontalPositionZScore(real val);

  real getHorizontalVelocityZScoreMin() const;

  bool setHorizontalVelocityZScoreMin(real val);

  real getHorizontalVelocityZScoreMax() const;

  bool setHorizontalVelocityZScoreMax(real val);

  real getHorizontalVelocityZDistance() const;

  bool setHorizontalVelocityZDistance(real val);

  real getHorizontalVelocityZDistance(const std::string& u) const;

  bool setHorizontalVelocityZDistance(real val, const std::string& u);

  real getVerticalPositionZScore() const;

  bool setVerticalPositionZScore(real val);

  real getVerticalSpeedZScore() const;

  bool setVerticalSpeedZScore(real val);

  /**
   * Get horizontal contour threshold, specified in internal units [rad] as an angle to
   * the left/right of current aircraft direction. A value of 0 means only conflict contours.
   * A value of pi means all contours.
   */
  real getHorizontalContourThreshold() const;

  /**
   * Get horizontal contour threshold, specified in given units [u] as an angle to
   * the left/right of current aircraft direction. A value of 0 means only conflict contours.
   * A value of pi means all contours.
   */
  real getHorizontalContourThreshold(const std::string& u) const;

  /**
   * Set horizontal contour threshold, specified in internal units [rad] [0 - pi] as an angle to
   * the left/right of current aircraft direction. A value of 0 means only conflict contours.
   * A value of pi means all contours.
   */
  bool setHorizontalContourThreshold(real val);

  /**
   * Set horizontal contour threshold, specified in given units [u] [0 - pi] as an angle to
   * the left/right of current aircraft direction. A value of 0 means only conflict contours.
   * A value of pi means all contours.
   */
  bool setHorizontalContourThreshold(real val, const std::string& u);

  /**
   * Get DTA Logic:
   * 0: Disabled
   * 1: Enabled special DTA maneuver guidance. Horizontal recovery is enabled,
   * but vertical recovery blocks down resolutions when alert is higher than corrective.
   * -1: Enabled special DTA maneuver guidance. Horizontal recovery is disabled,
   * vertical recovery blocks down resolutions when alert is higher than corrective.
   * NOTE:
   * When DTA logic is enabled, DAIDALUS automatically switches to DTA alerter and to
   * special maneuver guidance, when aircraft enters DTA volume (depending on ownship- vs
   * intruder-centric logic).
   */
  int getDTALogic() const;

  /**
   * Set DTA Logic:
   * 0: Disabled
   * 1: Enabled special DTA maneuver guidance. Horizontal recovery is enabled,
   * but vertical recovery blocks down resolutions when alert is higher than corrective.
   * -1: Enabled special DTA maneuver guidance. Horizontal recovery is disabled,
   * vertical recovery blocks down resolutions when alert is higher than corrective.
   * NOTE:
   * When DTA logic is enabled, DAIDALUS automatically switches to DTA alerter and to
   * special maneuver guidance, when aircraft enters DTA volume (depending on ownship- vs
   * intruder-centric logic).
   */
  void setDTALogic(int dta_logic);

  /**
   * Get DAA Terminal Area (DTA) position (lat/lon)
   */
  const Position& getDTAPosition() const;

  /**
   * Set DAA Terminal Area (DTA) latitude (internal units)
   */
  bool setDTALatitude(real lat);

  /**
   * Set DAA Terminal Area (DTA) latitude in given units
   */
  bool setDTALatitude(real lat, const std::string& ulat);

  /**
   * Set DAA Terminal Area (DTA) longitude (internal units)
   */
  bool setDTALongitude(real lon);

  /**
   * Set DAA Terminal Area (DTA) longitude in given units
   */
  bool setDTALongitude(real lon, const std::string& ulon);

  /**
   * Get DAA Terminal Area (DTA) radius (internal units)
   */
  real getDTARadius() const;

  /**
   * Get DAA Terminal Area (DTA) radius in given units
   */
  real getDTARadius(const std::string& u) const;

  /**
   * Set DAA Terminal Area (DTA) radius (internal units)
   */
  bool setDTARadius(real val);

  /**
   * Set DAA Terminal Area (DTA) radius in given units
   */
  bool setDTARadius(real val, const std::string& u);

  /**
   * Get DAA Terminal Area (DTA) height (internal units)
   */
  real getDTAHeight() const;

  /**
   * Get DAA Terminal Area (DTA) height in given units
   */
  real getDTAHeight(const std::string& u) const;

  /**
   * Set DAA Terminal Area (DTA) height (internal units)
   */
  bool setDTAHeight(real val);

  /**
   * Set DAA Terminal Area (DTA) height in given units
   */
  bool setDTAHeight(real val, const std::string& u);

  /**
   * Get DAA Terminal Area (DTA) alerter
   */
  int getDTAAlerter() const;

  /**
   * Set DAA Terminal Area (DTA) alerter
   */
  void setDTAAlerter(int alerter);

	/**
	 * Get Horizontal Direction Bands Logic When Below Airspeed Threshold: 
	 * 0: Special horizontal dir bands logic when airspeed is below airspeed_thr is disabled
	 * 1: Instantaneous horizontal direction bands computed assuming airspeed_thr 
	 * -1: Kinematic horizontal direction bands computed assumming airspeed_thr
	 */
  int getHorizontalDirBandsBelowAirspeedThr() const;

	/**
	 * Set Horizontal Direction Bands Logic When Below Airspeed Threshold: 
	 * 0: Special horizontal dir bands logic when airspeed is below airspeed_thr is disabled
	 * 1: Instantaneous horizontal direction bands computed assuming airspeed_thr 
	 * -1: Kinematic horizontal direction bands computed assumming airspeed_thr
	 */
  void setHorizontalDirBandsBelowAirspeedThr(int val);

  /**
   * Set alerting logic to the value indicated by ownship_centric.
   * If ownship_centric is true, alerting and guidance logic will use the alerter in ownship. Alerter
   * in every intruder will be disregarded.
   * If ownship_centric is false, alerting and guidance logic will use the alerter in every intruder. Alerter
   * in ownship will be disregarded.
   */
  void setAlertingLogic(bool ownship_centric);

  /**
   * Set alerting and guidance logic to ownship-centric. Alerting and guidance logic will use the alerter in ownship.
   * Alerter in every intruder will be disregarded.
   */
  void setOwnshipCentricAlertingLogic();

  /**
   * Set alerting and guidance logic to intruder-centric. Alerting and guidance logic will use the alerter in every intruder.
   * Alerter in ownship will be disregarded.
   */
  void setIntruderCentricAlertingLogic();

  /**
   * @return true if alerting/guidance logic is ownship centric.
   */
  bool isAlertingLogicOwnshipCentric() const;

  /**
   * Set logic of "bands add time to manevuer" to the value indicated by bands_add_time_to_maneuver
   * By default this parameter is set to false. When set to true, the time to maneuver is
   * considerd against alerting time when computing peripheral bands. The latter was the behavior in
   * v.2.0.1 and v1.
   */
  void setBandsAddTimeToManeuver(bool bands_add_time_to_maneuver);

  /**
   * Set logic of "bands add time to manevuer" to true
   */
  void enableBandsAddTimeToManeuver();

  /**
   * Set logic of "bands add time to manevuer" to false
   */
  void disableBandsAddTimeToManeuver();

  /**
   * @return true if peripheral bands logic adds time to maneuver (as used in v1 and v2.0.1)
   */
  bool isEnabledBandsAddTimeToManeuver() const;

  /**
   * Get corrective region for calculation of resolution maneuvers and bands saturation.
   */
  BandsRegion::Region getCorrectiveRegion() const;

  /**
   * Set corrective region for calculation of resolution maneuvers and bands saturation.
   */
  bool setCorrectiveRegion(BandsRegion::Region val);

  /**
   * @param alerter_idx Indice of an alerter (starting from 1)
   * @return corrective level of alerter at alerter_idx. The corrective level
   * is the first alert level that has a region equal to or more severe than corrective_region.
   * Return -1 if alerter_idx is out of range or if there is no corrective alert level
   * for this alerter.
   */
  int correctiveAlertLevel(int alerter_idx);

  /**
   * @return maximum number of alert levels for all alerters. Returns 0 if alerter list is empty.
   */
  int maxNumberOfAlertLevels() const;

  /**
   * Set instantaneous bands.
   */
  void setInstantaneousBands() ;

  /**
   * Set kinematic bands.
   * Set turn rate to 3 deg/s, when type is true; set turn rate to  1.5 deg/s
   * when type is false;
   */
  void setKinematicBands(bool type);

  void disableHysteresis();

#if DAIDALUS_FILE_IO
  /**
   *  Load parameters from file.
   */
  bool loadFromFile(const std::string& file);

  /**
   *  Write parameters to file.
   */
  bool saveToFile(const std::string& file);
#endif

  std::string toString() const;

  std::string toPVS() const;

  ParameterData getParameters() const;

  void updateParameterData(ParameterData& parameters) const;

  /**
    @Deprecated
   * Use maxNumberOfAlertLevels instead
   */
  virtual int maxAlertLevel() const {
    return maxNumberOfAlertLevels();
  }


private:

  /**
   * Return a ParameterData suitable to be read by readAlerterList() based on the supplied Alerter instances.
   * @param alerters list of Alerter instances
   * @param ordered true to modify identifiers to ensure they retain the input list's ordering when decoded, false will only modify identifiers if they are not unique
   * This is a cosmetic alteration to allow for the string representation to have parameters grouped together.
   */
  void writeAlerterList(ParameterData& p) const;

  /**
   * This parses parameters involved with a list of Alerter objects, if properly specified.  The list will be sorted by the instance identifiers.
   *
   * Specifically this looks for parameters:
   * alerterNameList : create Alerter objects for each item in the list
   * XXX_* : parameters associated with Alerter with id "XXX"
   *
   * The user still needs to assign these to the appropriate object(s).
   * If no alternates are loaded, return the empty list.
   *
   * If verbose is false (default true), suppress all status messages except exceptions
   */
  void readAlerterList(const std::vector<std::string>& alerter_list, const ParameterData& params);

  static bool contains(const ParameterData& p,const std::string& key);

  static std::string getUnit(const ParameterData& p,const std::string& key);

  static real getValue(const ParameterData& p,const std::string& key);

  static bool getBool(const ParameterData& p,const std::string& key);

  static int getInt(const ParameterData& p,const std::string& key);

  static std::string getString(const ParameterData& p,const std::string& key);

  static std::vector<std::string> getListString(const ParameterData& p,const std::string& key);

public:
  bool setParameterData(const ParameterData& p);

  void setParameters(const ParameterData& p);

  /** 
   * Get string units of parameters keys
   */
  std::string getUnitsOf(const std::string& key) const;

  /** 
   * Set string units of parameters keys. Return false if units are not compatible
   */
  bool setUnitsOf(const std::string& keystr, const std::string& units);

  bool hasError() const;

  bool hasMessage() const;

  std::string getMessage();

  std::string getMessageNoClear() const;

};

template <typename T>T parameter_data(const ParameterData& p, const std::string& key, Function<const std::string&,T> &f){
  T default_val = f.apply(key);
  if (p.contains(key)) {
    return default_val;
  }
  DaidalusParameters::aliasestype::const_iterator find_ptr = DaidalusParameters::getAliases().find(key);
  if (find_ptr == DaidalusParameters::getAliases().end()) {
    return default_val;
  }
  std::vector<std::string>::const_iterator alias_ptr;
  for (alias_ptr = find_ptr->second.begin(); alias_ptr != find_ptr->second.end(); ++alias_ptr) {
    if (p.contains(*alias_ptr)) {
      return f.apply(*alias_ptr);
    }
  }
  return default_val;
}

}
#endif
