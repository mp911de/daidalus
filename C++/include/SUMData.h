/*
 * Copyright (c) 2015-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */

#ifndef SUMDATA_H_
#define SUMDATA_H_
#include "Real.h"

#include <string>

namespace larcfm {

/** 
 * It is assumed that for all alert level i >= 1: detection(i+1) => detection(i) and that there is at least one
 * guidance level whose region is different from NONE. 
 */
class SUMData {

private:
  real s_EW_std_;     // East/West position standard deviation in internal units
  real s_NS_std_;     // North/South position standard deviation in internal units
  real s_EN_std_;     // East/North position standard deviation in internal units
  real sz_std_;       // Vertical position standard deviation in internal units
  real v_EW_std_;     // East/West velocity standard deviation in internal units
  real v_NS_std_;     // North/South velocity standard deviation in internal units
  real v_EN_std_;     // East/North velocity standard deviation in internal units
  real vz_std_;       // Vertical velocity standard deviation in internal units
  // The following errors are not multiplied by z-score yet
  real s_err_;        // Uncertainty error in the horizontal position
  real v_err_;        // Uncertainty error in the horizontal velocity

public:
  SUMData();

  SUMData(const SUMData& sum);

  static const SUMData& EMPTY();

  real get_s_EW_std() const;

  real get_s_EW_std(const std::string& u) const;

  real get_s_NS_std() const;

  real get_s_NS_std(const std::string& u) const;

  real get_s_EN_std() const;

  real get_s_EN_std(const std::string& u) const;

  real get_sz_std() const;

  real get_sz_std(const std::string& u) const;

  real get_v_EW_std() const;

  real get_v_EW_std(const std::string& u) const;

  real get_v_NS_std() const;

  real get_v_NS_std(const std::string& u) const;

  real get_v_EN_std() const;

  real get_v_EN_std(const std::string& u) const;

  real get_vz_std() const;

  real get_vz_std(const std::string& u) const;

  void set(const SUMData& sum);

  static real eigen_value_bound(real var1, real var2, real cov);

  /**
   * In PVS: covariance@h_pos_uncertainty and covariance@h_vel_uncertainty, but in Java they are not multiplied by the z-score yet.
   */
  static real horizontal_uncertainty(real x_std, real y_std, real xy_std);

  /**
   * s_EW_std: East/West position standard deviation in internal units
   * s_NS_std: North/South position standard deviation in internal units
   * s_EN_std: East/North position standard deviation in internal units
   */
  void setHorizontalPositionUncertainty(real s_EW_std, real s_NS_std, real s_EN_std);

  /**
   * sz_std : Vertical position standard deviation in internal units
   */
  void setVerticalPositionUncertainty(real sz_std);

  /**
   * v_EW_std: East/West velocity standard deviation in internal units
   * v_NS_std: North/South velocity standard deviation in internal units
   * v_EN_std: East/North velocity standard deviation in internal units
   */
  void setHorizontalVelocityUncertainty(real v_EW_std, real v_NS_std,  real v_EN_std);

  /**
   * vz_std : Vertical speed standard deviation in internal units
   */
  void setVerticalSpeedUncertainty(real vz_std);

  /**
   * Set all uncertainties to 0
   */
  void resetUncertainty();

  /**
   * @return Horizontal position error
   */
  real getHorizontalPositionError() const;

  /**
   * @return Vertical position error
   */
  real getVerticalPositionError() const;

  /**
     Set Horizontal speed error
   */
  real getHorizontalSpeedError() const;

  /**
   * @return Vertical speed error
   */
  real getVerticalSpeedError() const;

  /**
   * Check if aircraft is using sensor uncertainty mitigation
   */
  bool is_SUM() const;

};
}
#endif
