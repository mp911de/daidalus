/*
 * Copyright (c) 2015-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */

#include "Real.h"
#include "SUMData.h"
#include "Units.h"
#include "Util.h"
#include <string>

namespace larcfm {


SUMData::SUMData() {
    s_err_ = real(0.0);
    v_err_ = real(0.0);
    s_EW_std_ = real(0.0);
    s_NS_std_ = real(0.0);
    s_EN_std_ = real(0.0);
    sz_std_ = real(0.0);
    v_EW_std_ = real(0.0);
    v_NS_std_ = real(0.0);
    v_EN_std_ = real(0.0);
    vz_std_ = real(0.0);
}

SUMData::SUMData(const SUMData& sum) {
    set(sum);
}

const SUMData& SUMData::EMPTY() {
  static SUMData tmp;
  return tmp;
}

real SUMData::get_s_EW_std() const {
    return s_EW_std_;
}

real SUMData::get_s_EW_std(const std::string& u) const {
    return Units::to(u,s_EW_std_);
}

real SUMData::get_s_NS_std() const {
    return s_NS_std_;
}

real SUMData::get_s_NS_std(const std::string& u) const {
    return Units::to(u,s_NS_std_);
}

real SUMData::get_s_EN_std() const {
    return s_EN_std_;
}

real SUMData::get_s_EN_std(const std::string& u) const {
    return Units::to(u,s_EN_std_);
}

real SUMData::get_sz_std() const {
    return sz_std_;
}

real SUMData::get_sz_std(const std::string& u) const {
    return Units::to(u,sz_std_);
}

real SUMData::get_v_EW_std() const {
    return v_EW_std_;
}

real SUMData::get_v_EW_std(const std::string& u) const {
    return Units::to(u,v_EW_std_);
}

real SUMData::get_v_NS_std() const {
    return v_NS_std_;
}

real SUMData::get_v_NS_std(const std::string& u) const {
    return Units::to(u,v_NS_std_);
}

real SUMData::get_v_EN_std() const {
    return v_EN_std_;
}

real SUMData::get_v_EN_std(const std::string& u) const {
    return Units::to(u,v_EN_std_);
}

real SUMData::get_vz_std() const {
    return vz_std_;
}

real SUMData::get_vz_std(const std::string& u) const {
    return Units::to(u,vz_std_);
}

void SUMData::set(const SUMData& sum) {
    s_err_ = sum.s_err_;
    v_err_ = sum.v_err_;
    s_EW_std_ = sum.s_EW_std_;
    s_NS_std_ = sum.s_NS_std_;
    s_EN_std_ = sum.s_EN_std_;
    sz_std_ = sum.sz_std_;
    v_EW_std_ = sum.v_EW_std_;
    v_NS_std_ = sum.v_NS_std_;
    v_EN_std_ = sum.v_EN_std_;
    vz_std_ = sum.vz_std_;
}

real SUMData::eigen_value_bound(real var1, real var2, real cov) {
    real varAve = (var1 + var2)/real(2.0);
    real det = var1*var2 - Util::sq(cov);
    real temp1 = Util::sqrt_safe(varAve * varAve - det);
    return varAve + temp1;
}

/**
 * In PVS: covariance@h_pos_uncertainty and covariance@h_vel_uncertainty, but in Java they are not multiplied by the z-score yet.
 */
real SUMData::horizontal_uncertainty(real x_std, real y_std, real xy_std) {
    return Util::sqrt_safe(eigen_value_bound(Util::sq(x_std),Util::sq(y_std),Util::sign(xy_std)*Util::sq(xy_std)));
}

/**
 * s_EW_std: East/West position standard deviation in internal units
 * s_NS_std: North/South position standard deviation in internal units
 * s_EN_std: East/North position standard deviation in internal units
 */
void SUMData::setHorizontalPositionUncertainty(real s_EW_std, real s_NS_std, real s_EN_std) {
    s_EW_std_ = s_EW_std;
    s_NS_std_ = s_NS_std;
    s_EN_std_ = s_EN_std;
    s_err_ = horizontal_uncertainty(s_EW_std,s_NS_std,s_EN_std);
}

/**
 * sz_std : Vertical position standard deviation in internal units
 */
void SUMData::setVerticalPositionUncertainty(real sz_std) {
    sz_std_ = sz_std;
}

/**
 * v_EW_std: East/West velocity standard deviation in internal units
 * v_NS_std: North/South velocity standard deviation in internal units
 * v_EN_std: East/North velocity standard deviation in internal units
 */
void SUMData::setHorizontalVelocityUncertainty(real v_EW_std, real v_NS_std,  real v_EN_std) {
    v_EW_std_ = v_EW_std;
    v_NS_std_ = v_NS_std;
    v_EN_std_ = v_EN_std;
    v_err_ = horizontal_uncertainty(v_EW_std,v_NS_std,v_EN_std);
}

/**
 * vz_std : Vertical speed standard deviation in internal units
 */
void SUMData::setVerticalSpeedUncertainty(real vz_std) {
    vz_std_ = vz_std;
}

/**
 * Set all uncertainties to 0
 */
void SUMData::resetUncertainty() {
    s_err_ = real(0.0);
    v_err_ = real(0.0);
    s_EW_std_ = real(0.0);
    s_NS_std_ = real(0.0);
    s_EN_std_ = real(0.0);
    sz_std_ = real(0.0);
    v_EW_std_ = real(0.0);
    v_NS_std_ = real(0.0);
    v_EN_std_ = real(0.0);
    vz_std_ = real(0.0);
}

/**
 * @return Horizontal position rrror
 */
real SUMData::getHorizontalPositionError() const {
    return s_err_;
}

/**
 * @return vertical position error
 */
real SUMData::getVerticalPositionError() const {
    return sz_std_;
}

/**
 Set Horizontal speed error
 */
real SUMData::getHorizontalSpeedError() const {
    return v_err_;
}

/**
 * @return Vertical speed error
 */
real SUMData::getVerticalSpeedError() const {
    return vz_std_;
}

/**
 * Check if aircraft is using sensor uncertainty mitigation
 */
bool SUMData::is_SUM() const {
    return s_err_ != real(0.0) || sz_std_ != real(0.0) || v_err_ != real(0.0) || vz_std_ != real(0.0);
}



}
