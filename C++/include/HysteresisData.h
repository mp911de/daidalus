/*
 * Implementation of hysteresis logic that includes MofN and persistence.
 * Contact: Cesar A. Munoz
 * Organization: NASA/Langley Research Center
 *
 * Copyright (c) 2011-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */

#ifndef HYSTERESISDATA_H_
#define HYSTERESISDATA_H_
#include "Real.h"

#include "MofN.h"
#include <deque>

namespace larcfm {

class HysteresisData {
public:
  /*
   * Creates an empty object
   */
  HysteresisData();

  void init();

  bool isValid() const;

  /*
   * Set given values of hysteresis, persistence time, and M of N parameters
   */
  void setHysteresisData(real hysteresis_time, real persistence_time, int m, int n);

  void outdateIfCurrentTime(real current_time);

  bool isUpdatedAtCurrentTime(real current_time) const;

  real getInitTime() const;

  real getLastTime() const;

  int getLastValue() const;

  std::string toString() const;

  /*
   * In addition of m_of_n, this function applies persistence logic
   */
  int applyHysteresisLogic(int current_value, real current_time);

  virtual ~HysteresisData() {}

private:
  MofN   mofn_;
  real hysteresis_time_;
  real persistence_time_;
  real init_time_;
  real last_time_;
  int    last_value_;
  /* When this flag is true, setting a value at current time resets hysteresis values */
  bool outdated_;

  /*
   * Reset object with given M of N value
   */
  void reset(int val);

};

} /* namespace larcfm */

#endif /* HYSTERESISDATA_H_ */
