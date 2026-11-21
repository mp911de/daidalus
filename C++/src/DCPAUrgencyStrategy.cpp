/*
 * Copyright (c) 2015-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */
/*
 * DCPAUrgencyStrategy.cpp
 *
 * Most urgent strategy based on distance at closest point of approach. When this distance is less than the minimum
 * recovery separation given by D and H, time to closest point of approach is used.
 *
 */

#include "Real.h"
#include "ACCoRDConfig.h"
#include "CD3D.h"
#include "DCPAUrgencyStrategy.h"

namespace larcfm {

/**
 * @return most urgent traffic aircraft for given ownship, traffic and lookahead time T.
 * Return -1 if no aircraft is most urgent.
 */
int DCPAUrgencyStrategy::mostUrgentAircraft(const TrafficState& ownship, const std::vector<TrafficState>& traffic, real T) const {
  (void)T; //bypass unused parameter warning (needed for interface)
  int repac = -1;
  if (!ownship.isValid() || traffic.empty()) {
    return repac;
  }
  real mindcpa = 0;
  real mintcpa = 0;
  real D = ACCoRDConfig::NMAC_D;
  real H = ACCoRDConfig::NMAC_H;
  Vect3 so = ownship.get_s();
  Vect3 vo = ownship.get_v();
  for (TrafficState::nat ac = 0; ac < traffic.size(); ++ac) {
    Vect3 si = traffic[ac].get_s();
    Vect3 vi = traffic[ac].get_v();
    Vect3 s = so.Sub(si);
    Vect3 v = vo.Sub(vi);
    real tcpa = CD3D::tccpa(s,vo,vi,D,H);
    real dcpa = v.ScalAdd(tcpa,s).cyl_norm(D,H);
    // If aircraft have almost same tcpa, select the one with smallest dcpa
    // Otherwise,  select aircraft with smallest tcpa
    bool tcpa_strategy = Util::almost_equals(tcpa,mintcpa,PRECISION5) ? dcpa < mindcpa : tcpa < mintcpa;
    // If aircraft have almost same dcpa, select the one with smallest tcpa
    // Otherwise,  select aircraft with smallest dcpa
    bool dcpa_strategy = Util::almost_equals(dcpa,mindcpa,PRECISION5) ? tcpa < mintcpa : dcpa < mindcpa;
    // If aircraft are both in a min recovery trajectory, follows tcpa strategy. Otherwise follows dcpa strategy
    if (repac < 0 || // There are no candidates
        (dcpa <= 1 ? mindcpa > 1 || tcpa_strategy : dcpa_strategy)) {
      repac = ac;
      mindcpa = dcpa;
      mintcpa = tcpa;
    }
  }
  return repac;
}

UrgencyStrategy* DCPAUrgencyStrategy::copy() const {
  return new DCPAUrgencyStrategy();
}

}
