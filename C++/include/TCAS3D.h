/*
 * Copyright (c) 2012-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */
#ifndef TCAS3D_H_
#define TCAS3D_H_
#include "Real.h"

#include "Detection3D.h"
#include "TCAS2D.h"
#include "Vect2.h"
#include "TCASTable.h"
#include "ConflictData.h"
#include "Triple.h"

namespace larcfm {
class TCAS3D : public Detection3D {

private:
  TCASTable table_;
  std::string id;

  Triple<real,real,real> RA3D_interval(int sl, const Vect2& so2, real soz, const Vect2& vo2, real voz,
      const Vect2& si2, real siz, const Vect2& vi2, real viz, real B, real T, bool severity) const;
  LossData RA3D_times(const Vect3& so, const Vect3& vo, const Vect3& si, const Vect3& vi, real B, real T, real* critical_time) const;
protected:
  virtual LossData conflictInterval(const Vect3& so, const Vect3& vo, const Vect3& si, const Vect3& vi, real B, real T) const;
  virtual LossData conflictIntervalWithTrafficState(const TrafficState& ownship, const TrafficState& intruder, real B, real T) const;
public:

  /** Constructor that uses the default TCAS tables. */
  TCAS3D();

  /** Constructor that specifies a particular instance of the TCAS tables. */
  explicit TCAS3D(const TCASTable& tables);

  /**
   * @return one static TCAS3D
   */
  static const TCAS3D& A_TCAS3D();

  /**
   * @return one static TCASII_RA
   */
  static const TCAS3D& TCASII_RA();

  /**
   * @return one static TCASII_TA
   */
  static const TCAS3D& TCASII_TA();

  /** Make TCAS3D object with empty Table **/
  static TCAS3D make_Empty();

  /** Make TCAS3D object with an RA Table **/
  static TCAS3D make_TCASII_RA();

  /** Make TCAS3D objec with a TA Table **/
  static TCAS3D make_TCASII_TA();

  /** This returns a copy of the internal TCAS table */
  TCASTable& getTCASTable();

  /**
   * Set table to TCASII Thresholds (RA Table when ra is true, TA Table when ra is false)
   */
  void setDefaultTCASIIThresholds(bool ra);


  /** This sets the internal table to be a deep copy of the supplied one.  Any previous links will be broken. */
  void setTCASTable(const TCASTable& tables);

  // if true, then ownship has a TCAS resolution advisory at current time
  bool TCASII_RA(const Vect3& so, const Vect3& vo, const Vect3& si, const Vect3& vi) const;

  Triple<real,real,real> RA3D_interval(int sl, const Vect2& so2, real soz, const Vect2& vo2, real voz,
      const Vect2& si2, real siz, const Vect2& vi2, real viz, real B, real T) const;

  // if true, within lookahead time T, the ownship has a TCAS resolution advisory (at time time before T) (effectively conflict detection)
  ConflictData RA3D(const Vect3& so, const Vect3& vo,  const Vect3& si, const Vect3& vi, real B, real T) const;

  // The methods violation and conflict are inherited from Detection3DSum. This enable a uniform
  // treatment of border cases in the generic bands algorithms

  ConflictData conflictDetection(const Vect3& so, const Vect3& vo, const Vect3& si, const Vect3& vi, real B, real T) const;

  TCAS3D* copy() const;
  TCAS3D* make() const;

  static bool vertical_RA(real sz, real vz, real ZTHR, real TCOA);


  //  static double time_coalt(double sz, double voz, double viz, double h);


  static bool cd2d_TCAS_after(real HMD, Vect2 s, Vect2 vo, Vect2 vi, real t);

  static bool cd2d_TCAS(real HMD, Vect2 s, Vect2 vo, Vect2 vi);


  virtual std::string getSimpleClassName() const;


  /**
   * Returns TAU threshold for sensitivity level sl in seconds
   */
  real getTAU(int sl) const;

  /**
   * Returns TCOA threshold for sensitivity level sl in seconds
   */
  real getTCOA(int sl) const;

  /**
   * Returns DMOD for sensitivity level sl in internal units.
   */
  real getDMOD(int sl) const;

  /**
   * Returns DMOD for sensitivity level sl in u units.
   */
  real getDMOD(int sl, const std::string& u) const;

  /**
   * Returns Z threshold for sensitivity level sl in internal units.
   */
  real getZTHR(int sl) const;

  /**
   * Returns Z threshold for sensitivity level sl in u units.
   */
  real getZTHR(int sl,const std::string& u) const;

  /**
   * Returns HMD for sensitivity level sl in internal units.
   */
  real getHMD(int sl) const;

  /**
   * Returns HMD for sensitivity level sl in u units.
   */
  real getHMD(int sl, const std::string& u) const;

  /** Modify the value of Tau Threshold for a given sensitivity level (2-8)
   * Parameter val is given in seconds
   */
  void setTAU(int sl, real val);

  /** Modify the value of TCOA Threshold for a given sensitivity level (2-8)
   * Parameter val is given in seconds
   */
  void setTCOA(int sl, real val);

  /** Modify the value of DMOD for a given sensitivity level (2-8)
   * Parameter val is given in internal units
   */
  void setDMOD(int sl, real val);

  /** Modify the value of DMOD for a given sensitivity level (2-8)
   * Parameter val is given in u units
   */
  void setDMOD(int sl, real val, const std::string& u);

  /** Modify the value of ZTHR for a given sensitivity level (2-8)
   * Parameter val is given in internal units
   */
  void setZTHR(int sl, real val);

  /** Modify the value of ZTHR for a given sensitivity level (2-8)
   * Parameter val is given in u units
   */
  void setZTHR(int sl, real val, const std::string& u);

  /**
   * Modify the value of HMD for a given sensitivity level (2-8)
   * Parameter val is given in internal units
   */
  void setHMD(int sl, real val);

  /**
   * Modify the value of HMD for a given sensitivity level (2-8)
   * Parameter val is given in u units
   */
  void setHMD(int sl, real val, const std::string& u);

  void setHMDFilter(bool flag);

  bool getHMDFilter();

  virtual std::string toString() const;
  virtual std::string toPVS() const;

  virtual std::string getIdentifier() const;
  virtual void setIdentifier(const std::string& id);

  virtual bool equals(const Detection3D& d) const;

  ParameterData getParameters() const;
  void updateParameterData(ParameterData& p) const;
  void setParameters(const ParameterData& p);

  bool contains(const Detection3D& cd) const;

  virtual void horizontalHazardZone(std::vector<Position>& haz, const TrafficState& ownship, const TrafficState& intruder,
      real T) const;

};

}
#endif
