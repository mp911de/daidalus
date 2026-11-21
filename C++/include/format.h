/*
 * Copyright (c) 2011-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */

#ifndef FORMAT_H_
#define FORMAT_H_
#include "Real.h"

#include <string>
#include <vector>
#include <map>
#include "Vect2.h"
#include "Vect3.h"
#include "Triple.h"

namespace larcfm {
  /** Format an int with 0 decimal places */
  std::string Fmi(int v);
  /** Format an unsigned in with 0 decimal places */
  std::string Fmui(unsigned int v);
  /** Format an unsigned long with 0 decimal places */
  std::string Fmul(unsigned long v);
  /** Format a double with 0 decimal places */
  std::string Fm0(real v);
  /** Format a double with 1 decimal place */
  std::string Fm1(real v);
  /** Format a double with 2 decimal place */
  std::string Fm2(real v);
  /** Format a double with 2 decimal place, filled with zeros */
  std::string Fm2z(real v);
  /** Format a double with 3 decimal place */
  std::string Fm3(real v);
  /** Format a double with 3 decimal place, filled with zeros */
  std::string Fm3z(real v);
  /** Format a double with 4 decimal place */
  std::string Fm4(real v);
  /** Format a double with 6 decimal place */
  std::string Fm6(real v);
  /** Format a double with 8 decimal place */
  std::string Fm8(real v);
  /** Format a double with 12 decimal place */
  std::string Fm12(real v);
  /** Format a double with 16 decimal place */
  std::string Fm16(real v);

  /** Format an integer to have at least minLength digits, prepending zeros as necessary */
  std::string FmLead(int i, int minLength);

  std::string FmPrecision(real v);
  std::string FmPrecision(real v, int precision);
  std::string FmPrecision(real v, int precision, bool includeTrailingZeros);

  std::string Fmb(bool b);

  /** Format a vector */
  std::string FmVec(const Vect2& v);
  /** Format a vector */
  std::string FmVec(const Vect3& v);

    /** Return a string 'n' long with 's' left-justified */
  std::string padLeft(const std::string& s, int n);
  /** Return a string 'n' long with 's' right-justified */
  std::string padRight(const std::string& s, int n);

  /** Returns true/false string */
  std::string bool2str(bool b);

  /** send string to the console with a "carriage return" */
  void fpln(const std::string& str);
  /** send a string to the console without a "carriage return" */
  void fp(const std::string& str);
  /** send string to the console with a "carriage return" */
  void fpln(std::ostream* os, const std::string& str);
  /** send a string to the console without a "carriage return" */
  void fp(std::ostream* os, const std::string& str);
  /** send string to the error console with a "carriage return" */
  void fdln(const std::string& str);

  std::string FmPair(const std::pair<int,int>& p);
  std::string FmPair(const std::pair<real,real>& p);
  std::string FmTriple(const Triple<int,int,int>& p);
  std::string FmTriple(const Triple<real,real,real>& p);

  std::string Fobj(const std::vector<int>& v);
  std::string Fobj(const std::vector<real>& v);
  std::string Fobj(const std::vector<std::string>& v);
  std::string Fobj(const std::vector<std::pair<int,int> >& v);
  std::string Fobj(const std::vector<std::pair<real,real> >& v);
  std::string Fobj(const std::vector<Triple<int,int,int> >& v);
  std::string Fobj(const std::vector<Triple<real, real, real> >& v);

  std::string Fobj(const std::map<std::string,std::string>& v);
  std::string Fobj(const std::map<int,int>& v);
  std::string Fobj(const std::map<int,std::string>& v);
  std::string Fobj(const std::map<std::string,int>& v);
  std::string Fobj(const std::map<real,real>& v);
  std::string Fobj(const std::map<real,std::string>& v);
  std::string Fobj(const std::map<std::string,real>& v);

  std::string list2str(const std::vector<std::string>& l, const std::string& delimiter);

  std::string Farray(int const v[], int sz);
  std::string Farray(real const v[], int sz);
  std::string Farray(std::string const v[], int sz);

  std::string double2PVS(real val);
  std::string double2PVS(real val, int precision);

}
#endif /* FORMAT_H_ */
