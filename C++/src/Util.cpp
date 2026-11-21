/*
 * Util.cpp
 *
 * Contact: Cesar Munoz (cesar.a.munoz@nasa.gov), Jeff Maddalon (j.m.maddalon@nasa.gov)
 *
 * Utility functions.
 *
 * Copyright (c) 2011-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */

#include "Real.h"
#include "Util.h"
#include "Units.h"
#include "Constants.h"
#include <cmath>
#include <limits>
#include <string.h>
#include "format.h"
#include "string_util.h"
#include <stdexcept>
#include <type_traits>
#if defined(_MSC_VER)
#include <regex>
#else
#include <regex.h>
#endif
//#include <dirent.h>



namespace larcfm {

using std::string;
using std::vector;

static const real twopi = 2*Pi;


bool Util::almost_less(real a, real b) {
	if (almost_equals(a, b)) {
		return false;
	}
	return a < b;
}

bool Util::almost_less(real a, real b, INT64FM maxUlps) {
	if (almost_equals(a, b, maxUlps)) {
		return false;
	}
	return a < b;
}

bool Util::almost_greater(real a, real b) {
	if (almost_equals(a, b)) {
		return false;
	}
	return a > b;
}

bool Util::almost_greater(real a, real b, INT64FM maxUlps) {
	if (almost_equals(a, b, maxUlps)) {
		return false;
	}
	return a > b;
}


bool Util::almost_geq(real a, real b, INT64FM maxUlps) {
	return a >= b || almost_equals(a, b, maxUlps);
}

bool Util::almost_geq(real a, real b) {
	return almost_geq(a, b, PRECISION_DEFAULT);
}

bool Util::almost_leq(real a, real b, INT64FM maxUlps) {
	return a <= b || almost_equals(a, b, maxUlps);
}

bool Util::almost_leq(real a, real b) {
	return almost_leq(a, b, PRECISION_DEFAULT);
}

// This assumes that a INT64FM has the bits of a real, see Real.h
bool Util::almost_equals(real a, real b, INT64FM maxUlps) {
	if (a == b) { // if two numbers are equal, then the are almost_equal
		return true;
	}

	// special case of comparing to zero.
	if (a == real(0.0) || b == real(0.0)) {
		real comp = REAL_ZERO13;  // should correspond to PRECISION_DEFAULT;
		if (maxUlps == PRECISION5) comp = REAL_ZERO5;
		if (maxUlps == PRECISION7) comp = REAL_ZERO7;
		if (maxUlps == PRECISION9) comp = REAL_ZERO9;
		if (maxUlps == PRECISION13) comp = REAL_ZERO13;
		if (std::abs(a) < comp && std::abs(b) < comp) {
			return true;
		}
	}

	//if (ISNAN(a) || ISNAN(b)) { // this operation is sloooooow
	//	return false;
	//  }
	if (!(a < b || b < a)) { // idiom to filter out NaN's
		return false;
	}
	if (ISINF(a) || ISINF(b)) {
		return false;
	}

	//long aInt = Double.doubleToLongBits(a);
	INT64FM aInt;
	memcpy(&aInt,&a,sizeof(a)); // same as:  INT64FM aInt = *(INT64FM*)&a;

	// Make aInt lexicographically ordered as a twos-complement long
	if (aInt < 0)
		aInt = HIGH_BIT - aInt;

	// Make bInt lexicographically ordered as a twos-complement long
	//long bInt = Double.doubleToLongBits(b);
	INT64FM bInt;
	memcpy(&bInt,&b,sizeof(b)); // same as: INT64FM bInt = *(INT64FM*)&b;

	if (bInt < 0)
		bInt = HIGH_BIT - bInt;

	// Opposite signs can put the distance outside the signed integer range.
	typedef std::make_unsigned<INT64FM>::type bits;
	const bits distance = aInt >= bInt ? bits(aInt) - bits(bInt) : bits(bInt) - bits(aInt);
	return maxUlps >= 0 && distance <= bits(maxUlps);
}


// This assumes that a INT64FM has the bits of a real, see Real.h
bool Util::almost_equals(real a, real b) {
	return almost_equals(a,b,PRECISION_DEFAULT);
}


bool Util::within_epsilon(real a, real b, real epsilon) {
	return std::abs(a-b) < epsilon;
}

bool Util::within_epsilon(real a, real epsilon) {
	return std::abs(a) < epsilon;
}

INT64FM Util::llabs(const INT64FM x) {
	return x >= 0 ? x : -x;
}



real  Util::discretizeDir(real voz, real nvoz, real discreteUnits) {
	int sgn = -1;
	if (nvoz >= 0) sgn = 1;
	//double rtn =  sgn*ceil(std::abs(nvoz/discreteUnits))*discreteUnits;
	real rtn;
	if (sgn*nvoz >= sgn*voz)
		rtn =  sgn*std::ceil(std::abs(nvoz/discreteUnits))*discreteUnits;
	else
		rtn =  sgn*std::floor(std::abs(nvoz/discreteUnits))*discreteUnits;
	//f.pln(" #### vsDiscretize: voz = "+Units.str("fpm",voz)+" nvoz = "+Units.str("fpm",nvoz)+" rtn = "+Units.str("fpm",rtn));
	return rtn;
}



real Util::sq(const real x) {
	return x*x;
}

real Util::sqrt_safe(const real x) {
	return std::sqrt(Util::max(x,real(0.0)));
}

///**
// * Fast inverse square root approximation, as presented in https://en.wikipedia.org/wiki/Fast_inverse_square_root
// * (Initial testing did not show noteworthy performance improvements for daidalus.)
// */
//double Util::Q_rsqrt(const double number) {
//	double x2 = number * 0.5;
//	double y = number;
//	std::int64_t i  = * (std::int64_t *) &y;	// ensure right number of bits
//	i = 0x5fe6eb50c7b537a9 - ( i >> 1 );		// magic number to approximate sqrt
//	y = * ( double * ) &i;
//	y = y * ( 1.5 - ( x2 * y * y ) );   // 1st iteration of Newton's Method to refine number
////	y = y * ( 1.5 - ( x2 * y * y ) );   // 2nd iteration, this can be removed
//
//	return y;
//}
//
//double invsqrtQuake( double number )
//  {
//      double y = number;
//      double x2 = y * 0.5;
//      std::int64_t i = *(std::int64_t *) &y;
//      // The magic number is for doubles is from https://cs.uwaterloo.ca/~m32rober/rsqrt.pdf
//      i = 0x5fe6eb50c7b537a9 - (i >> 1);
//      y = *(double *) &i;
//      y = y * (1.5 - (x2 * y * y));   // 1st iteration
//      //      y  = y * ( 1.5 - ( x2 * y * y ) );   // 2nd iteration, this can be removed
//      return y;
//  }

real Util::atan2_safe(const real y, const real x) {
	if (y == 0 && x == 0)
		return 0;
	return std::atan2(y,x);
}

real Util::asin_safe(real x) {
	return std::asin(Util::max(-real(1.0),Util::min(x,real(1.0))));
}

real Util::acos_safe(real x) {
	return std::acos(Util::max(-real(1.0),Util::min(x,real(1.0))));
}

real Util::discr(const real a, const real b, const real c){
	return sq(b) - 4*a*c;
}

real Util::root(const real a, const real b, const real c, const int eps) {
	if (a == real(0.0) && b == real(0.0))
		return NaN;
	else if (a == real(0.0))
		return -c/b;
	else {
		real sqb = sq(b);
		real ac  = 4*a*c;
		if (almost_equals(sqb,ac) || sqb > ac)
			return (-b + eps*sqrt_safe(sqb-ac))/(2*a);
		return NaN;
	}
}

real Util::root2b(const real a, const real b, const real c, const int eps) {
	if (a == real(0.0) && b == real(0.0))
		return NaN;
	else if (a == real(0.0))
		return -c/(2*b);
	else {
		real sqb = sq(b);
		real ac  = a*c;
		if (almost_equals(sqb,ac) || sqb > ac)
			return (-b + eps*sqrt_safe(sqb-ac))/a;
		return NaN;
	}
}

real Util::rootNegC(const real a, const real b, const real c) {
	if (a == 0) return -c / b;
	real sqb = sq(b);
	real ac = 4 * a * c;
	if (almost_equals(sqb, ac) || sqb > ac)
		return (-b + sqrt_safe(sqb - ac)) / (2 * a);
	return -1;
}

int Util::sign(const bool b) {
	if (b)
		return 1;
	return -1;
}

int Util::sign(const int x) {
	if (x >= 0)
		return 1;
	return -1;
}

int Util::sign(const real x) {
	if (x >= real(0.0))
		return 1;
	return -1;
}

int Util::signTriple(const real x) {
	if (x > real(0.0)) {
		return 1;
	} else if (x < real(0.0)) {
		return -1;
	} else {
		return 0;
	}
}


const double& Util::min(const double& x, const double& y) {
	if (std::isnan(y)) return y;
	if (x == 0.0 && y == 0.0 && std::signbit(y)) return y; // if y is -0, return y
	return std::min(x,y);
}

const float& Util::min(const float& x, const float& y) {
	if (std::isnan(y)) return y;
	if (x == 0.0f && y == 0.0f && std::signbit(y)) return y; // if y is -0, return y
	return std::min(x,y);
}

const int& Util::min(const int& x, const int& y) {
	return (x < y)? x : y;
}

const long& Util::min(const long& x, const long& y) {
	return (x < y)? x : y;
}

const double& Util::max(const double& x, const double& y) {
	if (std::isnan(y)) return y;
	if (x == 0.0 && y == 0.0 && std::signbit(x)) return y; // if x is -0, return y
	return std::max(x,y);
}

const float& Util::max(const float& x, const float& y) {
	if (std::isnan(y)) return y;
	if (x == 0.0f && y == 0.0f && std::signbit(x)) return y; // if x is -0, return y
	return std::max(x,y);
}

const int& Util::max(const int& x, const int& y) {
	return (x > y)? x : y;
}

const long& Util::max(const long& x, const long& y) {
	return (x > y)? x : y;
}

/**
 * Computes the modulo of val and mod. The returned value is in the range [0,mod)
 *
 * @param val numerator
 * @param mod denominator
 * @return modulo value
 */
real Util::modulo(real val, real mod) {
	real n = std::floor(val / mod);
	real r = val - n * mod;
	return Util::almost_equals(r,mod) ? real(0.0) : r;
}

/**
 * Computes the modulo of val and mod. If mod > 0, the returned value is in the r
 * ange [0,mod). Otherwise, returns val.
 *
 * @param val numerator
 * @param mod denominator
 * @return modulo value
 */
real Util::safe_modulo(real val, real mod) {
	return mod > 0 ? modulo(val,mod) : val;
}


// To range [0,2*pi)
real Util::to_2pi(real rad) {
	return modulo(rad,twopi);
}

// To range [0,360)
real Util::to_360(real deg) {
	return modulo(deg,360);
}

// To range (-pi,pi]
real Util::to_pi(real rad) {
	real r = to_2pi(rad);
	if (r > Pi)
		return r-twopi;
	else {
		return r;
	}
}

// to range [-pi/2, +pi/2]; continuous, thus to_pi2_cont(PI/2+eps) = PI/2-eps
real Util::to_pi2_cont(real rad) {
	real r = to_pi(rad);
	if (r < -Pi / 2) {
		return -Pi - r;
	} else if (r < Pi / 2) {
		return r;
	} else {
		return Pi - r;
	}
}


// To range (-180,180]
real Util::to_180(real deg) {
	real d = to_360(deg);
	if (d > real(180.0))
		return d-real(360.0);
	else {
		return d;
	}
}


bool Util::less_or_equal(const std::string& s1, const std::string& s2) {
	if (s1.compare(s2) >= 0) return true;
	return false;
}


bool Util::clockwise(real alpha, real beta) {
	real a = to_2pi(alpha);
	real b = to_2pi(beta);
	if (std::abs(a-b) <= real(M_PI)) {
		return b >= a;
	}
	return a > b;
}

int Util::turnDir(real initTrack, real goalTrack) {
	if (Util::clockwise(initTrack,goalTrack)) return 1;
	else return -1;
}


real Util::turnDelta(real alpha, real beta) {
	real a = to_2pi(alpha);
	real b = to_2pi(beta);
	real delta = std::abs(a-b);
	if (delta <= real(M_PI)) {
		return delta;
	}
	return real(2.0)*real(M_PI) - delta;
}

real Util::signedTurnDelta(real alpha, real beta) {
	return turnDir(alpha,beta)*turnDelta(alpha,beta);
}


real Util::turnDelta(real alpha, real beta, bool turnRight) {
	if (Constants::almost_equals_radian(alpha,beta)) return 0;                 // do not want 2 PI returned
	//if (std::abs(alpha-beta) < 0.00001) return 0;                              // do not want 2 PI returned
	bool clkWise = Util::clockwise(alpha,beta);
	real rtn = Util::turnDelta(alpha,beta);
	if (turnRight != clkWise)  { // go the long way around
		rtn = real(2.0) * real(M_PI) - rtn;
	}
	return rtn;
}

real Util::turnDelta(real alpha, real beta, int dir) {
	return turnDelta(alpha, beta, dir > 0);
}


#if defined(_MSC_VER)
bool Util::is_double(const string& str) {
	std::string sb(str);
	trim(sb," \t");
	std::regex numre("^-?[0-9]*(\\.[0-9]*)?$");
	return std::regex_match(sb, numre);
}

#else
bool Util::is_double(const string& str) {

	string sb(str);
	regex_t regex;
	int reti;

	trim(sb," \t");

	reti = regcomp(&regex, "^-?[0-9]*(\\.[0-9]*)?$", REG_EXTENDED);
	if (reti != 0) {
		fdln("Could not compile regex\n");
	}

	/* Execute regular expression */
	reti = regexec(&regex, sb.c_str(), 0, NULL, 0);
	regfree(&regex);
	return !reti;
}
#endif

real Util::parse_double(const string& str) {
	std::istringstream stream;
	stream.str(str);
	real d;
	stream >> d;
	if (stream.fail()) {
		return real(0.0);
	}
	return d;
}

/**
 * Returns true if the stored value for key is likely a boolean
 * @param s name
 * @return true if string value is true/false/t/f, false otherwise
 */
bool Util::is_boolean(const std::string& s) {
	return (equalsIgnoreCase(s, "true") || equalsIgnoreCase(s, "T") || equalsIgnoreCase(s, "false") || equalsIgnoreCase(s, "F"));
}


/** @param degMinSec  Lat/Lon string of the form "46:55:00"  or "-111:57:00"
      @return   numbers of degrees decimal
 */
real Util::decimalDegrees(const std::string& degMinSec) {
	vector<std::string> dms = split(degMinSec,":");
	real degrees = parse_double(dms[0]);
	int sgn = Util::sign(degrees);
	real minutes = parse_double(dms[1]);
	real seconds = parse_double(dms[2]);
	return sgn*(std::abs(degrees) + minutes/60+ seconds/3600);
}


/** Reads in a clock string and converts it to seconds.
 * Accepts hh:mm:ss, mm:ss, and ss.
 */
real Util::parse_time(const string& s) {
	try {
		real tm = -real(1.0);
		string patternStr = "[:]";
		vector<std::string> fields2 = split(s, patternStr);
		if (fields2.size() >= 3) {
			tm = parse_double(fields2[2]) + 60 * parse_double(fields2[1]) + 3600 * parse_double(fields2[0]); // hrs:min:sec
		} else if (fields2.size() == 2) {
			tm = parse_double(fields2[1]) + 60 * parse_double(fields2[0]); // min:sec
		} else if (fields2.size() == 1){
			tm = parse_double(fields2[0]); //getColumn(_sec, head[TM_CLK]);
		}
		return tm;
	} catch (std::runtime_error const& e) {
		return -real(1.0);
	}
}


/**
 * @param t time in seconds
 * @return String of hours:mins:secs
 */
string Util::hoursMinutesSeconds(real t) {
	int hours = (int) t/3600;
	int rem = (int) t - hours*3600;
	int mins = rem / 60;
	rem = rem - mins*60;
	int secs = rem;
	return ""+FmLead(hours,1)+":"+FmLead(mins,2)+":"+FmLead(secs,2);
}

///**
// * The behavior of the x%y operator is different between Java and C++ if either x or y is negative.  Use this to always return a value between 0 and y.
// * @param x value
// * @param y range
// * @return x mod y, having the same sign as y (Java behavior)
// */
//int Util::mod(int x, int y) {
//	int r = std::abs(x) % std::abs(y);
//	return r*sign(y);
//}



//*******************************************
// deprecated functions:


INT64FM llabs(const INT64FM x) {
	return Util::llabs(x);
}

real sq(const real x) {
	return Util::sq(x);
}

real sqrt_safe(const real x) {
	return std::sqrt(Util::max(x,real(0.0)));
}

real atan2_safe(const real y, const real x) {
	return Util::atan2_safe(y,x);
}

real asin_safe(real x) {
	return Util::asin_safe(x);
}

real acos_safe(real x) {
	return Util::acos_safe(x);
}

real discr(const real a, const real b, const real c){
	return Util::discr(a,b,c);
}

real root(const real a, const real b, const real c, const int eps) {
	return Util::root(a,b,c,eps);
}

real root2b(const real a, const real b, const real c, const int eps) {
	return Util::root2b(a,b,c,eps);
}

int sign(const real x) {
	return Util::sign(x);
}


// To range [0,2*pi)
real to_2pi(real rad) {
	return Util::to_2pi(rad);
}


// To range [0,360)
real to_360(real deg) {
	return Util::to_360(deg);
}

// To range (-pi,pi]
real to_pi(real rad) {
	return Util::to_pi(rad);
}

// to range [-pi/2, +pi/2]; continuous, thus to_pi2_cont(PI/2+eps) = PI/2-eps
real to_pi2_cont(real rad) {
	return Util::to_pi2_cont(rad);
}


// To range (-180,180]
real to_180(real deg) {
	return Util::to_180(deg);
}

//#ifndef _MSC_VER
//// will (probably) not work in Visual C++
//// modified from code on net (http://www.linuxquestions.org/quecgstions/programming-9/c-list-files-in-directory-379323/)
//
//int filesInDir (string path, vector<string>& fv) {
//        DIR *directory;
//        struct dirent *dir_entry;
//        directory = opendir(path.c_str());
//        if(directory == NULL) {
//                return -1; // error
//        }
//        dir_entry = readdir(directory);
//        while (dir_entry != NULL) {
//                string s = string(dir_entry->d_name);
//                int x = s.find(".txt", 0); //doesn't work in the if test, for some reason
//                if (x >= 0) {
//                        fv.push_back(path+"/"+s);
//                }
//                dir_entry = readdir(directory);
//        }
//        closedir(directory);
//        return 1;  // ok
//}
//#endif
//
//bool fileExists(const char *filename) {
//        std::ifstream tmp(filename);
//        return tmp.good();
//}


}




