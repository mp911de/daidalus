/*
 * Units.cpp
 *
 * Contact: Jeff Maddalon
 * NASA LaRC
 *
 * Conversion to internal units: meters, kilogrames, seconds, radians, ...
 *
 * Copyright (c) 2011-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */

#include "Real.h"
#include "Units.h"
#include "Util.h"
#include "Constants.h"
#include "format.h"
#include "string_util.h"
#if defined(_MSC_VER)
#include <regex>
#else
#include <regex.h>
#include <memory>
#endif
#include <cstdio>
#include <sstream>
#include <string.h>



using namespace larcfm;
using namespace std;

// The functions below use the "construct-on-first-read" pattern, which avoids
// the static initialization fiasco problem.  This problem manifests itself differently
// between g++ and Visual C++, so be careful if you want to "optimize" this section.

static real _FormalATM_unspecified() {
	static real ans = real(1.0);
	return ans;
}
static real _FormalATM_unitless() {
	static real ans = real(1.0);
	return ans;
}
static real _FormalATM_internal() {
	static real ans = real(1.0);
	return ans;
}

static real _FormalATM_m() {
	static real ans = real(1.0);
	return ans;
}
static real _FormalATM_km() {
	static real ans = real(1000.0) * _FormalATM_m();
	return ans;
}
static real _FormalATM_NM() {
	static real ans = real(1852.0) * _FormalATM_m();
	return ans;
}
static real _FormalATM_ft() {
	static real ans = real(0.3048) * _FormalATM_m();
	return ans;
}
static real _FormalATM_mile() {
	static real ans = real(5280.0) * _FormalATM_ft();
	return ans;
}
static real _FormalATM_inch() {
	static real ans = _FormalATM_ft() / real(12.0);
	return ans;
}
static real _FormalATM_yard() {
	static real ans = real(3.0) * _FormalATM_ft();
	return ans;
}
static real _FormalATM_mm() {
	static real ans = real(0.001) * _FormalATM_m();
	return ans;
}

static real _FormalATM_m2() {
	static real ans = _FormalATM_m() * _FormalATM_m();
	return ans;
}
static real _FormalATM_ft2() {
	static real ans = _FormalATM_ft() * _FormalATM_ft();
	return ans;
}

static real _FormalATM_kg() {
	static real ans = real(1.0);
	return ans;
}
static real _FormalATM_pound_mass() {
	static real ans = real(0.45359237) * _FormalATM_kg();
	return ans;
}

static real _FormalATM_s() {
	static real ans = real(1.0);
	return ans;
}
static real _FormalATM_min() {
	static real ans = real(60.0) * _FormalATM_s();
	return ans;
}
static real _FormalATM_hour() {
	static real ans = real(3600.0) * _FormalATM_s();
	return ans;
}
static real _FormalATM_day() {
	static real ans = real(8640000.0) * _FormalATM_s();
	return ans;
}
static real _FormalATM_ms() {
	static real ans = real(0.001) * _FormalATM_s();
	return ans;
}

static real _FormalATM_rad() {
	static real ans = real(1.0);
	return ans;
}
static real _FormalATM_deg() {
	static real ans = real(M_PI) / real(180.0);
	return ans;
}

static real _FormalATM_radian_per_second() {
	static real ans = _FormalATM_rad() / _FormalATM_s();
	return ans;
}

static real _FormalATM_degree_per_second() {
	static real ans = _FormalATM_deg() / _FormalATM_s();
	return ans;
}

static real _FormalATM_mps() {
	static real ans = real(1.0);
	return ans;
}
static real _FormalATM_kph() {
	static real ans = _FormalATM_km() / _FormalATM_hour();
	return ans;
}
static real _FormalATM_kn() {
	static real ans = _FormalATM_NM() / _FormalATM_hour();
	return ans;
}
static real _FormalATM_fpm() {
	static real ans = _FormalATM_ft() / _FormalATM_min();
	return ans;
}
static real _FormalATM_mph() {
	static real ans = _FormalATM_mile() / _FormalATM_hour();
	return ans;
}
static real _FormalATM_foot_per_second() {
	static real ans = _FormalATM_ft() / _FormalATM_s();
	return ans;
}
static real _FormalATM_NMps() {
	static real ans = _FormalATM_NM() / _FormalATM_s();
	return ans;
}


static real _FormalATM_meter_per_second2() {
	static real ans = _FormalATM_m() / (_FormalATM_s() * _FormalATM_s());
	return ans;
}
real larcfm::_FormalATM_gn() {
	static real ans = real(9.80665) * _FormalATM_meter_per_second2();
	return ans;
}
static real _FormalATM_foot_per_second2() {
	static real ans = _FormalATM_ft() / (_FormalATM_s() * _FormalATM_s());
	return ans;
}

static real _FormalATM_newton() {
	static real ans = _FormalATM_kg() * _FormalATM_meter_per_second2();
	return ans;
}
static real _FormalATM_slug() {
	static real ans = _FormalATM_pound_mass() * _FormalATM_gn()
							/ _FormalATM_foot_per_second2();
	return ans;
}
static real _FormalATM_pound_force() {
	static real ans = _FormalATM_slug() * _FormalATM_foot_per_second2();
	return ans;
}

static real _FormalATM_joule() {
	static real ans = _FormalATM_newton() * _FormalATM_m();
	return ans;
}
static real _FormalATM_foot_pound_force() {
	static real ans = _FormalATM_ft() * _FormalATM_pound_force();
	return ans;
}


static real _FormalATM_pascal() {
	static real ans = _FormalATM_newton() / (_FormalATM_m() * _FormalATM_m());
	return ans;
}
real larcfm::_FormalATM_P0() {
	static real ans = real(101325.0) * _FormalATM_pascal();
	return ans;
}
static real _FormalATM_atm() {
	static real ans = real(1.0) * _FormalATM_P0();
	return ans;
}

static real _FormalATM_degreeC() {
  static real res = -real(1.0);
  return res;
}
static real _FormalATM_degreeF() {
  static real res = -real(2.0);
  return res;
}
static real _FormalATM_K() {
  static real res = real(1.0);
  return res;
}
static real _FormalATM_degreeR() {
  static real res = real(1.0) / real(1.8);
  return res;
}

const real Units::unspecified = Units::getFactor("unspecified");
const real Units::unitless = Units::getFactor("unitless");
const real Units::internal = Units::getFactor("internal");

const real Units::m = Units::getFactor("m");
const real Units::km = Units::getFactor("km");
const real Units::nmi = Units::getFactor("NM");
const real Units::NM = Units::getFactor("NM");
const real Units::ft = Units::getFactor("ft");
const real Units::mm = Units::getFactor("mm");

const real Units::s = Units::getFactor("s");
const real Units::sec = Units::getFactor("sec");
const real Units::min = Units::getFactor("min");
const real Units::hour = Units::getFactor("hour");
const real Units::day = Units::getFactor("day");
const real Units::ms = Units::getFactor("ms");

const real Units::rad = Units::getFactor("rad");
const real Units::deg = Units::getFactor("deg");

const real Units::kg = Units::getFactor("kg");
const real Units::pound_mass = Units::getFactor("pound_mass");

const real Units::mps = Units::getFactor("mps");
const real Units::meter_per_second = Units::getFactor("meter_per_second");
const real Units::kph = Units::getFactor("kph");
const real Units::kilometer_per_hour = Units::getFactor("kilometer_per_hour");
const real Units::knot = Units::getFactor("knot");
const real Units::kn = Units::getFactor("kn");
const real Units::kts = Units::getFactor("kts");
const real Units::fpm = Units::getFactor("fpm");
const real Units::foot_per_second = Units::getFactor("foot_per_second");
const real Units::nautical_mile_per_second = Units::getFactor("nautical_mile_per_second");
const real Units::mph = Units::getFactor("mph");

const real Units::meter_per_second2 = Units::getFactor("meter_per_second2");
const real Units::gn = _FormalATM_gn();
const real Units::G = Units::getFactor("G");

const real Units::newton = Units::getFactor("newton");
const real Units::pascal = Units::getFactor("pascal");

const real Units::P0 = _FormalATM_P0();
const real Units::atm = Units::getFactor("atm");

const real Units::degreeC = Units::getFactor("degreeC");
const real Units::degreeF = Units::getFactor("degreeF");
const real Units::K = Units::getFactor("K");
const real Units::degreeK = Units::getFactor("K");
const real Units::degreeR = Units::getFactor("degreeR");

// Converts value canonical units to [symbol] units
real Units::getFactor(const std::string& symbolp) {

	const std::string symbol = Units::canonical(symbolp);

	if (symbol == "m") {
		return _FormalATM_m();
	} else if (symbol == "ft") {
		return _FormalATM_ft();
	} else if (symbol == "yard") {
		return _FormalATM_yard();
	} else if (symbol == "in") {
		return _FormalATM_inch();
	} else if (symbol == "km") {
		return _FormalATM_km();
	} else if (symbol == "mi") {
		return _FormalATM_mile();
	} else if (symbol == "NM") {
		return _FormalATM_NM();
	} else if (symbol == "mm") {
		return _FormalATM_mm();

	} else if (symbol == "m^2") {
		return _FormalATM_m2();
	} else if (symbol == "ft^2") {
		return _FormalATM_ft2();

	} else if (symbol == "s") {
		return _FormalATM_s();
	} else if (symbol == "min") {
		return _FormalATM_min();
	} else if (symbol == "hour") {
		return _FormalATM_hour();
	} else if (symbol == "day") {
		return _FormalATM_day();
	} else if (symbol == "ms") {
		return _FormalATM_ms();

	} else if (symbol == "rad") {
		return _FormalATM_rad();
	} else if (symbol == "deg") {
		return _FormalATM_deg();
	} else if (symbol == "deg/s") {
		return _FormalATM_degree_per_second();
	} else if (symbol == "rad/s") {
		return _FormalATM_radian_per_second();

	} else if (symbol == "kg") {
		return _FormalATM_kg();
	} else if (symbol == "lbm") {
		return _FormalATM_pound_mass();
	} else if (symbol == "slug") {
		return _FormalATM_slug();

	} else if (symbol == "km/h") {
		return _FormalATM_kph();
	} else if (symbol == "kn") {
		return _FormalATM_kn();
	} else if (symbol == "fpm") {
		return _FormalATM_fpm();
	} else if (symbol == "ft/s") {
		return _FormalATM_foot_per_second();
	} else if (symbol == "mph") {
		return _FormalATM_mph();
	} else if (symbol == "NM/s") {
		return _FormalATM_NMps();
	} else if (symbol == "m/s") {
		return _FormalATM_mps();

	} else if (symbol == "m/s^2") {
		return _FormalATM_meter_per_second2();
	} else if (symbol == "G") {
		return _FormalATM_gn();
	} else if (symbol == "ft/s^2") {
		return _FormalATM_foot_per_second2();

	} else if (symbol == "N") {
		return _FormalATM_newton();
	} else if (symbol == "lbf") {
		return _FormalATM_pound_force();

	} else if (symbol == "J") {
		return _FormalATM_joule();
	} else if (symbol == "ft-lbf") {
		return _FormalATM_foot_pound_force();

	} else if (symbol == "atm") {
		return _FormalATM_atm();
	} else if (symbol == "Pa") {
		return _FormalATM_pascal();

	} else if (symbol == "degreeC") {
		return _FormalATM_degreeC();
	} else if (symbol == "degreeF") {
		return _FormalATM_degreeF();
	} else if (symbol == "K") {
		return _FormalATM_K();
	} else if (symbol == "degreeR") {
		return _FormalATM_degreeR();

	} else if (symbol == "internal") {
		return _FormalATM_internal();
	} else if (symbol == "unspecified") {
		return _FormalATM_unspecified();
	} else if (symbol == "unitless") {
		return _FormalATM_unitless();
	}
	return real(0.0);  // this is a special value that indicates an invalid unit
}

bool Units::isUnit(const std::string& unit) {
	return getFactor(unit) != real(0.0);
}

const std::string Units::canonical(const std::string& unit) {
	if (unit == "unspecified")
		return unit;
	if (unit == "unitless" || unit == "none")
		return "unitless";
	if (unit == "internal")
		return unit;
	if (unit == "m" || unit == "meter" || unit == "metre")
		return "m";
	if (unit == "ft" || unit == "feet" || unit == "foot")
		return "ft";
	if (unit == "km" || unit == "kilometer")
		return "km";
	if (unit == "NM" || unit == "nmi" || unit == "nm"  // Do not add nm, nm means nanometers, not nautical miles
			|| unit == "nautical_mile")
		return "NM";
	if (unit == "mile" || unit == "mi")
		return "mi";
	if (unit == "inch" || unit == "in")
		return "in";
	if (unit == "yard")
		return unit;
	if (unit == "mm" || unit == "millimeter")
		return "mm";

	if (unit == "m^2")
		return unit;
	if (unit == "ft^2")
		return unit;

	if (unit == "s" || unit == "sec")
		return "s";
	if (unit == "min" || unit == "minute" || unit == "minutes")
		return "min";
	if (unit == "hour" || unit == "h" || unit == "hr")
		return "hour";
	if (unit == "day")
		return "day";
	if (unit == "ms")
		return "ms";

	if (unit == "rad" || unit == "radian")
		return "rad";
	if (unit == "deg" || unit == "degree")
		return "deg";

	if (unit == "deg/s")
		return unit;
	if (unit == "rad/s" || unit == "radian_per_second")
		return "rad/s";

	if (unit == "kg" || unit == "kilogram")
		return "kg";
	if (unit == "lbm" || unit == "pound_mass")
		return "lbm";
	if (unit == "slug")
		return "slug";

	if (unit == "m/s" || unit == "mps" || unit == "meter_per_second")
		return "m/s";
	if (unit == "kph" || unit == "km/h" || unit == "kilometer_per_hour" || unit == "kilometre_per_hour")
		return "km/h";
	if (unit == "knot" || unit == "kn" || unit == "kts")
		return "kn";
	if (unit == "fpm" || unit == "ft/min" || unit == "foot/min"	|| unit == "feet/min")
		return "fpm";
	if (unit == "ft/s" || unit == "fps" || unit == "foot_per_second" || unit == "feet_per_second")
		return "ft/s";
	if (unit == "nautical_mile_per_second" || unit == "NM/s" || unit == "nmi/s")
		return "NM/s";
	if (unit == "mph")
		return "mph";

	if (unit == "m/s^2" || unit == "meter_per_second2")
		return "m/s^2";
	if (unit == "G")
		return "G";
	if (unit == "ft/s^2" || unit == "foot_per_second2")
		return "ft/s^2";

	if (unit == "newton" || unit == "N")
		return "N";
	if (unit == "lbf" || unit == "pound_force")
		return "lbf";

	if (unit == "joule" || unit == "J")
		return "J";
	if (unit == "foot_pound_force" || unit == "ft-lbf")
		return "ft-lbf";

	if (unit == "atm")
		return "atm";
	if (unit == "pascal" || unit == "Pa")
		return "Pa";

	if (unit == "degreeF") 
	  return "degreeF";
	if (unit == "degreeC")
		return "degreeC";
	if (unit == "degreeR" || unit == "Rankine")
		return "degreeR";
	if (unit == "K" || unit == "degreeK")
		return "K";

	return "UNKNOWN_UNIT_" + unit;

}

bool Units::isCompatible(const std::string& unit1p, const std::string& unit2p) {
	const std::string unit1 = Units::canonical(unit1p);
	const std::string unit2 = Units::canonical(unit2p);
	if (unit1 == "unspecified" || unit2 == "unspecified") {
		return true;
	}

	if (unit1 == "m" || unit1 == "ft" || unit1 == "km" || unit1 == "NM"
			|| unit1 == "mi" || unit1 == "in" || unit1 == "yard" || unit1 == "mm") {
		return unit2 == "m" || unit2 == "ft" || unit2 == "km" || unit2 == "NM"
				|| unit2 == "mi" || unit2 == "in" || unit2 == "yard" || unit2 == "mm";
	} else if (unit1 == "ft^2" || unit1 == "m^2") {
		return unit2 == "ft^2" || unit2 == "m^2";
	} else if (unit1 == "s" || unit1 == "sec" || unit1 == "min" || unit1 == "hour" || unit1 == "ms" || unit1 == "day") {
		return unit2 == "s" || unit2 == "sec" || unit2 == "min" || unit2 == "hour" || unit2 == "ms" || unit2 == "day";
	} else if (unit1 == "rad" || unit1 == "deg") {
		return unit2 == "rad" || unit2 == "deg";
	} else if (unit1 == "rad/s" || unit1 == "deg/s") {
		return unit2 == "rad/s" || unit2 == "deg/s";
	} else if (unit1 == "kg" || unit1 == "slug" || unit1 == "lbm") {
		return unit2 == "kg" || unit2 == "slug" || unit2 == "lbm";
	} else if (unit1 == "kph" || unit1 == "kn" || unit1 == "fpm"
			|| unit1 == "m/s" || unit1 == "mph" || unit1 == "ft/s" || unit1 == "NM/s") {
		return unit2 == "kph" || unit2 == "kn" || unit2 == "fpm"
				|| unit2 == "m/s" || unit2 == "mph" || unit2 == "ft/s" || unit2 == "NM/s";
	} else if (unit1 == "m/s^2" || unit1 == "G" || unit1 == "ft/s^2") {
		return unit2 == "m/s^2" || unit2 == "G" || unit2 == "ft/s^2";
	} else if (unit1 == "N" || unit1 == "lbf") {
		return unit2 == "N" || unit2 == "lbf";
	} else if (unit1 == "atm" || unit1 == "Pa") {
		return unit2 == "atm" || unit2 == "Pa";
	} else if (unit1 == "unitless") {
		return unit2 == "unitless";
	}
	return false;
}

/**
 * Return the internal (factor 1.0) unit compatible with this unit
 * @param unit unit to check
 * @return unit with factor 1.0 that is compatible with this unit.  In the event the unit is recognized but has no factor 1.0 equivalent, this will return the empty string.
 */
std::string Units::getCompatibleInternalUnit(const std::string& unitp) {
	const std::string unit = Units::canonical(unitp);
	if (unit == "unspecified") {
		return "";
	}
	if (unit == "m" || unit == "ft" || unit == "km" || unit == "NM"
			|| unit == "mi" || unit == "in" || unit == "yard" || unit == "mm") {
		return  "m";
	} else if (unit == "ft^2" || unit == "m^2") {
		return "m^2";
	} else if (unit == "s" || unit == "sec" || unit == "min" || unit == "hour" || unit == "ms" || unit == "day") {
		return "s";
	} else if (unit == "rad" || unit == "deg") {
		return "rad";
	} else if (unit == "rad/s" || unit == "deg/s") {
		return "rad/s";
	} else if (unit == "kg" || unit == "slug" || unit == "lbm") {
		return "kg";
	} else if (unit == "kph" || unit == "kn" || unit == "fpm"
			|| unit == "m/s" || unit == "mph" || unit == "ft/s" || unit == "NM/s") {
		return "m/s";
	} else if (unit == "m/s^2" || unit == "G" || unit == "ft/s^2") {
		return "m/s^2";
	} else if (unit == "N" || unit == "lbf") {
		return "N";
	} else if (unit == "atm" || unit == "Pa") {
		return "Pa";
	} else if (unit == "unitless") {
		return "unitless";
	}
	return "";
}


// Converts value canonical units to [symbol] units
real Units::to(const std::string& symbol, real value) {
  return to(Units::getFactor(symbol), value);
}

// Converts value canonical units to [symbol] units
real Units::to(const real symbol, const real value) {
  if (symbol == degreeC) {
    return value/Units::K - real(273.15);
  } else if (symbol == degreeF) {
    return value/Units::degreeR - (real(273.15) * real(1.8) - real(32.0));
  } else {
    return value / symbol;
  }
}

real Units::from(const std::string& symbol, real value) {
  return from(Units::getFactor(symbol), value);
}

real Units::from(const real symbol, const real value) {
  if (symbol == degreeC) {
    return (value+real(273.15))*Units::K;
  } else if (symbol == degreeF) {
    return (value + (real(273.15) * real(1.8) - real(32.0))) * Units::degreeR;
  } else {
    return value * symbol;
  }
}

real Units::fromInternal(const std::string& defaultUnits, const std::string& units, real value) {
  if (units == "unspecified") {
    return from(getFactor(defaultUnits), value);
  } else {
    return value;
  }
}

std::string Units::str(const std::string& symbol, real value) {
	return str(symbol,value,Constants::get_output_precision());
}

std::string Units::strX(const std::string& symbol, real value) {
	std::ostringstream st;
	st << Units::to(symbol, value) << " [" << symbol << "]";
	return st.str();
}

std::string Units::str(const std::string& symbol, real value, int precision) {
	return FmPrecision(to(symbol, value),precision) + " [" + symbol + "]";
}

static void trimBuilder(std::string& sb) {
	char const* delims = " \t\r\n";

	// trim leading whitespace
	std::string::size_type notwhite = sb.find_first_not_of(delims);
	sb.erase(0, notwhite);

	// trim trailing whitespace
	notwhite = sb.find_last_not_of(delims);
	sb.erase(notwhite + 1);
}

std::string Units::cleanOnly(const std::string& unit) {

	std::string sb(unit);
	trimBuilder(sb);

	  size_t start_idx = sb.find_first_of("[");
	  if (start_idx != string::npos) {    //if (start_idx >= 0) {
		  sb.erase(0,start_idx+1);
	  }
	  size_t end_idx = sb.find_first_of("]");
	  if (end_idx != string::npos) {   //if (end_idx >= 0) {
		  sb.erase(end_idx);
	  }
	  trimBuilder(sb);

	return sb;
}

std::string Units::clean(const std::string& unit) {
	std::string sb(cleanOnly(unit));

	if (isUnit(sb)) {
		return sb;
	} else {
		return "unspecified";
	}
}

real Units::parse(const std::string& s) {
	return parse(s,real(0.0));
}

real Units::parse(const string& str, real default_value) {
	  	return parse("internal", str, default_value);
}

real Units::parse(const std::string& defaultUnitsFrom, const std::string& str) {
	  return parse(defaultUnitsFrom, str, real(0.0));
  }


#if defined(_MSC_VER)

real getd(string str, real def) {
	std::istringstream stream;
	stream.str(str);
	real d;
	stream >> d;
	if (stream.fail()) {
		return def;
	}
	return d;
}


real Units::parse(const string& defaultUnitsFrom, const std::string& s, real default_value) {
	real ret = real(0.0);
	std::smatch m;
	std::regex numre("\\s*([-+0-9\\.]+)\\s*\\[?\\s*([-/^_a-zA-Z0-9]*).*"); //(.*)");   We want to add this unicode character \u00B0 (the degree symbol) to the units part  //TODO: does not recognize e-notation
//	std::regex numre("\\s*([-+0-9\\.]+)\\s*\\[?\\s*([/^_a-zA-Z0-9]*)\\s*\\]?\\s*$"); //(.*)");   We want to add this unicode character \u00B0 (the degree symbol) to the units part
	//Java: Pattern.compile("\\s*([-+0-9\\.]+)\\s*\\[?\\s*([-\\/^_a-zA-Z0-9\\u00B0]*).*");

	std::regex_match(s, m, numre);
	std::smatch group(m);
	if (group.size() > 0) { // should always be two groups if the regex was matched at all!
		std::string m1(group[1]);
		std::string unit = Units::cleanOnly(group[2]);

		// The logic here can be debated.  What should be returned when
		// an invalid value or invalid unit is provided?  If parse("ft", "10 fjkdsj", 5)
		// is called, what should be returned? Units.from("ft", 10) or Units.from("ft", 5)?
		// I chose Units.from("ft", 10) because in the degenerate case of
		// parse("10 jfkdjks", 5)--that is, with an implied default unit of "internal"--
		// returning 10 seems more correct than returning 5.
		//
		//   supplied                returned (aka converted)
		//   value        unit       value          unit
		//   -----        ----       -----          ----
		//   valid        valid      suppliedvalue  suppliedunit
		//   valid        invalid    suppliedvalue  defunit
		//   invalid      valid      defvalue       defunit
		//   invalid      invalid    defvalue       defunit
		//   illformed               defvalue       defunit

		if (Util::is_double(m1)) {
			ret = Util::parse_double(m1);
			if (Units::isUnit(unit)) {
				ret = Units::from(unit, ret);
			} else {
				ret = Units::from(defaultUnitsFrom, ret);
			}
		} else {  // no value
			ret = Units::from(defaultUnitsFrom, default_value);
		}

	} else { // no match
		ret = Units::from(defaultUnitsFrom, default_value);
	}

	return ret;
}

std::string Units::parseUnits(const std::string& s) {
	std::string unit = "unspecified";
	std::smatch m;
	std::regex numre("\\s*([-+0-9\\.]+)\\s*\\[?\\s*([-/^_a-zA-Z0-9]*).*"); //(.*)");   We want to add this unicode character \u00B0 (the degree symbol) to the units part
//	std::regex numre("\\s*([-+0-9\\.]+)\\s*\\[?\\s*([/^_a-zA-Z0-9]*)\\s*\\]?\\s*$"); //(.*)");   We want to add this unicode character \u00B0 (the degree symbol) to the units part
	//Java: Pattern.compile("\\s*([-+0-9\\.]+)\\s*\\[?\\s*([-\\/^_a-zA-Z0-9\\u00B0]*).*");
	std::regex_match(s, m, numre);
	std::smatch group(m);
	if (group.size() > 0) {
		unit = Units::clean(group[2]);
	}
	return unit;
}

#else

real Units::parse(const string& defaultUnitsFrom, const std::string& s, real default_value) {
	real ret = Units::from(defaultUnitsFrom,default_value);
	regex_t regex;
	int reti;
	char msgbuf[100];

	/* Compile regular expression, use C regular expression library */
	reti = regcomp(&regex, "^[[:blank:]]*([-+0-9\\.]+)[[:blank:]]*\\[?[[:blank:]]*([-\\/^_a-zA-Z0-9]*).*", REG_EXTENDED); //[[:blank:]]*\\]?(.*)", REG_EXTENDED); //TODO: does not recognize e-notation
	//Java: Pattern.compile("\\s*([-+0-9\\.]+)\\s*\\[?\\s*([-\\/^_a-zA-Z0-9\\u00B0]*).*");
	if (reti != 0) {
		fdln("Could not compile regex in Units::parse\n");
		//		DebugSupport::halt();
		return ret;
	}
	std::unique_ptr<regex_t, decltype(&regfree)> compiled(&regex, &regfree);

	/* Execute regular expression */
	regmatch_t matchptr[4];
	reti = regexec(&regex, s.c_str(), 4, matchptr, 0);
	if (!reti) {
		char match[100];
		int numchars;

		// Match #1
		numchars = (int)matchptr[1].rm_eo - (int)matchptr[1].rm_so;
		strncpy(match,s.c_str()+matchptr[1].rm_so,numchars);
		match[numchars] = '\0';
		string m1(match);
		//fpln("match="+m1);

		// Match #2
		numchars = (int)matchptr[2].rm_eo - (int)matchptr[2].rm_so;
		strncpy(match,s.c_str()+matchptr[2].rm_so,numchars);
		match[numchars] = '\0';
		string m2(match);
		//fpln("match2="+m2);
		std::string unit = Units::cleanOnly(m2);
		//fpln("match22="+unit);

		// The logic here can be debated.  What should be returned when
		// an invalid value or invalid unit is provided?  If parse("ft", "10 fjkdsj", 5)
		// is called, what should be returned? Units.from("ft", 10) or Units.from("ft", 5)?
		// I chose Units.from("ft", 10) because in the degenerate case of
		// parse("10 jfkdjks", 5)--that is, with an implied default unit of "internal"--
		// returning 10 seems more correct than returning 5.
		//
		//   supplied                returned (aka converted)
		//   value        unit       value          unit
		//   -----        ----       -----          ----
		//   valid        valid      suppliedvalue  suppliedunit
		//   valid        invalid    suppliedvalue  defunit
		//   invalid      valid      defvalue       defunit
		//   invalid      invalid    defvalue       defunit
		//   illformed               defvalue       defunit

		if (Util::is_double(m1)) {
			real dbl = Util::parse_double(m1);
			if (Units::isUnit(unit)) {
				ret = Units::from(unit, dbl);
			} else {
				ret = Units::from(defaultUnitsFrom, dbl);
			}
		} else {
			//return the default value
		}


		//fpln("ret "+Fm12(ret));
	} else if (reti == REG_NOMATCH) {
		// no match, return default value
	} else {
		regerror(reti, &regex, msgbuf, sizeof(msgbuf));
		string m1(msgbuf);
		fdln("Regex match failed: "+m1);
	}

	return ret;
}

std::string Units::parseUnits(const std::string& s) {
	string ret("unspecified");
	regex_t regex;
	int reti;
	char msgbuf[100];

	/* Compile regular expression */
	//reti = regcomp(&regex, "^[[:blank:]]*([-+0-9\\.]+)[[:blank:]]*\\[?[[:blank:]]*([/^_a-zA-Z0-9]*)[[:blank:]]*\\]?(.*)", REG_EXTENDED);
	reti = regcomp(&regex, "^[[:blank:]]*([-+0-9\\.]+)[[:blank:]]*\\[?[[:blank:]]*([-\\/^_a-zA-Z0-9]*).*", REG_EXTENDED); //[[:blank:]]*\\]?(.*)", REG_EXTENDED); //TODO: does not recognize e-notation
	if (reti != 0) {
		fdln("$$$ERROR$$$: Could not compile regex in Units::parseUnits\n");
//		DebugSupport::halt();
		return ret;
	}
	std::unique_ptr<regex_t, decltype(&regfree)> compiled(&regex, &regfree);

	/* Execute regular expression */
	regmatch_t matchptr[4];
	char match[100];

	//fpln("input X"+s+"X");
	reti = regexec(&regex, s.c_str(), 4, matchptr, 0);
	if (!reti) {
		int numchars;

		// Match #2
		numchars = (int)matchptr[2].rm_eo - (int)matchptr[2].rm_so;
		strncpy(match,s.c_str()+matchptr[2].rm_so,numchars);
		match[numchars] = '\0';
		ret = match;
		ret = Units::clean(ret);
		//fpln("match2 "+ret);
	} else if (reti == REG_NOMATCH) {
		// no match, return default value
	} else {
		regerror(reti, &regex, msgbuf, sizeof(msgbuf));
		fprintf(stderr, "$$$ERROR$$$ Regex match failed in Units::parseUnits: %s\n", msgbuf);
//		DebugSupport::halt();
	}

	return ret;
}
#endif
