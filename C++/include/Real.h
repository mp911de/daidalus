/*
 * Real.h
 *
 * Precision types and comparison tolerances adapted from the AIRGREP port.
 */

#ifndef REAL_H_
#define REAL_H_

#include <cstdint>
#include <limits>

namespace larcfm {

#ifdef DAIDALUS_SINGLE_PRECISION

typedef float real;
typedef std::int32_t real_bits;

// Units in the last place are 2^29 times coarser than those of double. PRECISION5 and PRECISION7
// keep their relative width. PRECISION9 and PRECISION13 would be less than one unit and share the
// floor of REAL_ULPS_FLOOR, and their thresholds near zero the floor of REAL_ZERO_FLOOR.
static const real_bits REAL_ULPS_FLOOR = 4;
static const real REAL_ZERO_FLOOR = 1.0e-6f;
static const real_bits REAL_ULPS5 = 1 << 11;
static const real_bits REAL_ULPS7 = 1 << 5;
static const real_bits REAL_ULPS9 = REAL_ULPS_FLOOR;
static const real_bits REAL_ULPS13 = REAL_ULPS_FLOOR;
static const real REAL_ZERO5 = 1.0e-5f;
static const real REAL_ZERO7 = REAL_ZERO_FLOOR;
static const real REAL_ZERO9 = REAL_ZERO_FLOOR;
static const real REAL_ZERO13 = REAL_ZERO_FLOOR;

#else

typedef double real;
typedef std::int64_t real_bits;

static const real_bits REAL_ULPS5 = 1LL << 40;
static const real_bits REAL_ULPS7 = 1LL << 34;
static const real_bits REAL_ULPS9 = 1LL << 27;
static const real_bits REAL_ULPS13 = 16348;
static const real REAL_ZERO5 = 1.0e-5;
static const real REAL_ZERO7 = 1.0e-7;
static const real REAL_ZERO9 = 1.0e-9;
static const real REAL_ZERO13 = 1.0e-13;

#endif

static const real_bits REAL_HIGH_BIT = std::numeric_limits<real_bits>::min();

static_assert(sizeof(real) == sizeof(real_bits), "real_bits holds the bits of a real");

}

#endif /* REAL_H_ */
