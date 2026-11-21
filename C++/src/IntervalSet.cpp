/* 
 * IntervalSet
 *
 * Contact: Jeff Maddalon
 * Organization: NASA/Langley Research Center
 *
 * Copyright (c) 2011-2021 United States Government as represented by
 * the National Aeronautics and Space Administration.  No copyright
 * is claimed in the United States under Title 17, U.S.Code. All Other
 * Rights Reserved.
 */

#include "Real.h"
#include "IntervalSet.h"
#include "Units.h"
#include "format.h"
#include "ErrorLog.h"
#include <iostream>
#include <stdexcept>
#include <cstring>
#include <vector>

using namespace std;
using namespace larcfm;

IntervalSet::IntervalSet() {
	length = 0;
}

IntervalSet::IntervalSet(const IntervalSet& l) {
	length = l.length;

	memcpy(r,l.r,sizeof(Interval)*length);
}

IntervalSet::IntervalSet(const std::vector<Interval>& v) {//: error("IntervalSet") {
	if (v.size() > static_cast<std::size_t>(max_intervals)) {
		throw std::length_error("IntervalSet capacity exceeded");
	}
	length = static_cast<int>(v.size());
	for (int i = 0; i < length; ++i) {
		r[i] = v[i];
	}
}

std::vector<Interval> IntervalSet::toVector() const {
	std::vector<Interval> v = std::vector<Interval>();
	v.reserve(length); // note: reserve() is appropriate because we exclusively use push_back to populate the vector
	for (int i = 0; i < length; i++) {
		v.push_back(r[i]);
	}
	return v;
}

void IntervalSet::clear() {
	length = 0;
}

const Interval& IntervalSet::getInterval(int i) const {
	if (i >= length || i < 0) {
		return Interval::EMPTY;
	}

	return r[i];
}

int IntervalSet::size() const {
	return length;
}

bool IntervalSet::isEmpty() const {
	return size()==0;
}

bool IntervalSet::in(real x) const {
	return order(x) >= 0;
}

void IntervalSet::unions(const Interval& rn) {
	if (rn.isEmpty()) {
		return;   //nothing to add
	}

	int iLow = order(rn.low);
	int iHigh = order(rn.up);

	real low, high;
	int start, end;

	if (iLow < 0) {
		low = rn.low;
		start = -(iLow+1);
	} else {
		low = r[iLow].low;
		start = iLow;
	}

	if (iHigh < 0) {
		high = rn.up;
		end = -(iHigh+1)-1;
	} else {
		high = r[iHigh].up;
		end = iHigh;
	}

	remove(start,end-start+1);  // start to end inclusive
	insert(start, Interval(low, high));
} // union

void IntervalSet::unions(const IntervalSet& n) {
	IntervalSet result(*this);
	for (int i = 0; i < n.length; i++) {
		result.unions(n.r[i]);
	}
	*this = result;
}

/**
 * Union the given IntervalSet into the current IntervalSet. Set n is
 * unmodified. This method uses "almost" inequalities to compute the intersection.
 * 
 * @param n set
 */
void IntervalSet::almost_unions(const IntervalSet& n, INT64FM maxUlps) {
	IntervalSet result(*this);
	for (int i=0; i < n.size(); ++i) {
		result.almost_add(n.getInterval(i).low,n.getInterval(i).up,maxUlps);
	}
	*this = result;
}

/**
 * Add the given interval into this set. If this interval overlaps any
 * interval in the set, then the intervals are merged.
 * This method uses "almost" inequalities to compute the addition.
 */
void IntervalSet::almost_add(real l, real u) {
	almost_add(l,u,PRECISION_DEFAULT);
}

/**
 * Add the given interval into this set. If this interval overlaps any
 * interval in the set, then the intervals are merged.
 * This method uses "almost" inequalities to compute the addition.
 */
void IntervalSet::almost_add(real l, real u, INT64FM maxUlps) {
	if (Util::almost_less(l,u,maxUlps)) {
		IntervalSet result;
		bool go = false;
		for (int i=0; i < size(); ++i) {
			Interval ii = getInterval(i);
			if (go) {
				result.unions(ii);
			} else if ((Util::almost_leq(ii.low,l,maxUlps) && Util::almost_leq(l,ii.up,maxUlps)) ||
					(Util::almost_leq(l,ii.low,maxUlps) && Util::almost_leq(ii.low,u,maxUlps))) {
				l = Util::min(ii.low,l);
				u = Util::max(ii.up,u);
			} else if (Util::almost_less(u,ii.low,maxUlps)) {
				result.unions(Interval(l,u));
				result.unions(ii);
				go = true;
			} else {
				result.unions(ii);
			}
		}
		if (!go) {
			result.unions(Interval(l,u));
		}
		*this = result;
	}
}

/**
 * Intersect the given IntervalSet into the current IntervalSet. Set n is
 * unmodified. This method uses "almost" inequalities to compute the intersection.
 */
void IntervalSet::almost_intersect(const IntervalSet& n) {
	almost_intersect(n,PRECISION_DEFAULT);
}

/**
 * Intersect the given IntervalSet into the current IntervalSet. Set n is
 * unmodified. This method uses "almost" inequalities to compute the intersection.
 */
void IntervalSet::almost_intersect(const IntervalSet& n, INT64FM maxUlps) {
	IntervalSet result;
	if (!isEmpty() && !n.isEmpty()) {
		int i=0;
		int j=0;
		while (i < size() && j < n.size()) {
			Interval ii = getInterval(i);
			Interval jj = n.getInterval(j);
			if (Util::almost_leq(jj.low,ii.low,maxUlps) &&
					Util::almost_less(ii.low,jj.up,maxUlps)) {
				if (Util::almost_leq(ii.up,jj.up,maxUlps)) {
					result.unions(ii);
					++i;
				} else {
					result.unions(Interval(ii.low,jj.up));
					++j;
				}
			} else if (Util::almost_leq(ii.low,jj.low,maxUlps) &&
					Util::almost_less(jj.low,ii.up,maxUlps)) {
				if (Util::almost_leq(jj.up,ii.up,maxUlps)) {
					result.unions(jj);
					++j;
				} else {
					result.unions(Interval(jj.low,ii.up));
					++i;
				}
			} else if (Util::almost_leq(ii.up,jj.low,maxUlps)){
				++i;
			} else if (Util::almost_leq(jj.up,ii.low,maxUlps)){
				++j;
			}
		}
	}
	*this = result;
}

void IntervalSet::diff(const Interval& rn) {
	if (rn.isEmpty()) {
		return;   //nothing for set difference
	}
	if (rn.isSingle()) {
		return;   // If rn is a single-valued open interval
		// then there is nothing to remove
	}

	int iLow = order(rn.low);
	int iHigh = order(rn.up);

	if (iLow >= 0 && iLow == iHigh) {
		if (length >= max_intervals) {
			throw std::length_error("IntervalSet capacity exceeded");
		}
		real r_iHigh_up = r[iHigh].up;
		r[iLow] = Interval(r[iLow].low, rn.low);
		insert(iLow+1, Interval(rn.up, r_iHigh_up));

		return;
	}

	int start, end;

	if (iLow < 0) {
		start = -(iLow+1);
	} else {
		r[iLow] = Interval(r[iLow].low, rn.low);
		start = iLow+1;
	}

	if (iHigh < 0) {
		end = -(iHigh+1)-1;
	} else {
		r[iHigh] = Interval(rn.up, r[iHigh].up);
		end = iHigh-1;
	}

	remove(start,end-start+1);  // start to end inclusive
} // diff

void IntervalSet::diff(const IntervalSet& n) {
	IntervalSet result(*this);
	for (int i = 0; i < n.length; i++) {
		result.diff(n.r[i]);
	}
	*this = result;
}

void IntervalSet::removeSingle(real x, real width) {
	int i = order(x);
	if (i >= 0 && r[i].isSingle(width)) {
		remove(i);
	}
}

void IntervalSet::removeSingle(real x) {
	removeSingle(x, real(0.0));
}

void IntervalSet::sweepSingle(real width) {
	int i = 0;
	while (i < length) {
		if (r[i].isSingle(width)) {
			remove(i);
		} else {
			i++;
		}
	}
}

void IntervalSet::sweepSingle() {
	sweepSingle(real(0.0));
}

void IntervalSet::sweepBreaks(real width) {
	int i = 0;
	while (i < length-1) {
		if (r[i].up+width > r[i+1].low) {
			unions(Interval(r[i].low, r[i+1].up));
		} else {
			i++;
		}
	}
}


/* 
 * Insert the given interval at point i.  If the point i is greater than
 * the number of Intervals in the list, then add this interval to the
 * end.  
 */
void IntervalSet::insert(int i, const Interval& region) {
	if (region.isEmpty()) {
		return;
	}

	if (i < 0) {
		i = 0;
	}

	if (i > length) {
		i = length;
	}

	if (length >= max_intervals) {
		throw std::length_error("IntervalSet capacity exceeded");
	}
	for (int c = length; c > i; --c) {
		r[c] = r[c-1];
	}
	r[i] = region;
	++length;
} // insert

/* 
 * Remove Interval i and return it
 */
void IntervalSet::remove(int i) {
	if (i < 0 || i >= length ) {
		return;
		//return Interval();
	}

	//	Interval t = r[i];

	while (i + 1 < length) {
		r[i] = r[i+1];
		i++;
	}

	length--;

	//return t;
}

/* 
 * Remove the len number of intervals starting at i.
 */
void IntervalSet::remove(int i, int len) {
	for( int j = 0; j < len; j++) {
		remove(i);
	}
}

/* 
 * Find the point where x fits into the list.  If x is a
 * member, then the returned value is the index of the interval.
 * If x is not a member, then -index-1 is returned, where index is
 * the index of the interval after x.
 *
 */
int IntervalSet::order(real x) const {
	for (int i = 0; i < length; i++) {
		if (r[i].in(x)) {
			return i;
		}
		if (x < r[i].low) {
			return -i-1;
		}

	}

	return -length-1;
}

/** Print the contents of this IntervalSet */
std::string  IntervalSet::toString() const {
	std::string s = "";
	for (int i = 0; i < length; i++) {
		if (i != 0) {
			s += ", ";
		}
		s += "Interval ["+Fm0(i)+"]: ";
		s += r[i].toString();
	}
	return s;
}


/** Print the contents of this IntervalSet */
std::string IntervalSet::toString(const std::string& unit) const {
	std::string s = "";
	for (int i = 0; i < length; i++) {
		if (i != 0) {
			s += ", ";
		}
		s += "Interval ["+Fm0(i)+"]: ";
		s += r[i].toStringUnits(unit);
	}
	return s;
}



