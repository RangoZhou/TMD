#pragma once

#include <cassert>
#include <cmath>
#include "common.h"

namespace tmd {

inline void normalize_angle(Float& x) { // subtract or add enough 2*k_pi's to make x be in [-k_pi, k_pi]
	if(x >  3*k_pi){ // very large
		Float n = (x-k_pi) / (2*k_pi); // how many 2*k_pi's do you want to subtract?
		x -= 2*k_pi*std::ceil(n); // ceil can be very slow, but this should not be called often
		normalize_angle(x);
	}else if(x < -3*k_pi){ // very small
		Float n = (-x-k_pi) / (2*k_pi); // how many 2*k_pi's do you want to add?
		x += 2*k_pi*std::ceil(n); // ceil can be very slow, but this should not be called often
		normalize_angle(x);
	}else if(x > k_pi){ // in (k_pi, 3*k_pi]
		x -= 2*k_pi;
	}else if(x < -k_pi){ // in [-3*k_pi,  -k_pi)
		x += 2*k_pi;
	}
	assert(x >= -k_pi && x <= k_pi);
	// in [-k_pi, k_pi]
}

inline const Float normalized_angle(Float x) {
	normalize_angle(x);
	return x;
}

}