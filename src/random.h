#pragma once

#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>
#include <cassert>
#include <ctime> // for time (for seeding)
#include "vec3d.h"
#include "process.h"

namespace tmd {

using RNGType = boost::mt19937;

inline const Float random_double(Float a, Float b, RNGType& generator) { // expects a < b, returns rand in [a, b]
	assert(a < b); // BOOST also asserts a < b
	using distr = boost::uniform_real<Float>;
	boost::variate_generator<RNGType&, distr> r(generator, distr(a, b));
	Float tmp = r();
	assert(tmp >= a);
	assert(tmp <= b);
	return tmp;
}

inline const Float random_normal_distribution(Float mean, Float sigma, RNGType& generator) { // expects sigma >= 0
	assert(sigma >= 0); // BOOST asserts this as well
	using distr = boost::normal_distribution<Float>;
	boost::variate_generator<RNGType&, distr> r(generator, distr(mean, sigma));
	return r();
}

inline const int random_int(int a, int b, RNGType& generator) { // expects a <= b, returns rand in [a, b]
	assert(a <= b); // BOOST asserts this as well
	using distr = boost::uniform_int<int>;
	boost::variate_generator<RNGType&, distr> r(generator, distr(a, b));
	int tmp = r();
	assert(tmp >= a);
	assert(tmp <= b);
	return tmp;
}

// inline const Size_Type random_size_t(Size_Type a, Size_Type b, RNGType& generator) { // expects a <= b, returns rand in [a, b]
// 	assert(a <= b);
// 	assert(int(a) >= 0);
// 	assert(int(b) >= 0);
// 	int i = random_int(int(a), int(b), generator);
// 	assert(i >= 0);
// 	assert(i >= int(a));
// 	assert(i <= int(b));
// 	return static_cast<Size_Type>(i);
// }

inline const Vec3d random_inside_sphere(RNGType& generator) {
	while(true) { // on average, this will have to be run about twice
		Float r1 = random_double(-1, 1, generator);
		Float r2 = random_double(-1, 1, generator);
		Float r3 = random_double(-1, 1, generator);
		Vec3d tmp(r1, r2, r3);
		if(tmp*tmp < 1)
			return tmp;
	}
}

inline const Vec3d random_unit_vec3d(RNGType& generator) {
    Vec3d v(random_normal_distribution(0, 1, generator),
			random_normal_distribution(0, 1, generator),
			random_normal_distribution(0, 1, generator));
    const Float nrm = v.norm();
    if(nrm > k_epsilon) {
        v *= 1/nrm;
        assert(v.is_normalized());
        return v;
    } else {
        return random_unit_vec3d(generator); // this call should almost never happen
    }
}

inline const Vec3d random_in_box(const Vec3d& corner1, const Vec3d& corner2, RNGType& generator) { // expects corner1[i] < corner2[i]
	Vec3d tmp;
	for(int i = 0; i < tmp.size(); ++i){
		tmp[i] = random_double(corner1[i], corner2[i], generator);
	}
	return tmp;
}

inline const int auto_seed() { // make seed from PID and time
	return my_pid() * int(std::time(NULL));
}

}