#pragma once

#include <cassert>
#include "common.h"
#include "vec3d.h"

namespace tmd {

class Vec3d;//forwar declaration

class Mat3x3 {// column-major
	Float data[9];
public:
	Mat3x3() {//default initialization with nan
		// assert(!(
		// 	std::isnan(data[0]) || std::isnan(data[1]) || std::isnan(data[2]) ||
		// 	std::isnan(data[3]) || std::isnan(data[4]) || std::isnan(data[5]) ||
		// 	std::isnan(data[6]) || std::isnan(data[7]) || std::isnan(data[8])
		// ));
		data[0] = data[1] = data[2] =
		data[3] = data[4] = data[5] =
		data[6] = data[7] = data[8] = k_not_a_num;
	}
	//initialize in column-major
	Mat3x3(Float xx, Float xy, Float xz,
		   Float yx, Float yy, Float yz,
		   Float zx, Float zy, Float zz) {
		   data[0] = xx; data[3] = xy; data[6] = xz;
		   data[1] = yx; data[4] = yy; data[7] = yz;
		   data[2] = zx; data[5] = zy; data[8] = zz;
	}
	const Float& operator()(int i, int j) const {
		#ifdef DEBUG
		assert(i < 3 && i >= 0);
		assert(j < 3 && j >= 0);
		#endif
		return data[i + 3*j];
	}
	Float& operator()(int i, int j) {
		#ifdef DEBUG
		assert(i < 3 && i >= 0);
		assert(j < 3 && j >= 0);
		#endif
		return data[i + 3*j];
	}
	Vec3d operator*(const Vec3d& v) const {
		return Vec3d(data[0]*v[0] + data[3]*v[1] + data[6]*v[2],
					data[1]*v[0] + data[4]*v[1] + data[7]*v[2],
					data[2]*v[0] + data[5]*v[1] + data[8]*v[2]);
	}
	const Mat3x3& operator*=(Float s) {
		for(int i = 0; i < 9; ++i) {
			data[i] *= s;
		}
		return *this;
	}

};

inline const Vec3d operator*(const Vec3d& lhs, const Mat3x3& rhs) {
	return rhs * lhs;
}

}