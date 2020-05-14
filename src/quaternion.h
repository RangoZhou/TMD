#pragma once

#include <iostream>
#include "random.h"
#include "common.h"
#include "angle.h"
#include "vec3d.h"
#include "mat3x3.h"

namespace tmd {

// class Vec3d; // forward declaration
// class Mat3x3; // forward declaration

class Quaternion {
    Float data[4];
public:
    Quaternion(const Float r1, const Float r2, const Float r3, const Float r4) {
        data[0] = r1;
        data[1] = r2;
        data[2] = r3;
        data[3] = r4;
    }
    Quaternion(const Vec3d& ax, Float ang) { // angle in radian, axis is assumed to be a unit vector
        assert(eq(ax.norm(), 1));
        normalize_angle(ang); // this is probably only necessary if angles can be very big
        const Float c = std::cos(ang/2);
        const Float s = std::sin(ang/2);
        data[0] = c;
        data[1] = s*ax[0];
        data[2] = s*ax[1];
        data[3] = s*ax[2];
    }

    Quaternion(const Vec3d& rotation) {//angle in radian
        Float angle = rotation.norm();
        if(angle > k_epsilon) {
            Vec3d axis = rotation * (1/angle);
            assert(eq(axis.norm(), 1));
            normalize_angle(angle); // this is probably only necessary if angles can be very big
            const Float c = std::cos(angle/2);
            const Float s = std::sin(angle/2);
            data[0] = c;
            data[1] = s*axis[0];
            data[2] = s*axis[1];
            data[3] = s*axis[2];
        }
        assert(false);
    }

    const Float r1() const {
        return data[0];
    }
    const Float r2() const {
        return data[1];
    }
    const Float r3() const {
        return data[2];
    }
    const Float r4() const {
        return data[3];
    }
    const Float norm() const {
        return std::sqrt(data[0]*data[0]+data[1]*data[1]+data[2]*data[2]+data[3]*data[3]);
    }
    const Float norm_square() const {
        return data[0]*data[0]+data[1]*data[1]+data[2]*data[2]+data[3]*data[3];
    }

    void normalize() {
        const Float a = this->norm();
        assert(a > k_epsilon);
        data[0] = data[0]/a; data[1] = data[1]/a;
        data[2] = data[2]/a; data[3] = data[3]/a;
        assert(this->is_normalized());
    }

    void normalize_approx() {
        const Float s = this->norm_square();
        if(std::abs(s - 1) < k_tolerance)
            ; // most likely scenario
        else {
            this->normalize();
        }
    }

    const bool is_normalized() const {
        return eq(this->norm(),1.0);
    }

    const Quaternion conjugate() const {
        return Quaternion(data[0],-data[1],-data[2],-data[3]);
    }

    const Quaternion inverse() const {
        const Float s = this->norm_square();
        if(s > k_epsilon)
        return this->conjugate() * (1/s);
        assert(false);
    }

    //////////////////////////////////////////////////////

    const bool operator==(const Quaternion& rhs) {
        return eq(this->r1(), rhs.r1()) && \
               eq(this->r2(), rhs.r2()) && \
               eq(this->r3(), rhs.r3()) && \
               eq(this->r4(), rhs.r4());
    }

    const Quaternion operator*(const Float s) const {
        return Quaternion(data[0]*s,data[1]*s,data[2]*s,data[3]*s);
    }

    const Quaternion& operator*=(const Float s) {
        data[0] = data[0]*s;
        data[1] = data[1]*s;
        data[2] = data[2]*s;
        data[3] = data[3]*s;
        return *this;
    }

    const Quaternion operator*(const Quaternion& rhs) const{
        return Quaternion(
            data[0]*rhs.data[0] - data[1]*rhs.data[1] - data[2]*rhs.data[2] - data[3]*rhs.data[3],
            data[0]*rhs.data[1] + data[1]*rhs.data[0] + data[2]*rhs.data[3] - data[3]*rhs.data[2],
            data[0]*rhs.data[2] - data[1]*rhs.data[3] + data[2]*rhs.data[0] + data[3]*rhs.data[1],
            data[0]*rhs.data[3] + data[1]*rhs.data[2] - data[2]*rhs.data[1] + data[3]*rhs.data[0]);
    }

    const Quaternion& operator*=(const Quaternion& rhs) {
        data[0] = data[0]*rhs.data[0] - data[1]*rhs.data[1] - data[2]*rhs.data[2] - data[3]*rhs.data[3];
        data[1] = data[0]*rhs.data[1] + data[1]*rhs.data[0] + data[2]*rhs.data[3] - data[3]*rhs.data[2];
        data[2] = data[0]*rhs.data[2] - data[1]*rhs.data[3] + data[2]*rhs.data[0] + data[3]*rhs.data[1];
        data[3] = data[0]*rhs.data[3] + data[1]*rhs.data[2] - data[2]*rhs.data[1] + data[3]*rhs.data[0];
        return *this;
    }

    const Quaternion operator/(const Float rhs) const {
        assert(rhs > k_epsilon);
        return (*this) * (1/rhs);
    }

    const Quaternion& operator/=(const Float rhs) {
        assert(rhs > k_epsilon);
        (*this) *= (1/rhs);
        return *this;
    }

    const Quaternion operator/(const Quaternion& rhs) const {
        return (*this) * rhs.inverse();//norm size check is already in inverse
    }

    const Quaternion& operator/=(const Quaternion& rhs) {
        (*this) *= rhs.inverse();//norm size check is already in inverse
        return *this;
    }



    //////////////////////////////////////////////////////



    const Vec3d to_vec3d() const {//the direction of this Vec3d is rotation axis and the norm of this Vec3d is rotation angle in radian
        assert(this->is_normalized());
        const Float c = this->r1();
        if(c > -1 && c < 1) { // c may in theory be outside [-1, 1] even with approximately normalized q, due to rounding errors
            Float angle = 2*std::acos(c); // acos is in [0, pi]
            if(angle > k_pi) {
                angle -= 2*k_pi; // now angle is in [-pi, pi]
            }
            Vec3d axis(this->r2(), this->r3(), this->r4());
            const Float s = std::sin(angle/2); // perhaps not very efficient to calculate sin of acos
            if(std::abs(s) < k_epsilon) {
                return zero_vec3d;
            }
            axis *= (angle / s);
            return axis;
        } else { // when c = -1 or 1, angle/2 = 0 or pi, therefore angle = 0
            return zero_vec3d;
        }
    }

    Mat3x3 to_mat3x3() {
        assert(this->is_normalized());
        const Float a = this->r1();
        const Float b = this->r2();
        const Float c = this->r3();
        const Float d = this->r4();
        const Float aa = a*a;
        const Float ab = a*b;
        const Float ac = a*c;
        const Float ad = a*d;
        const Float bb = b*b;
        const Float bc = b*c;
        const Float bd = b*d;
        const Float cc = c*c;
        const Float cd = c*d;
        const Float dd = d*d;
        // std::cout << "a b c d ---> " << a << " " << b << " " << c << " " << d << std::endl;
        assert(eq(aa+bb+cc+dd, 1));
        Mat3x3 tmp;
        // from http://www.boost.org/doc/libs/1_35_0/libs/math/quaternion/TQE.pdf
        tmp(0, 0) = (aa+bb-cc-dd);
        tmp(0, 1) = 2*(-ad+bc);
        tmp(0, 2) = 2*(ac+bd);
        tmp(1, 0) = 2*(ad+bc);
        tmp(1, 1) = (aa-bb+cc-dd);
        tmp(1, 2) = 2*(-ab+cd);
        tmp(2, 0) = 2*(-ac+bd);
        tmp(2, 1) = 2*(ab+cd);
        tmp(2, 2) = (aa-bb-cc+dd);
        return tmp;
    }
    // rotation = angle * axis, where angle in radian, axis is assumed to be a unit vector
    void increment(const Vec3d& rotation) {
        assert(this->is_normalized());
        (*this) *= Quaternion(rotation);
        this->normalize_approx();
    }
    // angle in radian, axis is assumed to be a unit vector
    void increment(const Vec3d& ax, Float ang) {
        assert(eq(ax.norm(), 1));
        normalize_angle(ang); // this is probably only necessary if angles can be very big
        assert(this->is_normalized());
        (*this) *= Quaternion(ax,ang);
        this->normalize_approx();
    }
};

const Quaternion quaternion_identity(1, 0, 0, 0);

inline const Quaternion random_unit_quaternion(RNGType& generator) {
    Quaternion q(random_normal_distribution(0, 1, generator),//from random.h
                 random_normal_distribution(0, 1, generator),
                 random_normal_distribution(0, 1, generator),
                 random_normal_distribution(0, 1, generator));
    const Float nrm = q.norm();
    if(nrm > k_epsilon) {
        q *= 1/nrm;
        assert(q.is_normalized());
        return q;
    } else {
        return random_unit_quaternion(generator); // this call should almost never happen
    }
}

inline const Quaternion unit_quaternion_difference(const Quaternion& b, const Quaternion& a) { // rotation that needs to be applied to convert a to b
	assert(a.is_normalized());
	assert(b.is_normalized());
    return b / a;// b = tmp * a    =>   b * inv(a) = tmp
	// return quaternion_to_rotation(tmp); // already assert normalization
}

}
