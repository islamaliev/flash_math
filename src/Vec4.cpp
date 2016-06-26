#include <math.h>
#include "Mat4.h"

using namespace flash::math;

namespace {
    float _squareRootOfSquareSums(float a, float b, float c) {
        return sqrtf(a * a + b * b + c * c);
    }
}

Vec4::Vec4(float x, float y, float z, float w): x(x), y(y), z(z), w(w) {
};

float Vec4::distanceBetween(const Vec4& vec1, const Vec4& vec2) {
    float x(vec2.x - vec1.x);
    float y(vec2.y - vec1.y);
    float z(vec2.z - vec1.z);
    return _squareRootOfSquareSums(x, y, z);
}

float Vec4::dotProduct(const Vec4& vec1, const Vec4& vec2) {
    return vec1 | vec2;
}

float Vec4::angleBetween(const Vec4& vec1, const Vec4& vec2) {
    float dotProd = dotProduct(vec1, vec2);
    float val = dotProd / vec1.length() / vec2.length();
    return (float) (acosf(val) * 180 / M_PI);
}

Vec4 Vec4::crossProduct(const Vec4& vec1, const Vec4& vec2) {
    return vec1 ^ vec2;
}

float Vec4::length() const {
    return _squareRootOfSquareSums(x, y, z);
}

void Vec4::setLength(float value) {
    normalize();
    *this *= value;
}

bool Vec4::equals(const Vec4& vec) const {
    return vec.x == x && vec.y == y && vec.z == z && vec.w == w;
}

void Vec4::normalize() {
    *this *= 1 / length();
}

float Vec4::operator|(const Vec4 &v) const {
    return x * v.x + y * v.y + z * v.z;
}

Vec4 Vec4::operator^(const Vec4 &v) const {
    Vec4 resultVec;
    resultVec.x = y * v.z - z * v.y;
    resultVec.y = z * v.x - x * v.z;
    resultVec.z = x * v.y - y * v.x;
    return resultVec;
}

Vec4 Vec4::operator+(const Vec4 &v) const {
    Vec4 resultVec = *this;
    resultVec += v;
    return resultVec;
}

Vec4& Vec4::operator+=(const Vec4 &v) {
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
}

Vec4 Vec4::operator/(const Vec4 &v) const {
    return Vec4::crossProduct(*this, v);
}

Vec4 Vec4::operator*(float scalar) const {
    Vec4 resultVector = *this;
    resultVector *= scalar;
    return resultVector;
}

Vec4& Vec4::operator*=(float scalar) {
    x *= scalar;
    y *= scalar;
    z *= scalar;
    return *this;
}

Vec4 Vec4::operator-(const Vec4 &v) const {
    Vec4 resultVec = *this;
    resultVec -= v;
    return resultVec;
}

Vec4& Vec4::operator-=(const Vec4 &v) {
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
}

bool Vec4::operator==(const Vec4 &v) const {
    return equals(v);
}

Vec4 Vec4::operator*(const Mat4& m) const {
    Vec4 result = *this;
    result *= m;
    return result;
}

Vec4& Vec4::operator*=(const Mat4& m) {
    float newX = x * m.v1.x + y * m.v2.x + z * m.v3.x + w * m.vt.x;
    float newY = x * m.v1.y + y * m.v2.y + z * m.v3.y + w * m.vt.y;
    z = x * m.v1.z + y * m.v2.z + z * m.v3.z + w * m.vt.z;
    x = newX;
    y = newY;
    return *this;
}