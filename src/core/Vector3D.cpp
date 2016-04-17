#import <math.h>
#include "Matrix3D.h"

using namespace flash::math;

Vec4::Vec4(float x, float y, float z, float w): x(x), y(y), z(z), w(w) {
};

float Vec4::distanceBetween(Vec4 const &vec1, Vec4 const &vec2) {
    float x(vec2.x - vec1.x);
    float y(vec2.y - vec1.y);
    float z(vec2.z - vec1.z);
    return _squareRootOfSquareSums(x, y, z);
}

float Vec4::dotProduct(Vec4 const &vec1, Vec4 const &vec2) {
    return vec2.x * vec1.x + vec2.y * vec1.y + vec2.z * vec1.z;
}

float Vec4::angleBetween(Vec4 const &vec1, Vec4 const &vec2) {
    float dotProd = dotProduct(vec1, vec2);
    float val = dotProd / vec1.length() / vec2.length();
    return (float) (acosf(val) * 180 / M_PI);
}

Vec4 Vec4::crossProduct(Vec4 const &vec1, Vec4 const &vec2) {
    Vec4 resultVec;
    resultVec.x = vec1.y * vec2.z - vec1.z * vec2.y;
    resultVec.y = vec1.z * vec2.x - vec1.x * vec2.z;
    resultVec.z = vec1.x * vec2.y - vec1.y * vec2.x;
    return resultVec;
}

float Vec4::_squareRootOfSquareSums(float a, float b, float c) {
    return sqrtf(a * a + b * b + c * c);
}

float Vec4::length() const {
    return _squareRootOfSquareSums(x, y, z);
}

void Vec4::setLength(float value) {
    normalize();
    multiplyByScalar(value);
}

void Vec4::add(Vec4 const &vec) {
    x += vec.x;
    y += vec.y;
    z += vec.z;
}

void Vec4::subtract(Vec4 const &vec) {
    x -= vec.x;
    y -= vec.y;
    z -= vec.z;
}

void Vec4::multiplyByMatrix(Matrix3D const &matrix) {
	Matrix3D::multiplyVectorByMatrix(*this, matrix);
}

Vec4 Vec4::clone() const {
    return Vec4(x, y, z, w);
}

bool Vec4::isEqualTo(Vec4 const &vec) const {
    return vec.x == x && vec.y == y && vec.z == z && vec.w == w;
}

void Vec4::multiplyByScalar(float scalar) {
    x *= scalar;
    y *= scalar;
    z *= scalar;
}

void Vec4::normalize() {
    multiplyByScalar(1 / length());
}

float Vec4::operator*(Vec4 &v) const {
    return dotProduct(*this, v);
}

Vec4 Vec4::operator+(Vec4 &v) const {
    Vec4 resultVec = *this;
    resultVec.add(v);
    return resultVec;
}

Vec4 Vec4::operator/(Vec4 &v) const {
    return Vec4::crossProduct(*this, v);
}

Vec4 Vec4::operator*(float scalar) const {
    Vec4 resultVector = *this;
    resultVector.multiplyByScalar(scalar);
    return resultVector;
}

Vec4 Vec4::operator-(Vec4 &v) const {
    Vec4 resultVector = *this;
    resultVector.subtract(v);
    return resultVector;
}

bool Vec4::operator==(Vec4 &v) const {
    return isEqualTo(v);
}

Vec4 Vec4::operator*(Matrix3D &m) const {
    Vec4 resultVector = *this;
    resultVector.multiplyByMatrix(m);
    return resultVector;
}
