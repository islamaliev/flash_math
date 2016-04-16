#import <math.h>
#include "Matrix3D.h"

using namespace flash::math;

Vector3D::Vector3D(float x, float y, float z, float w): x(x), y(y), z(z), w(w) {
};

float Vector3D::distanceBetween(Vector3D const &vec1, Vector3D const &vec2) {
    float x(vec2.x - vec1.x);
    float y(vec2.y - vec1.y);
    float z(vec2.z - vec1.z);
    return _squareRootOfSquareSums(x, y, z);
}

float Vector3D::dotProduct(Vector3D const &vec1, Vector3D const &vec2) {
    return vec2.x * vec1.x + vec2.y * vec1.y + vec2.z * vec1.z;
}

float Vector3D::angleBetween(Vector3D const &vec1, Vector3D const &vec2) {
    float dotProd = dotProduct(vec1, vec2);
    float val = dotProd / vec1.length() / vec2.length();
    return (float) (acosf(val) * 180 / M_PI);
}

Vector3D Vector3D::crossProduct(Vector3D const &vec1, Vector3D const &vec2) {
    Vector3D resultVec;
    resultVec.x = vec1.y * vec2.z - vec1.z * vec2.y;
    resultVec.y = vec1.z * vec2.x - vec1.x * vec2.z;
    resultVec.z = vec1.x * vec2.y - vec1.y * vec2.x;
    return resultVec;
}

float Vector3D::_squareRootOfSquareSums(float a, float b, float c) {
    return sqrtf(a * a + b * b + c * c);
}

void Vector3D::x(float value) {
    x = value;
}

void Vector3D::y(float value) {
    y = value;
}

void Vector3D::z(float value) {
    z = value;
}

void Vector3D::w(float value) {
    w = value;
}

float Vector3D::length() const {
    return _squareRootOfSquareSums(x, y, z);
}

void Vector3D::setLength(float value) {
    normalize();
    multiplyByScalar(value);
}

void Vector3D::add(Vector3D const &vec) {
    x += vec.x;
    y += vec.y;
    z += vec.z;
}

void Vector3D::subtract(Vector3D const &vec) {
    x -= vec.x;
    y -= vec.y;
    z -= vec.z;
}

void Vector3D::multiplyByMatrix(Matrix3D const &matrix) {
	Matrix3D::multiplyVectorByMatrix(*this, matrix);
}

Vector3D Vector3D::clone() const {
    return Vector3D(x, y, z, w);
}

bool Vector3D::isEqualTo(Vector3D const &vec) const {
    return vec.x == x && vec.y == y && vec.z == z && vec.w == w;
}

void Vector3D::multiplyByScalar(float scalar) {
    x *= scalar;
    y *= scalar;
    z *= scalar;
}

void Vector3D::normalize() {
    multiplyByScalar(1 / length());
}

float Vector3D::operator*(Vector3D &v) const {
    return dotProduct(*this, v);
}

Vector3D Vector3D::operator+(Vector3D &v) const {
    Vector3D resultVec = *this;
    resultVec.add(v);
    return resultVec;
}

Vector3D Vector3D::operator/(Vector3D &v) const {
    return Vector3D::crossProduct(*this, v);
}

Vector3D Vector3D::operator*(float scalar) const {
    Vector3D resultVector = *this;
    resultVector.multiplyByScalar(scalar);
    return resultVector;
}

Vector3D Vector3D::operator-(Vector3D &v) const {
    Vector3D resultVector = *this;
    resultVector.subtract(v);
    return resultVector;
}

bool Vector3D::operator==(Vector3D &v) const {
    return isEqualTo(v);
}

Vector3D Vector3D::operator*(Matrix3D &m) const {
    Vector3D resultVector = *this;
    resultVector.multiplyByMatrix(m);
    return resultVector;
}
