#include <math.h>
#include "Mat4.h"

using namespace flash::math;

static float ORTHOGONALIZE_FRACTION = 0.25;

const Mat4 Mat4::IDENTITY = Mat4();

namespace {

    void _performGrandSchmidtOrthogonalizingAlgorithm(Mat4& m) {
        Vec4 tempV1 = Vec4(m.v1.x, m.v1.y, m.v1.z);
        Vec4 tempV2 = Vec4(m.v2.x, m.v2.y, m.v2.z);
        tempV1.normalize();
        tempV2.normalize();
        m.v1.x = tempV1.x;
        m.v1.y = tempV1.y;
        m.v1.z = tempV1.z;

        float dp = Vec4::dotProduct(tempV2, tempV1);
        Vec4 auxVec = tempV2;
        auxVec *= dp;
        tempV2 -= auxVec;
        m.v2.x = tempV2.x;
        m.v2.y = tempV2.y;
        m.v2.z = tempV2.z;

        Vec4 tempV3 = Vec4::crossProduct(tempV1, tempV2);
        m.v3.x = tempV3.x;
        m.v3.y = tempV3.y;
        m.v3.z = tempV3.z;
    }

    Vec4 _deriveBasisVector(Vec4 const &vi, Vec4 const &vj, Vec4 const &vk) {
        float dp12 = Vec4::dotProduct(vi, vj);
        float dp22 = Vec4::dotProduct(vj, vj);
        float dp13 = Vec4::dotProduct(vi, vk);
        float dp33 = Vec4::dotProduct(vk, vk);
        Vec4 basis = vi;
        Vec4 v2Clone = vj;
        Vec4 v3Clone = vk;
        v2Clone *= ORTHOGONALIZE_FRACTION * dp12 / dp22;
        v3Clone *= ORTHOGONALIZE_FRACTION * dp13 / dp33;
        basis -= v2Clone;
        basis -= v2Clone;
        return basis;
    }

    void _performOrthogonalizingAlgorithm(Mat4& m) {
        Vec4 tempV1 = Vec4(m.v1.x, m.v1.y, m.v1.z);
        Vec4 tempV2 = Vec4(m.v2.x, m.v2.y, m.v2.z);
        Vec4 tempV3 = Vec4(m.v3.x, m.v3.y, m.v3.z);
        Vec4 basis1 = _deriveBasisVector(tempV1, tempV2, tempV3);
        Vec4 basis2 = _deriveBasisVector(tempV2, tempV1, tempV3);
        Vec4 basis3 = _deriveBasisVector(tempV3, tempV1, tempV2);
        m.v1.x = basis1.x;
        m.v1.y = basis1.y;
        m.v1.z = basis1.z;
        m.v2.x = basis2.x;
        m.v2.y = basis2.y;
        m.v2.z = basis2.z;
        m.v3.x = basis3.x;
        m.v3.y = basis3.y;
        m.v3.z = basis3.z;
    }

    bool _areClose(float value1, float value2, int factor) {
        const float d1 = roundf(value1 * factor);
        const float d2 = roundf(value2 * factor);
        return d1 == d2;
    }

}

Mat4 Mat4::operator*(float scalar) const {
    Mat4 result(*this);
    result *= scalar;
    return result;
}

Mat4& Mat4::operator*=(float scalar) {
    v1.x *= scalar;
    v1.y *= scalar;
    v1.z *= scalar;
    v2.x *= scalar;
    v2.y *= scalar;
    v2.z *= scalar;
    v3.x *= scalar;
    v3.y *= scalar;
    v3.z *= scalar;
    return *this;
}

Mat4::Mat4(float x1, float y1, float z1, float x2, float y2,
		float z2, float x3, float y3, float z3, float xt, float yt,
		float zt)
    : v1(x1, y1, z1, 0)
    , v2(x2, y2, z2, 0)
    , v3(x3, y3, z3, 0)
    , vt(xt, yt, zt, 1)
{}

Mat4 Mat4::perspectiveProjection(float fovy, float aspectRatio, float near, float far) {
	Mat4 matrix;
	float zoomY = 1 / tanf(fovy * (float) M_PI / 360);
	float zoomX = zoomY / aspectRatio;
	float z = (near + far) / (near - far);
	float zt = (2 * near * far) / (near - far);
	matrix.x1(zoomX);
	matrix.y2(zoomY);
	matrix.z3(z);
	matrix.zt(zt);
	matrix.w3(-1);
	matrix.wt(0);
	return matrix;
}

Mat4 Mat4::orthographicProjection(float left, float right, float bottom, float top, float near, float far) {
	Mat4 matrix;
    matrix.x1(2.0f / (right - left));
    matrix.y2(2.0f / (top - bottom));
    matrix.z3(2.0f / (near - far));
    matrix.xt((left + right) / (left - right));
    matrix.yt((bottom + top) / (bottom - top));
    matrix.zt((near + far) / (far - near));
    matrix.wt(1.0f);
	return matrix;
}

float Mat4::determinant() const {
    return v1.x * v2.y * v3.z + v1.y * v2.z * v3.x + v1.z * v2.x * v3.y -
            v1.x * v2.z * v3.y - v1.y * v2.x * v3.z - v1.z * v2.y * v3.x;
}

void Mat4::transpose() {
    float temp = v2.x;
    v2.x = v1.y;
    v1.y = temp;
    temp = v3.x;
    v3.x = v1.z;
    v1.z = temp;
    temp = v2.z;
    v2.z = v3.y;
    v3.y = temp;
}

void Mat4::identity() {
    v1.x = 1;
    v1.y = 0;
    v1.z = 0;
    v2.x = 0;
    v2.y = 1;
    v2.z = 0;
    v3.x = 0;
    v3.y = 0;
    v3.z = 1;
}

void Mat4::multiplyByScalar(float scalar) {
    v1.x *= scalar;
    v1.y *= scalar;
    v1.z *= scalar;
    v2.x *= scalar;
    v2.y *= scalar;
    v2.z *= scalar;
    v3.x *= scalar;
    v3.y *= scalar;
    v3.z *= scalar;
}

void Mat4::multiplyByMatrix(Mat4 const & m) {
    Mat4 result;
    for (int j = 0; j < 4; j++)
    {
        for (int i = 0; i < 4; i++)
        {
            float sum(0);

            for (int n = 0; n < 4; n++)
            {
                sum += operator[](n)[i] * m[j][n];
            }

            result[j][i] = sum;
        }
    }

    v1.x = result[0][0];
    v2.x = result[1][0];
    v3.x = result[2][0];
    vt.x = result[3][0];
    v1.y = result[0][1];
    v2.y = result[1][1];
    v3.y = result[2][1];
    vt.y = result[3][1];
    v1.z = result[0][2];
    v2.z = result[1][2];
    v3.z = result[2][2];
    vt.z = result[3][2];
    v1.w = result[0][3];
    v2.w = result[1][3];
    v3.w = result[2][3];
    vt.w = result[3][3];
}

void Mat4::rotateAboutX(float degrees) {
    float radians = (float) (degrees * M_PI / 180);
    float cosVal = cosf(radians);
    float sinVal = sinf(radians);
    v2.y = cosVal;
    v2.z = sinVal;
    v3.y = -sinVal;
    v3.z = cosVal;
}

void Mat4::rotateAboutY(float degrees) {
    float radians = (float) (degrees * M_PI / 180);
    float cosVal = cosf(radians);
    float sinVal = sinf(radians);
    v1.x = cosVal;
    v1.z = -sinVal;
    v3.x = sinVal;
    v3.z = cosVal;
}

void Mat4::rotateAboutZ(float degrees) {
    float radians = (float) (degrees * M_PI / 180);
    float cosVal = cosf(radians);
    float sinVal = sinf(radians);
    v1.x = cosVal;
    v1.y = sinVal;
    v2.x = -sinVal;
    v2.y = cosVal;
}

void Mat4::rotateAbout(Vec4 const &vector, float degrees) {
    float radians = (float) (degrees * M_PI / 180);
    float  vx = vector.x;
    float  vy = vector.y;
    float  vz = vector.z;
    float cosVal = cosf(radians);
    float sinVal = sinf(radians);
    float cosOp = 1 - cosVal;
    float vzSin = vz * sinVal;
    float vySin = vy * sinVal;
    float vxSin = vx * sinVal;
    v1.x = vx * vx * cosOp + cosVal;
    v1.y = vx * vy * cosOp + vzSin;
    v1.z = vx * vz * cosOp - vySin;
    v2.x = vx * vy * cosOp - vzSin;
    v2.y = vy * vy * cosOp + cosVal;
    v2.z = vy * vz * cosOp + vxSin;
    v3.x = vx * vz * cosOp + vySin;
    v3.y = vy * vz * cosOp - vxSin;
    v3.z = vz * vz * cosOp + cosVal;
}

void Mat4::scaleAlong(Vec4 const &vector, float factor) {
	scaleAlong(vector.x, vector.y, vector.z, factor);
}

void Mat4::scaleAlong(float x, float y, float z, float factor) {
	float facOp = factor - 1;
	float y1x2 = facOp * x * y;
	float z1x3 = facOp * x * z;
	float z2y3 = facOp * y * z;
	v1.x *= 1 + facOp * x * x;
	v1.y *= y1x2;
	v1.z *= z1x3;
	v2.x *= y1x2;
	v2.y *= 1 + facOp * y * y;
	v2.z *= z2y3;
	v3.x *= z1x3;
	v3.y *= z2y3;
	v3.z *= 1 + facOp * z * z;
}

void Mat4::scale(float scaleX, float scaleY, float scaleZ) {
	v1.x *= scaleX;
	v1.y *= scaleX;
	v1.z *= scaleX;
	v2.x *= scaleY;
	v2.y *= scaleY;
	v2.z *= scaleY;
	v3.x *= scaleZ;
	v3.y *= scaleZ;
	v3.z *= scaleZ;
}

void Mat4::inverse() {
    float det = determinant();
    float newX1 = v2.y * v3.z - v2.z * v3.y; // cofactor 0, 0 +
    float newY1 = v2.z * v3.x - v2.x * v3.z; // cofactor 0, 1 -
    float newZ1 = v2.x * v3.y - v2.y * v3.x; // cofactor 0, 2 +
    float newX2 = v1.z * v3.y - v1.y * v3.z; // cofactor 1, 0 -
    float newY2 = v1.x * v3.z - v1.z * v3.x; // cofactor 1, 1 +
    float newZ2 = v1.y * v3.x - v1.x * v3.y; // cofactor 1, 2 -
    float newX3 = v1.y * v2.z - v1.z * v2.y; // cofactor 2, 0 +
    float newY3 = v1.z * v2.x - v1.x * v2.z; // cofactor 2, 1 -
    float newZ3 = v1.x * v2.y - v1.y * v2.x; // cofactor 2, 2 +
    v1.x = newX1;
    v1.y = newY1;
    v1.z = newZ1;
    v2.x = newX2;
    v2.y = newY2;
    v2.z = newZ2;
    v3.x = newX3;
    v3.y = newY3;
    v3.z = newZ3;
    transpose();
    multiplyByScalar(1 / det);
}

bool Mat4::isEqual(Mat4 const &matrix) const {
    return matrix.x1() == v1.x && matrix.x2() == v2.x && matrix.x3() == v3.x && matrix.w1() == v1.w
           && matrix.y1() == v1.y && matrix.y2() == v2.y && matrix.y3() == v3.y && matrix.w2() == v2.w
           && matrix.z1() == v1.z && matrix.z2() == v2.z && matrix.z3() == v3.z && matrix.w3() == v3.w
           && matrix.xt() == vt.x && matrix.yt() == vt.y && matrix.zt() == vt.z && matrix.wt() == vt.w;
}

bool Mat4::isClose(Mat4 const &matrix, unsigned int precision) const {
    int factor = (int) powf(10, precision);
    return _areClose(matrix.x1(), v1.x, factor) && _areClose(matrix.y1(), v1.y, factor) && _areClose(matrix.z1(), v1.z,
    factor) && _areClose(matrix.x2(), v2.x, factor) && _areClose(matrix.y2(), v2.y,
    factor) && _areClose(matrix.z2(), v2.z, factor) && _areClose(matrix.x3(), v3.x,
    factor) && _areClose(matrix.y3(), v3.y, factor) && _areClose(matrix.z3(), v3.z, factor);
}

void Mat4::translate(float xt, float yt, float zt) {
	vt.x = xt;
	vt.y = yt;
	vt.z = zt;
}

void Mat4::transform(Vec4 &vector) const {
    vector *= *this;
}

Mat4 Mat4::clone() const {
    Mat4 matrix = Mat4(v1.x, v1.y, v1.z, v2.x, v2.y, v2.z,
            v3.x, v3.y, v3.z, vt.x, vt.y, vt.z);
    matrix.w1(v1.w);
    matrix.w2(v2.w);
    matrix.w3(v3.w);
    matrix.wt(vt.w);
    return matrix;
}

bool Mat4::isOrthogonal() const {
	Mat4 transposedMatrix = clone();
	transposedMatrix.transpose();
	transposedMatrix.multiplyByMatrix(*this);
	return transposedMatrix.isClose(Mat4(), 5);
}

void Mat4::orthogonalize() {
	Vec4 basis = Vec4(v1.x, v1.y, v1.z);
	basis.normalize();
	v1.x = basis.x;
	v1.y = basis.y;
	v1.z = basis.z;
	basis = Vec4(v2.x, v2.y, v2.z);
	basis.normalize();
	v2.x = basis.x;
	v2.y = basis.y;
	v2.z = basis.z;
	basis = Vec4(v3.x, v3.y, v3.z);
	basis.normalize();
	v3.x = basis.x;
	v3.y = basis.y;
	v3.z = basis.z;
	for (int i = 0; i < 5; i++) {
		_performOrthogonalizingAlgorithm(*this);
	}
	_performGrandSchmidtOrthogonalizingAlgorithm(*this);
}

Mat4 Mat4::operator+(const Mat4& m) const {
    Mat4 result(*this);
    result += m;
    return result;
}

Mat4& Mat4::operator+=(const Mat4& m) {
    v1.x = v1.x + m.v1.x;
    v1.y = v1.y + m.v1.y;
    v1.z = v1.z + m.v1.z;
    v1.w = v1.w + m.v1.w;
    v2.x = v2.x + m.v2.x;
    v2.y = v2.y + m.v2.y;
    v2.z = v2.z + m.v2.z;
    v2.w = v2.w + m.v2.w;
    v3.x = v3.x + m.v3.x;
    v3.y = v3.y + m.v3.y;
    v3.z = v3.z + m.v3.z;
    v3.w = v3.w + m.v3.w;
    vt.x = vt.x + m.vt.x;
    vt.y = vt.y + m.vt.y;
    vt.z = vt.z + m.vt.z;
    vt.w = vt.w + m.vt.w;
    return *this;
}

Mat4 Mat4::operator*(const Mat4& m) const {
    Mat4 result(*this);
    result *= m;
    return result;
}

Mat4& Mat4::operator*=(const Mat4& m) {
    Mat4 result;
    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 4; ++i) {
            float sum = 0;

            for (int n = 0; n < 4; ++n) {
                sum += operator[](n)[i] * m[j][n];
            }

            result[j][i] = sum;
        }
    }

    v1.x = result[0][0];
    v2.x = result[1][0];
    v3.x = result[2][0];
    vt.x = result[3][0];
    v1.y = result[0][1];
    v2.y = result[1][1];
    v3.y = result[2][1];
    vt.y = result[3][1];
    v1.z = result[0][2];
    v2.z = result[1][2];
    v3.z = result[2][2];
    vt.z = result[3][2];
    v1.w = result[0][3];
    v2.w = result[1][3];
    v3.w = result[2][3];
    vt.w = result[3][3];
    return *this;
}
