#include <math.h>
#include <stdexcept>
#include "Matrix3D.h"

using namespace flash::math;

static float ORTHOGONALIZE_FRACTION = 0.25;

Vector3D _deriveBasisVector(Vector3D const &vi, Vector3D const &vj, Vector3D const &vk) ;

Matrix3D::Matrix3D(float x1, float y1, float z1, float x2, float y2,
		float z2, float x3, float y3, float z3, float xt, float yt,
		float zt)
    : v1(x1, y1, z1, 0)
    , v2(x2, y2, z2, 0)
    , v3(x3, y3, z3, 0)
    , vt(xt, yt, zt, 1)
{}

Matrix3D Matrix3D::perspectiveProjection(float fovy, float aspectRatio, float near, float far) {
	Matrix3D matrix;
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

Matrix3D Matrix3D::orthographicProjection(float left, float right, float bottom, float top, float near, float far) {
	Matrix3D matrix;
    matrix.x1(2.0f / (right - left));
    matrix.y2(2.0f / (top - bottom));
    matrix.z3(2.0f / (near - far));
    matrix.xt((left + right) / (left - right));
    matrix.yt((bottom + top) / (bottom - top));
    matrix.zt((near + far) / (far - near));
    matrix.wt(1.0f);
	return matrix;
}

float Matrix3D::determinant() const {
    return v1.x * v2.y * v3.z + v1.y * v2.z * v3.x + v1.z * v2.x * v3.y -
            v1.x * v2.z * v3.y - v1.y * v2.x * v3.z - v1.z * v2.y * v3.x;
}

void Matrix3D::transpose() {
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

void Matrix3D::identity() {
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

void Matrix3D::multiplyByScalar(float scalar) {
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

void Matrix3D::multiplyByMatrix(Matrix3D const & m) {
    Matrix3D result;
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

void Matrix3D::rotateAboutX(float degrees) {
    float radians = (float) (degrees * M_PI / 180);
    float cosVal = cosf(radians);
    float sinVal = sinf(radians);
    v2.y = cosVal;
    v2.z = sinVal;
    v3.y = -sinVal;
    v3.z = cosVal;
}

void Matrix3D::rotateAboutY(float degrees) {
    float radians = (float) (degrees * M_PI / 180);
    float cosVal = cosf(radians);
    float sinVal = sinf(radians);
    v1.x = cosVal;
    v1.z = -sinVal;
    v3.x = sinVal;
    v3.z = cosVal;
}

void Matrix3D::rotateAboutZ(float degrees) {
    float radians = (float) (degrees * M_PI / 180);
    float cosVal = cosf(radians);
    float sinVal = sinf(radians);
    v1.x = cosVal;
    v1.y = sinVal;
    v2.x = -sinVal;
    v2.y = cosVal;
}

void Matrix3D::rotateAbout(Vector3D const &vector, float degrees) {
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

void Matrix3D::scaleAlong(Vector3D const &vector, float factor) {
	scaleAlong(vector.x, vector.y, vector.z, factor);
}

void Matrix3D::scaleAlong(float x, float y, float z, float factor) {
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

void Matrix3D::scale(float scaleX, float scaleY, float scaleZ) {
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

void Matrix3D::inverse() {
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

bool Matrix3D::isEqual(Matrix3D const &matrix) const {
    return matrix.x1() == v1.x && matrix.x2() == v2.x && matrix.x3() == v3.x && matrix.y1() == v1.y &&
            matrix.y2() == v2.y && matrix.y3() == v3.y && matrix.z1() == v1.z && matrix.z2() == v2.z &&
            matrix.z3() == v3.z;
}

bool Matrix3D::isClose(Matrix3D const &matrix, unsigned int precision) const {
    int factor = (int) powf(10, precision);
    return _areClose(matrix.x1(), v1.x, factor) && _areClose(matrix.y1(), v1.y, factor) && _areClose(matrix.z1(), v1.z,
    factor) && _areClose(matrix.x2(), v2.x, factor) && _areClose(matrix.y2(), v2.y,
    factor) && _areClose(matrix.z2(), v2.z, factor) && _areClose(matrix.x3(), v3.x,
    factor) && _areClose(matrix.y3(), v3.y, factor) && _areClose(matrix.z3(), v3.z, factor);
}

void Matrix3D::translate(float xt, float yt, float zt) {
	vt.x = xt;
	vt.y = yt;
	vt.z = zt;
}

void Matrix3D::transform(Vector3D &vector) const {
	multiplyVectorByMatrix(vector, *this);
}

Matrix3D Matrix3D::clone() const {
    Matrix3D matrix = Matrix3D(v1.x, v1.y, v1.z, v2.x, v2.y, v2.z,
            v3.x, v3.y, v3.z, vt.x, vt.y, vt.z);
    matrix.w1(v1.w);
    matrix.w2(v2.w);
    matrix.w3(v3.w);
    matrix.wt(vt.w);
    return matrix;
}

bool Matrix3D::isOrthogonal() const {
	Matrix3D transposedMatrix = clone();
	transposedMatrix.transpose();
	transposedMatrix.multiplyByMatrix(*this);
	return transposedMatrix.isClose(Matrix3D(), 5);
}

void Matrix3D::orthogonalize() {
	Vector3D basis = Vector3D(v1.x, v1.y, v1.z);
	basis.normalize();
	v1.x = basis.x;
	v1.y = basis.y;
	v1.z = basis.z;
	basis = Vector3D(v2.x, v2.y, v2.z);
	basis.normalize();
	v2.x = basis.x;
	v2.y = basis.y;
	v2.z = basis.z;
	basis = Vector3D(v3.x, v3.y, v3.z);
	basis.normalize();
	v3.x = basis.x;
	v3.y = basis.y;
	v3.z = basis.z;
	for (int i = 0; i < 5; i++) {
		_performOrthogonalizingAlgorithm();
	}
	_performGrandSchmidtOrthogonalizingAlgorithm();
}

void Matrix3D::_performGrandSchmidtOrthogonalizingAlgorithm() {
	Vector3D tempV1 = Vector3D(v1.x, v1.y, v1.z);
	Vector3D tempV2 = Vector3D(v2.x, v2.y, v2.z);
	tempV1.normalize();
	tempV2.normalize();
	v1.x = tempV1.x;
	v1.y = tempV1.y;
	v1.z = tempV1.z;

	float dp = Vector3D::dotProduct(tempV2, tempV1);
	Vector3D auxVec = tempV2.clone();
	auxVec.multiplyByScalar(dp);
	tempV2.subtract(auxVec);
	v2.x = tempV2.x;
	v2.y = tempV2.y;
	v2.z = tempV2.z;

	Vector3D tempV3 = Vector3D::crossProduct(tempV1, tempV2);
	v3.x = tempV3.x;
	v3.y = tempV3.y;
	v3.z = tempV3.z;
}

void Matrix3D::_performOrthogonalizingAlgorithm() {
	Vector3D tempV1 = Vector3D(v1.x, v1.y, v1.z);
	Vector3D tempV2 = Vector3D(v2.x, v2.y, v2.z);
	Vector3D tempV3 = Vector3D(v3.x, v3.y, v3.z);
	Vector3D basis1 = _deriveBasisVector(tempV1, tempV2, tempV3);
	Vector3D basis2 = _deriveBasisVector(tempV2, tempV1, tempV3);
	Vector3D basis3 = _deriveBasisVector(tempV3, tempV1, tempV2);
	v1.x = basis1.x;
	v1.y = basis1.y;
	v1.z = basis1.z;
	v2.x = basis2.x;
	v2.y = basis2.y;
	v2.z = basis2.z;
	v3.x = basis3.x;
	v3.y = basis3.y;
	v3.z = basis3.z;
}

Vector3D _deriveBasisVector(Vector3D const &vi, Vector3D const &vj, Vector3D const &vk) {
	float dp12 = Vector3D::dotProduct(vi, vj);
	float dp22 = Vector3D::dotProduct(vj, vj);
	float dp13 = Vector3D::dotProduct(vi, vk);
	float dp33 = Vector3D::dotProduct(vk, vk);
	Vector3D basis = vi.clone();
	Vector3D v2Clone = vj.clone();
	Vector3D v3Clone = vk.clone();
	v2Clone.multiplyByScalar(ORTHOGONALIZE_FRACTION * dp12 / dp22);
	v3Clone.multiplyByScalar(ORTHOGONALIZE_FRACTION * dp13 / dp33);
	basis.subtract(v2Clone);
	basis.subtract(v2Clone);
	return basis;
}

bool Matrix3D::_areClose(float value1, float value2, int factor) const {
    const float d1 = roundf(value1 * factor);
    const float d2 = roundf(value2 * factor);
    return d1 == d2;
}

void Matrix3D::multiplyVectorByMatrix(Vector3D &vector, const Matrix3D &matrix) {
	float x = vector.x;
	float y = vector.y;
	float z = vector.z;
	float newX = x * matrix.x1() + y * matrix.x2() + z * matrix.x3() + vector.w * matrix.xt();
	float newY = x * matrix.y1() + y * matrix.y2() + z * matrix.y3() + vector.w * matrix.yt();
    float zt = matrix.zt();
    vector.z = x * matrix.z1() + y * matrix.z2() + z * matrix.z3() + vector.w * zt;
	vector.x = newX;
	vector.y = newY;
}
