#pragma once


#include "Vector3D.h"

namespace flash {
	namespace math {

	class Vector3D;

	class Matrix3D {

	public:
		explicit Matrix3D(float x1, float y1 = 0, float z1 = 0, float x2 = 0,
				float y2 = 1, float z2 = 0, float x3 = 0, float y3 = 0,
				float z3 = 1, float xt = 0, float yt = 0, float zt = 0);

		Matrix3D() = default;

		static void multiplyVectorByMatrix(Vector3D& vector, const Matrix3D& matrix);

		static Matrix3D perspectiveProjection(float fovy, float aspectRatio, float near, float far);

		static Matrix3D orthographicProjection(float left, float right, float bottom, float top, float near, float far);

		float x1() const { return v1.x; };
		void x1(float value) { v1.x = value; }

		float x2() const { return v2.x; }
		void x2(float value) { v2.x = value; }

		float x3() const { return v3.x; }
		void x3(float value) { v3.x = value; }

		float y1() const { return v1.y; }
		void y1(float value) { v1.y = value; }

		float y2() const { return v2.y; }
		void y2(float value) { v2.y = value; }

		float y3() const { return v3.y; }
		void y3(float value) { v3.y = value; }

		float z1() const { return v1.z; }
		void z1(float value) { v1.z = value; }

		float z2() const { return v2.z; }
		void z2(float value) { v2.z = value; }

		float z3() const { return v3.z; }
		void z3(float value) { v3.z = value; }

		float xt() const { return vt.x; }
		void xt(float value) { vt.x = value; }

		float yt() const { return vt.y; }
		void yt(float value) { vt.y = value; }

		float zt() const { return vt.z; }
		void zt(float value) { vt.z = value; }

		float w1() const { return v1.w; }
		void w1(float value) { v1.w = value; }

		float w2() const { return v2.w; }
		void w2(float value) { v2.w = value; }

		float w3() const { return v3.w; }
		void w3(float value) { v3.w = value; }

		float wt() const { return vt.w; }
		void wt(float value) { vt.w = value; }

		float determinant() const;

		void transpose();

		void identity();

		void multiplyByScalar(float scalar);

		void multiplyByMatrix(const Matrix3D& m);

		void rotateAboutX(float degrees);

		void rotateAboutY(float degrees);

		void rotateAboutZ(float degrees);

		void rotateAbout(const Vector3D& vector, float degrees);

		void scaleAlong(const Vector3D& vector, float factor);
		void scaleAlong(float x, float y, float z, float factor);

		void scale(float scaleX, float scaleY, float scaleZ);

		void inverse();

		Matrix3D clone() const;

		bool isClose(const Matrix3D& matrix, unsigned precision) const;

		bool isEqual(const Matrix3D& matrix) const;

		bool isOrthogonal() const;

		void orthogonalize();

		void translate(float xt, float yt, float zt);

		void transform(Vector3D& vector) const;

		inline operator float*() { return &v1.x; }
		inline operator const float*() const { return &v1.x; }

		Vector3D v1 {1, 0, 0, 0};
		Vector3D v2 {0, 1, 0, 0};
		Vector3D v3 {0, 0, 1, 0};
		Vector3D vt {0, 0, 0, 1};

	private:
		Vector3D& operator[](int index) { return *(&v1 + index); }
		const Vector3D& operator[](int index) const { return *(&v1 + index); }

		bool _areClose(float value1, float value2, int factor) const;

		void _performGrandSchmidtOrthogonalizingAlgorithm();

		void _performOrthogonalizingAlgorithm();
	};

	}
}
