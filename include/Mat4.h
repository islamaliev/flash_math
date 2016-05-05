#pragma once


#include "Vec4.h"

namespace flash {
	namespace math {

	class Vec4;

	class Mat4 {

	public:
		static const Mat4 IDENTITY;

		explicit Mat4(float x1, float y1 = 0, float z1 = 0
				, float x2 = 0, float y2 = 1, float z2 = 0
				, float x3 = 0, float y3 = 0, float z3 = 1
				, float xt = 0, float yt = 0, float zt = 0);

		Mat4() = default;

		static Mat4 perspectiveProjection(float fovy, float aspectRatio, float near, float far);

		static Mat4 orthographicProjection(float left, float right, float bottom, float top, float near, float far);

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

		void multiplyByMatrix(const Mat4& m);

		void rotateAboutX(float degrees);

		void rotateAboutY(float degrees);

		void rotateAboutZ(float degrees);

		void rotateAbout(const Vec4& vector, float degrees);

		void scaleAlong(const Vec4& vector, float factor);
		void scaleAlong(float x, float y, float z, float factor);

		void scale(float scaleX, float scaleY, float scaleZ);

		void inverse();

		Mat4 clone() const;

		bool isClose(const Mat4& matrix, unsigned precision) const;

		bool isEqual(const Mat4& matrix) const;

		bool isOrthogonal() const;

		void orthogonalize();

		void translate(float xt, float yt, float zt);

		void transform(Vec4& vector) const;

		inline operator float*() { return &v1.x; }
		inline operator const float*() const { return &v1.x; }

		Mat4 operator+(const Mat4& m) const;
		Mat4& operator+=(const Mat4& m);

		Mat4 operator*(const Mat4& m) const;
		Mat4& operator*=(const Mat4& m);

		Mat4 operator*(float scalar) const;
		Mat4& operator*=(float scalar);

		Vec4 operator*(const Vec4& v) const;
		friend Vec4& operator*=(Vec4& v, const Mat4& m);

		Vec4 v1 {1, 0, 0, 0};
		Vec4 v2 {0, 1, 0, 0};
		Vec4 v3 {0, 0, 1, 0};
		Vec4 vt {0, 0, 0, 1};

	private:
		Vec4& operator[](int index) { return *(&v1 + index); }
		const Vec4& operator[](int index) const { return *(&v1 + index); }

	};

	Vec4& operator*=(Vec4& v, const Mat4& m);

}
}

