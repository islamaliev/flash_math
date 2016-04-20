#pragma once

namespace flash {
    namespace math {

    class Mat4;

    class Vec4 {
    public:
        Vec4(float x = 0, float y = 0, float z = 0, float w = 0);

        static float distanceBetween(Vec4 const &vector1, Vec4 const &vector2);

        static float dotProduct(Vec4 const &vector1, Vec4 const &vector2);

        static float angleBetween(Vec4 const &vector1, Vec4 const &vector2);

        static Vec4 crossProduct(Vec4 const &vector1, Vec4 const &vector2);

        float length() const;
        void setLength(float value);

        void multiplyByScalar(float scalar);

        void normalize();

        void add(Vec4 const &vector);

        void subtract(Vec4 const &vector);

        bool isEqualTo(Vec4 const &vector) const;

        void multiplyByMatrix(Mat4 const &matrix);

        Vec4 clone() const;

        float operator*(Vec4& v) const;

        Vec4 operator+(Vec4& v) const;

        Vec4 operator-(Vec4& v) const;

        Vec4 operator/(Vec4& v) const;

        Vec4 operator*(float scalar) const;

        Vec4 operator*(Mat4& m) const;

        bool operator==(Vec4& v) const;

        float& operator[](int index) { return *(&x + index); }
        const float& operator[](int index) const { return *(&x + index);  }

        float x, y, z, w;

    private:
        static float _squareRootOfSquareSums(float a, float b, float c);
    };

    }
}
