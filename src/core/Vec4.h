#pragma once

namespace flash {
    namespace math {

    class Matrix3D;

    class Vec4 {
    public:
        Vec4(float x = 0, float y = 0, float z = 0, float w = 0) noexcept;

        static float distanceBetween(Vec4 const &vector1, Vec4 const &vector2) noexcept;

        static float dotProduct(Vec4 const &vector1, Vec4 const &vector2) noexcept;

        static float angleBetween(Vec4 const &vector1, Vec4 const &vector2) noexcept;

        static Vec4 crossProduct(Vec4 const &vector1, Vec4 const &vector2) noexcept;

        float length() const  noexcept;
        void setLength(float value) noexcept;

        void multiplyByScalar(float scalar) noexcept;

        void normalize() noexcept;

        void add(Vec4 const &vector) noexcept;

        void subtract(Vec4 const &vector) noexcept;

        bool isEqualTo(Vec4 const &vector) const noexcept;

        void multiplyByMatrix(Matrix3D const &matrix) noexcept;

        Vec4 clone() const;

        float operator*(Vec4& v) const noexcept;

        Vec4 operator+(Vec4& v) const noexcept;

        Vec4 operator-(Vec4& v) const noexcept;

        Vec4 operator/(Vec4& v) const noexcept;

        Vec4 operator*(float scalar) const noexcept;

        Vec4 operator*(Matrix3D& m) const noexcept;

        bool operator==(Vec4& v) const noexcept;

        float& operator[](int index) { return *(&x + index); }
        const float& operator[](int index) const { return *(&x + index);  }

        float x, y, z, w;

    private:
        static float _squareRootOfSquareSums(float a, float b, float c) noexcept;
    };

    }
}
