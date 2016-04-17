#pragma once

namespace flash {
    namespace math {

    class Matrix3D;

    class Vector3D {
    public:
        Vector3D(float x = 0, float y = 0, float z = 0, float w = 0) noexcept;

        static float distanceBetween(Vector3D const &vector1, Vector3D const &vector2) noexcept;

        static float dotProduct(Vector3D const &vector1, Vector3D const &vector2) noexcept;

        static float angleBetween(Vector3D const &vector1, Vector3D const &vector2) noexcept;

        static Vector3D crossProduct(Vector3D const &vector1, Vector3D const &vector2) noexcept;

        float length() const  noexcept;
        void setLength(float value) noexcept;

        void multiplyByScalar(float scalar) noexcept;

        void normalize() noexcept;

        void add(Vector3D const &vector) noexcept;

        void subtract(Vector3D const &vector) noexcept;

        bool isEqualTo(Vector3D const &vector) const noexcept;

        void multiplyByMatrix(Matrix3D const &matrix) noexcept;

        Vector3D clone() const;

        float operator*(Vector3D& v) const noexcept;

        Vector3D operator+(Vector3D& v) const noexcept;

        Vector3D operator-(Vector3D& v) const noexcept;

        Vector3D operator/(Vector3D& v) const noexcept;

        Vector3D operator*(float scalar) const noexcept;

        Vector3D operator*(Matrix3D& m) const noexcept;

        bool operator==(Vector3D& v) const noexcept;

        float& operator[](int index) { return *(&x + index); }
        const float& operator[](int index) const { return *(&x + index);  }

        float x, y, z, w;

    private:
        static float _squareRootOfSquareSums(float a, float b, float c) noexcept;
    };

    }
}
