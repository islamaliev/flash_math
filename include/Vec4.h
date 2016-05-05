#pragma once

namespace flash {

namespace math {

    class Mat4;

    class Vec4 {
    public:
        Vec4(float x = 0, float y = 0, float z = 0, float w = 0);

        static float distanceBetween(const Vec4&, const Vec4&);

        static float dotProduct(const Vec4&, const Vec4&);

        static float angleBetween(const Vec4&, const Vec4&);

        static Vec4 crossProduct(const Vec4&, const Vec4&);

        float length() const;
        void setLength(float);

        void normalize();

        bool equals(const Vec4&) const;

        /**
         * Dot product
         */
        float operator|(const Vec4&) const;

        /**
         * Cross product
         */
        Vec4 operator^(const Vec4&) const;

        Vec4 operator+(const Vec4&) const;
        Vec4& operator+=(const Vec4&);

        Vec4 operator-(const Vec4&) const;
        Vec4& operator-=(const Vec4&);

        Vec4 operator/(const Vec4&) const;

        Vec4 operator*(float) const;
        Vec4& operator*=(float);

        Vec4 operator*(const Mat4&) const;
        Vec4& operator*=(const Mat4&);

        bool operator==(const Vec4&) const;

        bool operator!=(const Vec4& v) const { return !(*this == v); };

        float& operator[](int index) { return *(&x + index); }
        const float& operator[](int index) const { return *(&x + index);  }

        float x, y, z, w;
    };

}
}
