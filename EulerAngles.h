#pragma once


#include "Mat4.h"
#include "Quaternion.h"

namespace flash {
	namespace math
	{

	class Quaternion;

	class EulerAngles
	{
	public:
		static EulerAngles fromUprightMatrix(const Mat4& matrix);

		static EulerAngles fromObjectMatrix(const Mat4& matrix);

		EulerAngles(float heading = 0, float pitch = 0, float bank = 0);

		float heading() const
		{
			return _heading;
		}

		void heading(float value);

		float pitch() const
		{
			return _pitch;
		}

		void pitch(float value);

		float bank() const
		{
			return _bank;
		}

		void bank(float value);

		bool isCanonical() const;

		void canonize();

		Mat4 toUprightMatrix() const;

		Mat4 toObjectMatrix() const;

		Quaternion toUprightQuaternion() const;

		Quaternion toObjectQuaternion() const;

	private:
		const static int HALF_CIRCLE = 180;
		const static int FULL_CIRCLE = 360;

		float _heading;
		float _pitch;
		float _bank;

		float _wrap180(float theta);
	};

	}
}