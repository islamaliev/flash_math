#pragma once


namespace flash {
	namespace math {

	class Orientation {
	public:
		float static getShortestDifference(float angle1, float angle2);
	private:
		float static _wrapPi(float theta);
	};

	}
}
