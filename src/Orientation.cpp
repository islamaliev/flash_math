#include <math.h>
#include "Orientation.h"

using namespace flash::math;

const static double pi = M_PI;
const static double pi2 = 2 * M_PI;

float Orientation::getShortestDifference(float angle1, float angle2) {
	return _wrapPi(angle2 - angle1);
}

float Orientation::_wrapPi(float theta) {
	// Check if already in range
	if (fabsf(theta) > pi) {
		// out of range. Determine how many "revolutions" we need to add
		double revolutions = floor((theta + pi) / pi2);
		theta -= revolutions * pi2;
	}
	return theta;
}
