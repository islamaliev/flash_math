#include <iostream>
#include "Mat4.h"

using namespace flash::math;

int main() {
    Mat4 m;
	Vec4 v;
	m.rotateAbout(v, 45);
    std::cout << "Test" << std::endl;
}

