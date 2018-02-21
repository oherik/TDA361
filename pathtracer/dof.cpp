#include "dof.h"

using namespace std;
using namespace glm;

namespace pathtracer
{
	vec2 dofConcentricSampleDisk(const vec2 &u) {
		vec2 uOffset = 2.0f * u - vec2(1.0f);
		if (uOffset.x == 0 && uOffset.y == 0) {
			return vec2(0.0f);
		}
		float theta, r;
		if (std::abs(uOffset.x) > std::abs(uOffset.y)) {
			r = uOffset.x;
			theta = M_PI / 4.0f * (uOffset.y / uOffset.x);
		}
		else {
			r = uOffset.y;
			theta = M_PI / 2.0f - M_PI / 4.0f * (uOffset.x / uOffset.y);
		}
		return r * vec2(std::cos(theta), std::sin(theta));
	}
}