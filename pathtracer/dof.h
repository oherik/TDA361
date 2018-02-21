#pragma once
#include <glm/glm.hpp>
#include <vector>

#ifdef M_PI
#undef M_PI
#endif
#define M_PI 3.14159265359f
#define EPSILON 0.0001f

using namespace glm;

namespace pathtracer
{
	extern struct DepthOfField {
		float lensRadius;
		float focusDistance;
	} depthOfField;

	//PBRT
	vec2 dofConcentricSampleDisk(const vec2 &u);
}