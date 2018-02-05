#pragma once
#include <glm/glm.hpp>
#include "Pathtracer.h"
#include "sampling.h"

using namespace glm; 

namespace pathtracer
{
	///////////////////////////////////////////////////////////////////////////
	// The interface for any BRDF. 
	///////////////////////////////////////////////////////////////////////////
	class BRDF
	{
	public: 
		// Return the value of the brdf for specific directions
		virtual vec3 f(const vec3 & wi, const vec3 & wo, const vec3 & n) = 0; 
		// Sample a suitable direction and return the brdf in that direction as
		// well as the pdf (~probability) that the direction was chosen. 
		virtual vec3 sample_wi(vec3 & wi, const vec3 & wo, const vec3 & n, float & p) = 0;
	};


	class Transparent : public BRDF
	{
	public:

		float transparency;
		vec3 color;
		Transparent(float _transparency, vec3 _color) :
			transparency(_transparency), color(_color) {}
		virtual vec3 f(const vec3 & wi, const vec3 & wo, const vec3 & n) override;
		virtual vec3 sample_wi(vec3 & wi, const vec3 & wo, const vec3 & n, float & p) override;
	};

	///////////////////////////////////////////////////////////////////////////
	// A Lambertian (diffuse) material
	///////////////////////////////////////////////////////////////////////////
	class Diffuse : public BRDF
	{
	public: 
		vec3 color;
		Diffuse(vec3 c) : color(c) {}
		virtual vec3 f(const vec3 & wi, const vec3 & wo, const vec3 & n) override;
		virtual vec3 sample_wi(vec3 & wi, const vec3 & wo, const vec3 & n, float & p) override;
	};

	///////////////////////////////////////////////////////////////////////////
	// A Blinn Phong Dielectric Microfacet BRFD
	///////////////////////////////////////////////////////////////////////////
	class BlinnPhong : public BRDF
	{
	public: 
		float shininess; 
		float R0; 
		BRDF * refraction_layer;
		BlinnPhong(float _shininess, float _R0, BRDF * _refraction_layer = NULL) :
			shininess(_shininess), R0(_R0), refraction_layer(_refraction_layer) {}
		virtual vec3 refraction_brdf(const vec3 & wi, const vec3 & wo, const vec3 & n);
		virtual vec3 reflection_brdf(const vec3 & wi, const vec3 & wo, const vec3 & n);
		virtual vec3 f(const vec3 & wi, const vec3 & wo, const vec3 & n) override;
		virtual vec3 sample_wi(vec3 & wi, const vec3 & wo, const vec3 & n, float & p) override;
	};

	///////////////////////////////////////////////////////////////////////////
	// A Blinn Phong Metal Microfacet BRFD (extends the BlinnPhong class)
	///////////////////////////////////////////////////////////////////////////
	class BlinnPhongMetal : public BlinnPhong
	{
	public: 
		vec3 color, m_r, m_g; 
		BlinnPhongMetal(vec3 c, vec3 r, vec3 g, float _shininess, float _R0) : color(c), m_r(r), m_g(g), BlinnPhong(_shininess, _R0) {}
		virtual vec3 refraction_brdf(const vec3 & wi, const vec3 & wo, const vec3 & n);
		virtual vec3 reflection_brdf(const vec3 & wi, const vec3 & wo, const vec3 & n);
	};

	///////////////////////////////////////////////////////////////////////////
	// A Linear Blend between two BRDFs
	///////////////////////////////////////////////////////////////////////////
	class LinearBlend : public BRDF
	{
	public: 
		float w; 
		BRDF * bsdf0;
		BRDF * bsdf1;
		LinearBlend(float _w, BRDF * a, BRDF * b) : w(_w), bsdf0(a), bsdf1(b) {};
		virtual vec3 f(const vec3 & wi, const vec3 & wo, const vec3 & n) override; 
		virtual vec3 sample_wi(vec3 & wi, const vec3 & wo, const vec3 & n, float & p) override; 
	};


	///////////////////////////////////////////////////////////////////////////
	// A Transparency Blend 
	///////////////////////////////////////////////////////////////////////////

	class TransparencyBlend : public BRDF
	{
	public:
		float a;
		BRDF * transparency;
		vec3 color;
		TransparencyBlend(float _a, BRDF * _transparency, vec3 _color) : a(_a), transparency(_transparency), color(_color) {};
		virtual vec3 f(const vec3 & wi, const vec3 & wo, const vec3 & n) override;
		virtual vec3 sample_wi(vec3 & wi, const vec3 & wo, const vec3 & n, float & p) override;
	};

}