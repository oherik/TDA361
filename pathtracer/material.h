#pragma once
#include <glm/glm.hpp>
#include "Pathtracer.h"
#include "sampling.h"
#include "spectrum.h"

using namespace glm; 

namespace pathtracer
{
	///////////////////////////////////////////////////////////////////////////
	// Changes we've made
	///////////////////////////////////////////////////////////////////////////
	extern struct CustomSettings {
		bool bumpmap;
		bool diffusemap;
		bool spectrum;
	} customSettings;

	///////////////////////////////////////////////////////////////////////////
	// for any BRDF. 
	///////////////////////////////////////////////////////////////////////////
    extern struct Brdf {
        int fresnelCurrent;
        int diffuseCurrent;
        int geometricCurrent;
    } brdf;


	///////////////////////////////////////////////////////////////////////////
	// The interface for any BRDF. 
	///////////////////////////////////////////////////////////////////////////
	class BRDF
	{
	public: 
		// Return the value of the brdf for specific directions
		virtual Spectrum f(const vec3 & wi, const vec3 & wo, const vec3 & n) = 0; 
		// Sample a suitable direction and return the brdf in that direction as
		// well as the pdf (~probability) that the direction was chosen. 
		virtual Spectrum sample_wi(vec3 & wi, const vec3 & wo, const vec3 & n, float & p) = 0;
	};


	class Transparent : public BRDF
	{
	public:

		float transparency;
		vec3 color;
		Transparent(float _transparency, vec3 _color) :
			transparency(_transparency), color(_color) {}
		virtual Spectrum f(const vec3 & wi, const vec3 & wo, const vec3 & n) override;
		virtual Spectrum sample_wi(vec3 & wi, const vec3 & wo, const vec3 & n, float & p) override;
	};

	///////////////////////////////////////////////////////////////////////////
	// A Lambertian (diffuse) material
	///////////////////////////////////////////////////////////////////////////
	class Diffuse : public BRDF
	{
	public: 
		vec3 color;
		Diffuse(vec3 c) : color(c) {}
		virtual Spectrum f(const vec3 & wi, const vec3 & wo, const vec3 & n) override;
		virtual Spectrum sample_wi(vec3 & wi, const vec3 & wo, const vec3 & n, float & p) override;
	};

	///////////////////////////////////////////////////////////////////////////
	// A Blinn Phong Dielectric Microfacet BRFD
	///////////////////////////////////////////////////////////////////////////
	class CustomDefined: public BRDF
	{
	public: 
		float shininess; 
		float R0; 
		BRDF * refraction_layer;
		CustomDefined(float _shininess, float _R0, BRDF * _refraction_layer = NULL) :
			shininess(_shininess), R0(_R0), refraction_layer(_refraction_layer) {}
		virtual Spectrum refraction_brdf(const vec3 & wi, const vec3 & wo, const vec3 & n);
		virtual Spectrum reflection_brdf(const vec3 & wi, const vec3 & wo, const vec3 & n);
		virtual Spectrum f(const vec3 & wi, const vec3 & wo, const vec3 & n) override;
		virtual Spectrum sample_wi(vec3 & wi, const vec3 & wo, const vec3 & n, float & p) override;
	};

	///////////////////////////////////////////////////////////////////////////
	// A Blinn Phong Metal Microfacet BRFD (extends the BlinnPhong class)
	///////////////////////////////////////////////////////////////////////////
	class CustomDefinedMetal: public CustomDefined
	{
	public: 
		vec3 m_n, m_k, color;
		float *lambda; //The wavelengths associated with the n and k values
		CustomDefinedMetal(vec3 color, vec3 m_n, vec3 m_k, float *lambda, float _shininess, float _R0) : color(color), m_n(m_n), m_k(m_k), lambda(lambda), CustomDefined(_shininess, _R0) {}
		virtual Spectrum refraction_brdf(const vec3 & wi, const vec3 & wo, const vec3 & n);
		virtual Spectrum reflection_brdf(const vec3 & wi, const vec3 & wo, const vec3 & n);
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
		virtual Spectrum f(const vec3 & wi, const vec3 & wo, const vec3 & n) override;
		virtual Spectrum sample_wi(vec3 & wi, const vec3 & wo, const vec3 & n, float & p) override;
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
		virtual Spectrum f(const vec3 & wi, const vec3 & wo, const vec3 & n) override;
		virtual Spectrum sample_wi(vec3 & wi, const vec3 & wo, const vec3 & n, float & p) override;
	};

	

}
