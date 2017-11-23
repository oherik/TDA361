#include "material.h"
#include "sampling.h"

namespace pathtracer
{
	///////////////////////////////////////////////////////////////////////////
	// A Lambertian (diffuse) material
	///////////////////////////////////////////////////////////////////////////
	vec3 Diffuse::f(const vec3 & wi, const vec3 & wo, const vec3 & n) {
		if (dot(wi, n) <= 0.0f) return vec3(0.0f);
		if (!sameHemisphere(wi, wo, n)) return vec3(0.0f);
		return (1.0f / M_PI) * color;
	}

	vec3 Diffuse::sample_wi(vec3 & wi, const vec3 & wo, const vec3 & n, float & p) {
		vec3 tangent = normalize(perpendicular(n));
		vec3 bitangent = normalize(cross(tangent, n));
		vec3 sample = cosineSampleHemisphere();
		wi = normalize(sample.x * tangent + sample.y * bitangent + sample.z * n);
		if (dot(wi, n) <= 0.0f) p = 0.0f;
		else p = max(0.0f, dot(n, wi)) / M_PI;
		return f(wi, wo, n);
	}

	///////////////////////////////////////////////////////////////////////////
	// A Blinn Phong Dielectric Microfacet BRFD
	///////////////////////////////////////////////////////////////////////////
	vec3 BlinnPhong::refraction_brdf(const vec3 & wi, const vec3 & wo, const vec3 & n) {
		vec3 wh = normalize(wi + wo);
		float F_wi = R0 + (1 - R0)*pow(1 - dot(wh, wi), 5);
		if (refraction_layer == NULL){
			return vec3(0.0f);
		}
		return vec3((1-F_wi)*refraction_layer->f(wi, wo, n));
	}
	vec3 BlinnPhong::reflection_brdf(const vec3 & wi, const vec3 & wo, const vec3 & n) {
		vec3 wh = normalize(wi + wo);
		if (length(wo) == 0 || dot(wi, n) < 0 || dot(wo, n) < 0){
			return vec3(0.0f);
		}

		float F_wi = R0 + (1 - R0)*pow(1 - dot(wh, wi), 5);
		float D_wh = (shininess + 2) / (2 * M_PI) * pow(dot(n, wh), shininess);
		float G_wiwo = min(1.0f, min(2 * dot(n, wh)*dot(n, wo) / dot(wo, wh), 2 * dot(n, wh)*dot(n, wi) / dot(wo, wh)));

		return vec3(F_wi*D_wh*G_wiwo / (4 * dot(n, wo)*dot(n, wi)));

	}

	vec3 BlinnPhong::f(const vec3 & wi, const vec3 & wo, const vec3 & n) {
		return reflection_brdf(wi, wo, n) + refraction_brdf(wi, wo, n); 
	}

	vec3 BlinnPhong::sample_wi(vec3 & wi, const vec3 & wo, const vec3 & n, float & p) {
		vec3 tangent = normalize(perpendicular(n));
		vec3 bitangent = normalize(cross(tangent, n));
		vec3 sample = cosineSampleHemisphere();
		wi = normalize(sample.x * tangent + sample.y * bitangent + sample.z * n);
		if (dot(wi, n) <= 0.0f) p = 0.0f;
		else p = max(0.0f, dot(n, wi)) / M_PI;
		return f(wi, wo, n); 
	}

	///////////////////////////////////////////////////////////////////////////
	// A Blinn Phong Metal Microfacet BRFD (extends the BlinnPhong class)
	///////////////////////////////////////////////////////////////////////////
	vec3 BlinnPhongMetal::refraction_brdf(const vec3 & wi, const vec3 & wo, const vec3 & n) { 
		return vec3(0.0f); 
	}
	vec3 BlinnPhongMetal::reflection_brdf(const vec3 & wi, const vec3 & wo, const vec3 & n) { 
		return BlinnPhong::reflection_brdf(wi, wo, n) * color; 
	};

	///////////////////////////////////////////////////////////////////////////
	// A Linear Blend between two BRDFs
	///////////////////////////////////////////////////////////////////////////
	vec3 LinearBlend::f(const vec3 & wi, const vec3 & wo, const vec3 & n) {
		return vec3(0.0); 
	}

	vec3 LinearBlend::sample_wi(vec3 & wi, const vec3 & wo, const vec3 & n, float & p) {
		p = 0.0f; 
		return vec3(0.0f);
	}

	///////////////////////////////////////////////////////////////////////////
	// A perfect specular refraction.
	///////////////////////////////////////////////////////////////////////////
}