#include "material.h"
#include "sampling.h"
#include <complex>

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
	// Refract yo
	///////////////////////////////////////////////////////////////////////////
	vec3 Transparent::f(const vec3 & wi, const vec3 & wo, const vec3 & n) {
		return vec3(0.0f);
	}

	vec3 Transparent::sample_wi(vec3 & wi, const vec3 & wo, const vec3 & n, float & p) {
		vec3 L = wo;
		
		
		if (dot(L,n) > 0.0f){
			// Strålen på väg in i materialet
			float n1 = 1.0;
			float n2 = 1.52f;
			float index = n1/n2; //TODO: this is for windows glass
			float w = index*dot(L, n);
			
			float k = sqrt(1.0f + (w - index) * (w + index));

			wi =  (w - k)*n - index * L;
			p = 1.0f;
			float cosineTerm = abs(dot(wi, n)); // För att bli av med den i senare uträkning och slippa mörk rand
			return vec3(1.0) / cosineTerm; 
		}
		else {
			//Strålen lämnar materialet
			vec3 n_new = -n;
			float n1 = 1.52f;
			float n2 = 1.0;
			float index = n1/n2; //TODO: this is for windows glass
			float w = index*dot(L, n_new);
			float k_t = 1.0f + (w - index) * (w + index);
			if (k_t < 0) { // we've got TIR, baby
				wi = reflect(L,n_new);
			} else {
				float k = sqrt(k_t);
				wi =  (w - k)*n_new - index * L;
			}
			
			p = 1.0f;

			float cosineTerm = abs(dot(wi, n)); // För att bli av med den i senare uträkning och slippa mörk rand
			return vec3(1.0) / cosineTerm; 
		}
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
	
		float whdotwi = max(0.0f, dot(wh, wi));
		float ndotwh = max(0.0f, dot(n, wh));
		float wodotwh = max(0.0f, dot(wo, wh));
		float ndotwi = max(0.0f, dot(n, wi));
		float ndotwo = max(0.0f, dot(n, wo));

		float F_wi = R0 + (1.0f - R0)*pow(1.0f - whdotwi, 5.0f);
		float D_wh = (shininess + 2.0f) / (2.0f * M_PI) * pow(ndotwh, shininess);
		float G_wiwo = min(1.0f, min(2.0f * ndotwh*ndotwo / wodotwh, 2.0f * ndotwh*ndotwi / wodotwh));

		float den = (4.0f * ndotwo*ndotwi);

		if (den < EPSILON) return vec3(0.0f);

		return vec3(F_wi*D_wh*G_wiwo / den);

	}

	vec3 BlinnPhong::f(const vec3 & wi, const vec3 & wo, const vec3 & n) {
		return reflection_brdf(wi, wo, n) + refraction_brdf(wi, wo, n); 
	}

	vec3 BlinnPhong::sample_wi(vec3 & wi, const vec3 & wo, const vec3 & n, float & p) {
		
		if (randf() < 0.5){
			vec3 tangent = normalize(perpendicular(n));
			vec3 bitangent = normalize(cross(tangent, n));
			float phi = 2.0f * M_PI * randf();
			float cos_theta = pow(randf(), 1.0f / (shininess + 1.0f));
			float sin_theta = sqrt(max(0.0f, 1.0f - cos_theta * cos_theta));
			vec3 wh = normalize(sin_theta * cos(phi) * tangent +
				sin_theta * sin(phi) * bitangent +
				cos_theta * n);
			if (dot(wo, n) <= 0.0f){
				return vec3(0.0f);
			}
			float p_wh = (shininess + 1.0f)*(pow(dot(n, wh), shininess)) / (2.0f * M_PI);
			float den = 4.0f * dot(wo, wh);
			if (den < EPSILON){
				return vec3(0.0f);
			}
			p = p_wh / den;
			p *= 0.5f;
			wi = normalize(reflect(-normalize(wo),normalize(wh)));	
			return reflection_brdf(wi, wo, n);
		}
		else {
			if (refraction_layer == NULL){
				return vec3(0.0f);
			}
			vec3 brdf = refraction_layer->sample_wi(wi, wo, n, p);	
			p *= 0.5f;
			vec3 wh = normalize(wo + wi);
			float F = R0 + (1.0f - R0) * pow(1.0f - abs(dot(wh, wi)), 5.0f);	
			return (1 - F) * brdf;
		}

		/*


		vec3 tangent = normalize(perpendicular(n));
		vec3 bitangent = normalize(cross(tangent, n));
		vec3 sample = cosineSampleHemisphere();
		wi = normalize(sample.x * tangent + sample.y * bitangent + sample.z * n);
		if (dot(wi, n) <= 0.0f) p = 0.0f;
		else p = max(0.0f, dot(n, wi)) / M_PI;
		return f(wi, wo, n); 
		*/
		
	}

	///////////////////////////////////////////////////////////////////////////
	// A Blinn Phong Metal Microfacet BRFD (extends the BlinnPhong class)
	///////////////////////////////////////////////////////////////////////////
	vec3 BlinnPhongMetal::refraction_brdf(const vec3 & wi, const vec3 & wo, const vec3 & n) { 
		return vec3(0.0f); 
	}

	float exactReflection(float n, float k, float cost) {
		return( 0.5*(

				(pow(n,2.0f) + pow(k,2.0f) - 2 * n * cost + pow(cost,2) )/
				(pow(n, 2.0f) + pow(k, 2.0f) + 2 * n * cost + pow(cost, 2)) +
			
			((pow(n, 2.0f) + pow(k, 2.0f))*pow(cost,2.0f)-2*n*cost+1) /
			((pow(n, 2.0f) + pow(k, 2.0f))*pow(cost, 2.0f) + 2 * n*cost + 1)
			 
			));
	}

	float rgToN(float r, float g) {
		return g * (1.0f - r) / (1.0f - r) + (1.0f - g)  * (1.0f - sqrt(r)) / (1.0f - sqrt(r));
	}
	float rgToK(float r, float g) {
		return sqrt(1.0 / (1.0 - r) * r * (pow((rgToN(r, g) + 1), 2.0f) - pow(rgToN(r, g) - 1, 2.0f)));
	}

	vec3 BlinnPhongMetal::reflection_brdf(const vec3 & wi, const vec3 & wo, const vec3 & n) { 
		//Koppar
		//float n_m []= { 0.294f, 1.0697f, 1.2404f };
	//float k_m [] = { 3.2456f, 2.6866f, 2.3929f };
		
		
		//Guld
		//float n_m []= { 0.15557f, 0.42415f,1.3831f };
		//float k_m [] = { 3.6024f, 2.4721f,1.9155f };

		//Aluminium
		//float n_m[] = { 1.5580f, 1.0152f, 0.63324f };
		//float k_m[] = { 7.7124f, 6.6273f, 5.4544f };

		//Platina
		//float n_m[] = { 0.040000f, 0.059582f, 0.052225f };
		//float k_m[] = { 2.6484f,  3.5974f, 4.4094f };

		//jadu
		float n_m []= { 1.0f, 1.0f,1.0f };
		float k_m [] = { 1.0f, 1.0f,0.1f };

		std::complex<float> c1(n_m[0], k_m[0]);
		std::complex<float> c2(n_m[1], k_m[1]);
		std::complex<float> c3(n_m[2], k_m[2]);

		std::complex<float> c1_c = conj(c1);
		std::complex<float> c2_c = conj(c2);
		std::complex<float> c3_c = conj(c3);

		float r0_1 = ((c1 - 1.0f)*(c1_c - 1.0f) / ((c1 + 1.0f)*(c1_c + 1.0f))).real();
		float r0_2 = ((c2 - 1.0f)*(c2_c - 1.0f) / ((c2 + 1.0f)*(c2_c + 1.0f))).real();
		float r0_3 = ((c3 - 1.0f)*(c3_c - 1.0f) / ((c3 + 1.0f)*(c3_c + 1.0f))).real();

		// COpy
		vec3 wh = normalize(wi + wo);

		float whdotwi = max(0.0f, dot(wh, wi));
		float ndotwh = max(0.0f, dot(n, wh));
		float wodotwh = max(0.0f, dot(wo, wh));
		float ndotwi = max(0.0f, dot(n, wi));
		float ndotwo = max(0.0f, dot(n, wo));

		float cost = abs(dot(wi, n) / (n.length()*wi.length()));


		float rValues[] = { 0.0f, 0.4f,0.4f };
		float gValues[] = { 0.0f, 1.0f,0.0f };



		float F_wi_1 = exactReflection(rgToN(rValues[0], gValues[0]), rgToK(rValues[0], gValues[0]), cost); //n_m[0], k_m[0], cost);// r0_1 + (1.0f - r0_1)*pow(1.0f - whdotwi, 5.0f);
		float F_wi_2 = exactReflection(rgToN(rValues[1], gValues[1]), rgToK(rValues[1], gValues[1]), cost);//r0_2 + (1.0f - r0_2)*pow(1.0f - whdotwi, 5.0f);
		float F_wi_3 = exactReflection(rgToN(rValues[2], gValues[2]), rgToK(rValues[2], gValues[2]), cost);//r0_3 + (1.0f - r0_3)*pow(1.0f - whdotwi, 5.0f);


		printf("%f     ", rgToN(rValues[0], gValues[0]));


		float F_wi_tot = (F_wi_1 + F_wi_2 + F_wi_3);
		float F_to_use = 0.0f;

		vec3 color;
		float r = randf();
		if (r*F_wi_tot < F_wi_1) {
			F_to_use = F_wi_1;
			color = vec3(1.0f, 0, 0);
		}
		else if (r*F_wi_tot < F_wi_1+F_wi_2) {
			F_to_use = F_wi_2;
			color = vec3(0, 1.0f, 0);
		}
		else {
			F_to_use = F_wi_3;
			color = vec3(0, 0, 1.0f);
		}


		float D_wh = (shininess + 2.0f) / (2.0f * M_PI) * pow(ndotwh, shininess);
		float G_wiwo = min(1.0f, min(2.0f * ndotwh*ndotwo / wodotwh, 2.0f * ndotwh*ndotwi / wodotwh));

		float den = (4.0f * ndotwo*ndotwi);

		if (den < EPSILON) return vec3(0.0f);



		return vec3(vec3(F_wi_1, F_wi_2, F_wi_3));// *D_wh*G_wiwo / den);



		//return BlinnPhong::reflection_brdf(wi, wo, n) * color; 
	};

	///////////////////////////////////////////////////////////////////////////
	// A Transparency Blend between two BRDFs
	///////////////////////////////////////////////////////////////////////////
	vec3 TransparencyBlend::f(const vec3 & wi, const vec3 & wo, const vec3 & n) {
		vec3 sample = transparency->f(wi, wo, n);
		return sample * a + color * (1.0f - a);
	}

	vec3 TransparencyBlend::sample_wi(vec3 & wi, const vec3 & wo, const vec3 & n, float & p) {
		if (randf() < a){
			vec3 brdf = transparency->sample_wi(wi, wo, n, p);
			p = p * a;
			return brdf;
		}
		else {
			vec3 brdf = transparency->sample_wi(wi, wo, n, p);
			p = p * (1 - a);
			float cosineTerm = abs(dot(wi, n)); // För att bli av med den i senare uträkning och slippa mörk rand
			return color / cosineTerm;
		}
	}

	///////////////////////////////////////////////////////////////////////////
	// A Linear Blend between two BRDFs
	///////////////////////////////////////////////////////////////////////////
	vec3 LinearBlend::f(const vec3 & wi, const vec3 & wo, const vec3 & n) {
		vec3 sampleA = bsdf0->f(wi, wo, n);
		vec3 sampleB = bsdf1->f(wi, wo, n);
		return sampleA * w + sampleB * (1.0f - w);
	}

	vec3 LinearBlend::sample_wi(vec3 & wi, const vec3 & wo, const vec3 & n, float & p) {
		if (randf() < w){
			vec3 brdf = bsdf0->sample_wi(wi, wo, n, p);
			p = p * w;
			return brdf;
		}
		else {
			vec3 brdf = bsdf1->sample_wi(wi, wo, n, p);
			p = p * (1- w);
			return brdf;
		}
	}

	///////////////////////////////////////////////////////////////////////////
	// A perfect specular refraction.
	///////////////////////////////////////////////////////////////////////////
}
