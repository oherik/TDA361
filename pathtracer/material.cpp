#include "material.h"
#include "sampling.h"
#include <complex>

namespace pathtracer
{

    //Fresnell terms

    //Schlik linear approx
    float schlickFres(float R0, vec3 wh, vec3 wi){
		return R0 + (1 - R0)*pow(1 - dot(wh, wi), 5);
    }

    //Diffuse terms

    //Blinn Phong
    float blinnPhongDiff(float shininess, float ndotwh){
		return (shininess + 2.0f) / (2.0f * M_PI) * pow(ndotwh, shininess);
    }


    //Beckmann
    //m = sqrt(2/(shininess+2))
    float beckDiff(float ndotwh, float m){
        float ndotwh2 = ndotwh * ndotwh;
        float m2 = m * m;
        return exp((ndotwh2 - 1)/(m2 * ndotwh2))/(M_PI * m2 * ndotwh2 * ndotwh2);
    }


    //Geometric terms

    //Cook-torrance
    float cookGeom(float ndotwh, float ndotwo, float wodotwh, float ndotwi){
	    return min(1.0f, min(2.0f * ndotwh*ndotwo / wodotwh, 2.0f * ndotwh*ndotwi / wodotwh));
    }

    //Smith-Schlick
    //k = roughness * sqrt(2/M_PI)
    float smithSchlickGeom(float ndotwo, float ndotwi, float k){
        return ((ndotwi) / (ndotwi * (1-k) + k)) * ((ndotwo) / (ndotwo * (1-k) + k));
    }

    //Smith-Walter
    float smithWalterGeom(float ndotwo, float ndotwi, float wo, float wi, float wh, float n, float m){
        if (dot(wo,wh)/dot(wo,n) <= 0 || dot(wi,wh)/dot(wi,n) <= 0){
            return 0;
        }

        float awo = 1 / (m * tan(acos(ndotwo)));
        float first, second;
        if(awo < 1.6){
            float awo2 = awo * awo;
            first = (3.535f * awo + 2.181f * awo2) / (1.0f + 2.276f * awo + 2.577f * awo2);
        }else{
            first = 1;
        }

        float awi = 1 / (m * tan(acos(ndotwi)));
        if(awi < 1.6){
            float awi2 = awi2 * awi2;
            second = (3.535f * awi + 2.181f * awi2) / (1.0f + 2.276f * awi + 2.577f * awi2);
        }else{
            return first;
        }

        return first * second;
    }

    //

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
		float ndotwi = max(0.0f , dot(n, wi));
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

	vec3 BlinnPhongMetal::reflection_brdf(const vec3 & wi, const vec3 & wo, const vec3 & n) { 
		//Koppar
        //float m_n []= { 0.294f, 1.0697f, 1.2404f };
		//float m_k [] = { 3.2456f, 2.6866f, 2.3929f };
		
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
		//float n_m []= { 9.0f, 7.0f, 4.5f};
		//float k_m [] = { 5.0f, 4.3f, 3.4f };

		vec3 wh = normalize(wi + wo);

		float whdotwi = max(0.0f, dot(wh, wi));
		float ndotwh = max(0.0f, dot(n, wh));
		float wodotwh = max(0.0f, dot(wo, wh));
		float ndotwi = max(0.0f, dot(n, wi));
		float ndotwo = max(0.0f, dot(n, wo));

		vec3 n_normalized = normalize(n);
		vec3 wi_normalized = normalize(wi);

		float cost = dot(n_normalized, wi_normalized);

		float F_wi_1 = exactReflection(m_n[0], m_k[0], cost);
		float F_wi_2 = exactReflection(m_n[1], m_k[1], cost);
		float F_wi_3 = exactReflection(m_n[2], m_k[2], cost);

		float D_wh = (shininess + 2.0f) / (2.0f * M_PI) * pow(ndotwh, shininess);

		float G_wiwo = min(1.0f, min(2.0f * ndotwh*ndotwo / wodotwh, 2.0f * ndotwh*ndotwi / wodotwh));

		float den = (4.0f * ndotwo*ndotwi);

		if (den < EPSILON) return vec3(0.0f);

		return vec3(F_wi_1, F_wi_2, F_wi_3) * D_wh * G_wiwo / den;
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
