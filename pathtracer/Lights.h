#ifndef LIGHTS_H
#define LIGHTS_H
#pragma once
#include <glm/glm.hpp>
#include "sampling.h"
#include "spectrum.h"
#include "embree.h"
#include <iostream>
#include <memory>
#include <algorithm>
#ifndef M_PI
#define M_PI 3.14159265359f
#endif


using namespace std;
using namespace glm; 

namespace pathtracer{

    class Shape{
        public:
			Shape(mat4 *_worldToObject, mat4 *_objectToWorld)
				: worldToObject(_worldToObject), objectToWorld(_objectToWorld) {};
            virtual bool Intersect(Ray &ray, float *tHit, Intersection *hit) = 0;
            virtual float area() = 0;
            virtual Intersection Sample( vec2 &u) = 0;
            //virtual vec3 UniformSampleSphere( vec2 &u);
            float Pdf( Intersection &ref,  vec3 &wi);
            mat4 *worldToObject, *objectToWorld;
    };

    class Sphere : public Shape {
        private:
             float radius, zMin, zMax;
             float thetaMin, thetaMax, phiMax;
        public:
            Sphere(mat4 *_objectToWorld, mat4 *_worldToObject, float _radius, float _zMin, float _zMax, float _phiMax) : 
                Shape(_objectToWorld, _worldToObject),
                radius(_radius), 
                zMin(clamp(std::min(_zMin,_zMax), -radius, radius)), 
                zMax(clamp(std::max(_zMin, _zMax), -radius, radius)),
                thetaMin(acos(clamp(zMin / radius, -1.0f, 1.0f))),
                thetaMax(acos(clamp(zMax / radius, -1.0f, 1.0f))),
                phiMax(radians(clamp(_phiMax, 0.0f, 360.0f))) {}

            bool Intersect( Ray &_ray, float *tHit, Intersection *hit)  override;
            float area() override;
            Intersection Sample( vec2 &u)  override;
            vec3 UniformSampleSphere( vec2 &u);
			//float Pdf( Intersection &ref,  vec3 &wi)  override;
	};
	
	class Disk : public Shape {
	public:
		Disk(mat4 *_objectToWorld, mat4 *_worldToObject,
			float _height, float _radius, float _innerRadius, float _phiMax) :
			Shape(_objectToWorld, _worldToObject),
			height(_height), radius(_radius), innerRadius(_innerRadius), phiMax(M_PI / 180 * (Clamp(_phiMax, 0, 360))) {}
		bool Intersect( Ray &_ray, float * tHit, Intersection * hit) override;
		float area()  override;
		vec2 UniformSampleSphere( vec2 &u);
		Intersection Sample( vec2 &u) override;
		//float Pdf( Intersection &ref,  vec3 &wi)  override;
	private:
		 float height, radius, innerRadius, phiMax;
	};
	
    
    enum class LightFlags : int {
            DeltaPosition = 1, DeltaDirection = 2, Area = 4, Infinite = 8
    };
        
	class Light{
        public: 
            Light(int flags, mat4 * _LightToWorld, int nSamples = 1)
            : flags(flags), nSamples(std::max(1, nSamples)), LightToWorld(*_LightToWorld),
            WorldToLight(inverse(*_LightToWorld)) {
            };
            int flags;
            int nSamples;
            virtual Spectrum Sample_Li( Intersection &ref, Intersection *lightHit,  vec2 &u, vec3 *wi, float *pdf) = 0;
            virtual Spectrum Power() = 0;
			Spectrum Le(const Ray &ray) const {
				return Spectrum(0.f);
			}

        protected:
             mat4 LightToWorld, WorldToLight;
    };


    class AreaLight : public Light {
        public:
			Spectrum lEmit;
            AreaLight(mat4 * LightToWorld, Spectrum _lEmit, int nSamples = 1)
                : Light(4, LightToWorld, nSamples), lEmit(_lEmit){};
            virtual Spectrum L( vec3 &coordinate,  vec3 &w) = 0; 
            float pdf(vec3 lightPos, Intersection hit, vec3 n, vec3 wi);
			virtual float getArea() = 0;
			Spectrum getLEmit();
    };


    class DiffuseAreaLight : public AreaLight {
        public:
            float area;
            Shape * shape;
            DiffuseAreaLight(mat4 * LightToWorld, Spectrum _lEmit, int nSamples, Shape *shape)
            : AreaLight(LightToWorld, _lEmit, nSamples),
            shape(shape), area(shape->area()) {};
     
            float Pdf_Li( Intersection &ref,  vec3 &wi);
            Spectrum L( vec3 &n,  vec3 &w) override;
            Spectrum Sample_Li( Intersection &ref, Intersection *lightHit,  vec2 &u, vec3 *wi, float *pdf) override;
            Spectrum Power() override;
			float getArea() override;
	};
	
	bool quadratic(float a, float b, float c, float *t0, float *t1);
	inline constexpr float gamma(int n);

	inline float BitsToFloat(uint32_t ui) {
		float f;
		memcpy(&f, &ui, sizeof(uint32_t));
		return f;
	}

	inline uint32_t FloatToBits(float f) {
		uint32_t ui;
		memcpy(&ui, &f, sizeof(float));
		return ui;
	}

	inline float NextFloatUp(float v) {
		if (std::isinf(v) && v > 0.) {
			return v;
		}
		if (v == -0.f) {
			v = 0.f;
		}

		uint32_t ui = FloatToBits(v);
		if (v >= 0) ++ui;
		else --ui;
		return BitsToFloat(ui);

	}
	
	inline float NextFloatDown(float v) {
		if (std::isinf(v) && v > 0.) {
			return v;
		}
		if (v == 0.f) {
			v = -0.f;
		}
		
		uint32_t ui = FloatToBits(v);
		if (v >= 0) { --ui; }
		else { ++ui; }
		return BitsToFloat(ui);
	}

	class EFloat {
	public:
		EFloat() { }
		EFloat(float v, float err = 0.f) : v(v), err(err) {
		#ifndef NDEBUG
		ld = v;
		#endif // NDEBUG
		}
		float v;
		float err;
		#ifndef NDEBUG
		long double ld;
		#endif // NDEBUG
		
		EFloat operator+(EFloat f) {
			EFloat r;
			r.v = v + f.v;
			#ifndef NDEBUG
			r.ld = ld + f.ld;
			#endif // DEBUG
			r.err = err + f.err +
				gamma(1) * (std::abs(v + f.v) + err + f.err);
			return r;
		}

		EFloat operator-(EFloat f) {
			EFloat r;
			r.v = v - f.v;
			#ifndef NDEBUG
			r.ld = ld - f.ld;
			#endif // DEBUG
			r.err = err +	 f.err +
				gamma(1) * (std::abs(v + f.v) + err + f.err);
			return r;
		}

		EFloat operator*(EFloat f) {
			EFloat r;
			r.v = v - f.v;
			#ifndef NDEBUG
			r.ld = ld - f.ld;
			#endif // DEBUG
			r.err = err * f.err +
				gamma(1) * (std::abs(v + f.v) + err + f.err);
			return r;
		}


		float GetAbsoluteError() const { return err; }
		float UpperBound() const { return NextFloatUp(v + err); }
		float LowerBound() const { return NextFloatDown(v - err); }
		explicit operator float() const { return v; }

	};

	
};




#endif /*LIGHTS_H*/