#pragma once
#include <glm/glm.hpp>
#include "Pathtracer.h"
#include "sampling.h"
#include "spectrum.h"
#include "embree.h"
#include <iostream>
#include <memory>

using namespace glm; 

namespace pathtracer{

    class Shape{
        public:
            Shape(const mat4 *worldToObject, const mat4 *objectToWorld);
            virtual bool Intersect(const Ray &ray, float *tHit, Intersection *hit) const = 0;
            virtual float area() const = 0;
            virtual Intersection Sample(const vec2 &u) const = 0;
            virtual vec3 UniformSampleSphere(const vec2 &u) const = 0;
            virtual float Pdf(const Intersection &ref, const vec3 &wi) const = 0;
            const mat4 *worldToObject, *objectToWorld;
    };

    class Sphere : public Shape {
        private:
            const float radius, zMin, zMax;
            const float thetaMin, thetaMax, phiMax;
        public:
            Sphere(const mat4 *_objectToWorld, const mat4 *_worldToObject, float _radius, float _zMin, float _zMax, float _phiMax) : 
                Shape(objectToWorld, worldToObject),
                radius(_radius), 
                zMin(clamp(min(_zMin,_zMax), -radius, radius)), 
                zMax(clamp(max(_zMin, _zMax), -radius, radius)),
                thetaMin(acos(clamp(zMin / radius, -1.0f, 1.0f))),
                thetaMax(acos(clamp(zMax / radius, -1.0f, 1.0f))),
                phiMax(radians(clamp(_phiMax, 0.0f, 360.0f))) {}

            virtual bool Intersect(const Ray &ray, float *tHit, Intersection *hit) const override;
            virtual float area() const override;
            virtual Intersection Sample(const vec2 &u) const override;
            virtual vec3 UniformSampleSphere(const vec2 &u) const override;
    };


    
    enum class LightFlags : int {
            DeltaPosition = 1, DeltaDirection = 2, Area = 4, Infinite = 8
    };
        
	class Light{
        public: 
            Light();
            Light(int flags, const mat4 &LightToWorld, int nSamples = 1)
            : flags(flags), nSamples(max(1, nSamples)), LightToWorld(LightToWorld),
            WorldToLight(inverse(LightToWorld)) {
            };
            const int flags;
            const int nSamples;
            virtual Spectrum Sample_Li(const Intersection &ref, const vec2 &u, vec3 *wi, float *pdf) {};
            virtual Spectrum Power() const = 0;

        protected:
            const mat4 LightToWorld, WorldToLight;
    };


    class AreaLight : public Light {
        public:
            AreaLight();
            AreaLight(mat4 LightToWorld, int nSamples = 1)
                : Light(4, LightToWorld, nSamples){};
            virtual Spectrum L(const vec3 &coordinate, const vec3 &w) {}; 
            virtual float area() const {};
            float pdf(vec3 lightPos, Intersection hit, vec3 n, vec3 wi);
    };


    class DiffuseAreaLight : public AreaLight {
        public:
            DiffuseAreaLight();
            const Spectrum Lemit;
            const float area;
            const Shape * shape;
            DiffuseAreaLight(const mat4 &LightToWorld,
            const Spectrum &Lemit,
            int nSamples, const Shape *shape)
            : AreaLight(LightToWorld, nSamples), Lemit(Lemit),
            shape(shape), area(shape->area()) {};
     
            float Pdf_Li(const Intersection &ref, const vec3 &wi);
            Spectrum L(const vec3 &n, const vec3 &w) override;
            Spectrum Sample_Li(const Intersection &ref, const vec2 &u, vec3 *wi, float *pdf) override;
            Spectrum Power();
    };

    bool quadratic(float a, float b, float c, float *t0, float *t1);
};
