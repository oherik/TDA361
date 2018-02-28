#pragma once
#include <glm/glm.hpp>
#include "Pathtracer.h"
#include "sampling.h"
#include "spectrum.h"
#include "embree.h"

using namespace glm; 

namespace pathtracer
{

    
    enum class LightFlags : int {
            DeltaPosition = 1, DeltaDirection = 2, Area = 4, Infinite = 8
    };
        
	class Light{
        public: 
            Light(int flags, const mat4 &LightToWorld, int nSamples = 1);
            const int flags;
            const int nSamples;
            virtual Spectrum Sample_Li(const Intersection &ref, const float &u, const float &v, vec3 *wi, Float *pdf) const = 0;
            virtual Spectrum Power() const = 0;

        protected:
            const mat4 LightToWorld, WorldToLight;
    };


    class AreaLight : public Light {
        public:
           virtual Spectrum L(vec3 &coordinate, vec3 &wo, const vec3 &w) const = 0; 
    };


    class DiffuseAreaLight : public AreaLight {
        public:
            const Spectrum Lemit;
            const Float area;
        protected:
    };
}
