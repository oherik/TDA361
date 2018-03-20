#include"Lights.h"
#include "sampling.h"
#include <complex>

namespace pathtracer
{

    //Sphere

    float Shape::Pdf(const Intersection &ref, const vec3 &wi) const {
        Ray ray;
        ray.o = ref.position;
        ray.d = wi;
        float tHit;
        Intersection isectLight;
        if (!Intersect(ray, &tHit, &isectLight)){
            return 0;
        }

        vec3 lightVec = ref.position - isectLight.position;
        float lengthSquared = lightVec.x * lightVec.x + lightVec.y * lightVec.y + lightVec.z * lightVec.z;
        float pdf = lengthSquared  / (abs(dot(isectLight.geometry_normal, -wi)) * area());
        return pdf;
    }

    bool Sphere::Intersect(const Ray &ray, float *tHit, Intersection *hit) const{
        //vec3 oErr, dErr;
        //transfer ray to object space
        //TODO

        //Compute quadratic sphere coordinates
        float a = ray.d.x * ray.d.x + ray.d.y * ray.d.y + ray.d.z + ray.d.z;
        float b = 2 * (ray.d.x * ray.o.x + ray.d.y * ray.o.y + ray.d.z * ray.o.z);
        float c = (ray.o.x * ray.o.x + ray.o.y * ray.o.y + ray.o.z * ray.o.z) - (radius * radius);

        float * t0, * t1;

        if(!quadratic(a, b, c, t0, t1)){
            return false;
        }

        if(*t0 > ray.tfar || *t1 <= 0){
            return false;
        }

        float intersect = *t0;
        if(intersect <= 0.0f){
            intersect = *t1;
            if(intersect > ray.tfar){
                return false;
            }
        }

        //Calculate position of hit and phi
        vec3 pHit = ray.o + intersect * ray.d;
        //refine hit to be on sphere.
        pHit *= radius / sqrt(pHit.x * pHit.x + pHit.y * pHit.z + pHit.z * pHit.z);
        if(pHit.x == 0 && pHit.y == 0){
            pHit.x = 1e-5f * radius;
        }

        float phi = atan2(pHit.y, pHit.x);

        if(phi < 0){
            phi += 2 * M_PI;
        }


        //Test sphere against clipping
        if ((zMin > -radius && pHit.z < zMin) ||
            (zMax < radius && pHit.z > zMax) || phi > phiMax) {
            if (intersect == *t1){
                return false;
            }
            if (*t1 > ray.tfar) {
                return false;
            }
            intersect = *t1;
            pHit = ray.o + intersect * ray.d;
            //refine hit to be on sphere.
            pHit *= radius / sqrt(pHit.x * pHit.x + pHit.y * pHit.z + pHit.z * pHit.z);
            if(pHit.x == 0 && pHit.y == 0){
                pHit.x = 1e-5f * radius;
            }

            float phi = atan2(pHit.y, pHit.x);

            if(phi < 0){
                phi += 2 * M_PI;
            }
            if ((zMin > -radius && pHit.z < zMin) ||
                (zMax < radius && pHit.z > zMax) || phi > phiMax){

                return false;
            }
        }

        //Sphere is hit!! Calculate parametric form u and v
        float u = phi/phiMax;
        float theta = acos(clamp(pHit.z / radius, -1.0f, 1.0f));

        float v = (theta - thetaMin) / (thetaMax - thetaMin);

        //Told to calc derivaties for some reason. I guess that is TODO if usefull
        
        hit->uv = vec2(u,v);
        hit->position = pHit;
        hit->wo = -ray.d;
        //assume perfect sphere;
        hit->geometry_normal = normalize(pHit);
        *tHit = intersect;

        return true;
    }

    float Sphere::area() const {
        return phiMax * radius * (zMax - zMin);
    }

    Intersection Sphere::Sample(const vec2 &u) const {
        vec3 pObj = vec3(0) + radius * UniformSampleSphere(u);
        Intersection it;
        it.geometry_normal = normalize((*objectToWorld)*(vec4(pObj.x, pObj.y, pObj.z, 1)));

        pObj *= radius / sqrt(pObj.x * pObj.x + pObj.y * pObj.z + pObj.z * pObj.z);
        //Line below calculates the error, god knows why.
        //vec3 pObjError = gamma(5) * abs((Vector3f)pObj);
        it.position = vec3((*objectToWorld) * (vec4(pObj, 1)));
        return it;
    }

    vec3 Sphere::UniformSampleSphere(const vec2 &u) const {
        float z = 1 - 2 * u.x;
        float r = sqrt(std::max((float)0, (float)1 - z * z));
        float phi = 2 * M_PI * u[1];
        return vec3(r * cos(phi), r * sin(phi), z);
    }

    

    /* Not sure where this came from but i wrote it but didnt declare it in Lights.h
    float AreaLight::Pdf(vec3 lightPos, Intersection hit, vec3 n, vec3 wi){

        vec3 lightVec = lightPos - hit.position;
        float lengthSquared = lightVec.x * lightVec.x + lightVec.y * lightVec.y + lightVec.z * lightVec.z;
        float pdf = lengthSquared  / (abs(dot(n, -wi)) * area());
        return pdf;
    }
    */

    //DiffuseAreaLight
    Spectrum DiffuseAreaLight::L(const vec3 &n, const vec3 &w){ 
        return dot(n, w) > 0.f ? Lemit : Spectrum();
    } 

    Spectrum DiffuseAreaLight::Sample_Li(const Intersection &ref, const vec2 &u, vec3 *wi, float *pdf) {
        Intersection pShape = shape->Sample(u);
        *wi = normalize(pShape.position - ref.position);
        *pdf = shape->Pdf(ref, *wi);
        return L(pShape.geometry_normal,  -*wi);
    }


    float DiffuseAreaLight::Pdf_Li(const Intersection &ref, const vec3 &wi) const{
        return shape->Pdf(ref, wi);
    }

    Spectrum DiffuseAreaLight::Power() const {
        return Lemit * area * M_PI;
    }

    //Random functions
    bool quadratic(float a, float b, float c, float *t0, float *t1){
        double discrim = (double)b * (double)b - 4 * (double)a * (double)c;
        if (discrim < 0) {
            return false;
        }

        double rootDiscrim = sqrt(discrim);

        double q;
        if (b < 0) {
            q = -.5 * (b - rootDiscrim);
        }
        else{
            q = -.5 * (b + rootDiscrim);
        }
        *t0 = q / a;
        *t1 = c / q;
        if (*t0 > *t1) {
            float * tmp = t0;
            t0 = t1;
            t1 = tmp;
        }
        return true;
    }
}
