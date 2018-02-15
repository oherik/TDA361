#include <complex>
#include <glm/glm.hpp>
#include "Model.h"
#include <embree2/rtcore.h>
#include <embree2/rtcore_ray.h>
#include <algorithm>
#include <iostream>
#include <memory>
#include <map>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <random>
#include <omp.h>

using namespace glm;



//Values
//Koppar
float m_n []= { 0.294f, 1.0697f, 1.2404f };
float m_k [] = { 3.2456f, 2.6866f, 2.3929f };

//Guld
//float n_m []= { 0.15557f, 0.42415f,1.3831f };
//float k_m [] = { 3.6024f, 2.4721f,1.9155f };

//Aluminium
//float n_m[] = { 1.5580f, 1.0152f, 0.63324f };
//float k_m[] = { 7.7124f, 6.6273f, 5.4544f };

//Platina
//float n_m[] = { 0.040000f, 0.059582f, 0.052225f };
//float k_m[] = { 2.6484f,  3.5974f, 4.4094f };


float shininess = 25000.0f;
float EPSILON = 0.00001f;
float pi = 3.14159265359f;

///////////////////////////////////////////////////////////////////////////////
// Get a random float. Note that we need one "generator" per thread, or we 
// would need to lock everytime someone called randf(). 
///////////////////////////////////////////////////////////////////////////////
std::mt19937 generators[24]; // Assuming no more than 24 cores
float randf() {
    return float(generators[omp_get_thread_num()]() /
        double(generators[omp_get_thread_num()].max()));
}


float exactReflection(float n, float k, float cost) {
    return( 0.5*(

            (pow(n,2.0f) + pow(k,2.0f) - 2 * n * cost + pow(cost,2) )/
            (pow(n, 2.0f) + pow(k, 2.0f) + 2 * n * cost + pow(cost, 2)) +
        
        ((pow(n, 2.0f) + pow(k, 2.0f))*pow(cost,2.0f)-2*n*cost+1) /
        ((pow(n, 2.0f) + pow(k, 2.0f))*pow(cost, 2.0f) + 2 * n*cost + 1)
        ));
}

glm::vec3 perpendicular(const glm::vec3 &v)
{
    if (fabsf(v.x) < fabsf(v.y))
    {
        return glm::vec3(0.0f, -v.z, v.y);
    }
    return glm::vec3(-v.z, 0.0f, v.x);
}



vec3 reflection_brdf(const vec3 & wi, const vec3 & wo, const vec3 & n) { 

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

    float D_wh = (shininess + 2.0f) / (2.0f * pi) * pow(ndotwh, shininess);
    float G_wiwo = min(1.0f, min(2.0f * ndotwh*ndotwo / wodotwh, 2.0f * ndotwh*ndotwi / wodotwh));

    float den = (4.0f * ndotwo*ndotwi);

    if (den < EPSILON) return vec3(0.0f);
    
    printf("%f %f %f\n", F_wi_1, F_wi_2, F_wi_3);
    return vec3(F_wi_1, F_wi_2, F_wi_3) * D_wh * G_wiwo / den;
};





vec3 sample_wi(vec3 & wi, const vec3 & wo, const vec3 & n, float & p) {
    
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
    wi = normalize(reflect(-normalize(wo),normalize(wh)));	
    return reflection_brdf(wi, wo, n);
}

int main() { 

    float theta = 0.02f;
    float angle = 0.0f;
    vec3 n = vec3(1.0f, 0.0f, 0.0f);
    vec3 wo = vec3(1.0f, 0.0f, 0.0f);
    vec3 wotar = vec3(0.0f, -1.0f, 0.0f);
    vec3 wi;
    float pdf = 0.0f;
    vec3 brdf;

    while(angle<(3.14f/2.0f)){


        printf("angle: %f\n", degrees(angle));
        brdf = sample_wi(wi, wo, n, pdf);
        float cosineTerm = abs(dot(wi, n));

        printf("Red: %f Green: %f Blue: %f\n", brdf.x/pdf, brdf.y/pdf, brdf.z/pdf);



        //calc new wo value
        angle += theta;
        vec3 C = normalize(cross(wo, wotar));
        vec3 F = cross(C, wo);
        wo = cos(theta) * wo + sin(theta) * F;
    }
    //return BlinnPhong::reflection_brdf(wi, wo, n) * color; 
};


