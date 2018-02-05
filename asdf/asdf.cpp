#include <complex>
#include <glm/glm.hpp>
#include "Model.h"
#include <embree2/rtcore.h>
#include <embree2/rtcore_ray.h>
#include <algorithm>
#include <iostream>
#include <memory>
#include <map>
#include<stdio.h>
#include<math.h>
#include<time.h>
#include<stdlib.h>
#include<random>
#include<omp.h>

using namespace glm;

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

float rgToN(float r, float g) {
    return g * (1.0f - r) / (1.0f - r) + (1.0f - g)  * (1.0f - sqrt(r)) / (1.0f - sqrt(r));
}
float rgToK(float r, float g) {
    return sqrt(1.0 / (1.0 - r) * r * (pow((rgToN(r, g) + 1), 2.0f) - pow(rgToN(r, g) - 1, 2.0f)));
}


int main() { 


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

    float rValues[] = { 0.0f, 0.4f,0.4f };
    float gValues[] = { 0.0f, 1.0f,0.0f };


    std::complex<float> c1(n_m[0], k_m[0]);
    std::complex<float> c2(n_m[1], k_m[1]);
    std::complex<float> c3(n_m[2], k_m[2]);

    std::complex<float> c1_c = conj(c1);
    std::complex<float> c2_c = conj(c2);
    std::complex<float> c3_c = conj(c3);

    float r0_1 = ((c1 - 1.0f)*(c1_c - 1.0f) / ((c1 + 1.0f)*(c1_c + 1.0f))).real();
    float r0_2 = ((c2 - 1.0f)*(c2_c - 1.0f) / ((c2 + 1.0f)*(c2_c + 1.0f))).real();
    float r0_3 = ((c3 - 1.0f)*(c3_c - 1.0f) / ((c3 + 1.0f)*(c3_c + 1.0f))).real();

    // Copy
    float theta = 0.08f;
    float angle = 0.0f;
    float shininess = 203.0f;
    float EPSILON = 0.00001f;
    vec3 n = vec3(1.0f, 0.0f, 0.0f);
    vec3 wi = vec3(1.0f, 0.0f, 0.0f);
    vec3 witar = vec3(0.0f, 1.0f, 0.0f);
    vec3 wo = vec3(1.0f, 0.0f, 0.0f);
    vec3 wotar = vec3(0.0f, -1.0f, 0.0f);








    while(angle<(3.14f/2.0f)){



        // COpy
        vec3 wh = normalize(wi + wo);

        float whdotwi = max(0.0f, dot(wh, wi));
        float ndotwh = max(0.0f, dot(n, wh));
        float wodotwh = max(0.0f, dot(wo, wh));
        float ndotwi = max(0.0f, dot(n, wi));
        float ndotwo = max(0.0f, dot(n, wo));

        float cost = abs(dot(wi, n) / (n.length()*wi.length()));


        vec3 rValues = vec3(1.0f, 0.0f, 0.0f);
        vec3 gValues = vec3(0.0f, 1.0f, 0.0f);

        //printf("%f      \n", rValues[0]);

        float F_wi_1 = exactReflection(rgToN(rValues[0], gValues[0]), rgToK(rValues[0], gValues[0]), cost); //n_m[0], k_m[0], cost);// r0_1 + (1.0f - r0_1)*pow(1.0f - whdotwi, 5.0f);
        float F_wi_2 = exactReflection(rgToN(rValues[1], gValues[1]), rgToK(rValues[1], gValues[1]), cost);//r0_2 + (1.0f - r0_2)*pow(1.0f - whdotwi, 5.0f);
        float F_wi_3 = exactReflection(rgToN(rValues[2], gValues[2]), rgToK(rValues[2], gValues[2]), cost);//r0_3 + (1.0f - r0_3)*pow(1.0f - whdotwi, 5.0f);


        //printf("%f     ", rgToN(rValues[0], gValues[0]));


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

        if (den < EPSILON) return 0;

        vec3 colors = (vec3(F_wi_1, F_wi_2, F_wi_3) *D_wh*G_wiwo) / den;
        //print res
        printf("angle: %f Red: %f Green: %f Blue: %f\n", angle, colors.x, colors.y, colors.z);

        //calc new wi value
        angle += theta;

        vec3 C = normalize(cross(wi, witar));
        vec3 F = cross(C, wi);
        wi = cos(theta) * wi + sin(theta) * F;

        //calc new wo value
        C = normalize(cross(wo, wotar));
        F = cross(C, wo);
        wi = cos(theta) * wo + sin(theta) * F;
    }
    //return BlinnPhong::reflection_brdf(wi, wo, n) * color; 
};


