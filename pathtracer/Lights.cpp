#include"Lights.h"
#include "sampling.h"
#include <complex>
#include "Pathtracer.h"
#include <limits>


namespace pathtracer
{

	//Sphere
	float pathtracer::Shape::Pdf(Intersection &ref, vec3 &wi) {
		Ray ray;
		ray.o = ref.position;
		ray.d = wi;
		float tHit;
		Intersection isectLight;
		if (!Intersect(ray, &tHit, &isectLight)) {
			return 0;
		}

		vec3 lightVec = ref.position - isectLight.position;
		float lengthSquared = lightVec.x * lightVec.x + lightVec.y * lightVec.y + lightVec.z * lightVec.z;
		float pdf = lengthSquared / (abs(dot(isectLight.geometry_normal, -wi)) * area());
		return pdf;
	}


	bool pathtracer::Sphere::Intersect(Ray & _ray, float * tHit, Intersection * hit) {
		float phi;
		vec3 pHit;

		//transform ray into objectspace

		Ray ray((*(this->worldToObject) * vec4(_ray.o, 1)), normalize(*(this->worldToObject) * vec4(_ray.d, 1)), 0.0f, _ray.tfar);

		//Compute quadric sphere coordinates

		vec3 d = ray.d;
		vec3 o = ray.o;

		float a = d.x * d.x + d.y * d.y + d.z * d.z;
		float b = 2 * (d.x * o.x + d.y * o.y + d.z * o.z);
		float c = o.x * o.x + o.y * o.y + o.z * o.z - radius * radius;

		//solve quadratic for t values
		float t0, t1;
		if (!quadratic(a, b, c, &t0, &t1)) {
			return false;
		}

		if (t1 <= 0) {
			return false;
		}
		float tShapeHit = t0;

		if (tShapeHit <= 0) {
			tShapeHit = t1;
		}

		//Compute sphere hit pos and phi
		pHit = ray.o + ray.d * tShapeHit;

		//Refine sphere intersection point
		pHit *= radius / distance(pHit, vec3(0, 0, 0));

		if (pHit.x == 0 && pHit.y == 0) {
			pHit.x = 1e-5f * radius;
		}
		phi = std::atan2(pHit.y, pHit.x);
		if (phi < 0) {
			phi += 2 * M_PI;
		}


		//Find parametric representation of sphere hit
		float u = phi / phiMax;
		float theta = std::acos(Clamp(pHit.z / radius, -1, 1));
		float v = (theta - thetaMin) / (thetaMax - thetaMin);

		//initialize intersection from parametric information
		//Init hit
		hit->uv = vec2(u, v);
		hit->position = vec3((*(this->objectToWorld)) * vec4(pHit, 1.0f));
		hit->wo = vec3(*this->objectToWorld * vec4(-ray.d, 1.0f));
		//assume perfect sphere;
		hit->geometry_normal = normalize(vec3(*this->objectToWorld * vec4(pHit, 1.0f)));
		//update tHit
		*tHit = tShapeHit;		//update tHit for quadric intersection

		return true;
	}
	float pathtracer::Sphere::area() {
		return phiMax * radius * (zMax - zMin);
	}

	Intersection pathtracer::Sphere::Sample(vec2 & u) {
		vec3 pObj = vec3(0) + radius * UniformSampleSphere(u);
		Intersection it;
		it.geometry_normal = normalize(*(this->objectToWorld) * (vec4(pObj.x, pObj.y, pObj.z, 1)));
		pObj *= radius / sqrt(pObj.x * pObj.x + pObj.y * pObj.y + pObj.z * pObj.z);

		//Line below calculates the error, god knows why.
		//vec3 pObjError = gamma(5) * abs((Vector3f)pObj);

		it.position = vec3(*(this->objectToWorld) * (vec4(pObj, 1)));
		return it;
	}


	vec3 pathtracer::Sphere::UniformSampleSphere(vec2 & u) {
		float z = 1 - 2 * u.x;
		float r = sqrt(std::max((float)0, (float)1 - z * z));
		float phi = 2 * M_PI * u[1];
		return vec3(r * cos(phi), r * sin(phi), z);
	}

	//Cube

	bool checkSide(float e, float f, float h, float &tmin, float &tmax) {
		if (abs(f) > EPSILON) {
			float t1 = (e + h) / f;
			float t2 = (e - h) / f;

			if (t1 > t2) {
				swap(t1, t2);
			}
			if (t1 > tmin) {
				tmin = t1;
			}
			if (t2 < tmax) {
				tmax = t2;
			}
			if (tmax < 0) {
				return false;
			}
		}
		else if (-e - h > 0 || -e + h < 0) {
			return false;
		}
		return true;
	}

	bool Cube::Intersect(Ray & ray, float * tHit, Intersection * hit) {

		//transform ray into objectspace
		//Ray ray((*(this->worldToObject) * vec4(_ray.o, 1)), normalize(*(this->worldToObject) * vec4(_ray.d, 1)), 0.0f, _ray.tfar);

		vec3 hw = vec3(0.0f, 0.0f, side / 2.0f);
		vec3 hv = vec3(0.0f, side / 2.0f, 0.0f);
		vec3 hu = vec3(side / 2.0f, 0.0f, 0.0f);

		vec3 aw = normalize(hw);
		vec3 av = normalize(hv);
		vec3 au = normalize(hu);
		float tmin = -std::numeric_limits<float>::max();
		float tmax = std::numeric_limits<float>::max();

		vec3 p = vec3(0.0f) - ray.o;

		float h = side / 2.0f;

		//for each side

		//au
		float e = dot(au, p);
		float f =  dot(au, ray.d);
		bool ret = !checkSide(e, f, h, tmin, tmax);
		if (ret) {
			return false;
		}

		//av
		e = dot(av, p);
		f = dot(av, ray.d);
		ret = !checkSide(e, f, h, tmin, tmax);
		if (ret) {
			return false;
		}

		//aw
		e = dot(aw, p);
		f = dot(aw, ray.d);

		ret = !checkSide(e, f, h, tmin, tmax);
		if (ret) {
			return false;
		}

		//get hit position
		if (tmin > 0) {
			hit->position = ray.o + tmin * ray.d;
			*tHit = tmin;
		}
		else {
			hit->position = ray.o + tmin * ray.d;
			*tHit = tmax;
		}

		//figure out normal
		if (abs(hit->position.x) > abs(hit->position.y) && abs(hit->position.x) > abs(hit->position.z)) {
			hit->geometry_normal = normalize(vec3(hit->position.x, 0.0f, 0.0f));
		}
		else if (abs(hit->position.y) > abs(hit->position.z)) {
			hit->geometry_normal = normalize(vec3(0.0f, hit->position.y, 0.0f));
		}
		else {
			hit->geometry_normal = normalize(vec3(0.0f, 0.0f, hit->position.z));
		}

		//translate back to worldspace
		//hit->geometry_normal = normalize(vec3((*objectToWorld) * vec4(hit->geometry_normal, 1.0f)));
		//hit->position = vec3((*objectToWorld) * vec4(hit->position, 1.0f));


		return true;
	}

	float Cube::area(){
		return side * side * 6;
	}

	Intersection Cube::Sample(vec2 & u){

		int side = randf() * 6;
		float rand1 = (randf() - 0.5f) * side;
		float rand2 = (randf() - 0.5f) * side;


		Intersection it = Intersection();

		switch (side) {
		case 0: 
			it.position = vec3(rand1, rand2, side);
			it.geometry_normal = vec3(0.0f, 0.0f, 1.0f);
			break;
		case 1: 
			it.position = vec3(rand1, rand2, -side);
			it.geometry_normal = vec3(0.0f, 0.0f, -1.0f);

			break;
		case 2: 
			it.position = vec3(side, rand1, rand2);
			it.geometry_normal = vec3(1.0f, 0.0f, 0.0f);

			break;
		case 3: 
			it.position = vec3(-side, rand1, rand2);
			it.geometry_normal = vec3(-1.0f, 0.0f, 0.0f);

			break;
		case 4: 
			it.position = vec3(rand1, side, rand2);
			it.geometry_normal = vec3(0.0f, 1.0f, 0.0f);

			break;
		case 5: 
			it.position = vec3(rand1, -side, rand2);
			it.geometry_normal = vec3(0.0f, -1.0f, 0.0f);

			break;
		}

		it.position = vec3(*(this->objectToWorld) * vec4(it.position, 1.0f));
		it.geometry_normal = vec3(*(this->objectToWorld) * vec4(it.geometry_normal, 1.0f));

		return Intersection();
	}

	//obsolete for cube.
	vec3 Cube::UniformSampleSphere(vec2 & u)
	{
		return vec3();
	}

	//Disc

	bool pathtracer::Disk::Intersect(Ray & _ray, float * tHit, Intersection * hit){
		//Transform ray to object space
		Ray ray;
		ray.o = (*(this->worldToObject) * vec4(_ray.o, 1));
		ray.d = normalize(*(this->worldToObject) * vec4(_ray.d, 1));
		ray.tfar = _ray.tfar;

		//Compute plane intersection for disk
		float tShapeHit = (this->height - ray.o.z) / ray.d.z;
		if (tShapeHit <= 0.0f /*|| tShapeHit >= ray.tfar*/) {
			return false;
		}

		if (ray.d.z == 0) {
			return false;
		}

		vec3 pHit = ray.o + ray.d * tShapeHit;

		//See if hit point is inside disk radii and phimax
		float dist2 = pHit.x * pHit.x + pHit.y * pHit.y;
		if (dist2 > radius*radius || dist2 < innerRadius * innerRadius) {
			return false;
		}
		//test disc phi value against phimax, i don't care about that. TODO if you want

		//Find parametric representation
		float phi = std::atan2(pHit.y, pHit.x);
		float u = phi / phiMax;
		float rHit = std::sqrt(dist2);
		float oneMinusV = ((rHit - innerRadius) / (radius - innerRadius));
		float v = 1 - oneMinusV;
		/*
		vec3 dpdu (-phiMax * pHit.y, phiMax * pHit.x, 0);
		vec3 dpdv = vec3(pHit.x, pHit.y, 0) * (innerRadius - radius) / rHit;
		vec3 dndu(0, 0, 0), dndv(0, 0, 0);
		*/

		//Refine disk intersection point
		pHit *= radius / distance(pHit, vec3(0, 0, 0));
		//vec3 pError = gamma(5) * abs((vec3)pHit);
		pHit.z = height;
		//compute error bounds for disk intersection
		vec3 pError(0, 0, 0);

		//Init hit
		hit->uv = vec2(u, v);
		hit->position = vec3((*(this->objectToWorld)) * vec4(pHit, 1.0f));
		hit->wo = vec3(*this->objectToWorld * vec4(-ray.d, 1.0f));
		//assume perfect sphere;
		hit->geometry_normal = normalize(vec3(*this->objectToWorld * vec4(pHit, 1.0f)));
		//update tHit
		*tHit = tShapeHit;

		return true;
	}

	float pathtracer::Disk::area(){
		return phiMax * 0.5 * (radius * radius - innerRadius * innerRadius);
	}

	vec2 pathtracer::Disk::UniformSampleSphere(vec2 & u){
		float r = std::sqrt(u.x);
		float theta = 2 * M_PI * u.y;
		return vec2(r*std::cos(theta), r * std::sin(theta));
	}

	   
	Intersection pathtracer::Disk::Sample(vec2 & u){
		vec2 pd = UniformSampleSphere(u);
		vec3 pObj(pd.x * radius, pd.y * radius, height);
		Intersection it;
		it.geometry_normal = normalize(vec3((*objectToWorld) * vec4(0.0f, 0.0f, 1.0f, 1.0f)));
		it.position = vec3((*objectToWorld) * vec4(pObj, 1.0f));
		return Intersection();
	}

	Spectrum pathtracer::AreaLight::Le(vec3 &n, vec3 &w){
		return L(n, w);
	}

	float pathtracer::AreaLight::pdf(vec3 lightPos, Intersection hit, vec3 n, vec3 wi){
		vec3 lightVec = lightPos - hit.position;
		float lengthSquared = lightVec.x * lightVec.x + lightVec.y * lightVec.y + lightVec.z * lightVec.z;
		float pdf = lengthSquared / (abs(dot(n, -wi)) * getArea());
		return pdf;
	}
 
	float pathtracer::AreaLight::getArea(){
		return 0.0f;
	}

	Spectrum pathtracer::AreaLight::getLEmit(){
		return this->lEmit;
	}

	float pathtracer::AreaLight::Pdf_Li(Intersection & it, vec3 & wi){
		return shape->Pdf(it, wi);
	}

	
    //DiffuseAreaLight
    Spectrum DiffuseAreaLight::L(vec3 &n, vec3 &w){ 
        //return dot(n, w) > 0.f ? lEmit : Spectrum(0.0f);
		return lEmit;
    } 

    Spectrum DiffuseAreaLight::Sample_Li(Intersection &ref, Intersection *lightHit, vec2 &u, vec3 *wi, float *pdf) {
        *lightHit = AreaLight::shape->Sample(u);
		*wi = normalize(lightHit->position - ref.position);
        *pdf = this->shape->Pdf(ref, *wi);
		Spectrum lEmit = L(lightHit->geometry_normal,  -*wi);
		return lEmit;
    }

	Spectrum pathtracer::DiffuseAreaLight::Power(){
		return lEmit * this->area * M_PI;
	}

	float pathtracer::DiffuseAreaLight::getArea(){
		return this->area;
	}

	

    //Random functions
    bool quadratic(float a, float b, float c, float *t0, float *t1){
		double discrim = (double)b * (double)b - 4 * (double)a * (double)c;
		if (discrim < 0) return false;
		double rootDiscrim = std::sqrt(discrim);

		double q;
		if (b < 0) q = -.5 * (b - rootDiscrim);
		else q = -.5 * (b + rootDiscrim);
		*t0 = q / a;
		*t1 = c / q;
		if (*t0 > *t1) std::swap(*t0, *t1);
		return true;
    }

	static constexpr Float MachineEpsilon =
		std::numeric_limits<Float>::epsilon() * 0.5;

	inline constexpr float gamma(int n) {
		return (n * MachineEpsilon) / (1 - n * MachineEpsilon);
	}

	
	
	
}
