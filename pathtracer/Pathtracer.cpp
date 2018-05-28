#include "Pathtracer.h"
#include <memory>
#include <iostream>
#include <map>
#include <algorithm>
#include "material.h"
#include "embree.h"
#include "sampling.h"
#include "spectrum.h"
#include "Lights.h"
#include <math.h>  
#include <glm/gtx/rotate_vector.hpp>

using namespace std; 
using namespace glm; 



namespace pathtracer
{ 

	Spectrum EstimateDirect(Intersection &it,
		vec2 &uScattering, AreaLight &light,
		vec2 &uLight, BRDF &mat, bool isSurfaceHit, bool handleMedia, bool specular); 

	float PowerHeuristic(int nf, float fPdf, int ng, float gPdf);


	///////////////////////////////////////////////////////////////////////////////
	// Global variables
	///////////////////////////////////////////////////////////////////////////////
	Settings settings;
	Environment environment; 
	Image rendered_image; 
	Image corners_image;
	PointLight point_light; 
	
    Brdf brdf;
	CustomSettings customSettings;
	DepthOfField depthOfField;
	static float epsilon = 0.1f; // for AASS
	static int maxLevel = 2; //for AASS, should be >=1

	///////////////////////////////////////////////////////////////////////////
	// Defining AreaLights here, hope you dont mind.. 
	///////////////////////////////////////////////////////////////////////////
	
	Shape * shape;
	AreaLight * areaLight;

	///////////////////////////////////////////////////////////////////////////
	// Restart rendering of image
	///////////////////////////////////////////////////////////////////////////
	void restart()
	{
		// No need to clear image, 
		rendered_image.number_of_samples = 0; 
		corners_image.number_of_samples = 0;
	}

	///////////////////////////////////////////////////////////////////////////
	// On window resize, window size is passed in, actual size of pathtraced
	// image may be smaller (if we're subsampling for speed)
	///////////////////////////////////////////////////////////////////////////
	void resize(int w, int h)
	{
		rendered_image.width = w / settings.subsampling; 
		rendered_image.height = h / settings.subsampling; 
		rendered_image.data.resize(rendered_image.width * rendered_image.height);

		corners_image.width = rendered_image.width + 1;
		corners_image.height = rendered_image.height + 1;
		corners_image.data.resize(corners_image.width * corners_image.height);

		restart(); 
	}

///////////////////////////////////////////////////////////////////////////
// Return the radiance from a certain direction wi from the environment
// map. 
///////////////////////////////////////////////////////////////////////////
Spectrum Lenvironment(const vec3 & wi) {
	const float theta = acos(std::max(-1.0f, std::min(1.0f, wi.y)));
	float phi = atan(wi.z, wi.x);
	if (phi < 0.0f) phi = phi + 2.0f * M_PI;
	vec2 lookup = vec2(phi / (2.0 * M_PI), theta / M_PI);
	return environment.multiplier * Spectrum::FromRGB(environment.map.sample(lookup.x, lookup.y), SpectrumType::Illuminant);
}


//sampleLi for pointlight
Spectrum sampleLi(Ray lightRay, Intersection hit, vec3 *wi, float *pdf) {

	if (!occluded(lightRay)) {
		vec3 lightVec = point_light.position - hit.position;
		*wi = normalize(lightVec);
		*pdf = 1.0f;

		Spectrum lightSpectrum = Spectrum::FromRGB(point_light.color, SpectrumType::Illuminant);
		float lengthSquared = lightVec.x * lightVec.x + lightVec.y * lightVec.y + lightVec.z * lightVec.z;

		return lightSpectrum / lengthSquared;
	}
	else {
		return Spectrum();
	}
}

///////////////////////////////////////////////////////////////////////////
// Calculate the radiance going from one point (r.hitPosition()) in one 
// direction (-r.d), through path tracing.  
///////////////////////////////////////////////////////////////////////////
Spectrum Li(Ray & primary_ray) {
	Spectrum spectrumSample;
	const vec3 fullLight = vec3(1.0f, 1.0f, 1.0f);
	Spectrum throughput = Spectrum(1.f);
	Ray current_ray = primary_ray;
	vec3 last_position = vec3(0.0f);


	//Task 5: bounce it up
	for (int i = 0; i < settings.max_bounces; i++) {

		float thit;
		Intersection lightIntersection;
		bool hitLight = shape->Intersect(current_ray, &thit, &lightIntersection);

		// Get the intersection information from the ray
		Intersection hit = getIntersection(current_ray);

		//if light is in the path of the ray check which one is closer
		if (hitLight) {
			//If the ray hit light first return light.
			float lightray = distance(current_ray.o, lightIntersection.position);
			float hitray = distance(current_ray.o, hit.position);
			float lightHitDistance = lightray - hitray;
			if (lightHitDistance <= 0.0f) {
				Spectrum areaLH = areaLight->getLEmit();
				return spectrumSample + areaLH * throughput;
			}
		}


			//Bump it a bit
			labhelper::Texture bumpmap = hit.material->m_bumpmap_texture;
			float u = hit.texture_coordinates.x;
			float v = hit.texture_coordinates.y;
			//vec3 originalShadingNormal = hit.shading_normal;
			if (bumpmap.valid && pathtracer::customSettings.bumpmap) {
				int i = round(u * (float)bumpmap.width);
				int j = round(v * (float)bumpmap.height);




				unsigned bytePerPixel = 3;
				unsigned char* pixelOffset = bumpmap.data + (i + bumpmap.height * j) * bytePerPixel;
				float r = (float)pixelOffset[0] / 255.f;
				float g = (float)pixelOffset[1] / 255.f;
				float b = (float)pixelOffset[2] / 255.f;
				vec3 offset = 2.f*vec3(r, g, b) - vec3(1.f);


				vec3 tangent = normalize(perpendicular(hit.shading_normal));
				vec3 bitangent = normalize(cross(tangent, hit.shading_normal));
				mat3 tbn = mat3(tangent, bitangent, hit.shading_normal);
				vec3 new_normal = normalize(tbn * offset);


				//Clamp it so it doesn't point underneath the surface
				if (dot(hit.geometry_normal, hit.shading_normal) < 0) {
					return Spectrum::FromRGB(vec3(1.f, 0.f, 0.f), SpectrumType::Illuminant);
				}

				//Clamp it to 90 degrees to avoid nasty artifacts.
				if (dot(new_normal, hit.wo) <= 0) {
					//return Spectrum::FromRGB(vec3(1.f, 1.f, 1.f), SpectrumType::Illuminant);
					vec3 rotation_vector = normalize(cross(cross(hit.wo, new_normal), hit.wo));
					hit.shading_normal = normalize(0.02f*hit.wo + rotation_vector); // 89 degrees
				}
				else {
					hit.shading_normal = new_normal;
				}


			}



			// Create a material tree

			//Use color map if available
			labhelper::Texture color_texture = hit.material->m_color_texture;

			vec3 diffuseColor;
			if (color_texture.valid && pathtracer::customSettings.diffusemap) {
				int i = round(u * (float)color_texture.width);
				int j = round(v * (float)color_texture.height);
				unsigned bytePerPixel = 4;
				unsigned char* pixelOffset = color_texture.data + (i + color_texture.height * j) * bytePerPixel;
				float r = (float)pixelOffset[0] / 255.f;
				float g = (float)pixelOffset[1] / 255.f;
				float b = (float)pixelOffset[2] / 255.f;
				diffuseColor = vec3(r, g, b);
			}
			else {
				diffuseColor = hit.material->m_color;
			}
			Diffuse diffuse(diffuseColor);

			/*
			Transparent transparent(hit.material->m_transparency, hit.material->m_color);

			//Calculate opacity

			float dist = (length(hit.position - last_position));
			float a_c = -log(hit.material->m_transparency);
			float exp = -a_c*dist;
			float T = pow(2.72, exp);
			float a = 0.0f;
			if (dot(hit.shading_normal, current_ray.d) < 0.0f){
				a = 1.0f;
			}
			else {
				if (hit.material->m_transparency < EPSILON){
					a = 0.0f;
				}
				else {
					a = T;
				}
			}


			TransparencyBlend transparency_blend(a, &transparent, hit.material->m_color);
			*/



			float wavelengths[3] = { hit.material->m_RGB_wavelengths.x, hit.material->m_RGB_wavelengths.y, hit.material->m_RGB_wavelengths.z };

			/*
						labhelper::Texture texture = hit.material->m_color_texture;// m_bumpmap_texture;



						int _i = round(u * (float)texture.width);
						int _j = round(v * (float)texture.height);


						unsigned bytePerPixel = 4;
						unsigned char* pixelOffset =texture.data  + (_i + texture.height * _j) * bytePerPixel; //static_cast<int>(round(asd * bytePerPixel));
						unsigned char r = pixelOffset[0];
						unsigned char g = pixelOffset[1];
						unsigned char b = pixelOffset[2];
						//printf("%i, %i, \n", _i, _j);
						if (bytePerPixel == 4){ // alpha channel
							unsigned char a = (bytePerPixel >= 4) ? pixelOffset[3] : 0xff;
						}

						vec3 coords = vec3(hit.texture_coordinates.x, hit.texture_coordinates.y, 0.f);
						vec3 textureSample = vec3((float)r/255.f, (float)g / 255.f, (float)b / 255.f);
						//return Spectrum::FromRGB(textureSample, SpectrumType::Illuminant);

						*/



			CustomDefined dielectric(hit.material->m_shininess, hit.material->m_fresnel, &diffuse);
			CustomDefinedMetal metal(hit.material->m_color, hit.material->m_n, hit.material->m_k, &wavelengths[0], hit.material->m_shininess,
				hit.material->m_fresnel);
			LinearBlend metal_blend(hit.material->m_metalness, &metal, &dielectric);
			LinearBlend reflectivity_blend(hit.material->m_reflectivity, &metal_blend, &diffuse);

			//LinearBlend transparency_blend_final(hit.material->m_transparency, &transparency_blend, &reflectivity_blend);


			BRDF & mat = reflectivity_blend;

			// Update last hit position for the next distance calculation
			last_position = vec3(hit.position);


			// Calculate Direct Illumination from light.

			//Sample random position from Area light
			vec2 lightPos = vec2(randf(), randf());

			spectrumSample += EstimateDirect(hit, vec2(randf(), randf()), *areaLight, lightPos, mat, true, false, false);

			/*
			vec3 lightWi;
			float lightPdf;
			Intersection * lightHit = new Intersection();
			Spectrum areaSample = areaLight->Sample_Li(hit, lightHit, lightPos, &lightWi, &lightPdf);
			

			const float distance_to_light = length(lightHit->position - hit.position);

			const float falloff_factor = 1.0f / (distance_to_light*distance_to_light);
			Ray lightRay(hit.position + EPSILON * hit.geometry_normal, normalize(lightHit->position - hit.position), 0.0f, distance_to_light);
			

			if (!occluded(lightRay)) {
				//vec3 wi = normalize(point_light.position - hit.position);
				Spectrum reflectance = mat.f(lightWi, hit.wo, hit.shading_normal);
				Spectrum Li = areaSample * reflectance * falloff_factor * throughput * std::max(0.0f, dot(lightWi, hit.shading_normal))*point_light.intensity_multiplier;
				spectrumSample = spectrumSample + Li;
			}
			*/

			// Emitted radiance from intersection
			spectrumSample += throughput * hit.material->m_emission;

			// Sample incoming direction
			vec3 wi;
			float pdf = 0.0f;
			Spectrum brdf = mat.sample_wi(wi, hit.wo, hit.shading_normal, pdf);

			if (pdf < EPSILON) {
				return spectrumSample;
			}

			float cosineTerm = abs(dot(wi, hit.shading_normal));
			throughput = throughput * (brdf * cosineTerm) / pdf;

			//Break if the throughput is 0
			if (throughput.IsBlack())
			{
				return spectrumSample;
			}

			//Next ray
			if (dot(hit.shading_normal, wi) > 0.0f) {
				current_ray = Ray(hit.position + EPSILON * hit.geometry_normal, wi);
			}
			else {
				current_ray = Ray(hit.position - EPSILON * hit.geometry_normal, wi);
			}

			//Intersect the new ray
			if (!intersect(current_ray)) {
				return spectrumSample + throughput * Lenvironment(current_ray.d);
			}

		}

		/*


		///////////////////////////////////////////////////////////////////
		// Get the intersection information from the ray
		///////////////////////////////////////////////////////////////////
		Intersection hit = getIntersection(current_ray);

		///////////////////////////////////////////////////////////////////
		// Create a Material tree for evaluating brdfs and calculating
		// sample directions.
		///////////////////////////////////////////////////////////////////


		//Task 3: ping-pong
		Diffuse diffuse(hit.material->m_color);
		BlinnPhong dielectric(hit.material->m_shininess, hit.material->m_fresnel, &diffuse);

		//Task 4: Metal
		BlinnPhongMetal metal(hit.material->m_color, hit.material->m_shininess, hit.material->m_fresnel);
		LinearBlend metal_blend(hit.material->m_metalness, &metal, &dielectric);
		LinearBlend reflectivity_blend(hit.material->m_reflectivity, &metal_blend, &diffuse);


		BRDF & mat = reflectivity_blend;

		//Task 2: shadows
		Ray lightRay;
		lightRay.o = hit.position + EPSILON * hit.geometry_normal;
		lightRay.d = normalize(point_light.position - lightRay.o);

		if (occluded(lightRay)){
			return vec3(0.0);
		}

		///////////////////////////////////////////////////////////////////
		// Calculate Direct Illumination from light.
		///////////////////////////////////////////////////////////////////
		{
			const float distance_to_light = length(point_light.position - hit.position);
			const float falloff_factor = 1.0f / (distance_to_light*distance_to_light);
			vec3 Li = point_light.intensity_multiplier * point_light.color * falloff_factor;
			vec3 wi = normalize(point_light.position - hit.position);
			L = mat.f(wi, hit.wo, hit.shading_normal) * Li * std::max(0.0f, dot(wi, hit.shading_normal));
		}
		// Return the final outgoing radiance for the primary ray
		*/
		return spectrumSample;
		//return L;
	
	}

	float squaredError(vec3 a, vec3 b) {		
		return pow(a.x - b.x, 2) + pow(a.y - b.y, 2) + pow(a.z - b.z, 2);
	}

	// Let's get some AASS in this app
	vec3 adaptiveSupersampling(int* numberOfRays, vec3 camera_pos, vec3 camera_dir, vec3 camera_right, vec3 camera_up, float focusDistance, vec3 lower_right_corner, float x, float y, vec3 X, vec3 Y, int level, vec3 firstQp = vec3(0.f), vec3 secondQp = vec3(0.f), vec3 thirdQp = vec3(0.f), vec3 fourthQp = vec3(0.f)) {
		vec3 centerColor, firstQuadrant, secondQuadrant, thirdQuadrant, fourthQuadrant;
		Ray primaryRay;
		(*numberOfRays)++;
		primaryRay.o = camera_pos;
		vec2 screenCoord = vec2(x / float(rendered_image.width), y / float(rendered_image.height));

		primaryRay.d = normalize(lower_right_corner + screenCoord.x * X + screenCoord.y * Y);
		float sideStep = 0.5f / level;

		//Task 1: introduce some randomness and jittering
		screenCoord.x += ((randf() - 0.5)*sideStep) / rendered_image.width;
		screenCoord.y += ((randf() - 0.5)*sideStep) / rendered_image.height;

		primaryRay.d = normalize(lower_right_corner + screenCoord.x * X + screenCoord.y * Y);

		//Add some DoF
		if (depthOfField.lensRadius > 0) {
			vec2 ret;
			concentricSampleDisk(&ret.x, &ret.y);
			vec2 lensPoint = depthOfField.lensRadius *  ret;
			float t = focusDistance / (dot(primaryRay.d, camera_dir));
			vec3 focusPoint = primaryRay.o + primaryRay.d * t;
			primaryRay.o = primaryRay.o + camera_right * lensPoint.x + camera_up * lensPoint.y;
			primaryRay.d = normalize(focusPoint - primaryRay.o);
		}

		if (intersect(primaryRay)) {
			// If it hit something, evaluate the radiance from that point
			centerColor = Li(primaryRay).ToRGB();
		}else {
			// Otherwise evaluate environment
			centerColor = Lenvironment(primaryRay.d).ToRGB();
		}

		if (level > maxLevel) {
			return centerColor;
		}
			
		//	Shoot new rays if necessary
		if (length(firstQp) > EPSILON) {
			firstQuadrant = firstQp;
		}
		else {
			Ray firstQuadrantRay;
			(*numberOfRays)++;
			firstQuadrantRay.o = camera_pos;
			float newX = x + sideStep;
			float newY = y + sideStep;
			vec2 screenCoord = vec2(newX / float(rendered_image.width), newY / float(rendered_image.height));
			firstQuadrantRay.d = normalize(lower_right_corner + screenCoord.x * X + screenCoord.y * Y);
			//Add some DoF
			if (depthOfField.lensRadius > 0) {
				vec2 ret;
				concentricSampleDisk(&ret.x, &ret.y);
				vec2 lensPoint = depthOfField.lensRadius *  ret;
				float t = focusDistance / (dot(firstQuadrantRay.d, camera_dir));
				vec3 focusPoint = firstQuadrantRay.o + firstQuadrantRay.d * t;
				firstQuadrantRay.o = firstQuadrantRay.o + camera_right * lensPoint.x + camera_up * lensPoint.y;
				firstQuadrantRay.d = normalize(focusPoint - firstQuadrantRay.o);
			}
			if (intersect(firstQuadrantRay)) {
				firstQuadrant = Li(firstQuadrantRay).ToRGB();
			}
			else {
				float thit;
				Intersection lightIntersection;
				bool hitLight = shape->Intersect(firstQuadrantRay, &thit, &lightIntersection);
				if (hitLight) {
					firstQuadrant = areaLight->getLEmit().ToRGB();
				}
				else {
					firstQuadrant = Lenvironment(firstQuadrantRay.d).ToRGB();
				}
			}
		}

		if (length(secondQp) > EPSILON) {
			secondQuadrant = secondQp;
		}
		else {
			Ray secondQuadrantRay;
			(*numberOfRays)++;
			secondQuadrantRay.o = camera_pos;
			float newX = x + sideStep;
			float newY = y - sideStep;
			vec2 screenCoord = vec2(newX / float(rendered_image.width), newY / float(rendered_image.height));
			secondQuadrantRay.d = normalize(lower_right_corner + screenCoord.x * X + screenCoord.y * Y);
			//Add some DoF
			if (depthOfField.lensRadius > 0) {
				vec2 ret;
				concentricSampleDisk(&ret.x, &ret.y);
				vec2 lensPoint = depthOfField.lensRadius *  ret;
				float t = focusDistance / (dot(secondQuadrantRay.d, camera_dir));
				vec3 focusPoint = secondQuadrantRay.o + secondQuadrantRay.d * t;
				secondQuadrantRay.o = secondQuadrantRay.o + camera_right * lensPoint.x + camera_up * lensPoint.y;
				secondQuadrantRay.d = normalize(focusPoint - secondQuadrantRay.o);
			}
			if (intersect(secondQuadrantRay)) {
				secondQuadrant = Li(secondQuadrantRay).ToRGB();
			}
			else {
				float thit;
				Intersection lightIntersection;
				bool hitLight = shape->Intersect(secondQuadrantRay, &thit, &lightIntersection);
				if (hitLight) {
					secondQuadrant = areaLight->getLEmit().ToRGB();
				}
				else {
					secondQuadrant = Lenvironment(secondQuadrantRay.d).ToRGB();
				}
			}

		}

		if (length(thirdQp) > EPSILON) {
			thirdQuadrant = thirdQp;
		}
		else {
			Ray thirdQuadrantRay;
			(*numberOfRays)++;
			thirdQuadrantRay.o = camera_pos;
			float newX = x - sideStep;
			float newY = y - sideStep;
			vec2 screenCoord = vec2(newX / float(rendered_image.width), newY / float(rendered_image.height));
			thirdQuadrantRay.d = normalize(lower_right_corner + screenCoord.x * X + screenCoord.y * Y);
			//Add some DoF
			if (depthOfField.lensRadius > 0) {
				vec2 ret;
				concentricSampleDisk(&ret.x, &ret.y);
				vec2 lensPoint = depthOfField.lensRadius *  ret;
				float t = focusDistance / (dot(thirdQuadrantRay.d, camera_dir));
				vec3 focusPoint = thirdQuadrantRay.o + thirdQuadrantRay.d * t;
				thirdQuadrantRay.o = thirdQuadrantRay.o + camera_right * lensPoint.x + camera_up * lensPoint.y;
				thirdQuadrantRay.d = normalize(focusPoint - thirdQuadrantRay.o);
			}
			if (intersect(thirdQuadrantRay)) {
				thirdQuadrant = Li(thirdQuadrantRay).ToRGB();
			}
			else {
				float thit;
				Intersection lightIntersection;
				bool hitLight = shape->Intersect(thirdQuadrantRay, &thit, &lightIntersection);
				if (hitLight) {
					thirdQuadrantRay = areaLight->getLEmit().ToRGB();
				}
				else {
					thirdQuadrantRay = Lenvironment(thirdQuadrantRay.d).ToRGB();
				}
			}
		}

		if (length(fourthQp) > EPSILON) {
			fourthQuadrant = fourthQp;
		}
		else {
			Ray fourthQuadrantRay;
			(*numberOfRays)++;
			fourthQuadrantRay.o = camera_pos;
			float newX = x - sideStep;
			float newY = y + sideStep;
			vec2 screenCoord = vec2(newX / float(rendered_image.width), newY / float(rendered_image.height));
			fourthQuadrantRay.d = normalize(lower_right_corner + screenCoord.x * X + screenCoord.y * Y);
			//Add some DoF
			if (depthOfField.lensRadius > 0) {
				vec2 ret;
				concentricSampleDisk(&ret.x, &ret.y);
				vec2 lensPoint = depthOfField.lensRadius *  ret;
				float t = focusDistance / (dot(fourthQuadrantRay.d, camera_dir));
				vec3 focusPoint = fourthQuadrantRay.o + fourthQuadrantRay.d * t;
				fourthQuadrantRay.o = fourthQuadrantRay.o + camera_right * lensPoint.x + camera_up * lensPoint.y;
				fourthQuadrantRay.d = normalize(focusPoint - fourthQuadrantRay.o);
			}
			if (intersect(fourthQuadrantRay)) {
				fourthQuadrant = Li(fourthQuadrantRay).ToRGB();
			}
			else {
				float thit;
				Intersection lightIntersection;
				bool hitLight = shape->Intersect(fourthQuadrantRay, &thit, &lightIntersection);
				if (hitLight) {
					fourthQuadrantRay = areaLight->getLEmit().ToRGB();
				}
				else {
					fourthQuadrantRay = Lenvironment(fourthQuadrantRay.d).ToRGB();
				}
			}
		}



		vec3 color = vec3(0.f);
		if (squaredError(centerColor, firstQuadrant) < epsilon) {
			color += 0.125f * (centerColor + firstQuadrant);
		} else {
			color += 0.25f * adaptiveSupersampling(numberOfRays, camera_pos, camera_dir, camera_right, camera_up, focusDistance, lower_right_corner, x + 0.5f*sideStep, y + 0.5f*sideStep, X, Y, level + 1, firstQuadrant, vec3(0.f), centerColor, vec3(0.f));
		}
		if (squaredError(centerColor, secondQuadrant) < epsilon) {
			color += 0.125f * (centerColor + secondQuadrant);
		}
		else {
			color += 0.25f * adaptiveSupersampling(numberOfRays, camera_pos, camera_dir, camera_right, camera_up, focusDistance, lower_right_corner, x + 0.5f*sideStep, y - 0.5f*sideStep, X, Y, level + 1, vec3(0.f), secondQuadrant, vec3(0.f), centerColor);
		}
		if (squaredError(centerColor, thirdQuadrant) < epsilon) {
			color += 0.125f * (centerColor + thirdQuadrant);
		}
		else {
			color += 0.25f * adaptiveSupersampling(numberOfRays, camera_pos, camera_dir, camera_right, camera_up, focusDistance, lower_right_corner,  x - 0.5f*sideStep, y - 0.5f*sideStep, X, Y, level + 1, centerColor, vec3(0.f), thirdQuadrant, vec3(0.f));
		}
		if (squaredError(centerColor, fourthQuadrant) < epsilon) {
			color += 0.125f * (centerColor + fourthQuadrant);
		}
		else {
			color += 0.25f * adaptiveSupersampling(numberOfRays, camera_pos, camera_dir, camera_right, camera_up, focusDistance, lower_right_corner, x - 0.5f*sideStep, y + 0.5f*sideStep, X, Y, level + 1, vec3(0.f), centerColor, vec3(0.f), fourthQuadrant);
		}
		
		return color;// getRayColor(camera_pos, camera_dir, camera_right, camera_up, focusDistance, lower_right_corner, float(x), float(y), X, Y);
	}

	///////////////////////////////////////////////////////////////////////////
	// Trace one path per pixel and accumulate the result in an image
	///////////////////////////////////////////////////////////////////////////
	void tracePaths(vec3 camera_pos, vec3 camera_dir, vec3 camera_up)
	{
		// Calculate where to shoot rays from the camera
		vec3 camera_right = normalize(cross(camera_dir, camera_up));
		camera_up = normalize(cross(camera_right, camera_dir));
		float camera_fov = 45.0f;
		float camera_aspectRatio = float(rendered_image.width) / float(rendered_image.height);
		// Calculate the lower left corner of a virtual screen, a vector X
		// that points from there to the lower left corner, and a vector Y
		// that points to the upper left corner
		glm::vec3 A = camera_dir * cos(camera_fov / 2.0f * (M_PI / 180.0f));
		glm::vec3 B = camera_up * sin(camera_fov / 2.0f * (M_PI / 180.0f));
		glm::vec3 C = camera_right * sin(camera_fov / 2.0f * (M_PI / 180.0f)) * camera_aspectRatio;
		glm::vec3 lower_right_corner = A - C - B;
		glm::vec3 X = 2.0f * ((A - B) - lower_right_corner);
		glm::vec3 Y = 2.0f * ((A - C) - lower_right_corner);

		// Depth of field focus, in the center
		float focusDistance;
		Ray centerRay;
		centerRay.o = camera_pos;
		centerRay.d = normalize(lower_right_corner + 0.5f * X + 0.5f * Y);

		// Intersect center ray with scene
		if (intersect(centerRay)) {
			Intersection hit = getIntersection(centerRay);
			focusDistance = glm::length(hit.position - camera_pos);
		}
		else {
			focusDistance = 1000;
		}

		// Stop here if we have as many samples as we want
		if ((int(rendered_image.number_of_samples) > settings.max_paths_per_pixel) &&
			(settings.max_paths_per_pixel != 0)) return;

		// For AASS
		if (settings.supersampling_method != 0) {	
#pragma omp parallel for
			for (int y = 0; y < corners_image.height; y++) {
				for (int x = 0; x < corners_image.width; x++) {
					vec3 color;
					Ray primaryRay;

					primaryRay.o = camera_pos;
					vec2 screenCoord = vec2(x / float(corners_image.width), y / float(corners_image.height));

					//Task 1: introduce some randomness and jittering
					screenCoord.x += (randf() - 0.5) / corners_image.width;
					screenCoord.y += (randf() - 0.5) / corners_image.height;

					primaryRay.d = normalize(lower_right_corner + screenCoord.x * X + screenCoord.y * Y);

					//Add some DoF
					if (depthOfField.lensRadius > 0) {
						vec2 ret;
						concentricSampleDisk(&ret.x, &ret.y);
						vec2 lensPoint = depthOfField.lensRadius *  ret;
						float t = focusDistance / (dot(primaryRay.d, camera_dir));
						vec3 focusPoint = primaryRay.o + primaryRay.d * t;
						primaryRay.o = primaryRay.o + camera_right * lensPoint.x + camera_up * lensPoint.y;
						primaryRay.d = normalize(focusPoint - primaryRay.o);
					}

					if (intersect(primaryRay)) {
						color = Li(primaryRay).ToRGB();
					}
					else {
						float thit;
						Intersection lightIntersection;
						bool hitLight = shape->Intersect(primaryRay, &thit, &lightIntersection);
						if (hitLight) {
							color = areaLight->getLEmit().ToRGB();
						}else{
							color = Lenvironment(primaryRay.d).ToRGB();
						}
					}

					float n = float(corners_image.number_of_samples);
					corners_image.data[y * corners_image.width + x] =
					corners_image.data[y * corners_image.width + x] * (n / (n + 1.0f)) +
					(1.0f / (n + 1.0f)) * color;
				}
			}
			corners_image.number_of_samples += 1;

		}

		// Trace one path per pixel (the omp parallel stuf magically distributes the 
		// pathtracing on all cores of your CPU).
#pragma omp parallel for
		for (int y = 0; y < rendered_image.height; y++) {
			for (int x = 0; x < rendered_image.width; x++) {
				//if (settings.supersampling_method == 0) {
				vec3 color;
				Ray primaryRay;
		

				primaryRay.o = camera_pos;
				// Create a ray that starts in the camera position and points toward
				// the current pixel on a virtual screen. 
				vec2 screenCoord = vec2(x / float(rendered_image.width), y / float(rendered_image.height));

				//Task 1: introduce some randomness and jittering
				screenCoord.x += (randf() - 0.5) / rendered_image.width;
				screenCoord.y += (randf() - 0.5) / rendered_image.height;

				primaryRay.d = normalize(lower_right_corner + screenCoord.x * X + screenCoord.y * Y);

				//Add some DoF
				if (depthOfField.lensRadius > 0) {
					vec2 ret;
					concentricSampleDisk(&ret.x, &ret.y);
					vec2 lensPoint = depthOfField.lensRadius *  ret;
					float t = focusDistance / (dot(primaryRay.d, camera_dir));
					vec3 focusPoint = primaryRay.o + primaryRay.d * t;
					primaryRay.o = primaryRay.o + camera_right * lensPoint.x + camera_up * lensPoint.y;
					primaryRay.d = normalize(focusPoint - primaryRay.o);
				}
				
				if(settings.supersampling_method == 0){
					if (intersect(primaryRay)) {
						// If it hit something, evaluate the radiance from that point
						color = Li(primaryRay).ToRGB();
					}
					else {
						float thit;
						Intersection lightIntersection;
						bool hitLight = shape->Intersect(primaryRay, &thit, &lightIntersection);
						if (hitLight) {
							color = areaLight->getLEmit().ToRGB();
						}
						else {
							color = Lenvironment(primaryRay.d).ToRGB();
						}
					}
				}else {
					vec3 firstQuadrant = corners_image.data[(y+1) * corners_image.width + (x+1)];
					vec3 secondQuadrant = corners_image.data[(y) * corners_image.width + (x + 1)];
					vec3 thirdQuadrant = corners_image.data[(y) * corners_image.width + (x)];
					vec3 fourthQuadrant = corners_image.data[(y + 1) * corners_image.width + (x)];
					int numberOfRays = 1; // On average, one corner per pixel (since many are shared)
					color = adaptiveSupersampling(&numberOfRays, camera_pos, camera_dir, camera_right, camera_up, focusDistance, lower_right_corner, float(x), float(y), X, Y, 1, firstQuadrant, secondQuadrant, thirdQuadrant, fourthQuadrant);
					if(customSettings.aassDensity) {
						color = vec3((float)numberOfRays / 30.f, 0.f, 0.f);
					}
				}
					
				
				// Accumulate the obtained radiance to the pixels color
				float n = float(rendered_image.number_of_samples);

				if (x == rendered_image.width / 2 && y == rendered_image.height / 2)
					rendered_image.data[y * rendered_image.width + x] =
					rendered_image.data[y * rendered_image.width + x] * (n / (n + 1.0f)) +
					(1.0f / (n + 1.0f)) * vec3(0.f, 0.f, 1.0f);
				else

				rendered_image.data[y * rendered_image.width + x] =
					rendered_image.data[y * rendered_image.width + x] * (n / (n + 1.0f)) +
					(1.0f / (n + 1.0f)) * color;
			}
		}
		rendered_image.number_of_samples += 1;
	}

	Spectrum EstimateDirect(Intersection &it,
		vec2 &uScattering, AreaLight &light,
		vec2 &uLight, BRDF &mat, bool isSurfaceHit, bool handleMedia, bool specular) {
		
		Spectrum Ld = Spectrum();
		//Sample light source with multiple importance sampling 858
		vec3 wi;
		Intersection lightHit;
		float lightPdf = 0, scatteringPdf = 0;
		Spectrum Li = light.Sample_Li(it, &lightHit, uLight, &wi, &lightPdf);

		if (lightPdf > 0 && !Li.IsBlack()) {
			//Compute BSDF value for light sample 859
			Spectrum f;
			
			
			Intersection &isect = (Intersection &)it;
			Spectrum reflectance = mat.f(wi, isect.wo, isect.shading_normal);
			f = reflectance * abs(dot(wi, isect.shading_normal));
			scatteringPdf = mat.PDF(wi, it.shading_normal, it.wo);
			
			
			if (!f.IsBlack()) {

				//Compute effect of visibility for light source sample 859

				const float distance_to_light = length(lightHit.position - isect.position);
				const float falloff_factor = 1.0f / (distance_to_light*distance_to_light);
				Ray lightRay(isect.position + EPSILON * isect.geometry_normal, normalize(lightHit.position - isect.position), 0.0f, distance_to_light);

				if (!occluded(lightRay)) {
					//vec3 wi = normalize(point_light.position - hit.position);
					Li  *= reflectance * falloff_factor * std::max(0.0f, dot(wi, isect.shading_normal))*point_light.intensity_multiplier;
				}
				else {
					Li = Spectrum(0.0f);
				}
				//Add light’s contribution to reflected radiance 860
				if (!Li.IsBlack()) {
					float weight = PowerHeuristic(1, lightPdf, 1, scatteringPdf);
					Ld += f * Li * weight / lightPdf;
				}
			}
		}

		//Sample BSDF with multiple importance sampling 860
		Spectrum f;
		bool sampledSpecular = false;
		//Sample scattered direction for surface interactions 860
		const Intersection &isect = (const Intersection &)it;
		f = mat.sample_wi(wi, isect.wo, isect.shading_normal, scatteringPdf);
		f *= abs(dot(wi, isect.shading_normal));
		//sampledSpecular = sampledType & BSDF_SPECULAR;
		
		if (!f.IsBlack() && scatteringPdf > 0) {
			//Account for light contributions along sampled direction wi 861
			
			lightPdf = light.Pdf_Li(it, wi);
			if (lightPdf == 0) {
				return Ld;
			}

			
			float weight = PowerHeuristic(1, scatteringPdf, 1, lightPdf);

			//Find intersection and compute transmittance 861
			Intersection lightIsect;
			Ray ray;
			if (dot(it.shading_normal, wi) > 0.0f) {
				ray = Ray(it.position + EPSILON * it.geometry_normal, wi);
			}
			else {
				ray = Ray(it.position - EPSILON * it.geometry_normal, wi);
			}

			Spectrum Tr(1.f);
			float tHit;
			bool lightIntersect = shape->Intersect(ray, &tHit, &it);
			if (lightIntersect) {
				ray.tfar = tHit;
				lightIntersect = !occluded(ray);
			}
			//Add light contribution from material sampling 861
			
			
			Li = Spectrum();
			if (lightIntersect) {
				Li = light.Le(vec3(0.0f), -wi);
			}
			if (!Li.IsBlack()){
				Ld += f * Li * Tr * weight / scatteringPdf;
			}
		}

		return Ld;
	}

	float PowerHeuristic(int nf, float fPdf, int ng, float gPdf) {
		float f = nf * fPdf, g = ng * gPdf;
		return (f * f) / (f * f + g * g);
	}
};



