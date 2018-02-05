#include "Pathtracer.h"
#include <memory>
#include <iostream>
#include <map>
#include <algorithm>
#include "material.h"
#include "embree.h"
#include "sampling.h"

using namespace std; 
using namespace glm; 

namespace pathtracer
{
	///////////////////////////////////////////////////////////////////////////////
	// Global variables
	///////////////////////////////////////////////////////////////////////////////
	Settings settings;
	Environment environment; 
	Image rendered_image; 
	PointLight point_light; 

	///////////////////////////////////////////////////////////////////////////
	// Restart rendering of image
	///////////////////////////////////////////////////////////////////////////
	void restart()
	{
		// No need to clear image, 
		rendered_image.number_of_samples = 0; 
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
		restart(); 
	}

	///////////////////////////////////////////////////////////////////////////
	// Return the radiance from a certain direction wi from the environment
	// map. 
	///////////////////////////////////////////////////////////////////////////
	vec3 Lenvironment(const vec3 & wi) {
		const float theta = acos(std::max(-1.0f, std::min(1.0f, wi.y)));
		float phi = atan(wi.z, wi.x);
		if (phi < 0.0f) phi = phi + 2.0f * M_PI;
		vec2 lookup = vec2(phi / (2.0 * M_PI), theta / M_PI);
		return environment.multiplier * environment.map.sample(lookup.x, lookup.y);
	}

	///////////////////////////////////////////////////////////////////////////
	// Calculate the radiance going from one point (r.hitPosition()) in one 
	// direction (-r.d), through path tracing.  
	///////////////////////////////////////////////////////////////////////////
	vec3 Li(Ray & primary_ray) {
		vec3 L = vec3(0.0f);
		vec3 path_throughput = vec3(1.0f);
		Ray current_ray = primary_ray;
		vec3 last_position = vec3(0.0f);

		//Task 5: bounce it up
		for (int i = 0; i < settings.max_bounces; i++){
			
			// Get the intersection information from the ray
			Intersection hit = getIntersection(current_ray);

			// Create a material tree
			Diffuse diffuse(hit.material->m_color);
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
				//float a_new = 10*(hit.material->m_transparency)/dist; // /5 is a bit arbitrary
				//a = T;// fmin(1.0f, a_new);
			}

			TransparencyBlend transparency_blend(a, &transparent, hit.material->m_color);
			BlinnPhong dielectric(hit.material->m_shininess, hit.material->m_fresnel, &diffuse);
			BlinnPhongMetal metal(hit.material->m_color, hit.material->m_r, hit.material->m_g, hit.material->m_shininess,
				hit.material->m_fresnel);
			LinearBlend metal_blend(hit.material->m_metalness, &metal, &dielectric);
			LinearBlend reflectivity_blend(hit.material->m_reflectivity, &metal_blend, &diffuse);

			LinearBlend transparency_blend_final(hit.material->m_transparency, &transparency_blend, &reflectivity_blend);
			
			
			BRDF & mat = metal;

			// Update last hit position for the next distance calculation
			last_position = vec3(hit.position);


			// Calculate Direct Illumination from light.
			const float distance_to_light = length(point_light.position - hit.position);
			const float falloff_factor = 1.0f / (distance_to_light*distance_to_light);
			Ray lightRay(hit.position + EPSILON * hit.geometry_normal, normalize(point_light.position - hit.position), 0.0f, distance_to_light);

			if (!occluded(lightRay)){
				vec3 Li = point_light.intensity_multiplier * point_light.color * falloff_factor;
				vec3 wi = normalize(point_light.position - hit.position);
				L = L + path_throughput * (mat.f(wi, hit.wo, hit.shading_normal) * Li * std::max(0.0f, dot(wi, hit.shading_normal)));
			}

			// Emitted radiance from intersection
			L = L + path_throughput * hit.material->m_emission;

			// Sample incoming direction
			vec3 wi;
			float pdf = 0.0f;
			vec3 brdf = mat.sample_wi(wi, hit.wo, hit.shading_normal, pdf);

			if (pdf < EPSILON) return L;
			
			float cosineTerm = abs(dot(wi, hit.shading_normal));

			path_throughput = path_throughput * (brdf * cosineTerm) / pdf;

			//Break if the throughput is 0
			if (length(path_throughput) < EPSILON)
			{
				return L;
			}
			
			//Next ray
			if (dot(hit.shading_normal, wi) > 0.0f){
				current_ray = Ray(hit.position + EPSILON * hit.geometry_normal, wi);
			}
			else {
				current_ray = Ray(hit.position - EPSILON * hit.geometry_normal, wi);
			}

			//Intersect the new ray
			if (!intersect(current_ray)){
				return L + path_throughput * Lenvironment(current_ray.d);
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
		return L;
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
		// Stop here if we have as many samples as we want
		if ((int(rendered_image.number_of_samples) > settings.max_paths_per_pixel) &&
			(settings.max_paths_per_pixel != 0)) return;
		// Trace one path per pixel (the omp parallel stuf magically distributes the 
		// pathtracing on all cores of your CPU).
#pragma omp parallel for
		for (int y = 0; y < rendered_image.height; y++) {
			for (int x = 0; x < rendered_image.width; x++) {
				vec3 color;
				Ray primaryRay;
				primaryRay.o = camera_pos;
				// Create a ray that starts in the camera position and points toward
				// the current pixel on a virtual screen. 
				vec2 screenCoord = vec2(float(x) / float(rendered_image.width), float(y) / float(rendered_image.height));
				
				//Task 1: introduce some randomness and jittering
				screenCoord.x += (pathtracer::randf() -0.5) / rendered_image.width;
				screenCoord.y += (pathtracer::randf() - 0.5) / rendered_image.height;
		
				primaryRay.d = normalize(lower_right_corner + screenCoord.x * X + screenCoord.y * Y);
				
				// Intersect ray with scene
				if (intersect(primaryRay)) {
					// If it hit something, evaluate the radiance from that point
					color = Li(primaryRay);
				}
				else {
					// Otherwise evaluate environment
					color = Lenvironment(primaryRay.d);
				}
				// Accumulate the obtained radiance to the pixels color
				float n = float(rendered_image.number_of_samples);
				rendered_image.data[y * rendered_image.width + x] =
					rendered_image.data[y * rendered_image.width + x] * (n / (n + 1.0f)) +
					(1.0f / (n + 1.0f)) * color;
			}
		}
		rendered_image.number_of_samples += 1;
	}
};
