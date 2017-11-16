#version 420

// required by GLSL spec Sect 4.5.3 (though nvidia does not, amd does)
precision highp float;

///////////////////////////////////////////////////////////////////////////////
// Material
///////////////////////////////////////////////////////////////////////////////
uniform vec3 material_color;
uniform float material_reflectivity;
uniform float material_metalness;
uniform float material_fresnel;
uniform float material_shininess;
uniform float material_emission;

uniform int has_emission_texture;
layout(binding = 5) uniform sampler2D emissiveMap;

///////////////////////////////////////////////////////////////////////////////
// Environment
///////////////////////////////////////////////////////////////////////////////
layout(binding = 6) uniform sampler2D environmentMap;
layout(binding = 7) uniform sampler2D irradianceMap;
layout(binding = 8) uniform sampler2D reflectionMap;
uniform float environment_multiplier;

///////////////////////////////////////////////////////////////////////////////
// Light source
///////////////////////////////////////////////////////////////////////////////
uniform vec3 point_light_color = vec3(1.0, 1.0, 1.0);
uniform float point_light_intensity_multiplier = 50.0;

///////////////////////////////////////////////////////////////////////////////
// Constants
///////////////////////////////////////////////////////////////////////////////
#define PI 3.14159265359

///////////////////////////////////////////////////////////////////////////////
// Input varyings from vertex shader
///////////////////////////////////////////////////////////////////////////////
in vec2 texCoord;
in vec3 viewSpaceNormal;
in vec3 viewSpacePosition;

//Task 2
in vec4 shadowMapCoord;
layout(binding = 10) uniform sampler2D shadowMapTex;

///////////////////////////////////////////////////////////////////////////////
// Input uniform variables
///////////////////////////////////////////////////////////////////////////////
uniform mat4 viewInverse;
uniform vec3 viewSpaceLightPosition;

///////////////////////////////////////////////////////////////////////////////
// Output color
///////////////////////////////////////////////////////////////////////////////
layout(location = 0) out vec4 fragmentColor;




////////////////////////////
// Från Lab 4
vec3 calculateIndirectIllumination(vec3 wo, vec3 n)
{
	vec3 nws = vec3(viewInverse*vec4(n, 0));


	// Calculate the world-space direction from the camera to that position
	vec3 dir = nws;

	// Calculate the spherical coordinates of the direction
	float theta = acos(max(-1.0f, min(1.0f, dir.y)));
	float phi = atan(dir.z, dir.x);
	if (phi < 0.0f) phi = phi + 2.0f * PI;

	// Use these to lookup the color in the environment map
	vec2 lookup = vec2(phi / (2.0 * PI), theta / PI);
	vec3 irradiance = vec3(environment_multiplier * texture(irradianceMap, lookup));

	vec3 diffuse_term = material_color * (1.0 / PI) * irradiance;



	vec3 wi = reflect(-wo, n);
	vec3 wi_world = (viewInverse * vec4(wi, 0.0)).xyz;

	dir = wi_world;

	// Calculate the spherical coordinates of the direction
	theta = acos(max(-1.0f, min(1.0f, dir.y)));
	phi = atan(dir.z, dir.x);
	if (phi < 0.0f) phi = phi + 2.0f * PI;

	// Use these to lookup the color in the environment map
	lookup = vec2(phi / (2.0 * PI), theta / PI);



	float roughness = sqrt(sqrt(2.0 / (material_shininess + 2.0)));
	vec3 Li = environment_multiplier * textureLod(reflectionMap, lookup, roughness * 7.0).xyz;


	// Calculate the half-way vector
	vec3 wh = normalize(wi + wo);

	// Approcimation for the fresnel term
	float F = material_fresnel + (1.0 - material_fresnel)*pow((1.0 - dot(wh, wi)), 5.0);

	vec3 dielectric_term = F * Li + (1.0 - F) * diffuse_term;
	vec3 metal_term = F * material_color * Li;

	vec3 microfacet_term = material_metalness * metal_term + (1 - material_metalness)*dielectric_term;

	return material_reflectivity*microfacet_term + (1.0 - material_reflectivity) * diffuse_term;

}


vec3 calculateDirectIllumiunation(vec3 wo, vec3 n)
{

	float d = length(viewSpaceLightPosition - viewSpacePosition);

	// Inverse-square law: the intensity is proportional to the spherical area with radius d.
	vec3 Li = point_light_intensity_multiplier * point_light_color * (1.0 / pow(d, 2));

	vec3 wi = normalize(viewSpaceLightPosition - viewSpacePosition);

	if (dot(n, wi) <= 0){
		return vec3(0);
	}

	vec3 diffuse_term = material_color * (1.0 / PI) * abs(dot(n, wi)) * Li;


	// Calculate the half-way vector
	vec3 wh = normalize(wi + wo);

	// Approcimation for the fresnel term
	float F = material_fresnel + (1 - material_fresnel)*pow((1 - dot(wh, wi)), 5);

	// Calculate the microfacet distribution function
	// Normalized Blinn-Phong
	//Fattar andra delen, men varför (s+2)/(2pi)?
	float D = (material_shininess + 2) / (2 * PI)*pow(dot(n, wh), material_shininess);

	//Calculate the shadowing/masking function
	//Probability of being blocked???
	float G = min(1, min(2 * dot(n, wh)*dot(n, wo) / dot(wo, wh), 2 * dot(n, wh)*dot(n, wi) / dot(wo, wh)));

	//Put it all together
	float brdf = F*D*G / (4 * dot(n, wo)*dot(n, wi));

	vec3 dielectric_term = brdf * dot(n, wi)*Li + (1 - F)*diffuse_term;
	vec3 metal_term = brdf * material_color * dot(n, wi)*Li;
	vec3 microfacet_term = material_metalness * metal_term + (1 - material_metalness)*dielectric_term;

	return material_reflectivity*microfacet_term + (1 - material_reflectivity) * diffuse_term;
	//return brdf * dot(n,wi)*Li;
}

void main() 
{
	// Task 2: shadow map
	float depth = texture(shadowMapTex, shadowMapCoord.xy / shadowMapCoord.w).x; // Juste fö den kommer ha ett visst svartvitt x-värde kanske eller?
	float visibility = (depth >= (shadowMapCoord.z / shadowMapCoord.w)) ? 1.0 : 0.0;

	float attenuation = 1.0;
	
	vec3 wo = -normalize(viewSpacePosition);
	vec3 n = normalize(viewSpaceNormal);

	// Direct illumination
	vec3 direct_illumination_term = visibility * calculateDirectIllumiunation(wo, n);

	// Indirect illumination
	vec3 indirect_illumination_term = calculateIndirectIllumination(wo, n);

	

	///////////////////////////////////////////////////////////////////////////
	// Add emissive term. If emissive texture exists, sample this term.
	///////////////////////////////////////////////////////////////////////////
	vec3 emission_term = material_emission * material_color;
	if (has_emission_texture == 1) {
		emission_term = texture(emissiveMap, texCoord).xyz;
	}

	vec3 shading = 
		direct_illumination_term +
		indirect_illumination_term +
		emission_term;

	fragmentColor = vec4(shading, 1.0);
	return;

}
