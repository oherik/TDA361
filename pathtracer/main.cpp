#include <GL/glew.h>
#include <stb_image.h>
#include <chrono>
#include <iostream>
#include <labhelper.h>
#include <imgui.h>
#include <imgui_impl_sdl_gl3.h>
#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>
#include <Model.h>
#include <string>
#include "Pathtracer.h"
#include "embree.h"
#include "spectrum.h"
#include "material.h"
#include "Lights.h"

#define _CRTDBG_MAP_ALLOC  
#include <stdlib.h>  
#include <crtdbg.h>  


using namespace glm;
using namespace std; 

///////////////////////////////////////////////////////////////////////////////
// Various globals
///////////////////////////////////////////////////////////////////////////////
SDL_Window* g_window = nullptr;
int windowWidth = 0, windowHeight = 0;
vec3 tempR, tempG, tempN, tempK;

///////////////////////////////////////////////////////////////////////////////
// Shader programs
///////////////////////////////////////////////////////////////////////////////
GLuint shaderProgram;

///////////////////////////////////////////////////////////////////////////////
// GL texture to put pathtracing result into
///////////////////////////////////////////////////////////////////////////////
uint32_t pathtracer_result_txt_id; 

///////////////////////////////////////////////////////////////////////////////
// Camera parameters.
///////////////////////////////////////////////////////////////////////////////
vec3 cameraPosition(-30.0f, 10.0f, 30.0f);
vec3 cameraDirection = normalize(vec3(0.0f, 10.0f, 0.0f) - cameraPosition);
vec3 worldUp(0.0f, 1.0f, 0.0f);

///////////////////////////////////////////////////////////////////////////////
// Models
///////////////////////////////////////////////////////////////////////////////
vector<pair<labhelper::Model *, mat4>> models; 

///////////////////////////////////////////////////////////////////////////////
// Load shaders, environment maps, models and so on
///////////////////////////////////////////////////////////////////////////////
void initialize()
{

	///////////////////////////////////////////////////////////////////////////
	// Set initial custom settings
	///////////////////////////////////////////////////////////////////////////
	pathtracer::customSettings.bumpmap = true;
	pathtracer::customSettings.diffusemap = true;
	pathtracer::customSettings.spectrum = true;
	pathtracer::customSettings.aassDensity = false;

	///////////////////////////////////////////////////////////////////////////
	// Set initial value of BRDF 
	///////////////////////////////////////////////////////////////////////////
    pathtracer::brdf.fresnelCurrent = 0;
    pathtracer::brdf.diffuseCurrent = 0;
    pathtracer::brdf.geometricCurrent = 0;

	///////////////////////////////////////////////////////////////////////////
	// Load shader program
	///////////////////////////////////////////////////////////////////////////
	shaderProgram = labhelper::loadShaderProgram("../pathtracer/simple.vert", "../pathtracer/simple.frag");

	///////////////////////////////////////////////////////////////////////////
	// Initial path-tracer settings
	///////////////////////////////////////////////////////////////////////////
	pathtracer::settings.max_bounces = 8;
	pathtracer::settings.max_paths_per_pixel = 0; // 0 = Infinite
	#ifdef _DEBUG
	pathtracer::settings.subsampling = 16; 
	#else
	pathtracer::settings.subsampling = 4;
	#endif
	pathtracer::settings.supersampling_method = 0;

	///////////////////////////////////////////////////////////////////////////
	// Set up light
	///////////////////////////////////////////////////////////////////////////
	pathtracer::point_light.intensity_multiplier = 2500.0f; 
	pathtracer::point_light.color = vec3(1.f, 1.f, 1.f);
	pathtracer::point_light.position = vec3(-10.0f, 40.0f, 10.0f);

	///////////////////////////////////////////////////////////////////////////
	// Set up arealight
	///////////////////////////////////////////////////////////////////////////
	//Initialize position and transformation matrixes
	vec3 * aLightPos = new vec3(0.0f, 0.0f, 0.0f);
	mat4 * lightToWorld = new mat4(1.0f);
	*lightToWorld += translate(*aLightPos);
	mat4 * worldToLight = new mat4(1.0f);
	*worldToLight = inverse(*lightToWorld);
	
	//Create shape and arealight
	//pathtracer::shape = new pathtracer::Disk(lightToWorld, worldToLight, 0.0f, 20.0f, 0.0f, 360.0f);
	pathtracer::shape = new pathtracer::Sphere(lightToWorld, worldToLight, 5.0f, -5.0f, 5.0f, 360.0f);
	pathtracer::areaLight = new pathtracer::DiffuseAreaLight(lightToWorld, Spectrum(1.0f), pathtracer::shape, 1);

	///////////////////////////////////////////////////////////////////////////
	// Set up depth of field
	///////////////////////////////////////////////////////////////////////////
	pathtracer::depthOfField.lensRadius = 0.0f; // 0 = disable
	pathtracer::depthOfField.focusDistance = 10;

	///////////////////////////////////////////////////////////////////////////
	// Load environment map 
	///////////////////////////////////////////////////////////////////////////
	pathtracer::environment.map.load("../scenes/envmaps/001.hdr");
	pathtracer::environment.multiplier = 1.0f; 

	///////////////////////////////////////////////////////////////////////////
	// Load .obj models to scene
	///////////////////////////////////////////////////////////////////////////

	models.push_back(make_pair(labhelper::loadModelFromOBJ("../scenes/NewShip.obj"), translate(vec3(0.0f, -20.0f, 0.0f))));
	//models.push_back(make_pair(labhelper::loadModelFromOBJ("../scenes/NewShip.obj"), translate(vec3(0.0f, 10.0f, 10.0f))));
	models.push_back(make_pair(labhelper::loadModelFromOBJ("../scenes/landingpad2.obj"), translate(vec3(0.0f, -40.0f, 0.0f))));
	//models.push_back(make_pair(labhelper::loadModelFromOBJ("../scenes/cornellbottle2.obj"), translate(vec3(0.0f, 0.0f, -15.0f))));
	//models.push_back(make_pair(labhelper::loadModelFromOBJ("../scenes/Monk.obj"), translate(vec3(0.0f, 0.0f, -15.0f))));

	//models.push_back(make_pair(labhelper::loadModelFromOBJ("../scenes/bigsphere.obj"), scale(mat4(1.0f), vec3(0.05f))));

	//models.push_back(make_pair(labhelper::loadModelFromOBJ("../scenes/water.obj"), scale(mat4(1.0f),vec3(0.02f))));
	//mat4 asd = mat4(1.0f);
	//models.push_back(make_pair(labhelper::loadModelFromOBJ("../scenes/island.obj"), translate(
	//	glm::rotate(asd, (glm::mediump_float)1.56, glm::vec3(1.0f, 0.0f, 0.0f))
	//	, vec3(-100, -10, 0))));
	//models.push_back(make_pair(labhelper::loadModelFromOBJ("../scenes/island.obj"), translate(mat4(1.0f), vec3(0, -10, 0))));
	
	///////////////////////////////////////////////////////////////////////////
	// Add models to pathtracer scene
	///////////////////////////////////////////////////////////////////////////
	for (auto m : models) {
		pathtracer::addModel(m.first, m.second);
	}	
	pathtracer::buildBVH();

	///////////////////////////////////////////////////////////////////////////
	// Generate result texture
	///////////////////////////////////////////////////////////////////////////
	glGenTextures(1, &pathtracer_result_txt_id); 
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, pathtracer_result_txt_id);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);	

	///////////////////////////////////////////////////////////////////////////
	// This is INCORRECT! But an easy way to get us a brighter image that 
	// just looks a little better... 
	///////////////////////////////////////////////////////////////////////////
	//glEnable(GL_FRAMEBUFFER_SRGB);

	// Initialize spectrum sampling
	SampledSpectrum::Init();
}

void display(void)
{
	{	///////////////////////////////////////////////////////////////////////
		// If first frame, or window resized, or subsampling changes, 
		// inform the pathtracer
		///////////////////////////////////////////////////////////////////////
		int w, h; 
		SDL_GetWindowSize(g_window, &w, &h);
		static int old_subsampling; 
		if (windowWidth != w || windowHeight != h || old_subsampling != pathtracer::settings.subsampling) {
			pathtracer::resize(w, h);
			windowWidth = w; 
			windowWidth = h;
			old_subsampling = pathtracer::settings.subsampling; 
		}
	}

	///////////////////////////////////////////////////////////////////////////
	// Trace one path per pixel
	///////////////////////////////////////////////////////////////////////////
	vec3 cameraRight = normalize(cross(cameraDirection, worldUp));
	vec3 cameraUp = normalize(cross(cameraRight, cameraDirection));
	pathtracer::tracePaths(cameraPosition, cameraDirection, cameraUp);

	///////////////////////////////////////////////////////////////////////////
	// Copy pathtraced image to texture for display
	///////////////////////////////////////////////////////////////////////////
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, pathtracer::rendered_image.width, pathtracer::rendered_image.height,
		0, GL_RGB, GL_FLOAT, pathtracer::rendered_image.getPtr());

	///////////////////////////////////////////////////////////////////////////
	// Render a fullscreen quad, textured with our pathtraced image.
	///////////////////////////////////////////////////////////////////////////
	glViewport(0, 0, windowWidth, windowHeight);
	glClearColor(0.1,0.1,0.6,1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	SDL_GetWindowSize(g_window, &windowWidth, &windowHeight);
	glUseProgram(shaderProgram);
	labhelper::drawFullScreenQuad();
}

bool handleEvents(void)
{
	// check events (keyboard among other)
	SDL_Event event;
	bool quitEvent = false;

	// Allow ImGui to capture events.
	ImGuiIO& io = ImGui::GetIO();

	while (SDL_PollEvent(&event)) {
		ImGui_ImplSdlGL3_ProcessEvent(&event);

		if (event.type == SDL_QUIT || (event.type == SDL_KEYUP && event.key.keysym.sym == SDLK_ESCAPE)) {
			quitEvent = true;
		}
		
		if (!io.WantCaptureMouse) {
			if (event.type == SDL_MOUSEMOTION) {
				static int prev_xcoord = event.motion.x;
				static int prev_ycoord = event.motion.y;
				int delta_x = event.motion.x - prev_xcoord;
				int delta_y = event.motion.y - prev_ycoord;
				if (event.button.button & SDL_BUTTON(SDL_BUTTON_LEFT)) {
					float rotationSpeed = 0.005f;
					mat4 yaw = rotate(rotationSpeed * -delta_x, worldUp);
					mat4 pitch = rotate(rotationSpeed * -delta_y, normalize(cross(cameraDirection, worldUp)));
					cameraDirection = vec3(pitch * yaw * vec4(cameraDirection, 0.0f));
					pathtracer::restart(); 
				}
				prev_xcoord = event.motion.x;
				prev_ycoord = event.motion.y;
			}
		}
	}

	if (!io.WantCaptureKeyboard)
	{
		// check keyboard state (which keys are still pressed)
		const uint8_t *state = SDL_GetKeyboardState(nullptr);
		vec3 cameraRight = cross(cameraDirection, worldUp);
		const float speed = 0.5f; 
		if (state[SDL_SCANCODE_W]) {
			cameraPosition += speed * cameraDirection;
			pathtracer::restart();
		}
		if (state[SDL_SCANCODE_S]) {
			cameraPosition -= speed * cameraDirection;
			pathtracer::restart();
		}
		if (state[SDL_SCANCODE_A]) {
			cameraPosition -= speed * cameraRight;
			pathtracer::restart();
		}
		if (state[SDL_SCANCODE_D]) {
			cameraPosition += speed * cameraRight;
			pathtracer::restart();
		}
		if (state[SDL_SCANCODE_Q]) {
			cameraPosition -= speed * worldUp;
			pathtracer::restart();
		}
		if (state[SDL_SCANCODE_E]) {
			cameraPosition += speed * worldUp;
			pathtracer::restart();
		}
	}

	return quitEvent;
}


float rgToN(float r_val, float g_val) {
	float r_fix = std::min(r_val, 1.0f - EPSILON);
	return g_val * (1.0f - r_fix) / (1.0f + r_fix) + (1.0f - g_val)  * (1.0f + sqrt(r_fix)) / (1.0f - sqrt(r_fix));
}
float rgToK(float r_val, float g_val) {
	float r_fix = std::min(r_val, 1.0f - EPSILON);
	return sqrt(1.0f / (1.0f - r_fix) * r_fix * (pow((rgToN(r_fix, g_val) + 1.0f), 2.0f) - pow(rgToN(r_fix, g_val) - 1.0f, 2.0f)));
}
float nkToR(float n, float k) {
	return (pow(n - 1.0f, 2.0f) + pow(k, 2.0f)) / (pow(n + 1.0, 2.0f) + pow(k, 2.0f));
}
float nkToG(float n, float k) {
	return   ((1.0f + sqrt(nkToR(n, k))) / (1.0f - sqrt(nkToR(n, k))) - n) / ((1.0f + sqrt(nkToR(n, k))) / (1.0f - sqrt(nkToR(n, k))) - (1.0f - nkToR(n, k)) / (1.0f + nkToR(n, k)));
}

void gui() {
	// Inform imgui of new frame
	ImGui_ImplSdlGL3_NewFrame(g_window);

	///////////////////////////////////////////////////////////////////////////
	// Helpers for getting lists of materials and meshes into widgets
	///////////////////////////////////////////////////////////////////////////
	static auto model_getter = [](void * vec, int idx, const char ** text){
		auto& v = *(static_cast<std::vector<pair<labhelper::Model *, mat4>> *>(vec));
		if (idx < 0 || idx >= static_cast<int>(v.size())) { return false; }
		*text = v[idx].first->m_name.c_str();
		return true;
	};


	static auto mesh_getter = [](void * vec, int idx, const char ** text){
		auto& vector = *static_cast<std::vector<labhelper::Mesh>*>(vec);
		if (idx < 0 || idx >= static_cast<int>(vector.size())) { return false; }
		*text = vector[idx].m_name.c_str();
		return true;
	};

	static auto material_getter = [](void * vec, int idx, const char ** text){
		auto& vector = *static_cast<std::vector<labhelper::Material>*>(vec);
		if (idx < 0 || idx >= static_cast<int>(vector.size())) { return false; }
		*text = vector[idx].m_name.c_str();
		return true;
	};

	///////////////////////////////////////////////////////////////////////////
	// Choose a model to modify
	///////////////////////////////////////////////////////////////////////////
	static int model_index = 0;
	static labhelper::Model * model = models[0].first; 
	static int mesh_index = 0;
	static int material_index = model->m_meshes[mesh_index].m_material_idx;

	if (ImGui::Combo("Model", &model_index, model_getter,
		(void *)&models, models.size())) {
		model = models[model_index].first;
		mesh_index = 0; 
		material_index = model->m_meshes[mesh_index].m_material_idx;
	}

	///////////////////////////////////////////////////////////////////////////
	// Toggle our modifications
	///////////////////////////////////////////////////////////////////////////

	if (ImGui::CollapsingHeader("Enable modifications", "changes_ch", true, true)) {
		ImGui::Checkbox("Enable diffuse map", &pathtracer::customSettings.diffusemap); 
		ImGui::Checkbox("Enable bumpmap", &pathtracer::customSettings.bumpmap); 
		ImGui::Checkbox("Enable spectrum", &pathtracer::customSettings.spectrum); 
		ImGui::Checkbox("Show AASS ray density", &pathtracer::customSettings.aassDensity);
	}

	///////////////////////////////////////////////////////////////////////////
	// List all meshes in the model and show properties for the selected
	///////////////////////////////////////////////////////////////////////////

	if (ImGui::CollapsingHeader("Meshes", "meshes_ch", true, true))
	{
		if (ImGui::ListBox("Meshes", &mesh_index, mesh_getter,
			(void*)&model->m_meshes, model->m_meshes.size(), 8)) {
			material_index = model->m_meshes[mesh_index].m_material_idx;
		}

		labhelper::Mesh & mesh = model->m_meshes[mesh_index];
		char name[256];
		strcpy(name, mesh.m_name.c_str());
		if (ImGui::InputText("Mesh Name", name, 256)) { mesh.m_name = name; }
		labhelper::Material & selected_material = model->m_materials[material_index];
		if (ImGui::Combo("Material", &material_index, material_getter,
			(void *)&model->m_materials, model->m_materials.size())) {
			mesh.m_material_idx = material_index;
		}
	}

	///////////////////////////////////////////////////////////////////////////
	// List all materials in the model and show properties for the selected
	///////////////////////////////////////////////////////////////////////////
	bool rgChanged = false;
	if (ImGui::CollapsingHeader("Materials", "materials_ch", true, true))
	{

		
		
		ImGui::ListBox("Materials", &material_index, material_getter,
			(void*)&model->m_materials, model->m_materials.size(), 8);
		labhelper::Material & material = model->m_materials[material_index];
		char name[256];
		strcpy(name, material.m_name.c_str());

		tempG = material.m_g;
		tempR = material.m_r;
		tempN = material.m_n;
		tempK = material.m_k;

		if (ImGui::InputText("Material Name", name, 256)) { material.m_name = name; }
		ImGui::ColorEdit3("Color", &material.m_color.x);
		ImGui::ColorEdit3("Normal reflection color", &material.m_r.x);
		ImGui::ColorEdit3("Grazing reflection color", &material.m_g.x);
		ImGui::InputFloat("n, red", &material.m_n.x, 0.1f, 0.1f, -1, 0);
		ImGui::InputFloat("n, green", &material.m_n.y, 0.1f, 0.1f, -1, 0);
		ImGui::InputFloat("n, blue", &material.m_n.z, 0.1f, 0.1f, -1, 0);
		ImGui::InputFloat("k, red", &material.m_k.x, 0.1f, 0.1f, -1, 0);
		ImGui::InputFloat("k, green", &material.m_k.y, 0.1f, 0.1f, -1, 0);
		ImGui::InputFloat("k, blue", &material.m_k.z, 0.1f, 0.1f, -1, 0);
		ImGui::SliderFloat("Reflectivity", &material.m_reflectivity, 0.0f, 1.0f);
		ImGui::SliderFloat("Metalness", &material.m_metalness, 0.0f, 1.0f);
		ImGui::SliderFloat("Fresnel", &material.m_fresnel, 0.0f, 1.0f);
		ImGui::SliderFloat("shininess", &material.m_shininess, 0.0f, 25000.0f);
		ImGui::SliderFloat("Emission", &material.m_emission, 0.0f, 10.0f);
		ImGui::SliderFloat("Transparency", &material.m_transparency, 0.0f, 1.0f);

		if (material.m_g != tempG || material.m_r != tempR) {
			material.m_n.x = rgToN(material.m_r.x, material.m_g.x);
			material.m_n.y = rgToN(material.m_r.y, material.m_g.y);
			material.m_n.z = rgToN(material.m_r.z, material.m_g.z);

			material.m_k.x = rgToK(material.m_r.x, material.m_g.x);
			material.m_k.y = rgToK(material.m_r.y, material.m_g.y);
			material.m_k.z = rgToK(material.m_r.z, material.m_g.z);
		} else if (material.m_n != tempN || material.m_k != tempK) {
			material.m_r.x = nkToR(material.m_n.x, material.m_k.x);
			material.m_r.y = nkToR(material.m_n.y, material.m_k.y);
			material.m_r.z = nkToR(material.m_n.z, material.m_k.z);

			material.m_g.x = nkToG(material.m_n.x, material.m_k.x);
			material.m_g.y = nkToG(material.m_n.y, material.m_k.y);
			material.m_g.z = nkToG(material.m_n.z, material.m_k.z);
		}
		
	}

	
	///////////////////////////////////////////////////////////////////////////
	// BRDF modifiers 
	///////////////////////////////////////////////////////////////////////////

    const char* fresnelTerms[] { "Schlick's approximation"};
    const char* diffuseTerms[] = { "Blinn Phong", "Beckmann"};
    const char* geometricTerms[] = { "Cook-Torrance", "Smith-Schlick", "Smith-Walter"};

	if (ImGui::CollapsingHeader("BRDF modifiers", "brdf_ch", true, true))
    {
        ImGui::ListBox("Fresnell Term", &pathtracer::brdf.fresnelCurrent, fresnelTerms, sizeof(fresnelTerms)/sizeof(fresnelTerms[0]), 2);
        ImGui::ListBox("Diffuse Term", &pathtracer::brdf.diffuseCurrent, diffuseTerms, sizeof(diffuseTerms)/sizeof(diffuseTerms[0]), 2);
        ImGui::ListBox("Geometric Term", &pathtracer::brdf.geometricCurrent, geometricTerms, sizeof(geometricTerms)/sizeof(geometricTerms[0]), 2);
    }

	///////////////////////////////////////////////////////////////////////////
	// Depth of Field
	///////////////////////////////////////////////////////////////////////////

	if (ImGui::CollapsingHeader("Depth of field", "dof_ch", true, true))
	{
		ImGui::SliderFloat("Lens radius (0 to disable)", &pathtracer::depthOfField.lensRadius, 0.f, 1.f);
	}

	///////////////////////////////////////////////////////////////////////////
	// Light and environment map
	///////////////////////////////////////////////////////////////////////////
	if (ImGui::CollapsingHeader("Light sources", "lights_ch", true, true))
	{
		ImGui::SliderFloat("Environment multiplier", &pathtracer::environment.multiplier, 0.0f, 10.0f);
		ImGui::ColorEdit3("Point light color", &pathtracer::point_light.color.x);
		ImGui::SliderFloat("Point light intensity multiplier", &pathtracer::point_light.intensity_multiplier, 0.0f, 10000.0f);
	}

	///////////////////////////////////////////////////////////////////////////
	// Pathtracer settings
	///////////////////////////////////////////////////////////////////////////
	const char* superSamplingMethods[]{ "None", "Adaptive supersampling" };
	if (ImGui::CollapsingHeader("Pathtracer", "pathtracer_ch", true, true))
	{
		ImGui::SliderInt("Subsampling", &pathtracer::settings.subsampling, 1, 16);
		ImGui::ListBox("Supersampling method", &pathtracer::settings.supersampling_method, superSamplingMethods, sizeof(superSamplingMethods) / sizeof(superSamplingMethods[0]), 2);
		ImGui::SliderInt("Max Bounces", &pathtracer::settings.max_bounces, 0, 16);
		ImGui::SliderInt("Max Paths Per Pixel", &pathtracer::settings.max_paths_per_pixel, 0, 1024);

	}

	///////////////////////////////////////////////////////////////////////////
	// A button for saving your results
	///////////////////////////////////////////////////////////////////////////
	if (ImGui::Button("Save Materials")) {
		for (auto & m : models) {
			labhelper::saveModelToOBJ(m.first, m.first->m_filename);
		}
	}

	// Render the GUI.
	ImGui::Render();
}

int main(int argc, char *argv[])
{
	g_window = labhelper::init_window_SDL("Pathtracer", 1280, 720);

	initialize();

	bool stopRendering = false;
	auto startTime = std::chrono::system_clock::now();

	while (!stopRendering) {
		// render to window
		display();

		// Then render overlay GUI.
		gui();

		// Swap front and back buffer. This frame will now be displayed.
		SDL_GL_SwapWindow(g_window);  

		// check events (keyboard among other)
		stopRendering = handleEvents();
	}

	// Delete Models
	for (auto & m : models) {
		labhelper::freeModel(m.first);
	}
	// Shut down everything. This includes the window and all other subsystems.
	labhelper::shutDown(g_window);
	return 0;          
}
