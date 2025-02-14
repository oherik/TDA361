
/* WARNING: if you change one the following, you also have to change the
* fragment shader!
*/
#define SOLUTION_USE_BUILTIN_SHADOW_TEST 0

#include <GL/glew.h>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <chrono>

#include <labhelper.h>
#include <imgui.h>
#include <imgui_impl_sdl_gl3.h>

#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>
using namespace glm;

#include <Model.h>
#include "hdr.h"
#include "fbo.h"

using std::min;
using std::max;

///////////////////////////////////////////////////////////////////////////////
// Various globals
///////////////////////////////////////////////////////////////////////////////
SDL_Window* g_window = nullptr;
float currentTime = 0.0f;
bool showUI = false;

///////////////////////////////////////////////////////////////////////////////
// Shader programs
///////////////////////////////////////////////////////////////////////////////
GLuint shaderProgram; // Shader for rendering the final image
GLuint simpleShaderProgram; // Shader used to draw the shadow map
GLuint backgroundProgram;

///////////////////////////////////////////////////////////////////////////////
// Environment
///////////////////////////////////////////////////////////////////////////////
float environment_multiplier = 0.9f;
GLuint environmentMap, irradianceMap, reflectionMap;
const std::string envmap_base_name = "001";

///////////////////////////////////////////////////////////////////////////////
// Light source
///////////////////////////////////////////////////////////////////////////////
vec3 lightPosition;
vec3 point_light_color = vec3(1.f, 1.f, 1.f);
bool useSpotLight = false;
float innerSpotlightAngle = 21.0f;
float outerSpotlightAngle = 22.5f;
float point_light_intensity_multiplier = 10000.0f;


///////////////////////////////////////////////////////////////////////////////
// Shadow map
///////////////////////////////////////////////////////////////////////////////
enum ClampMode {
	Edge = 1,
	Border = 2
};

FboInfo shadowMapFB;
int shadowMapResolution = 128;
int shadowMapClampMode = ClampMode::Edge;
bool usePolygonOffest = true;
bool useSoftFalloff = false;
bool useHardwarePCF = false;
float polygonOffset_factor = 1.1f;
float polygonOffset_units = 1.0f;


///////////////////////////////////////////////////////////////////////////////
// Camera parameters.
///////////////////////////////////////////////////////////////////////////////
vec3 cameraPosition(-70.0f, 50.0f, 70.0f);
vec3 cameraDirection = normalize(vec3(0.0f) - cameraPosition);
float cameraSpeed = 1.0f;

vec3 worldUp(0.0f, 1.0f, 0.0f);

///////////////////////////////////////////////////////////////////////////////
// Models
///////////////////////////////////////////////////////////////////////////////
labhelper::Model *fighterModel = nullptr;
labhelper::Model *landingpadModel = nullptr;
labhelper::Model *sphereModel = nullptr;

mat4 roomModelMatrix;
mat4 fighterModelMatrix;

void initGL()
{
	///////////////////////////////////////////////////////////////////////
	//		Load Shaders
	///////////////////////////////////////////////////////////////////////
	backgroundProgram   = labhelper::loadShaderProgram("../lab6-shadowmaps/background.vert", "../lab6-shadowmaps/background.frag");
	shaderProgram       = labhelper::loadShaderProgram("../lab6-shadowmaps/shading.vert",    "../lab6-shadowmaps/shading.frag");
	simpleShaderProgram = labhelper::loadShaderProgram("../lab6-shadowmaps/simple.vert",     "../lab6-shadowmaps/simple.frag");

	///////////////////////////////////////////////////////////////////////
	// Load models and set up model matrices
	///////////////////////////////////////////////////////////////////////
	fighterModel = labhelper::loadModelFromOBJ("../scenes/NewShip.obj");
	landingpadModel = labhelper::loadModelFromOBJ("../scenes/landingpad.obj");
	sphereModel = labhelper::loadModelFromOBJ("../scenes/sphere.obj");
  
	roomModelMatrix =  mat4(1.0f);
	fighterModelMatrix = translate(15.0f * worldUp);

	///////////////////////////////////////////////////////////////////////
	// Load environment map
	///////////////////////////////////////////////////////////////////////
	const int roughnesses = 8;
	std::vector<std::string> filenames;
	for (int i = 0; i < roughnesses; i++)
		filenames.push_back("../scenes/envmaps/" + envmap_base_name + "_dl_" + std::to_string(i) + ".hdr");

	reflectionMap = labhelper::loadHdrMipmapTexture(filenames);
	environmentMap = labhelper::loadHdrTexture("../scenes/envmaps/" + envmap_base_name + ".hdr");
	irradianceMap = labhelper::loadHdrTexture("../scenes/envmaps/" + envmap_base_name + "_irradiance.hdr");

	///////////////////////////////////////////////////////////////////////
	// Setup Framebuffer for shadow map rendering
	///////////////////////////////////////////////////////////////////////

	shadowMapFB.resize(shadowMapResolution, shadowMapResolution);

	glEnable(GL_DEPTH_TEST);	// enable Z-buffering 
	glEnable(GL_CULL_FACE);		// enables backface culling

	// Task 5
	glBindTexture(GL_TEXTURE_2D, shadowMapFB.depthBuffer);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC, GL_LEQUAL);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_COMPARE_REF_TO_TEXTURE);
}

void debugDrawLight(const glm::mat4 &viewMatrix, const glm::mat4 &projectionMatrix, const glm::vec3 &worldSpaceLightPos)
{
	mat4 modelMatrix = glm::translate(worldSpaceLightPos);
	glUseProgram(shaderProgram);
	labhelper::setUniformSlow(shaderProgram, "modelViewProjectionMatrix", projectionMatrix * viewMatrix * modelMatrix);
	labhelper::render(sphereModel);
	labhelper::setUniformSlow(shaderProgram, "modelViewProjectionMatrix", projectionMatrix * viewMatrix);
	labhelper::debugDrawLine(viewMatrix, projectionMatrix, worldSpaceLightPos);
}


void drawBackground(const mat4 &viewMatrix, const mat4 &projectionMatrix)
{
	glUseProgram(backgroundProgram);
	labhelper::setUniformSlow(backgroundProgram, "environment_multiplier", environment_multiplier);
	labhelper::setUniformSlow(backgroundProgram, "inv_PV", inverse(projectionMatrix * viewMatrix));
	labhelper::setUniformSlow(backgroundProgram, "camera_pos", cameraPosition);
	labhelper::drawFullScreenQuad();
}

void drawScene(GLuint currentShaderProgram, const mat4 &viewMatrix, const mat4 &projectionMatrix, const mat4 &lightViewMatrix, const mat4 &lightProjectionMatrix, const mat4 &lightMatrix)
{
	glUseProgram(currentShaderProgram);
	// Light source
	vec4 viewSpaceLightPosition = viewMatrix * vec4(lightPosition, 1.0f);
	labhelper::setUniformSlow(currentShaderProgram, "point_light_color", point_light_color);
	labhelper::setUniformSlow(currentShaderProgram, "point_light_intensity_multiplier", point_light_intensity_multiplier);
	labhelper::setUniformSlow(currentShaderProgram, "viewSpaceLightPosition", vec3(viewSpaceLightPosition));
	labhelper::setUniformSlow(currentShaderProgram, "viewSpaceLightDir", normalize(vec3(viewMatrix * vec4(-lightPosition, 0.0f))));
	labhelper::setUniformSlow(currentShaderProgram, "spotOuterAngle", std::cos(radians(outerSpotlightAngle)));
	labhelper::setUniformSlow(currentShaderProgram, "spotInnerAngle", std::cos(radians(innerSpotlightAngle)));

	// Environment
	labhelper::setUniformSlow(currentShaderProgram, "environment_multiplier", environment_multiplier);

	// camera
	labhelper::setUniformSlow(currentShaderProgram, "viewInverse", inverse(viewMatrix));

	// Task 2: light camera
	labhelper::setUniformSlow(currentShaderProgram, "lightMatrix", lightMatrix);

	// landing pad 
	mat4 modelMatrix(1.0f);
	labhelper::setUniformSlow(currentShaderProgram, "modelViewProjectionMatrix", projectionMatrix * viewMatrix * modelMatrix);
	labhelper::setUniformSlow(currentShaderProgram, "modelViewMatrix", viewMatrix * modelMatrix);
	labhelper::setUniformSlow(currentShaderProgram, "normalMatrix", inverse(transpose(viewMatrix * modelMatrix)));

	labhelper::render(landingpadModel);

	// Fighter
	labhelper::setUniformSlow(currentShaderProgram, "modelViewProjectionMatrix", projectionMatrix * viewMatrix * fighterModelMatrix);
	labhelper::setUniformSlow(currentShaderProgram, "modelViewMatrix", viewMatrix * fighterModelMatrix);
	labhelper::setUniformSlow(currentShaderProgram, "normalMatrix", inverse(transpose(viewMatrix * fighterModelMatrix)));

	labhelper::render(fighterModel);
}


void display(void)
{
	int w, h;
	SDL_GetWindowSize(g_window, &w, &h);

	///////////////////////////////////////////////////////////////////////////
	// setup matrices
	///////////////////////////////////////////////////////////////////////////
	mat4 projMatrix = perspective(radians(45.0f), float(w) / float(h), 5.0f, 500.0f);
	mat4 viewMatrix = lookAt(cameraPosition, cameraPosition + cameraDirection, worldUp);

	vec4 lightStartPosition = vec4(40.0f, 40.0f, 0.0f, 1.0f);
	lightPosition = vec3(rotate(currentTime, worldUp) * lightStartPosition);
	mat4 lightViewMatrix = lookAt(lightPosition, vec3(0.0f), worldUp);
	mat4 lightProjMatrix = perspective(radians(45.0f), 1.0f, 25.0f, 100.0f);

	//Task 2
	mat4 lightMatrix = lightProjMatrix * lightViewMatrix * inverse(viewMatrix);

	//Task 5
	mat4 scaleMatrix = scale(vec3(0.5));
	mat4 translateMatrix = translate(vec3(0.5));
	lightMatrix = translateMatrix * scaleMatrix * lightMatrix;

	///////////////////////////////////////////////////////////////////////////
	// Bind the environment map(s) to unused texture units
	///////////////////////////////////////////////////////////////////////////
	glActiveTexture(GL_TEXTURE6);
	glBindTexture(GL_TEXTURE_2D, environmentMap);
	glActiveTexture(GL_TEXTURE7);
	glBindTexture(GL_TEXTURE_2D, irradianceMap);
	glActiveTexture(GL_TEXTURE8);
	glBindTexture(GL_TEXTURE_2D, reflectionMap);
	glActiveTexture(GL_TEXTURE0);

	///////////////////////////////////////////////////////////////////////////
	// Set up shadow map parameters
	///////////////////////////////////////////////////////////////////////////
	// >>> @task 1
	if (shadowMapFB.width != shadowMapResolution || shadowMapFB.height != shadowMapResolution){
		shadowMapFB.resize(shadowMapResolution, shadowMapResolution);
	}

	//Task 3.2 clampilamp

	if (shadowMapClampMode == ClampMode::Edge){
		glBindTexture(GL_TEXTURE_2D, shadowMapFB.depthBuffer);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	}

	if (shadowMapClampMode == ClampMode::Border){
		glBindTexture(GL_TEXTURE_2D, shadowMapFB.depthBuffer);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
		vec4 zeros = vec4(0.0);
		glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, &zeros.x);
	}

	///////////////////////////////////////////////////////////////////////////
	// Draw Shadow Map
	///////////////////////////////////////////////////////////////////////////

	glBindFramebuffer(GL_FRAMEBUFFER, shadowMapFB.framebufferId);
	glViewport(0, 0, shadowMapFB.width, shadowMapFB.height); // Set the correct viewport mapping. Hur stort f�nster vi har. N�r vi renderar till v�r framebuffern, v�r textur.
	glClearColor(0.2, 0.2, 0.8, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//Remove acne
	if (usePolygonOffest){
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(polygonOffset_factor, polygonOffset_factor);
	}

	drawScene(simpleShaderProgram, lightViewMatrix, lightProjMatrix, lightViewMatrix, lightProjMatrix, lightMatrix);

	//Visualize the render from the camera's point of view
	labhelper::Material &screen = landingpadModel->m_materials[8]; //[8] is the TV screens
	screen.m_emission_texture.gl_id = shadowMapFB.colorTextureTarget; //Emissiva det som renderats fr�n ljusets synvinkel

	//Remove acne, step 2
	if (usePolygonOffest){
		glDisable(GL_POLYGON_OFFSET_FILL);
	}

	///////////////////////////////////////////////////////////////////////////
	// Draw from camera
	///////////////////////////////////////////////////////////////////////////
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glViewport(0, 0, w, h);
	glClearColor(0.2, 0.2, 0.8, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	
	drawBackground(viewMatrix, projMatrix);

	//Task 2: bind the depth buffer from the shadow map frame buffer to texture 10
	glActiveTexture(GL_TEXTURE10);
	glBindTexture(GL_TEXTURE_2D, shadowMapFB.depthBuffer);


	drawScene(shaderProgram, viewMatrix, projMatrix, lightViewMatrix, lightProjMatrix, lightMatrix);
	debugDrawLight(viewMatrix, projMatrix, vec3(lightPosition));

	
	CHECK_GL_ERROR();
}

bool handleEvents(void)
{
	// check events (keyboard among other)
	SDL_Event event;
	bool quitEvent = false;
	while (SDL_PollEvent(&event)) {
		if (event.type == SDL_QUIT || (event.type == SDL_KEYUP && event.key.keysym.sym == SDLK_ESCAPE)) {
			quitEvent = true;
		}
		if (event.type == SDL_KEYUP && event.key.keysym.sym == SDLK_g) {
			showUI = !showUI;
		}
		if (event.type == SDL_MOUSEMOTION && !ImGui::IsMouseHoveringAnyWindow()) {
			static int prev_xcoord = event.motion.x;
			static int prev_ycoord = event.motion.y;
			int delta_x = event.motion.x - prev_xcoord;
			int delta_y = event.motion.y - prev_ycoord;

			if (event.button.button & SDL_BUTTON(SDL_BUTTON_LEFT)) {
				float rotationSpeed = 0.005f;
				mat4 yaw = rotate(rotationSpeed * -delta_x, worldUp);
				mat4 pitch = rotate(rotationSpeed * -delta_y, normalize(cross(cameraDirection, worldUp)));
				cameraDirection = vec3(pitch * yaw * vec4(cameraDirection, 0.0f));
			}
			prev_xcoord = event.motion.x;
			prev_ycoord = event.motion.y;
		}
	}

	// check keyboard state (which keys are still pressed)
	const uint8_t *state = SDL_GetKeyboardState(nullptr);
	vec3 cameraRight = cross(cameraDirection, worldUp);

	if (state[SDL_SCANCODE_W]) {
		cameraPosition += cameraSpeed* cameraDirection;
	}
	if (state[SDL_SCANCODE_S]) {
		cameraPosition -= cameraSpeed * cameraDirection;
	}
	if (state[SDL_SCANCODE_A]) {
		cameraPosition -= cameraSpeed * cameraRight;
	}
	if (state[SDL_SCANCODE_D]) {
		cameraPosition += cameraSpeed * cameraRight;
	}
	if (state[SDL_SCANCODE_Q]) {
		cameraPosition -= cameraSpeed * worldUp;
	}
	if (state[SDL_SCANCODE_E]) {
		cameraPosition += cameraSpeed * worldUp;
	}

	return quitEvent;
}

void gui()
{
	// Inform imgui of new frame
	ImGui_ImplSdlGL3_NewFrame(g_window);

	// ----------------- Set variables --------------------------
	ImGui::SliderInt("Shadow Map Resolution", &shadowMapResolution, 32, 2048);
	ImGui::Text("Polygon Offset");
	ImGui::Checkbox("Use polygon offset", &usePolygonOffest);
	ImGui::SliderFloat("Units", &polygonOffset_factor, 0.0f, 2.0f);
	ImGui::SliderFloat("Factor", &polygonOffset_units, 0.0f, 10.0f);
	ImGui::Text("Clamp Mode");
	ImGui::RadioButton("Clamp to edge", &shadowMapClampMode, ClampMode::Edge);
	ImGui::RadioButton("Clamp to border", &shadowMapClampMode, ClampMode::Border);
	ImGui::Checkbox("Use spot light", &useSpotLight);
	ImGui::Checkbox("Use soft falloff", &useSoftFalloff);
	ImGui::SliderFloat("Inner Deg.", &innerSpotlightAngle, 0.0f, 90.0f);
	ImGui::SliderFloat("Outer Deg.", &outerSpotlightAngle, 0.0f, 90.0f);
	ImGui::Checkbox("Use hardware PCF", &useHardwarePCF);
	ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
	// ----------------------------------------------------------

	// Render the GUI.
	ImGui::Render();
}

int main(int argc, char *argv[])
{
	g_window = labhelper::init_window_SDL("OpenGL Lab 6");

	initGL();

	bool stopRendering = false;
	auto startTime = std::chrono::system_clock::now();

	while (!stopRendering) {
		//update currentTime
		std::chrono::duration<float> timeSinceStart = std::chrono::system_clock::now() - startTime;
		currentTime = timeSinceStart.count();

		// render to window
		display();

		// Render overlay GUI.
		if (showUI) {
			gui();
		}
		
		// Swap front and back buffer. This frame will now been displayed.
		SDL_GL_SwapWindow(g_window);

		// check events (keyboard among other)
		stopRendering = handleEvents();
	}
	// Free Models
	labhelper::freeModel(fighterModel);
	labhelper::freeModel(landingpadModel);
	labhelper::freeModel(sphereModel);

	// Shut down everything. This includes the window and all other subsystems.
	labhelper::shutDown(g_window);
	return 0;
}
