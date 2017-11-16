#version 420
///////////////////////////////////////////////////////////////////////////////
// Input vertex attributes
///////////////////////////////////////////////////////////////////////////////
layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normalIn;
layout(location = 2) in vec2 texCoordIn;

///////////////////////////////////////////////////////////////////////////////
// Input uniform variables
///////////////////////////////////////////////////////////////////////////////
uniform mat4 normalMatrix;
uniform mat4 modelViewMatrix;
uniform mat4 modelViewProjectionMatrix;
uniform mat4 lightMatrix;

///////////////////////////////////////////////////////////////////////////////
// Output to fragment shader
///////////////////////////////////////////////////////////////////////////////
out vec2 texCoord;
out vec3 viewSpaceNormal;
out vec3 viewSpacePosition;

//Task 2
out vec4 shadowMapCoord;


void main() 
{
	gl_Position = modelViewProjectionMatrix * vec4(position, 1.0);
	texCoord = texCoordIn; 
	viewSpaceNormal = (normalMatrix * vec4(normalIn, 0.0)).xyz;
	viewSpacePosition = (modelViewMatrix * vec4(position, 1.0)).xyz;

	//Task 2
	shadowMapCoord = lightMatrix * vec4(viewSpacePosition, 1.0); // Så blir -1 till 1 på alla, i light clipping space

	shadowMapCoord.xyz *= vec3(0.5, 0.5, 0.5); //Skala alla med 0.5. * är **komponentvis** multiplikation
	shadowMapCoord.xyz  += (shadowMapCoord.w * vec3(0.5, 0.5, 0.5)); //Lägg till 0.5 på alla. ".w" for later division
}
