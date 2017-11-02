#version 420

// required by GLSL spec Sect 4.5.3 (though nvidia does not, amd does)
precision highp float;

in vec3 outColor;
in vec2 texCoord;

// >>> @task 3.4
//Sample the texture from colorTexture
layout(binding = 0) uniform sampler2D colorTexture;
layout(location = 0) out vec4 fragmentColor;

void main() 
{
	// Task 1
	//fragmentColor = vec4(1.0, texCoord.x, texCoord.y, 0.0); 

	// >>> @task 3.5
	fragmentColor = texture2D(colorTexture, texCoord.xy); 
}