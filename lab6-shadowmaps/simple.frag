#version 420

// required by GLSL spec Sect 4.5.3 (though nvidia does not, amd does)
precision highp float;

uniform vec3 material_color;

layout(location = 0) out vec4 fragmentColor;

void main() 
{
	fragmentColor = vec4(gl_FragCoord.z); //Let's B&W this image, yo :D
}
