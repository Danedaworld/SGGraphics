#version 130
// ^ Change this to version 130 if you have compatibility issues

uniform mat4 u_Model;
uniform mat4 u_ModelInvTr;
uniform mat4 u_ViewProj;
uniform vec4 u_LightPos;

in vec3 vs_Position;
in vec3 vs_Normal;
in vec3 vs_Color;

out vec4 fs_Normal;
out vec4 fs_LightVector;
out vec4 fs_Color;
out vec3 fs_Position;

void main()
{
    fs_Color = vec4(vs_Color, 1);
    fs_Normal = u_ModelInvTr * vec4(vs_Normal, 0);

    vec4 modelposition = u_Model * vec4(vs_Position, 1);

    // Set up our vector for the light
    fs_LightVector = u_LightPos - modelposition;

    //built-in things to pass down the pipeline
    gl_Position = u_ViewProj * modelposition;
	fs_Position = vs_Position;
}
