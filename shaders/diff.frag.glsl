#version 130
// ^ Change this to version 130 if you have compatibility issues

//these are the interpolated values out of the rasterizer, so you can't know
//their specific values without knowing the vertices that contributed to them
in vec4 fs_Normal;
in vec4 fs_LightVector;
in vec4 fs_Color;
in vec3 fs_Position;
out vec4 out_Color;

uniform vec3 u_LightColor;
uniform vec3 u_CameraPos;


void main()
{
    // Material base color (before shading)
    vec4 diffuseColor = fs_Color;

    // Calculate the diffuse term
    float diffuseTerm = dot(normalize(fs_Normal), normalize(fs_LightVector));
    // Avoid negative lighting values
    diffuseTerm = clamp(diffuseTerm, 0, 1);

    vec4 R = normalize(reflect(fs_LightVector, fs_Normal));
	vec3 V = normalize(fs_Position - u_CameraPos);
    float specularTerm = dot(R.rgb, V);
    specularTerm = clamp(specularTerm, 0, 1);
	specularTerm = pow(specularTerm, 10);

    float ka = 0.2, kd = 1, ks = 1;

    vec4 Color = vec4(diffuseColor.rgb * ka + u_LightColor * (kd * diffuseColor.rgb * diffuseTerm + ks * specularTerm * vec3(1,1,1) ) , diffuseColor.a);

    // Compute final shaded color

    out_Color = Color;
}
