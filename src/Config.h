#pragma once

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <fstream>
#include <string>


using namespace std;
using namespace glm;

class Config
{
public:
	Config(string filename);
	~Config(void);

	// Config
	vec3 RESO, EYEP, VDIR, UVEC;
	float FOVY;

	// light
	vec3 LPOS, LCOL;

	// node
	string NODE, PARENT, SHAPE;
	vec3 TRANSLATION, ROTATION, SCALE, CENTER, RGBA;
};

