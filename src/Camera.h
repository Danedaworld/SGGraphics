#pragma once


#include "Config.h"
#include <fstream>
#include <vector>

class Camera
{
public:
	Camera();
	Camera(string f);
	~Camera(void);
	void parseFile();
	vec3 RESO, EYEP, VDIR, UVEC;
	vec3 H, V, E, C, U, A, B, M;
	float theta, phi;
	float FOVY, aspRatio;
private:
	string filename;
};

