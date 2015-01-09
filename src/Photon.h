#pragma once

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <string>

struct Photon;

struct Photon {
	glm::vec3 loc;
	glm::vec3 dir;
	glm::vec3 pow;
	short flag;

	inline Photon()
	{
	}

	inline Photon(glm::vec3 l, glm::vec3 d, glm::vec3 p)
		: loc(l), dir(d), pow(p)
	{
	}
};