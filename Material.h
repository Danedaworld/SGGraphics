#pragma once

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <string>

using namespace std;
using namespace glm;

struct material;

struct material {

	material(string &name, vec3 &diff, vec3 &refl, float &expo, float &ior, float &mirr, float &tran, float emit) {
		NAME = name;
		DIFF = diff;
		REFL = refl;
		EXPO = expo;
		IOR = ior;
		MIRR = mirr;
		TRAN = tran;
		EMIT = emit;
		//cout << NAME<< EMIT << endl;
	}
	material() {
		NAME = "";
		REFL = DIFF = vec3(0);
		EXPO = IOR = MIRR = TRAN = EMIT = 0;
	}
	string NAME;
	vec3 DIFF, REFL;
	float EXPO, IOR, MIRR, TRAN, EMIT;
};
	