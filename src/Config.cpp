#include "Config.h"

using namespace std;

Config::Config(string filename)
{
	ifstream file(filename.c_str());
	while(!file.eof()) {
		string tmp;
		file >> tmp;
		if (tmp == "CAMERA") {
			for(int i = 0; i < 6; i++) {
				file >> tmp;
				if (tmp == "RESO") {
					for (int i = 0; i < 2; i++) {
						file >> tmp;
						RESO[i] = atof(tmp.c_str());
					}
				} else if (tmp == "EYEP") {
					for (int i = 0; i < 3; i++) {
						file >> tmp;
						EYEP[i] = atof(tmp.c_str());
					}
				} else if (tmp == "VDIR") {
					for (int i = 0; i < 3; i++) {
						file >> tmp;
						VDIR[i] = atof(tmp.c_str());
					}
				} else if (tmp == "UVEC") {
					for (int i = 0; i < 3; i++) {
						file >> tmp;
						UVEC[i] = atof(tmp.c_str());
					}
				} else if (tmp == "FOVY") {
					file >> tmp;
					FOVY= atof(tmp.c_str());
				}
			}
		} else if (tmp == "LIGHT") {
			for (int i = 0; i < 2; i++) {
				file >> tmp;
				if (tmp == "LPOS") {
					for (int i = 0; i < 3; i++) {
						file >> tmp;
						LPOS[i] = atof(tmp.c_str());
					}
				} else if (tmp == "LCOL") {
					for (int i = 0; i < 3; i++) {
						file >> tmp;
						LCOL[i] = atof(tmp.c_str());
					}
				}
			}
		} else if (tmp == "NODE") {
			file >> NODE;
			for (int i = 0; i < 7; i++) {
				file >> tmp;
				if (tmp == "TRANSLATION") {
					for (int i = 0; i < 3; i++) {
						file >> tmp;
						TRANSLATION[i] = atof(tmp.c_str());
					}
				} else if (tmp == "ROTATION") {
					for (int i = 0; i < 3; i++) {
						file >> tmp;
						ROTATION[i] = atof(tmp.c_str());
					}
				} else if (tmp == "SCALE") {
					for (int i = 0; i < 3; i++) {
						file >> tmp;
						SCALE[i] = atof(tmp.c_str());
					}
				} else if (tmp == "CENTER") {
					for (int i = 0; i < 3; i++) {
						file >> tmp;
						CENTER[i] = atof(tmp.c_str());
					}
				} else if (tmp == "RGBA") {
					for (int i = 0; i < 3; i++) {
						file >> tmp;
						RGBA[i] = atof(tmp.c_str());
					}
				} else if (tmp == "PARENT") {
					file >> PARENT;
				} else if (tmp == "SHAPE") {
					file >> SHAPE;
				}
			}
		}
	}
}


Config::~Config(void)
{
}
