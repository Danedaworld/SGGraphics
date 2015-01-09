#include "Camera.h"

Camera::Camera() {

}

Camera::Camera(string f)
{
	filename = f;
	parseFile();
	//FOVY    /= 2;
	E        = EYEP;
	C        = VDIR;
	U        = UVEC;
	aspRatio = RESO.y / (float)RESO.x;
	phi      = FOVY/2;
	A        = cross(C,U);
	B        = cross(A,C);
	M        = E + C;
	theta    = degrees(glm::atan(RESO.x*glm::tan(radians(phi))/RESO.y));
	V        = B * (length(C) * glm::tan(radians(phi)) / length(B));
	H        = A * (length(C) * glm::tan(radians(theta)) / length(A));
}


Camera::~Camera(void)
{
}

void Camera::parseFile()
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
		}
	}
	file.close();
}