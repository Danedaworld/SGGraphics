#include "Cube.h"

Cube::Cube() :
	Geometry(CUBE),
	center_(glm::vec3(0.f,0.f,0.f)),
	apothem_(0.5f),
	mat(new material())
{
	buildGeometry();
}
// Create a unit cube
Cube::Cube(material* mat) :
	Geometry(CUBE),
	center_(glm::vec3(0.f,0.f,0.f)),
	apothem_(0.5f),
	mat(mat)
{
	buildGeometry();
}

Cube::~Cube() {}

void Cube::buildGeometry()
{
	vertices_.clear();
	colors_.clear();
	normals_.clear();
	indices_.clear();
	color = mat->DIFF;
	//printf("%f %f %f\n", color.x, color.y, color.z);
		// face 1
		vertices_.push_back(glm::vec3(0.5, 0.5, 0.5));
		vertices_.push_back(glm::vec3(0.5, -0.5, 0.5));
		vertices_.push_back(glm::vec3(-0.5, -0.5, 0.5));
		vertices_.push_back(glm::vec3(-0.5, 0.5, 0.5));
		// face 2
		vertices_.push_back(glm::vec3(0.5, -0.5, 0.5));
		vertices_.push_back(glm::vec3(0.5, 0.5, 0.5));
		vertices_.push_back(glm::vec3(0.5, 0.5, -0.5));
		vertices_.push_back(glm::vec3(0.5, -0.5, -0.5));
		// face 3
		vertices_.push_back(glm::vec3(0.5, 0.5, -0.5));
		vertices_.push_back(glm::vec3(0.5, -0.5, -0.5));
		vertices_.push_back(glm::vec3(-0.5, -0.5, -0.5));
		vertices_.push_back(glm::vec3(-0.5, 0.5, -0.5));
		// face 4
		vertices_.push_back(glm::vec3(-0.5, -0.5, -0.5));
		vertices_.push_back(glm::vec3(-0.5, 0.5, -0.5));
		vertices_.push_back(glm::vec3(-0.5, 0.5, 0.5));
		vertices_.push_back(glm::vec3(-0.5, -0.5, 0.5));
		// face 5
		vertices_.push_back(glm::vec3(-0.5, -0.5, 0.5));
		vertices_.push_back(glm::vec3(0.5, -0.5, 0.5));
		vertices_.push_back(glm::vec3(0.5, -0.5, -0.5));
		vertices_.push_back(glm::vec3(-0.5, -0.5, -0.5));
		// face 6
		vertices_.push_back(glm::vec3(-0.5, 0.5, 0.5));
		vertices_.push_back(glm::vec3(0.5, 0.5, 0.5));
		vertices_.push_back(glm::vec3(0.5, 0.5, -0.5));
		vertices_.push_back(glm::vec3(-0.5, 0.5, -0.5));
		// create 6 normals
	
		// face 1
		for (int i = 0; i < 4; i++)
		normals_.push_back(glm::vec3(0,1,0));
		
		// face 2
		for (int i = 0; i < 4; i++)
		normals_.push_back(glm::vec3(1,0,0));

		// face 3
		for (int i = 0; i < 4; i++)
		normals_.push_back(glm::vec3(0,0,-1));

		//face 4
		for (int i = 0; i < 4; i++)
		normals_.push_back(glm::vec3(-1,0,0));

		//face 5
		for (int i = 0; i < 4; i++)
		normals_.push_back(glm::vec3(0,-1,0));

		// face 6
		for (int i = 0; i < 4; i++)
		normals_.push_back(glm::vec3(0,1,0));

		//// face 1
		//for (int i = 0; i < 4; i++)
		//colors_.push_back(glm::vec3(0,0,1));
		//
		//// face 2
		//for (int i = 0; i < 4; i++)
		//colors_.push_back(glm::vec3(0,1,0));

		//// face 3
		//for (int i = 0; i < 4; i++)
		//colors_.push_back(glm::vec3(0,1,1));

		////face 4
		//for (int i = 0; i < 4; i++)
		//colors_.push_back(glm::vec3(1,1,0));

		////face 5
		//for (int i = 0; i < 4; i++)
		//colors_.push_back(glm::vec3(1,0,0));

		//// face 6
		//for (int i = 0; i < 4; i++)
		//colors_.push_back(glm::vec3(1,0,1));

	// indices and colors

	for (unsigned int i = 0; i < vertices_.size(); i++) {
		colors_.push_back(color);
	}

	unsigned int indices[36] = {
		0, 1, 2, 0, 2, 3,
		4, 5, 6, 4, 6, 7,
		8, 9, 10, 8, 10, 11,
		12, 13, 14, 12, 14, 15,
		16, 17, 18, 16, 18, 19,
		20, 21, 22, 20, 22, 23,
	};
	indices_ = vector<unsigned int>(indices, indices+sizeof(indices)/sizeof(indices[0]));
}

Intersection Cube::intersectImpl(const Ray &ray) const{
	Intersection isx;
	glm::vec3 bMin = -vertices_[0], bMax = vertices_[0];
	float tNear = -999999, tFar = 999999;
	if(abs(ray.dir.x)  == 0) {
		if(ray.orig.x < bMin.x || ray.orig.x > bMax.x) {
			isx.t = -1;
			return isx;
		}
	}
	float T1 = (bMin.x - ray.orig.x)/ray.dir.x;
	float T2 = (bMax.x - ray.orig.x)/ray.dir.x;

	if(T1 > T2){
		swap(T1,T2);
	}
	if(T1 > tNear) {
		tNear = T1;
	}
	if(T2 < tFar) {
		tFar = T2;
	}
	if(tNear > tFar) {
		isx.t = -1;
		return isx;
	}
	if(tFar < 0) {
		isx.t = -1;
		return isx;
	}

	if(abs(ray.dir.y) == 0) {
		if(ray.orig.y < bMin.y || ray.orig.y > bMax.y) {
			isx.t = -1;
			return isx;
		}
	}

	T1 = (bMin.y - ray.orig.y)/ray.dir.y;
	T2 = (bMax.y - ray.orig.y)/ray.dir.y;

	if(T1 > T2){
		swap(T1,T2);
	}
	if(T1 > tNear) {
		tNear = T1;
	}
	if(T2 < tFar) {
		tFar = T2;
	}
	if(tNear > tFar) {
		isx.t = -1;
		return isx;
	}
	if(tFar < 0) {
		isx.t = -1;
		return isx;
	}

	if(abs(ray.dir.z) == 0) {
		if(ray.orig.z < bMin.z || ray.orig.z > bMax.z) {
			isx.t = -1;
			return isx;
		}
	}

	T1 = (bMin.z - ray.orig.z)/ray.dir.z;
	T2 = (bMax.z - ray.orig.z)/ray.dir.z;

	if(T1 > T2){
		swap(T1,T2);
	}
	if(T1 > tNear) {
		tNear = T1;
	}
	if(T2 < tFar) {
		tFar = T2;
	}
	if(tNear > tFar) {
		isx.t = -1;
		return isx;
	}
	if(tFar < 0) {
		isx.t = -1;
		return isx;
	}
	if (tNear != 0) {
		isx.t = tNear;
	}
	else {
		isx.t = tFar;
	}
	glm::vec3 intersection = ray.orig + (ray.dir * (float) isx.t);
	for(int i = 0; i < 3; i++) {
		if(0.5 - abs(intersection[i]) < 0.001f) {
			isx.normal[i] = intersection[i];
		} else {
			isx.normal[i] = 0;
		}
	}
	//cout<<"NORMAL :"<<isx.normal.x<<" "<<isx.normal.y<<" "<<isx.normal.z<<" "<<isx.t<<endl;
	//std::cout<<isx.t<<endl;
	isx.m = mat;
	return isx;
}