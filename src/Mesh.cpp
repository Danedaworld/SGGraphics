#include "Mesh.h"
#include <iostream>
// Create a unit cube

Mesh::Mesh(glm::vec3 p0, glm::vec3 p1, glm::vec3 p2) :
	Geometry(MESH),
    center_(glm::vec3(0.f,0.f,0.f)),
    color(glm::vec3(1.f,1.f,1.f))
	{
		vertices_.push_back(p0);
		vertices_.push_back(p1);
		vertices_.push_back(p2);
		glm::vec3 n = glm::cross(p0,p1);
		normals_.push_back(n);
		normals_.push_back(n);
		normals_.push_back(n);
		indices_.push_back(1);
		indices_.push_back(2);
		indices_.push_back(3);
		colors_.push_back(glm::vec3(1.f,1.f,1.f));
		colors_.push_back(glm::vec3(1.f,1.f,1.f));
		colors_.push_back(glm::vec3(1.f,1.f,1.f));
	}

Mesh::Mesh(material *mat, string filename) :
    Geometry(MESH),
    center_(glm::vec3(0.f,0.f,0.f)),
    mat(mat)
{
	color = mat->DIFF;
	obj_vertices.clear();
	//obj_vertices.push_back(glm::vec3(0,0,0)); //indices start at 1 so i'm sticking this in as a placeholder
	
	ifstream in(filename.c_str());
	string tmp, tx, ty, tz;
	in >> tmp;
	while(!in.eof()) {
		in >> tmp;
		if(tmp == "v") {
			in >> tx >> ty >> tz;
			obj_vertices.push_back(glm::vec3(atof(tx.c_str()),atof(ty.c_str()),atof(tz.c_str())));
		} else if (tmp == "vn") {
			in >> tx >> ty >> tz;
			obj_normals.push_back(glm::vec3(atof(tx.c_str()),atof(ty.c_str()),atof(tz.c_str())));
		} else if (tmp == "f") {
			in >> tx >> ty >> tz;
			int index;
			if(tx.find("/") != -1) {
				obj_vertexpos.push_back(glm::vec3(atof(tx.substr(0,tx.find("/")).c_str()), 
												  atof(ty.substr(0,ty.find("/")).c_str()), 
												  atof(tz.substr(0,tz.find("/")).c_str())));
				tx = tx.substr(tx.find("/")+1);
				ty = tx.substr(tx.find("/")+1);
				tz = tx.substr(tx.find("/")+1);
				if(tx.find("/") != -1) {
					obj_texturepos.push_back(glm::vec3(atof(tx.substr(0,tx.find("/")).c_str()), 
													   atof(ty.substr(0,ty.find("/")).c_str()), 
													   atof(tz.substr(0,tz.find("/")).c_str())));
					tx = tx.substr(tx.find("/")+1);
					ty = tx.substr(tx.find("/")+1);
					tz = tx.substr(tx.find("/")+1);
					if(tx.find("/") != -1) {
						obj_normalpos.push_back(glm::vec3(atof(tx.substr(0,tx.find("/")).c_str()), 
														  atof(ty.substr(0,ty.find("/")).c_str()), 
														  atof(tz.substr(0,tz.find("/")).c_str())));
					}
				}
			}
		}
	}
    buildGeometry();
}

Mesh::~Mesh() {}

void Mesh::buildGeometry()
{
	vertices_.clear();
	colors_.clear();
	normals_.clear();
	indices_.clear();

	//std::cout<<"Mesh "<<normals_.size()<<" "<<vertices_.size()<<" "<<obj_vertexpos.size()<<endl;
	if(!obj_normals.empty()) {
		normals_ = obj_normals;
	} else {
		for (unsigned int i = 0; i < obj_vertexpos.size(); i++) {

			//vertex
			vertices_.push_back(obj_vertices[obj_vertexpos[i].x-1]);
			vertices_.push_back(obj_vertices[obj_vertexpos[i].y-1]);
			vertices_.push_back(obj_vertices[obj_vertexpos[i].z-1]);

			glm::vec3 v1,v2;
			v1 = obj_vertices[obj_vertexpos[i].x-1] - obj_vertices[obj_vertexpos[i].y-1];
			v2 = obj_vertices[obj_vertexpos[i].y-1] - obj_vertices[obj_vertexpos[i].z-1];
			glm::vec3 cross = glm::normalize(glm::cross(v1,v2));

			normals_.push_back(cross);
			normals_.push_back(cross);
			normals_.push_back(cross);

			colors_.push_back(color);
			colors_.push_back(color);
			colors_.push_back(color);
		}
	}
	//std::cout<<"Mesh "<<normals_.size()<<" "<<vertices_.size()<<" "<<obj_vertexpos.size()<<endl;

	//indices
	for(unsigned int i = 0; i < vertices_.size(); i++) {
		indices_.push_back(i);
	}

	//bounding sphere
	radius = 0;
	for (int i = 0; i < vertices_.size(); i++) {
		float len = length(vertices_[i]);
		if (len > radius) {
			radius = len;
		}
	}
	radius *= 2;
	printf("RADIUS: %f\n", radius);
}

bool Mesh::boundCheck(const Ray &ray, const float &radius) {
	float a = glm::dot(ray.dir, ray.dir);
	float b = glm::dot(2.0f*ray.orig, ray.dir);
	float c = glm::dot(ray.orig, ray.orig) - radius;
	float d = b * b - 4 * a * c;
	if (d >= 0) {
		return true;
	}
	else {
		return false;
	}
}

Intersection Mesh::intersectImpl(const Ray &ray) const{
	Intersection isx; isx.t = -1;
	if (boundCheck(ray, radius)) {
		for (int i = 0; i < vertices_.size() - 2; i += 3) {
			Intersection isx_tmp = intersectTri(ray, vertices_[i], vertices_[i + 1], vertices_[i + 2]);
			if (isx.t == -1) {
				isx = isx_tmp;
			}
			if (isx_tmp.t != -1 && isx_tmp.t < isx.t) {
				isx = isx_tmp;
			}
			//std::cout<<i<<endl;
		}
		isx.m = mat;
	}
	return isx;
}
#define EPSILON 0.0001
Intersection Mesh::intersectTri(const Ray &ray, glm::vec3 p0, glm::vec3 p1, glm::vec3 p2) const{
	Intersection isx;
	//moller-trumbore black magic wtf-wizardry
	glm::vec3 e1, e2, P, Q, T;
	float det, inv_det, u, v, t;
	isx.t = -1;
	e1 = p1-p0;
	e2 = p2-p0;
	P = glm::cross(ray.dir, e2);
	det = glm::dot(e1, P);
	if(det > -EPSILON && det < EPSILON) {
		return isx;
	}
	inv_det = 1.f/det;
	T = ray.orig - p0;
	u = glm::dot(T,P) * inv_det;
	if(u < 0.f || u > 1.f) {
		return isx;
	}
	Q = glm::cross(T,e1);
	v = glm::dot(ray.dir,Q) * inv_det;
	if(v < 0.f || u + v > 1.f) {
		return isx;
	}
	t = glm::dot(e2,Q) * inv_det;
	if(t > EPSILON) {
		isx.t = t;
		isx.normal = glm::cross(e1,e2);
		return isx;
	}
	
	return isx;
}