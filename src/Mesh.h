#ifndef MESH_H
#define MESH_H

#include "Geometry.h"
#include <fstream>
#include <string>

class Mesh : public Geometry
{
public:
	Mesh(glm::vec3 p0, glm::vec3 p1, glm::vec3 p2);
    Mesh(material *mat, string filename);
    virtual ~Mesh();

    virtual void buildGeometry();
	static bool boundCheck(const Ray &ray, const float &radius);
	Intersection intersectTri(const Ray &ray, glm::vec3 p0, glm::vec3 p1, glm::vec3 p2) const;
private:
    glm::vec3 center_;
    glm::vec3 color;
	vector<glm::vec3> obj_vertices;
	vector<glm::vec3> obj_normals;
	vector<glm::vec3> obj_vertexpos;
	vector<glm::vec3> obj_texturepos;
	vector<glm::vec3> obj_normalpos;
	material *mat;
	//bounding sphere
	float radius;
	glm::vec3 center;

protected:
	virtual Intersection intersectImpl(const Ray &ray) const;
};

#endif