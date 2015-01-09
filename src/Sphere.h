#ifndef SPHERE_H
#define SPHERE_H

#include "Geometry.h"


class Sphere : public Geometry
{
public:
	Sphere();
	Sphere(float radius);
    Sphere(material *mat);
    virtual ~Sphere();

    virtual void buildGeometry();
private:
    glm::vec3 center_;
    float radius_;
	glm::vec3 color;
	material *mat;
protected:
	virtual Intersection intersectImpl(const Ray &ray) const;
};

#endif