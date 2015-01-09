#ifndef CYLINDER_H
#define CYLINDER_H

#include "Geometry.h"

class Cylinder : public Geometry
{
public:
	Cylinder();
    Cylinder(material *mat);
    virtual ~Cylinder();

    virtual void buildGeometry();
private:
    glm::vec3 center_;
    float radius_;
    float height_;
	glm::vec3 color;
	material *mat;
protected:
	virtual Intersection intersectImpl(const Ray &ray) const;
};

#endif