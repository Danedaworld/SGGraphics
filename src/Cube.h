#ifndef CUBE_H
#define CUBE_H

#include "Geometry.h"

class Cube : public Geometry
{
public:
	Cube();
    Cube(material* mat);
    virtual ~Cube();

    virtual void buildGeometry();
private:
    glm::vec3 center_;
    float apothem_;
	glm::vec3 color;
	material *mat;
protected:
	virtual Intersection intersectImpl(const Ray &ray) const;
};

#endif