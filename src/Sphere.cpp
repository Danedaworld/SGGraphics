#include "Sphere.h"

static const float PI = 3.141592653589f;
Sphere::Sphere() :
	Geometry(SPHERE),
    center_(glm::vec3(0.f, 0.f, 0.f)),
    radius_(1.f),
	color(glm::vec3(1.f,1.f,1.f))
{
    buildGeometry();
}
// Creates a unit sphere.
Sphere::Sphere(float radius) :
Geometry(SPHERE),
center_(glm::vec3(0.f, 0.f, 0.f)),
radius_(radius),
color(glm::vec3(1.f, 1.f, 1.f))
{
}

Sphere::Sphere(material *mat) :
    Geometry(SPHERE),
    center_(glm::vec3(0.f, 0.f, 0.f)),
    radius_(1.f),
	mat(mat)
{
    buildGeometry();
}

Sphere::~Sphere() {}

void Sphere::buildGeometry()
{
    vertices_.clear();
    colors_.clear();
    normals_.clear();
    indices_.clear();
	color = mat->DIFF;
    // Find vertex positions for the sphere.
    unsigned int subdiv_axis = 16;      // vertical slices
    unsigned int subdiv_height = 16;        // horizontal slices
    float dphi = PI / subdiv_height;
    float dtheta = 2.0f * PI / subdiv_axis;
    float epsilon = 0.0001f;

    // North pole
    glm::vec3 point (0.0f, 1.0f, 0.0f);
    normals_.push_back(point);
    // scale by radius_ and translate by center_
    vertices_.push_back(center_ + radius_ * point);

    for (float phi = dphi; phi < PI; phi += dphi) {
        for (float theta = dtheta; theta <= 2.0f * PI + epsilon; theta += dtheta) {
            float sin_phi = sin(phi);

            point[0] = sin_phi * sin(theta);
            point[1] = cos(phi);
            point[2] = sin_phi * cos(theta);

            normals_.push_back(point);
            vertices_.push_back(center_ + radius_ * point);
        }
    }
    // South pole
    point = glm::vec3(0.0f, -1.0f, 0.0f);
    normals_.push_back(point);
    vertices_.push_back(center_ + radius_ * point);

    // fill in index array.
    // top cap
    for (unsigned int i = 0; i < subdiv_axis - 1; ++i) {
        indices_.push_back(0);
        indices_.push_back(i + 1);
        indices_.push_back(i + 2);
    }
    indices_.push_back(0);
    indices_.push_back(subdiv_axis);
    indices_.push_back(1);

    // middle subdivs
    unsigned int index = 1;
    for (unsigned int i = 0; i < subdiv_height - 2; i++) {
        for (unsigned int j = 0; j < subdiv_axis - 1; j++) {
            // first triangle
            indices_.push_back(index);
            indices_.push_back(index + subdiv_axis);
            indices_.push_back(index + subdiv_axis + 1);

            // second triangle
            indices_.push_back(index);
            indices_.push_back(index + subdiv_axis + 1);
            indices_.push_back(index + 1);

            index++;
        }
        // reuse vertices from start and end point of subdiv_axis slice
        indices_.push_back(index);
        indices_.push_back(index + subdiv_axis);
        indices_.push_back(index + 1);

        indices_.push_back(index);
        indices_.push_back(index + 1);
        indices_.push_back(index + 1 - subdiv_axis);

        index++;
    }

    // end cap
    unsigned int bottom = (subdiv_height - 1) * subdiv_axis + 1;
    unsigned int offset = bottom - subdiv_axis;
    for (unsigned int i = 0; i < subdiv_axis - 1 ; ++i) {
        indices_.push_back(bottom);
        indices_.push_back(i + offset);
        indices_.push_back(i + offset + 1);
    }
    indices_.push_back(bottom);
    indices_.push_back(bottom - 1);
    indices_.push_back(offset);

    // colors
    for (unsigned int i = 0; i < vertices_.size(); ++i) {
        colors_.push_back(color);
    }
}

Intersection Sphere::intersectImpl(const Ray &ray) const{
	Intersection isx;

	float a = glm::dot(ray.dir, ray.dir);
	float b = glm::dot(2.0f*ray.orig, ray.dir);
	float c = glm::dot(ray.orig, ray.orig) - radius_;
	float d = b * b - 4 * a * c;
	if(d >= 0) {
		float t0 = (-b - glm::sqrt(d))/(2*a);
		float t1 = (-b + glm::sqrt(d))/(2*a);
		if(t0 <= 0) t0 = 1e37;
		if(t1 <= 0) t0 = 1e37;
		isx.t = glm::min(t0,t1);
	} else {
		isx.t = -1;
		return isx;
	}

	isx.normal = (ray.orig + (ray.dir * (float) isx.t));
	isx.m = mat;
	return isx;
}