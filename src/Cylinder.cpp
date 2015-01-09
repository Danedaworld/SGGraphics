#include "Cylinder.h"

static const float PI = 3.141592653589f;

Cylinder::Cylinder() :
	Geometry(CYLINDER),
	center_(glm::vec3(0.f,0.f,0.f)),
	radius_(0.5f),
    height_(1.0f),
	mat(new material())
{
	buildGeometry();

}
// Creates a unit cylinder centered at (0, 0, 0)
Cylinder::Cylinder(material *mat) :
    Geometry(CYLINDER),
    center_(glm::vec3(0.f, 0.f, 0.f)),
    radius_(0.5f),
    height_(1.0f),
	mat(mat)
{
    buildGeometry();
}

Cylinder::~Cylinder() {}

void Cylinder::buildGeometry()
{
    vertices_.clear();
    colors_.clear();
    normals_.clear();
    indices_.clear();
	color = mat->DIFF; //temp
    unsigned short subdiv = 20;
    float dtheta = 2 * PI / subdiv;

    glm::vec4 point_top(0.0f, 0.5f * height_, radius_, 1.0f),
        point_bottom (0.0f, -0.5f * height_, radius_, 1.0f);
    vector<glm::vec3> cap_top, cap_bottom;

    // top and bottom cap vertices
    for (int i = 0; i < subdiv + 1; ++i) {
        glm::mat4 rotate = glm::rotate(glm::mat4(1.0f), i * dtheta, glm::vec3(0.f, 1.f, 0.f));
        glm::mat4 translate = glm::translate(glm::mat4(1.0f), center_);

        cap_top.push_back(glm::vec3(translate * rotate * point_top));
        cap_bottom.push_back(glm::vec3(translate * rotate * point_bottom));
    }

    //Create top cap.
    for ( int i = 0; i < subdiv - 2; i++) {
        vertices_.push_back(cap_top[0]);
        vertices_.push_back(cap_top[i + 1]);
        vertices_.push_back(cap_top[i + 2]);
    }
    //Create bottom cap.
    for (int i = 0; i < subdiv - 2; i++) {
        vertices_.push_back(cap_bottom[0]);
        vertices_.push_back(cap_bottom[i + 1]);
        vertices_.push_back(cap_bottom[i + 2]);
    }
    //Create barrel
    for (int i = 0; i < subdiv; i++) {
        //Right-side up triangle
        vertices_.push_back(cap_top[i]);
        vertices_.push_back(cap_bottom[i + 1]);
        vertices_.push_back(cap_bottom[i]);
        //Upside-down triangle
        vertices_.push_back(cap_top[i]);
        vertices_.push_back(cap_top[i + 1]);
        vertices_.push_back(cap_bottom[i + 1]);
    }

    // create normals
    glm::vec3 top_centerpoint(0.0f , 0.5f * height_ , 0.0f),
        bottom_centerpoint(0.0f, -0.5f * height_, 0.0f);
    glm::vec3 normal(0, 1, 0);

    // Create top cap.
    for (int i = 0; i < subdiv - 2; i++) {
        normals_.push_back(normal);
        normals_.push_back(normal);
        normals_.push_back(normal);
    }
    // Create bottom cap.
    for (int i = 0; i < subdiv - 2; i++) {
        normals_.push_back(-normal);
        normals_.push_back(-normal);
        normals_.push_back(-normal);
    }

    // Create barrel
    for (int i = 0; i < subdiv; i++) {
        //Right-side up triangle
        normals_.push_back(glm::normalize(cap_top[i] - top_centerpoint));
        normals_.push_back(glm::normalize(cap_bottom[i + 1] - bottom_centerpoint));
        normals_.push_back(glm::normalize(cap_bottom[i] - bottom_centerpoint));
        //Upside-down triangle
        normals_.push_back(glm::normalize(cap_top[i] - top_centerpoint));
        normals_.push_back(glm::normalize(cap_top[i + 1] - top_centerpoint));
        normals_.push_back(glm::normalize(cap_bottom[i + 1] - bottom_centerpoint));
    }

    // indices and colors

    for (unsigned int i = 0; i < vertices_.size(); ++i) {
        colors_.push_back(color);
    }

    for (unsigned int i = 0; i < vertices_.size(); ++i) {
        indices_.push_back(i);
    }
}

Intersection Cylinder::intersectImpl(const Ray &ray) const{
	Intersection isx;
	float a = ray.dir.x * ray.dir.x + ray.dir.z * ray.dir.z;
	float b = 2 * ray.dir.x * ray.orig.x + 2 * ray.dir.z * ray.orig.z;
	float c = ray.orig.x * ray.orig.x + ray.orig.z * ray.orig.z - 0.25f;
	float d = b * b - 4 * a * c;
	//std::cout<<d<<endl;
	float F_MAX = 1e23;
	float t0 = F_MAX, t1 = F_MAX, t2 = F_MAX, t3 = F_MAX;
	if(d >= 0 && a != 0) {
		t0 = (-b - glm::sqrt(d))/(2*a);
		t1 = (-b + glm::sqrt(d))/(2*a);
	}
	//std::cout<<t0<<" "<<t1<<endl;
	//check if it's within caps
	glm::vec3 intersect0 = ray.orig + ray.dir * (float) t0;
	glm::vec3 intersect1 = ray.orig + ray.dir * (float) t1;
	float y0 = intersect0.y;
	float y1 = intersect1.y;
	float ymin = -0.5f, ymax = 0.5f;
	
	if(y0 < ymin || y0 > ymax) {
		t0 = F_MAX;
	}
	
	if(y1 < ymin || y0 > ymax) {
		t1 = F_MAX;
	}
	if(abs(ray.dir.y) > 0) {
		t2 = (ymin - ray.orig.y)/ray.dir.y;
		if(glm::distance(ray.orig + ray.dir * (float) t2, glm::vec3(0,-0.5,0)) > radius_) {
			t2 = F_MAX;
		}
		t3 = (ymax - ray.orig.y)/ray.dir.y;
		if(glm::distance(ray.orig + ray.dir * (float) t3, glm::vec3(0,0.5,0)) > radius_) {
			t3 = F_MAX;
		}
	}
	if(t0 <= 0) t0 = F_MAX;
	if(t1 <= 0) t1 = F_MAX;
	if(t2 <= 0) t2 = F_MAX;
	if(t3 <= 0) t3 = F_MAX;
	isx.t = glm::min(t0,glm::min(t1,glm::min(t2,t3)));
	if(isx.t == F_MAX) isx.t = -1;
	glm::vec3 intersect = ray.orig + ray.dir * (float) isx.t;
	if(0.5f - abs(intersect.y) < 0.01f) {
		isx.normal = -glm::vec3(0,0.5,0); //I dont have to check this here
	} else {
		isx.normal = intersect * glm::vec3(1,0,1);
	}
	//std::cout<<isx.t<<endl;
	isx.m = mat;
	return isx;
}