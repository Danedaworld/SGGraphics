#include "Geometry.h"

using namespace glm;
static const float PI = 3.141592653589f;

Geometry::Geometry(geometryType geomType) :
    type_(geomType)
{
}

Geometry::~Geometry()
{
    vertices_.clear();
    normals_.clear();
    colors_.clear();
    indices_.clear();
}

Intersection Geometry::intersect(const mat4 &T, Ray ray_world) const
{
    // The input ray here is in WORLD-space. It may not be normalized!

    // TODO: normalize ray_world.
	ray_world.dir = normalize(ray_world.dir);
    // Transform the ray into OBJECT-LOCAL-space, for intersection calculation.
    Ray ray_local;  // TODO: COMPUTE THIS AS FOLLOWS:
    // Transform the ray by the inverse transformation to get ray_local.
    // (Remember that position = vec4(vec3, 1) while direction = vec4(vec3, 0).)
	mat4 T_inv = inverse(T);
	ray_local.orig = vec3(T_inv * vec4(ray_world.orig,1));
	ray_local.dir = vec3(T_inv * vec4(ray_world.orig + ray_world.dir, 0) - T_inv*vec4(ray_world.orig,0));

    // Compute the intersection in LOCAL-space.
    Intersection isx = intersectImpl(ray_local);
	isx.g = (Geometry *) this;
    if (isx.t != -1) {
        // Transform the local-space intersection BACK into world-space.
        //     (Note that, as long as you didn't re-normalize the ray direction
        //     earlier, `t` doesn't need to change.)
        const vec3 normal_local = isx.normal;
		
        vec3 normal_world;  // TODO: COMPUTE THIS AS FOLLOWS:
        // Inverse-transpose-transform the normal to get it back from
        // local-space to world-space. (If you were transforming a position,
        // you would just use the unmodified transform T.)
        // http://www.arcsynthesis.org/gltut/Illumination/Tut09%20Normal%20Transformation.html
		normal_world = vec3(transpose(inverse(T))*vec4(normal_local,0));
		
        isx.normal = normal_world;
		
        // TODO: You might want to do this here: make sure here that your
        // normal is pointing the right way (toward, not away from, the ray
        // origin). Instead of doing this inside intersectImpl, you can do so
        // here by just negating normal_world depending on the sign of
        //     dot(normal_world, ray_world.dir).
		
		isx.normal *= dot(normal_world, ray_world.dir) < 0 ? 1 : -1;
		isx.normal = normalize(isx.normal);
    }

    // The final output intersection data is in WORLD-space.
	/*if(isx.t>0)
	std::cout<<"Final: "<<isx.t<<endl;*/
    return isx;
}

// Returns a random point on a sphere
glm::vec3 Geometry::getRandomPointOnSphere(glm::mat4 matrix) {
	// generate u, v, in the range (0, 1)
	float u = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
	float v = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

	float theta = 2.0f * PI * u;
	float phi = acos(2.0f * v - 1.0f);

	// find x, y, z coordinates assuming unit sphere in object space
	glm::vec3 point;
	point[0] = sin(phi) * cos(theta);
	point[1] = sin(phi) * sin(theta);
	point[2] = cos(phi);

	// TODO: transform point to world space
	point = vec3(matrix * vec4(point, 1));
	return point;
}

// Returns a random point on a cube. Adapted from CIS 565
glm::vec3 Geometry::getRandomPointOnCube(glm::mat4 matrix) {
	// TODO: get the dimensions of the transformed cube in world space
	//matrix = inverse(matrix);
	//matrix = transpose(matrix);
	glm::vec3 p1 = vec3(matrix * vec4(0, 0, 0, 1));
	glm::vec3 p2 = vec3(matrix * vec4(1, 1, 1, 1));
	vec3 dim = vec3(abs(p2.x-p1.x), abs(p2.y-p1.y),abs(p2.z-p1.z));
	float side1 = dim.x * dim.y;		// x-y
	float side2 = dim.y * dim.z;		// y-z
	float side3 = dim.x * dim.z;		// x-z
	float totalArea = 2.0f * (side1 + side2 + side3);
	//printf("SIDE SIZES: %f %f %f\n", side1, side2, side3);
	// pick random face weighted by surface area
	float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
	// pick 2 random components for the point in the range (-0.5, 0.5)
	float c1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX)-0.5f;
	float c2 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX)-0.5f;
	float c3 = 0.5f;
	glm::vec3 point;
	if (r < side1 / totalArea) {
		// x-y front
		point = glm::vec3(c1, c2, c3);
	}
	else if (r < (side1 * 2) / totalArea) {
		// x-y back
		point = glm::vec3(c1, c2, -c3);
	}
	else if (r < (side1 * 2 + side2) / totalArea) {
		// y-z front
		point = glm::vec3(c3, c1, c2);
	}
	else if (r < (side1 * 2 + side2 * 2) / totalArea) {
		// y-z back
		point = glm::vec3(-c3, c1, c2);
	}
	else if (r < (side1 * 2 + side2 * 2 + side3) / totalArea) {
		// x-z front 
		point = glm::vec3(c1, c3, c2);
	}
	else {
		// x-z back
		point = glm::vec3(c1, -c3, c2);
	}
	point = vec3(matrix * vec4(point, 1)); 

	// TODO: transform point to world space
	return point;
}

// Given a normal vector, find a cosine weighted random direction in a hemisphere
// Adapted from CIS 565
glm::vec3 Geometry::getCosineWeightedDirection(const glm::vec3& normal) {

	// Pick 2 random numbers in the range (0, 1)
	float xi1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
	float xi2 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

	float up = glm::sqrt(xi1); 			// cos(theta)
	float over = glm::sqrt(1 - up * up); // sin(theta)
	float around = xi2 * 2.0f * PI;

	// Find a direction that is not the normal based off of whether or not the normal's components 
	// are all equal to sqrt(1/3) or whether or not at least one component is less than sqrt(1/3).
	const float SQRT_OF_ONE_THIRD = glm::sqrt(1.0f / 3.0f);
	glm::vec3 directionNotNormal;
	if (abs(normal.x) < SQRT_OF_ONE_THIRD) {
		directionNotNormal = glm::vec3(1.f, 0.f, 0.f);
	}
	else if (abs(normal.y) < SQRT_OF_ONE_THIRD) {
		directionNotNormal = glm::vec3(0.f, 1.f, 0.f);
	}
	else {
		directionNotNormal = glm::vec3(0.f, 0.f, 1.f);
	}

	//Use not-normal direction to generate two perpendicular directions
	glm::vec3 perpendicularDirection1 = glm::normalize(glm::cross(normal, directionNotNormal));
	glm::vec3 perpendicularDirection2 = glm::normalize(glm::cross(normal, perpendicularDirection1));

	return (up * normal) + (cos(around) * over * perpendicularDirection1) + (sin(around) * over * perpendicularDirection2);
}