// Sampling fuction stubs. Please fill in the TODOs to fit your scenegraph implementation
#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
static const float PI = 3.141592653589f;

// Returns a random point on a sphere
glm::vec3 getRandomPointOnSphere(Geometry* sphere) {
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
	return point;
}

// Returns a random point on a cube. Adapted from CIS 565
glm::vec3 getRandomPointOnCube(Geometry* cube) {
	// TODO: get the dimensions of the transformed cube in world space
	glm::vec3 dim (0, 0, 0);

	// Get surface area of the cube
	float side1 = dim[0] * dim[1];		// x-y
	float side2 = dim[1] * dim[2];		// y-z
	float side3 = dim[0] * dim[2];		// x-z
	float totalArea = 2.0f * (side1 + side2 + side3);	

	// pick random face weighted by surface area
	float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
	// pick 2 random components for the point in the range (-0.5, 0.5)
	float c1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX) - 0.5f;
	float c2 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX) - 0.5f;

	glm::vec3 point;
	if (r < side1 / totalArea) {				
		// x-y front
		point = glm::vec3(c1, c2, 0.5f);
	} else if (r < (side1 * 2) / totalArea) {
		// x-y back
		point = glm::vec3(c1, c2, -0.5f);
	} else if (r < (side1 * 2 + side2) / totalArea) {
		// y-z front
		point = glm::vec3(0.5f, c1, c2);
	} else if (r < (side1 * 2 + side2 * 2) / totalArea) {
		// y-z back
		point = glm::vec3(-0.5f, c1, c2);
	} else if (r < (side1 * 2 + side2 * 2 + side3) / totalArea) {
		// x-z front 
		point = glm::vec3(c1, 0.5f, c2);
	} else {
		// x-z back
		point = glm::vec3(c1, -0.5f, c2);
	}

	// TODO: transform point to world space
	return point;
}

// Given a normal vector, find a cosine weighted random direction in a hemisphere
// Adapted from CIS 565
glm::vec3 getCosineWeightedDirection(const glm::vec3& normal) {

	// Pick 2 random numbers in the range (0, 1)
	float xi1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
	float xi2 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

    float up = glm::sqrt(xi1); 			// cos(theta)
    float over = glm::sqrt(1 - up * up); // sin(theta)
    float around = xi2 * 2.0f * PI;
    
    // Find a direction that is not the normal based off of whether or not the normal's components 
    // are all equal to sqrt(1/3) or whether or not at least one component is less than sqrt(1/3).
    const float SQRT_OF_ONE_THIRD = glm::sqrt(1.0f/3.0f);
    glm::vec3 directionNotNormal;
    if (abs(normal.x) < SQRT_OF_ONE_THIRD) {
      directionNotNormal = glm::vec3(1.f, 0.f, 0.f);
    } else if (abs(normal.y) < SQRT_OF_ONE_THIRD) {
      directionNotNormal = glm::vec3(0.f, 1.f, 0.f);
    } else {
      directionNotNormal = glm::vec3(0.f, 0.f, 1.f);
    }
    
    //Use not-normal direction to generate two perpendicular directions
    glm::vec3 perpendicularDirection1 = glm::normalize(glm::cross(normal, directionNotNormal));
    glm::vec3 perpendicularDirection2 = glm::normalize(glm::cross(normal, perpendicularDirection1));
    
    return (up * normal) + (cos(around) * over * perpendicularDirection1) + (sin(around) * over * perpendicularDirection2);
}