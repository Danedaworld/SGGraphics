// NOTE: This definition forces GLM to use radians (not degrees) for ALL of its
// angle arguments. The documentation may not always reflect this fact.
// YOU SHOULD USE THIS IN ALL FILES YOU CREATE WHICH INCLUDE GLM
#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include "glew/glew.h"
#include <GL/glut.h>

#include "EasyBMP/EasyBMP.h"

#include <fstream>
#include <iostream>
#include <string>
#include <time.h>
#include <stack>
#include <random>
#include <chrono>

#include "Geometry.h"
#include "Sphere.h"
#include "Cube.h"
#include "Cylinder.h"
#include "Mesh.h"
#include "SceneGraph.h"
#include "Camera.h"
#include "tests.h"
#include "Photon.h"
#include "KDTree.h"
#include <omp.h>
//static const float PI = 3.141592653589f;

#define PI 3.141592653589793

// Vertex arrays needed for drawing
unsigned int *vboIdx;
unsigned int *vboConsolidate;

// Attributes
unsigned int locationPos;
unsigned int locationCol;
unsigned int locationNor;

// Uniforms
unsigned int unifModel;
unsigned int unifModelInvTr;
unsigned int unifViewProj;
unsigned int unifLightPos;
unsigned int unifLightColor;
unsigned int unifCameraPos;
// Needed to compile and link and use the shaders
unsigned int shaderProgram;

// Window dimensions, change if you want a bigger or smaller window
unsigned int windowWidth = 1600;
unsigned int windowHeight = 900;

// Animation/transformation stuff
clock_t old_time;
float rotation = 0.0f;

// TESTING GEOMETRY
Geometry* geometry;
void sampleDrawSphere(glm::mat4 model);

// Helper function to read shader source and put it in a char array
// thanks to Swiftless
std::string textFileRead(const char*);

// Some other helper functions from CIS 565 and CIS 277
void printLinkInfoLog(int);
void printShaderInfoLog(int);
void printGLErrorLog();

// Standard glut-based program functions
void init(void);
void resize(int, int);
void display(void);
void keypress(unsigned char, int, int);
void mousepress(int button, int state, int x, int y);
void cleanup(void);

void initShader();
void cleanupShader();

// modified from CIS 277
void uploadPrimitive(node *n, Geometry *g);
void drawPrimitive(glm::mat4 model, node *n, Geometry *g);

//make the scene graph and camera
SceneGraph s;
Camera c;
void renderScene(SceneGraph &s);
void uploadScene(node* n);
node* currentNode;
stack<node*> nodeStack;
void nextNode();
void freeMem();
void getLight(node* n, glm::mat4 matrix);
//raycasting
void getRayIntersections(Ray &ray, node *n, glm::mat4 matrix, vector<Intersection> &isxs);
void generateImage(SceneGraph s);
Intersection RayIntersect(Ray& ray, node* n);
bool shadowRayUnblocked(vec3& p1, vec3& p2, Intersection& orig_isx);
void traceRay(Ray &ray, int depth, glm::vec3& color, bool inside, int i, int j);

node *lightNode;
//monte carlo
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);
std::uniform_real_distribution<double> distribution(-1.0, 1.0);
bool isAbsorbed(float absorb);
vec3 getCosineWeightedDirection(const glm::vec3& normal);
vec3 MonteCarlo(Ray &ray, int depth, glm::vec3& transmittance, bool inside, int hitcount, int canGenerate);
vector<vector<float> > weightmap;
int extern_i, extern_j;


using namespace std;

int main(int argc, char** argv)
{
	s = SceneGraph(argv[1]);
	c = Camera(argv[1]);
	windowHeight = c.RESO.y;
	windowWidth = c.RESO.x;
	glutInit(&argc, argv);
	// Use RGBA double buffered window
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowSize(windowWidth, windowHeight);
	glutCreateWindow("Scene Graph");

	glewInit();

	vboIdx = new unsigned int[s.numNodes];
	vboConsolidate = new unsigned int[s.numNodes];
	init();
	currentNode = s.root;
	nodeStack.push(s.root);
	//RunTests(); //previous assignment
	getLight(s.root, glm::mat4());
	uploadScene(s.root);

	glutDisplayFunc(display);
	glutReshapeFunc(resize);
	glutKeyboardFunc(keypress);
	glutMouseFunc(mousepress);
	glutIdleFunc(display);

	glutMainLoop();

	freeMem();
	return 0;
}

/*==========================================================
**********************MAINTENANCE***************************
============================================================*/

void freeMem() {
	s.release();
}

void init()
{

	// Set the color which clears the screen between frames
	glClearColor(0, 0, 0, 1);
	// Enable and clear the depth buffer
	glEnable(GL_DEPTH_TEST);
	glClearDepth(1.0);
	glDepthFunc(GL_LEQUAL);

	// Set up our shaders here
	initShader();

	// example for uploading data for drawing.
	//uploadPrimitive();

	resize(windowWidth, windowHeight);
	old_time = clock();
}

void initShader()
{
	// Read in the shader program source files
	std::string vertSourceS = textFileRead("shaders/diff.vert.glsl");
	const char *vertSource = vertSourceS.c_str();
	std::string fragSourceS = textFileRead("shaders/diff.frag.glsl");
	const char *fragSource = fragSourceS.c_str();

	// Tell the GPU to create new shaders and a shader program
	GLuint shadVert = glCreateShader(GL_VERTEX_SHADER);
	GLuint shadFrag = glCreateShader(GL_FRAGMENT_SHADER);
	shaderProgram = glCreateProgram();

	// Load and compiler each shader program
	// Then check to make sure the shaders complied correctly
	// - Vertex shader
	glShaderSource(shadVert, 1, &vertSource, NULL);
	glCompileShader(shadVert);
	printShaderInfoLog(shadVert);
	// - Diffuse fragment shader
	glShaderSource(shadFrag, 1, &fragSource, NULL);
	glCompileShader(shadFrag);
	printShaderInfoLog(shadFrag);

	// Link the shader programs together from compiled bits
	glAttachShader(shaderProgram, shadVert);
	glAttachShader(shaderProgram, shadFrag);
	glLinkProgram(shaderProgram);
	printLinkInfoLog(shaderProgram);

	// Clean up the shaders now that they are linked
	glDetachShader(shaderProgram, shadVert);
	glDetachShader(shaderProgram, shadFrag);
	glDeleteShader(shadVert);
	glDeleteShader(shadFrag);

	// Find out what the GLSL locations are, since we can't pre-define these
	locationPos = glGetAttribLocation(shaderProgram, "vs_Position");
	locationNor = glGetAttribLocation(shaderProgram, "vs_Normal");
	locationCol = glGetAttribLocation(shaderProgram, "vs_Color");
	unifViewProj = glGetUniformLocation(shaderProgram, "u_ViewProj");
	unifModel = glGetUniformLocation(shaderProgram, "u_Model");
	unifModelInvTr = glGetUniformLocation(shaderProgram, "u_ModelInvTr");
	unifLightPos = glGetUniformLocation(shaderProgram, "u_LightPos");
	unifLightColor = glGetUniformLocation(shaderProgram, "u_LightColor");
	unifCameraPos = glGetUniformLocation(shaderProgram, "u_CameraPos");

	printGLErrorLog();
}

void cleanup()
{
	//glDeleteBuffers(1, &vboIdx);
	//glDeleteBuffers(1, &vboConsolidate);

	glDeleteProgram(shaderProgram);

	delete geometry;
}

/*==========================================================
**********************RAY TRACING***************************
============================================================*/

Intersection RayIntersect(Ray& ray, node* n) {
	vector<Intersection> isxs;
	getRayIntersections(ray, s.root, glm::mat4(), isxs);
	Intersection isx; isx.normal = glm::vec3(0, 0, 0);
	isx.t = -1;
	if (isxs.size()) {
		isx = isxs[0];
		for (int i = 1; i < isxs.size(); i++) {
			if (isxs[i].t < isx.t) {
				isx = isxs[i];
			}
		}
	}

	return isx;
}

void getRayIntersections(Ray &ray, node *n, glm::mat4 matrix, vector<Intersection> &isxs) {
	matrix = glm::translate(matrix, n->TRANSLATION);
	matrix = glm::translate(matrix, n->CENTER);
	matrix = glm::rotate(matrix, radians(n->ROTATION.x), vec3(1, 0, 0));
	matrix = glm::rotate(matrix, radians(n->ROTATION.y), vec3(0, 1, 0));
	matrix = glm::rotate(matrix, radians(n->ROTATION.z), vec3(0, 0, 1));
	matrix = glm::scale(matrix, n->SCALE);
	matrix = glm::translate(matrix, vec3(0, 0, 0) - n->CENTER);
	Intersection isx; isx.t = -1;
	if (n->geometry) {
		isx = n->geometry->intersect(matrix, ray);
		if (isx.t != -1 && isx.t > 0) {
			isxs.push_back(isx);
		}
	}
	for (int i = 0; i < n->children.size(); i++) {
		getRayIntersections(ray, n->children[i], matrix, isxs);
	}
}

float shadowFeelers(vec3& p1, Intersection& orig_isx) {
	vector<Intersection> isxs;
	Intersection isx;
	vec3 p2, dir;
	int sampleSize = 30;
	int count = 0;

	//mat4 transformation = getTransformationMatrix(s.root, mat4(), s.lightnode->geometry);
	mat4 transformation = *lightNode->matrix;

	for (int i = 0; i < sampleSize; i++) {
		if (lightNode->geometry->getGeometryType() == Geometry::SPHERE) {
			p2 = Geometry::getRandomPointOnSphere(transformation);
		}
		else {
			p2 = Geometry::getRandomPointOnCube(transformation);
		}
		vec3 v = p2 - p1;
		vec3 dir = normalize(v);
		float light_t = length(v);
		float epsilon = 0.001f;
		Ray r(p1 + epsilon * dir, dir);
		getRayIntersections(r, s.root, glm::mat4(), isxs);

		if (isxs.size()) {
 			isx = isxs[0];
		}
		//printf("%d ", isxs.size());
		for (int j = 0; j < isxs.size(); j++) {
			if (isxs[j].t < isx.t && orig_isx.g != isxs[j].g) {
				isx = isxs[j];
			}
		}
		if (isx.g == lightNode->geometry) {
			count++;
		}

		isxs.clear();
	}
	

	float retVal = count / (float)sampleSize;
	//if (retVal < 0.15f) retVal = 0.15f;
	return retVal;
}

void getLight() {
	for (int i = 0; i < s.nodelist.size(); i++) {
		if (s.nodelist[i]->MAT)
			if (s.nodelist[i]->MAT->EMIT) {
				lightNode = s.nodelist[i];
			}
	}
}

void traceRay(Ray& ray, int depth, glm::vec3& color, bool inside, int i, int j) {
	Intersection isx;
	glm::vec3 spec, ReflectedColor, refr, RefractedColor, trans, TransmittedColor;
	float specularTerm = 0, diffuseTerm = 0;
	float epsilon = 0.001f;
	glm::vec3 ambient(0.1, 0.1, 0.1);
	if (depth > 10) { //we're done and found nothing
		return;
	}
	isx = RayIntersect(ray, s.root);

	if (isx.t == -1) {//no intersection
		color = glm::vec3(0, 0, 0);
		return;
	}
	if (isx.m->EMIT) {
		color = vec3(1, 1, 1);
		weightmap[i][j] = 1.f;
		return;
	}
	else { //there is an intersection
		float schlick = 1;
		vec3 lightVector = (s.LPOS - (ray.orig + ray.dir * (float)isx.t));
		if (isx.m->MIRR) {// reflect
			vec3 new_dir = normalize(glm::reflect(ray.dir, isx.normal));
			vec3 new_orig = ray.orig + ray.dir * (float)isx.t + (new_dir * epsilon);
			traceRay(Ray(new_orig, new_dir), depth + 1, ReflectedColor, inside, i, j);
			if (isx.m->TRAN) {
				float r0 = pow(((isx.m->IOR - 1) / (isx.m->IOR + 1)), 2);
				schlick = r0 + (1 - r0) * pow((1 - dot(isx.normal, new_dir)), 5);
			}
			color += isx.m->REFL * ReflectedColor * schlick;
		}
		if (isx.m->EXPO) {//calculating local phong model
			vec3 R = normalize(reflect(lightVector, isx.normal));
			vec3 V = normalize(ray.dir * (float)isx.t);
			specularTerm = dot(R, V);
			specularTerm = glm::max(glm::min(specularTerm, 1.f), 0.f);
			specularTerm = glm::pow(specularTerm, isx.m->EXPO);
		}
		else {
			specularTerm = 0;
		}
		if (isx.m->TRAN) {//transmittance
			vec3 new_dir;
			if (inside) {
				new_dir = normalize(refract(ray.dir, isx.normal, isx.m->IOR));
				vec3 new_orig = ray.orig + ray.dir * (float)isx.t + (new_dir * epsilon);
				traceRay(Ray(new_orig, new_dir), depth + 1, RefractedColor, false, i, j);
			}
			else {
				new_dir = normalize(refract(ray.dir, isx.normal, 1.f / isx.m->IOR));
				vec3 new_orig = ray.orig + ray.dir * (float)isx.t + (new_dir * epsilon);
				traceRay(Ray(new_orig, new_dir), depth + 1, RefractedColor, true, i, j);
			}

			color += RefractedColor * (schlick != 1 ? 1 - schlick : schlick); //fresnel calculation

		}
		else {
			diffuseTerm = dot(normalize(isx.normal), normalize(lightVector));
			diffuseTerm = glm::max(glm::min(diffuseTerm, 1.f), 0.f);
		}

		float ka = 0.1, kd = 0.8, ks = 0.8;
		//printf("%d %d %d\n", color[0], color[1], color[2]);
		float shadow = 1;
		//float shadow = shadowFeelers(ray.orig + ray.dir * (float)isx.t, isx);
		color += ka * isx.m->DIFF + s.LCOL * (kd * isx.m->DIFF * diffuseTerm + ks * specularTerm);
		color *= shadow;
		if (weightmap[i][j] == -1.f ) {
			weightmap[i][j] = shadow;
		}
		color.x = glm::max(glm::min(color.x, 1.f), 0.f);
		color.y = glm::max(glm::min(color.y, 1.f), 0.f);
		color.z = glm::max(glm::min(color.z, 1.f), 0.f);
	}
}
void generateImage(SceneGraph s) {
	BMP Indirect, Direct, Combination;
	Indirect.SetSize(c.RESO.x, c.RESO.y);
	Direct.SetSize(c.RESO.x, c.RESO.y);
	Combination.SetSize(c.RESO.x, c.RESO.y);
	clock_t t_start = clock();
	//generatePhotonMap();
	int width = c.RESO.x;
	int height = c.RESO.y;
	weightmap = vector<vector<float> >(width, vector<float>(height,-1));
	float sampleSize = 0.25f;
	getLight();
	bool antialiasing = false; //set anti-aliasing here
	cout << "Generating Direct Illumination" << endl;
#pragma omp parallel for
	for (int i = 0; i < width; i++) {
		if (i % 10 == 0) cout << "-";
		for (int j = 0; j < height; j++) {
			//if (i != 30 || j != 59) continue;
			glm::vec3 final_color(0, 0, 0);
			float PNDCx, PNDCy;
			//anti-aliasing with supersampling
			if (antialiasing) {
				for (int k = 0; k < 9; k++) {
					float offset_i = i - sampleSize + (k / 3) * sampleSize;
					float offset_j = j - sampleSize + (k % 3) * sampleSize;
					PNDCx = offset_i / (float)(c.RESO.x - 1);
					PNDCy = offset_j / (float)(c.RESO.y - 1);
					glm::vec3 PW = c.M + (2 * PNDCx - 1)*c.H + (-1 * (2 * PNDCy - 1)*c.V);
					glm::vec3 direction = normalize(PW - c.E);
					Ray ray(c.EYEP, direction);
					vec3 color(0, 0, 0);
					traceRay(ray, 0, color, false, i, j);
					final_color += color / 9.f;
				}
			}
			else { //regular sampling
				PNDCx = i / (float)(c.RESO.x - 1);
				PNDCy = j / (float)(c.RESO.y - 1);
				glm::vec3 PW = c.M + (2 * PNDCx - 1)*c.H + (-1 * (2 * PNDCy - 1)*c.V);
				glm::vec3 direction = normalize(PW - c.E);
				direction.x += (abs((float)glm::cos((float)i*PI / 20.f)) * 0.1f * vec3(1, 1, 1)).x;

				Ray ray(c.EYEP, direction);
				traceRay(ray, 0, final_color, false, i, j);
			}
			vec3 out;
			final_color *= 255;
			out = final_color;
			Direct(i, j)->Red = out.x;
			Direct(i, j)->Green = out.y;
			Direct(i, j)->Blue = out.z;
		}
	}
	Direct.WriteToFile("Direct.bmp");
	antialiasing = false;
	int MCsamples = 10000;
	glm::vec3 transmittance(1, 1, 1);

	vector<vector<vec3> > cache(width, vector<vec3>(height));

	cout << "\nGenerating Indirect Illumination" << endl;
	for (int h = 1; h <= MCsamples; h++) {
		cout << "Sample #" << h << endl;
#pragma omp parallel for
		for (int i = 0; i < width; i++) {
			//if (!(i % 10))printf("%d\n", i);
			for (int j = 0; j < height; j++) {
				glm::vec3 final_color(0, 0, 0);
				float PNDCx, PNDCy;
				//anti-aliasing with supersampling
				if (antialiasing) {
					for (int k = 0; k < 9; k++) {
						float offset_i = i - sampleSize + (k / 3) * sampleSize;
						float offset_j = j - sampleSize + (k % 3) * sampleSize;
						PNDCx = offset_i / (float)(c.RESO.x - 1);
						PNDCy = offset_j / (float)(c.RESO.y - 1);
						glm::vec3 PW = c.M + (2 * PNDCx - 1)*c.H + (-1 * (2 * PNDCy - 1)*c.V);
						glm::vec3 direction = normalize(PW - c.E);	
						Ray ray(c.EYEP, direction);
						final_color += MonteCarlo(ray, 0, transmittance, false, 0,1) / 9.f;
					}
				}
				else { //regular sampling
					PNDCx = i / (float)(c.RESO.x - 1);
					PNDCy = j / (float)(c.RESO.y - 1);
					glm::vec3 PW = c.M + (2 * PNDCx - 1)*c.H + (-1 * (2 * PNDCy - 1)*c.V);
					glm::vec3 direction = normalize(PW - c.E);
					Ray ray(c.EYEP, direction);
					final_color += MonteCarlo(ray, 0, transmittance, false, 0,1);
				}
				final_color *= 255;
				vec3 out((Indirect(i, j)->Red * (h - 1.f) + final_color.x) / (float)h,
							(Indirect(i, j)->Green * (h - 1.f) + final_color.y) / (float)h,
							(Indirect(i, j)->Blue * (h - 1.f) + final_color.z) / (float)h);
				
				cache[i][j].r = (cache[i][j].r * (h - 1.f) + final_color.x) / (float)h;
				cache[i][j].g = (cache[i][j].g * (h - 1.f) + final_color.y) / (float)h;
				cache[i][j].b = (cache[i][j].b * (h - 1.f) + final_color.z) / (float)h;


				out.x = max(0.f, min(out.x, 255.f));
				out.y = max(0.f, min(out.y, 255.f));
				out.z = max(0.f, min(out.z, 255.f));

				Indirect(i, j)->Red = max(0.f, min(cache[i][j].r, 255.f));
				Indirect(i, j)->Green = max(0.f, min(cache[i][j].g, 255.f));
				Indirect(i, j)->Blue = max(0.f, min(cache[i][j].b, 255.f));
				
				float weight = abs(weightmap[i][j]);
				weight = (9.f * weight + 0.5f) / 10.f;

				Combination(i, j)->Red = Direct(i, j)->Red * weight + Indirect(i, j)->Red * (1 - weight);
				Combination(i, j)->Green = Direct(i, j)->Green * weight + Indirect(i, j)->Green * (1 - weight);
				Combination(i, j)->Blue = Direct(i, j)->Blue * weight + Indirect(i, j)->Blue * (1 - weight);
			}
		}
		Indirect.WriteToFile("Indirect.bmp");
		Combination.WriteToFile("Combination.bmp");
	}
	cout << endl << endl;
	clock_t t_end = clock();
	cout << ((float)t_end - (float)t_start) / 1000 << endl;
	//cout<<"] Done!"<<endl;
}

/*==========================================================
**********************MONTE CARLO************************
============================================================*/

vec3 MonteCarlo(Ray& ray, int depth, glm::vec3& transmittance, bool inside, int hitcount, int canGenerate) {
	Intersection isx;
	float epsilon = 0.001f, absorb;
	vec3 color(0, 0, 0);
	glm::vec3 ambient(0.1, 0.1, 0.1);
	if (depth > 20) { //we're done and found nothing
		return vec3(0, 0, 0);
	}

	isx = RayIntersect(ray, s.root);

	if (isx.t == -1) {
		return vec3(0, 0, 0);
	}
	else if (isx.m->EMIT) { //it hit a light
		if (hitcount < -1) return vec3(0, 0, 0);//YOU CAN CHANGE THE TYPE OF MC PATH TRACING HERE 2 FOR PURE INDIRECT -1 FOR VANILLA MC TRACING
		return isx.m->EMIT * transmittance * isx.m->DIFF;
	}
	else { //there is an intersection
		float schlick = 0;
		vec3 new_dir, new_orig;
		vec3 lightVector = (s.LPOS - (ray.orig + ray.dir * (float)isx.t));
		
		if (isx.m->TRAN && isx.m->MIRR) {
			float r0 = pow(((isx.m->IOR - 1) / (isx.m->IOR + 1)), 2);
			new_dir = normalize(reflect(ray.dir, isx.normal));
			schlick = r0 + (1 - r0) * pow((1 - dot(isx.normal, new_dir)), 5);
		}
		else if (isx.m->MIRR) {
			schlick = 1;
		}
		float random = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		
		//if(schlick != 0 || schlick != 1) printf("%f\n", schlick);

		if (random < isx.m->TRAN && random > schlick) {//transmittance
			vec3 new_dir, new_orig;
			if (inside) {
				new_dir = normalize(refract(ray.dir, isx.normal, isx.m->IOR));
				new_orig = ray.orig + ray.dir * (float)isx.t + (new_dir * epsilon);
			}
			else {
				new_dir = normalize(refract(ray.dir, isx.normal, 1.f / isx.m->IOR));
				new_orig = ray.orig + ray.dir * (float)isx.t + (new_dir * epsilon);
			}
			return MonteCarlo(Ray(new_orig, new_dir), depth + 1, transmittance, !inside, hitcount + (hitcount ? 1 : 0),1);

		}
		else if (random < isx.m->MIRR && random <= schlick) {// reflect
			new_dir = normalize(reflect(ray.dir, isx.normal));
			new_orig = ray.orig + ray.dir * (float)isx.t + (new_dir * epsilon);
			vec3 new_transmittance = transmittance * isx.m->REFL;
			return MonteCarlo(Ray(new_orig, new_dir), depth + 1, new_transmittance, inside, hitcount + (hitcount ? 1 : 0), 1);
		}
		else {
			absorb = 1 - max(isx.m->DIFF.z, max(isx.m->DIFF.x, isx.m->DIFF.y));
			if (isAbsorbed(absorb)) {
				return vec3(0, 0, 0);
			}
			else {
				isx.normal = normalize(isx.normal);
				vec3 new_dir = getCosineWeightedDirection(isx.normal);
				new_dir = normalize(new_dir);
				vec3 new_orig = ray.orig + ray.dir * (float)isx.t + (new_dir * epsilon);
				return MonteCarlo(Ray(new_orig, new_dir), depth + 1, transmittance * isx.m->DIFF / (1 - absorb), inside, hitcount + 1,1);
			}
		}
	}
}

bool isAbsorbed(const float absorb) {
	return (static_cast <float> (rand()) / static_cast <float> (RAND_MAX) < absorb);
}

vec3 getCosineWeightedDirection(const glm::vec3& normal) {

	// Pick 2 random numbers in the range (0, 1)
	float xi1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
	float xi2 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

	float up = glm::sqrt(xi1);          // cos(theta)
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


/*==========================================================
********************USER INTERFACE**************************
============================================================*/

void removeChildren(node *n) {
	for (int i = 0; i < n->children.size(); i++) {
		removeChildren(n->children[i]);
	}
	n->children.clear();
	delete n;
}

void removeNode(node *n) {
	node *p = n->parent;
	if (p != NULL) {
		for (int i = 0; i < p->children.size(); i++) {
			if (p->children[i] == n) {
				p->children.erase(p->children.begin() + i);
				break;
			}
		}
		removeChildren(n);
	}
	else {
		for (int i = 0; i < n->children.size(); i++) {
			removeChildren(n->children[i]);
		}
		n->children.clear();
		delete n;
		s.root = new node();
	}

}

void keypress(unsigned char key, int x, int y)
{
	cout << key << endl;
	switch (key) {
	case 'q':
		cleanup();
		exit(0);
		break;
	case 'n':
		nextNode();
		break;

	case 'u':
	{
		removeNode(currentNode);
		currentNode = s.root;
		break;
	}
	case 'a':
		currentNode->TRANSLATION.x += -0.5f;
		break;
	case 'd':
		currentNode->TRANSLATION.x += 0.5f;
		break;
	case 'w':
		currentNode->TRANSLATION.y += 0.5f;
		break;
	case 's':
		currentNode->TRANSLATION.y += -0.5f;
		break;
	case 'e':
		currentNode->TRANSLATION.z += 0.5f;
		break;
	case 'r':
		currentNode->TRANSLATION.z += -0.5f;
		break;
	case 'x':
		currentNode->SCALE.x += 0.5;
		break;
	case 'X':
		currentNode->SCALE.x += -0.5;
		break;
	case 'y':
		currentNode->SCALE.y += 0.5;
		break;
	case 'Y':
		currentNode->SCALE.y += -0.5;
		break;
	case 'z':
		currentNode->SCALE.z += 0.5;
		break;
	case 'Z':
		currentNode->SCALE.z += -0.5;
		break;
	case 'j':
		currentNode->ROTATION.x += (10.0f);
		break;
	case 'J':
		currentNode->ROTATION.x += (-10.0f);
		break;
	case 'k':
		currentNode->ROTATION.y += (10.0f);
		break;
	case 'K':
		currentNode->ROTATION.y += (-10.0f);
		break;
	case 'l':
		currentNode->ROTATION.z += (10.0f);
		break;
	case 'L':
		currentNode->ROTATION.z += (-10.0f);
		break;
	case 'f':
		s.LPOS.x += 0.5f;
		break;
	case 'F':
		s.LPOS.x -= 0.5f;
		break;
	case 'g':
		s.LPOS.y += 0.5f;
		break;
	case 'G':
		s.LPOS.y -= 0.5f;
		break;
	case 'h':
		s.LPOS.z += 0.5f;
		break;
	case 'H':
		s.LPOS.z -= 0.5f;
		break;
	case 'p':
		generateImage(s);
		break;
	default:
		break;
	}
	glutPostRedisplay();
}

void mousepress(int button, int state, int x, int y)
{
	// Put any mouse events here
}

void nextNode() {
	if (nodeStack.empty()) {
		nodeStack.push(s.root);
	}
	for (int i = currentNode->children.size() - 1; i >= 0; i--) {
		nodeStack.push(currentNode->children[i]);
	}
	currentNode = nodeStack.top();
	nodeStack.pop();
}

/*==========================================================
********************OPENGL RENDERING************************
============================================================*/
void getLight(node* n, glm::mat4 matrix) {
	matrix = glm::translate(matrix, n->TRANSLATION);
	matrix = glm::translate(matrix, n->CENTER);
	matrix = glm::rotate(matrix, radians(n->ROTATION.x), vec3(1, 0, 0));
	matrix = glm::rotate(matrix, radians(n->ROTATION.y), vec3(0, 1, 0));
	matrix = glm::rotate(matrix, radians(n->ROTATION.z), vec3(0, 0, 1));
	matrix = glm::scale(matrix, n->SCALE);
	matrix = glm::translate(matrix, vec3(0, 0, 0) - n->CENTER);

	if (n->geometry) {
		if (n->MAT->EMIT) {
			s.LCOL = n->MAT->DIFF;
			s.LPOS = vec3(matrix * vec4(0, 0, 0, 1));
		}
	}
	for (int i = 0; i < n->children.size(); i++) {
		getLight(n->children[i], matrix);
	}
}

void traverse(node* n, glm::mat4 matrix) {
	matrix = glm::translate(matrix, n->TRANSLATION);
	matrix = glm::translate(matrix, n->CENTER);
	matrix = glm::rotate(matrix, radians(n->ROTATION.x), vec3(1, 0, 0));
	matrix = glm::rotate(matrix, radians(n->ROTATION.y), vec3(0, 1, 0));
	matrix = glm::rotate(matrix, radians(n->ROTATION.z), vec3(0, 0, 1));
	matrix = glm::scale(matrix, n->SCALE);
	matrix = glm::translate(matrix, vec3(0, 0, 0) - n->CENTER);

	if (n->geometry) {
		//uploadPrimitive(n, n->geometry);
		drawPrimitive(matrix, n, n->geometry);
		n->matrix = new mat4(matrix);
	}
	for (int i = 0; i < n->children.size(); i++) {
		traverse(n->children[i], matrix);
	}
}

void renderScene(SceneGraph &s) {
	glm::mat4 matrix = glm::mat4();
	traverse(s.root, matrix);
}

bool isSelected(node *n) {
	if (n == currentNode) {
		return true;
	}
	else if (n == s.root) {
		return false;
	}
	else {
		return isSelected(n->parent);
	}
}

void uploadScene(node* n) {
	if (n->geometry) {
		uploadPrimitive(n, n->geometry);
	}
	for (int i = 0; i < n->children.size(); i++) {
		uploadScene(n->children[i]);
	}
}

void uploadPrimitive(node *n, Geometry* g)
{
	// Create the VBOs and vboIdx we'll be using to render Indirects in OpenGL
	glGenBuffers(1, &vboIdx[n->NUM]);
	glGenBuffers(1, &vboConsolidate[n->NUM]);

	const int VERTICES = g->getVertexCount();
	const int TRIANGLES = g->getIndexCount() / 3;

	static const GLsizei SIZE_POS = sizeof(glm::vec3);
	static const GLsizei SIZE_NOR = sizeof(glm::vec3);
	static const GLsizei SIZE_COL = sizeof(glm::vec3);
	static const GLsizei SIZE_TRI = 3 * sizeof(GLuint);
	vector<vec3> Colors(g->getColors());

	//if (isSelected(n)) { //THIS WILL GIVE YOU HIGHLIGHTING COLOR ON THE ACTIVE NODE
	//  Colors = vector<vec3>(VERTICES, vec3(0.95f, 0.83f, 0.3f));
	//}

	vector<float> consolidate;
	for (int i = 0; i < g->getVertexCount(); i++) {
		consolidate.push_back(g->getVertices()[i].x);
		consolidate.push_back(g->getVertices()[i].y);
		consolidate.push_back(g->getVertices()[i].z);

		consolidate.push_back(g->getNormals()[i].x);
		consolidate.push_back(g->getNormals()[i].y);
		consolidate.push_back(g->getNormals()[i].z);

		consolidate.push_back(Colors[i].x);
		consolidate.push_back(Colors[i].y);
		consolidate.push_back(Colors[i].z);

	}
	//std::cout<<"CONSOLIDATE SIZE: "<<consolidate.size()<<endl;

	glBindBuffer(GL_ARRAY_BUFFER, vboConsolidate[n->NUM]);
	glBufferData(GL_ARRAY_BUFFER, consolidate.size()*SIZE_POS / 3, &consolidate[0], GL_STATIC_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vboIdx[n->NUM]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, TRIANGLES * SIZE_TRI, &g->getIndices()[0], GL_STATIC_DRAW);

}

void drawPrimitive(glm::mat4 model, node *n, Geometry *g)
{
	glUseProgram(shaderProgram);
	//model = glm::rotate(model, 0.f, glm::vec3(1, 0, 1));

	const int FACES = 1;
	const int TRIANGLES = g->getIndexCount() / 3;

	glEnableVertexAttribArray(locationPos);
	glEnableVertexAttribArray(locationCol);
	glEnableVertexAttribArray(locationNor);

	glUniformMatrix4fv(unifModel, 1, GL_FALSE, &model[0][0]);
	const glm::mat4 modelInvTranspose = glm::inverse(glm::transpose(model));
	glUniformMatrix4fv(unifModelInvTr, 1, GL_FALSE, &modelInvTranspose[0][0]);

	glUniform4fv(unifLightPos, 1, &glm::vec4(s.LPOS.x, s.LPOS.y, s.LPOS.z, 1)[0]);
	glUniform3fv(unifLightColor, 1, &glm::vec3(s.LCOL.x, s.LCOL.y, s.LCOL.z)[0]);
	glUniform3fv(unifCameraPos, 1, &glm::vec3(c.EYEP.x, c.EYEP.y, c.EYEP.z)[0]);

	glBindBuffer(GL_ARRAY_BUFFER, vboConsolidate[n->NUM]);
	glVertexAttribPointer(locationPos, 3, GL_FLOAT, false, 9 * sizeof(float), (void*)(0 * (sizeof(float))));

	glBindBuffer(GL_ARRAY_BUFFER, vboConsolidate[n->NUM]);
	glVertexAttribPointer(locationNor, 3, GL_FLOAT, false, 9 * sizeof(float), (void*)(3 * (sizeof(float))));

	glBindBuffer(GL_ARRAY_BUFFER, vboConsolidate[n->NUM]);
	glVertexAttribPointer(locationCol, 3, GL_FLOAT, false, 9 * sizeof(float), (void*)(6 * (sizeof(float))));

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vboIdx[n->NUM]);

	glDrawElements(GL_TRIANGLES, TRIANGLES * 3, GL_UNSIGNED_INT, 0);

	glDisableVertexAttribArray(locationPos);
	glDisableVertexAttribArray(locationCol);
	glDisableVertexAttribArray(locationNor);

	printGLErrorLog();
}

void display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	clock_t newTime = clock();
	rotation += 2.5f * (static_cast<float>(newTime - old_time) / static_cast<float>(CLOCKS_PER_SEC));
	old_time = newTime;

	glUseProgram(shaderProgram);

	renderScene(s);

	glutSwapBuffers();

	printGLErrorLog();
}

void resize(int width, int height)
{
	// Set viewport
	glViewport(0, 0, c.RESO.x, c.RESO.y);

	// Get camera information
	// Add code here if you want to play with camera settings/ make camera interactive.
	glm::mat4 projection = glm::perspective(glm::radians(c.FOVY), c.RESO.x / (float)c.RESO.y, 0.1f, 100.0f);
	glm::mat4 camera = glm::lookAt(c.EYEP, c.EYEP + c.VDIR, c.UVEC);
	projection = projection * camera;

	// Upload the projection matrix, which changes only when the screen or
	// camera changes
	glUseProgram(shaderProgram);
	glUniformMatrix4fv(unifViewProj, 1, GL_FALSE, &projection[0][0]);

	glutPostRedisplay();
}

/*==========================================================
********************OPENGL DEBUGGING************************
============================================================*/

std::string textFileRead(const char *filename)
{
	// http://insanecoding.blogspot.com/2011/11/how-to-read-in-file-in-c.html
	std::ifstream in(filename, std::ios::in);
	if (!in) {
		std::cerr << "Error reading file" << std::endl;
		throw (errno);
	}
	return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

void printGLErrorLog()
{
	GLenum error = glGetError();
	if (error != GL_NO_ERROR) {
		std::cerr << "OpenGL error " << error << ": ";
		const char *e =
			error == GL_INVALID_OPERATION ? "GL_INVALID_OPERATION" :
			error == GL_INVALID_ENUM ? "GL_INVALID_ENUM" :
			error == GL_INVALID_VALUE ? "GL_INVALID_VALUE" :
			error == GL_INVALID_INDEX ? "GL_INVALID_INDEX" :
			"unknown";
		std::cerr << e << std::endl;

		// Throwing here allows us to use the debugger stack trace to track
		// down the error.
#ifndef __APPLE__
		// But don't do this on OS X. It might cause a premature crash.
		// http://lists.apple.com/archives/mac-opengl/2012/Jul/msg00038.html
		throw;
#endif
	}
}

void printLinkInfoLog(int prog)
{
	GLint linked;
	glGetProgramiv(prog, GL_LINK_STATUS, &linked);
	if (linked == GL_TRUE) {
		return;
	}
	std::cerr << "GLSL LINK ERROR" << std::endl;

	int infoLogLen = 0;
	int charsWritten = 0;
	GLchar *infoLog;

	glGetProgramiv(prog, GL_INFO_LOG_LENGTH, &infoLogLen);

	if (infoLogLen > 0) {
		infoLog = new GLchar[infoLogLen];
		// error check for fail to allocate memory omitted
		glGetProgramInfoLog(prog, infoLogLen, &charsWritten, infoLog);
		std::cerr << "InfoLog:" << std::endl << infoLog << std::endl;
		delete[] infoLog;
	}
	// Throwing here allows us to use the debugger to track down the error.
	throw;
}

void printShaderInfoLog(int shader)
{
	GLint compiled;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &compiled);
	if (compiled == GL_TRUE) {
		return;
	}
	std::cerr << "GLSL COMPILE ERROR" << std::endl;

	int infoLogLen = 0;
	int charsWritten = 0;
	GLchar *infoLog;

	glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &infoLogLen);

	if (infoLogLen > 0) {
		infoLog = new GLchar[infoLogLen];
		// error check for fail to allocate memory omitted
		glGetShaderInfoLog(shader, infoLogLen, &charsWritten, infoLog);
		std::cerr << "InfoLog:" << std::endl << infoLog << std::endl;
		delete[] infoLog;
	}
	// Throwing here allows us to use the debugger to track down the error.
	throw;
}
