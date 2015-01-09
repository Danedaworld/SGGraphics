#pragma once

#include "Geometry.h"
#include "Cube.h"
#include "Cylinder.h"
#include "Sphere.h"
#include "Mesh.h"
#include "Config.h"
#include <fstream>
#include <vector>

using namespace glm;

struct node;

struct node
{
	node *parent;
	vector<node*> children;
	Geometry *geometry;
	material *MAT;
	string NODE, PARENT, SHAPE, FILE;
	vec3 TRANSLATION, ROTATION, SCALE, CENTER;
	mat4 *matrix;
	int NUM;
	node() {
		parent = NULL;
		geometry = NULL;
	}
	node(string n, string p, string s, vec3 t, vec3 r, vec3 sc, vec3 c, material *mat, int N, string f = "") {
		parent = NULL;
		NODE = n;
		PARENT = p;
		SHAPE = s;
		TRANSLATION = t;
		ROTATION = r;
		SCALE = sc;
		CENTER = c;
		FILE = f;
		NUM = N;
		MAT = mat;
		//get the correct material


		if(SHAPE == "cube") {
			geometry = new Cube(mat);
		} else if (SHAPE == "sphere") {
			geometry = new Sphere(mat);
		} else if (SHAPE == "cylinder") {
			geometry = new Cylinder(mat);
		} else if (SHAPE == "mesh") {
			geometry = new Mesh(mat, FILE);
		} else {
			geometry = NULL;
		}
	}
	~node() {
		delete geometry;
	}


};

class SceneGraph
{
public:
	SceneGraph();
	SceneGraph(string f);
	~SceneGraph(void);
	
	void parseFile();
	void addNode(node* n);
	node* getRoot();
	material* getMat(string MATNAME);
	void release();

	node* root;
	vec3 LPOS, LCOL;
	node* lightnode;
	void setLight() {
		for (int i = 0; i < nodelist.size(); i++) {
			if (nodelist[i]->MAT) {
				if (nodelist[i]->MAT->EMIT == 1) {
					lightnode = nodelist[i];
					cout << lightnode->SHAPE << endl;
					return;
				}
			}	
		}
	}
	int numNodes;
	string filename;
	vector<node*> nodelist;
	vector<material*> matlist;
};

