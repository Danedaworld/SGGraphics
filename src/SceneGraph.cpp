#include "SceneGraph.h"

SceneGraph::SceneGraph() {

}

SceneGraph::SceneGraph(string f) : filename(f)
{
	numNodes = 0;
	parseFile();
}


SceneGraph::~SceneGraph(void)
{

}

void SceneGraph::parseFile() {

	ifstream file(filename.c_str());

	while(!file.eof()) {
		string tmp;
		file >> tmp;
		if (tmp == "MAT") {
			string NAME;
			vec3 DIFF, REFL;
			float EXPO, IOR, MIRR, TRAN;
			int EMIT;
			file >> NAME;
			for(int i = 0; i < 7; i++) {
				file >> tmp;
				if(tmp == "DIFF") {
					for (int i = 0; i < 3; i++) {
						file >> tmp;
						DIFF[i] = atof(tmp.c_str());
						
					}
				} else if (tmp == "REFL") {
					for (int i = 0; i < 3; i++) {
						file >> tmp;
						REFL[i] = atof(tmp.c_str());
					}
				} else if (tmp == "EXPO") {
					file >> tmp;
					EXPO = atof(tmp.c_str());
				} else if (tmp == "IOR") {
					file >> tmp;
					IOR = atof(tmp.c_str());
				} else if (tmp == "MIRR") {
					file >> tmp;
					MIRR = atof(tmp.c_str());
				} else if (tmp == "TRAN") {
					file >> tmp;
					TRAN = atof(tmp.c_str());
				} else if (tmp == "EMIT") {
					file >> tmp;
					EMIT = atof(tmp.c_str());
					//cout << EMIT << endl;
				}
			}
			matlist.push_back(new material(NAME, DIFF, REFL, EXPO, IOR, MIRR, TRAN, EMIT));
		} else if (tmp == "NODE") {
			string NODE, PARENT, SHAPE, MATNAME, FILE = "";
			vec3 TRANSLATION, ROTATION, SCALE, CENTER;
			// parse and load everything
			file >> NODE;
			for (int i = 0; i < 7; i++) {
				file >> tmp;
				if (tmp == "TRANSLATION") {
					for (int i = 0; i < 3; i++) {
						file >> tmp;
						TRANSLATION[i] = atof(tmp.c_str());
					}
				} else if (tmp == "ROTATION") {
					for (int i = 0; i < 3; i++) {
						file >> tmp;
						ROTATION[i] = atof(tmp.c_str());
					}
				} else if (tmp == "SCALE") {
					for (int i = 0; i < 3; i++) {
						file >> tmp;
						SCALE[i] = atof(tmp.c_str());
					}
				} else if (tmp == "CENTER") {
					for (int i = 0; i < 3; i++) {
						file >> tmp;
						CENTER[i] = atof(tmp.c_str());
					}
				} else if (tmp == "MAT") {
					file >> MATNAME;
				} else if (tmp == "PARENT") {
					file >> PARENT;
				} else if (tmp == "SHAPE") {
					file >> SHAPE;
					if(SHAPE == "mesh") {
						file >> tmp;
						file >> FILE;
					}
				}
			}
			// create node structures
			if(PARENT == "null") {
				
				root = new node(NODE, PARENT, SHAPE, TRANSLATION, ROTATION, SCALE, CENTER, getMat(MATNAME), numNodes, FILE);
				addNode(root);
			} else {
				addNode(new node(NODE, PARENT, SHAPE, TRANSLATION, ROTATION, SCALE, CENTER, getMat(MATNAME), numNodes, FILE));
			}
		} else if(tmp == "LIGHT") {
			for(int i = 0; i < 2; i++) {
				file >> tmp;
				if(tmp == "LPOS") {
					for (int i = 0; i < 3; i++) {
						file >> tmp;
						LPOS[i] = atof(tmp.c_str());
					}
				} else if(tmp == "LCOL") {
					for (int i = 0; i < 3; i++) {
						file >> tmp;
						LCOL[i] = atof(tmp.c_str());
					}
				}
			}
		}
	}
	file.close();
}
// lazy tree insertion
void SceneGraph::addNode(node* n) {
	for(int i = 0; i < nodelist.size(); i++) {
		if (n->PARENT == nodelist[i]->NODE) {
			n->parent = nodelist[i];
			nodelist[i]->children.push_back(n);
		}
	}
	nodelist.push_back(n);
	numNodes++;
}

node* SceneGraph::getRoot() {
	return root;
}

material* SceneGraph::getMat(string MATNAME) {
	for(int i = 0; i < matlist.size(); i++) {
		if(MATNAME == matlist[i]->NAME) {
			return matlist[i];
		}
	}
	return NULL;
}

void SceneGraph::release() {
	for(int i = 0; i < nodelist.size(); i++) {
		free(nodelist[i]);
	}
	for(int i = 0; i < matlist.size() ;i++) {
		free(matlist[i]);
	}
}