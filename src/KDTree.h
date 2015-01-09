#pragma once
#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <vector>
#include "Photon.h"

using namespace std;
using namespace glm;
struct KDTree;

class KDTree {
private:
	struct node {
		node() { isLeaf = false; }
		node(Photon p) {
			data = p;
			isLeaf = true;
		}
		Photon data;
		bool isLeaf;
		node *left, *right;
	};
	node *root;
public:
	KDTree() {
	}
	KDTree(const vector<Photon>& p) {
		P = p;
		vector<int> A;
		for (int i = 0; i < P.size(); i++) {
			A.push_back(i);
		}
		root = buildTree(A, 0);
	}
	
	vector<Photon> P;
	float(*KDTree::axis)(Photon);

	node* buildTree(vector<int>& A, int depth) {
		if (A.size() == 1) {
			return &node(P[0]);
		}
		else if (!A.size()) {
			return 0;
		}
		int axisNum = depth % 3;
		if (axisNum == 0) {
			axis = getX;
		}
		else if (axisNum == 1) {
			axis = getY;
		}
		else if (axisNum == 2) {
			axis = getZ;
		}

		float pivot_index = getMedianIndex(A, (A.size() + 1) / 2);
		vector<int> less, greater;
		for (int i = 0; i < A.size(); i++) {
			//if (i == pivot_index) continue;
			float pivot_val = (*axis)(P[pivot_index]);
			if ((*axis)(P[A[i]]) <= pivot_val) {
				less.push_back(A[i]);
			}
			else {
				greater.push_back(A[i]);
			}
		}
		node n;
		n.left = buildTree(less, depth + 1);
		n.right = buildTree(greater, depth + 1);
		return &n;
	}
	int getMedianIndex(vector<int>& A, int k) {
		if (A.size() == 0) {
			return -1;
		}
		else {
			float pivot = (*axis)(P[A[0]]);
			vector<int> less, greater;
			for (int i = 1; i < A.size(); i++) {
				float pval = (*axis)(P[A[i]]);
				if (pval <= pivot) {
					less.push_back(A[i]);
				}
				else {
					greater.push_back(A[i]);
				}
			}
			if (less.size() == k - 1) {
				return A[0];
			}
			else if (less.size() > k - 1) {
				return getMedianIndex(less, k);
			}
			else {
				return getMedianIndex(greater, k - less.size() - 1);
			}
		}
	}
	//vector<vec3> query(node* n, Photon& p, float radius) {
	//	vector<Photon> neighbors;
	//	if (n->isLeaf) {
	//		if (length(n->data.loc - p.loc) < radius) {
	//			neighbors.push_back(p);
	//		}
	//	}
	//	if ()
	//}
	static float getX(Photon p) {
		return p.loc.x;
	}
	static float getY(Photon p) {
		return p.loc.y;
	}
	static float getZ(Photon p) {
		return p.loc.z;
	}


};