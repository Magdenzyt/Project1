#pragma once
#include<vector>
#include"Node.h"
class Element {

public:

	vector <Node*> elemNodes;
	double K;

	Element();
	Element(vector <Node*>, double);
	~Element();

	void show();

};
