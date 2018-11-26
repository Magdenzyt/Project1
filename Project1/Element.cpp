#include <stdio.h>
#include <ctime>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "Element.h"
using namespace std;

Element::Element() {}

Element::Element(vector <Node*> elemNodes, double K) {
	this->elemNodes = elemNodes;
	this->K = K;
}

Element::~Element() {}

void Element::show() {
	cout << "Element:" << endl;
	for (int i = 0; i < elemNodes.size(); i++) {
		elemNodes[i]->show();
	}
}