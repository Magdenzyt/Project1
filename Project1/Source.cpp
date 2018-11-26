#include<iostream>
#include"Grid.h"
#include"Node.h"
#include"json.hpp"
#include"Element.h"
#include"UniversalElement.h"
#include<fstream>

using namespace std;

Grid* generateGrid(double, double, double, double, double, double);

int main(){

	double H = 0.1;
	double L = 0.1;
	double nH = 4.0;
	double nL = 4.0;
	double K = 25.0;
	double t = 100.0;

	
	Grid *myGrid = generateGrid(H,L,nH,nL,t,K);
	myGrid->solve();
	

	system("pause");
	return 0;
}
Grid* generateGrid(double H, double L, double nH, double nL, double t, double K){

	vector <Node*> myNodes;
	double dH = H / (nH - 1);
	double dL = L / (nL - 1);

	for (int i = 0; i < nL; i ++) {
		for (int j = 0; j < nH; j ++) {
			Node *myNode = new Node(i*nL + j, i*dH, j*dL, t);
			if ((i == 0 || i == nL - 1) || (j == 0 || j == nL - 1))
				myNode->bc = true;
			myNodes.push_back(myNode);
		}
	}

	vector <Element*> myElements;
	for (int i = 0; i < nL - 1; i++) {
		for (int j = 0; j < nH - 1; j++) {
			vector <Node*> tmpN;
			tmpN.push_back(myNodes[i * nH + j]);
			tmpN.push_back(myNodes[(i + 1) * nH + j]);
			tmpN.push_back(myNodes[(i + 1) * nH + j + 1]);
			tmpN.push_back(myNodes[i * nH + j + 1]);

			Element *tmpE = new Element(tmpN, K);

			myElements.push_back(tmpE);
		}
	}
	return (new Grid(myNodes, myElements, K));
}