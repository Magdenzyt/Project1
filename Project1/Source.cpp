#include<iostream>
#include<string>
#include"Grid.h"
#include"Node.h"
#include"json.hpp"
#include"Element.h"
#include"UniversalElement.h"
#include<fstream>

using namespace std;

Grid* generateGrid(double, double, double, double, double, double, double, double, double, double, double, double);
void readFromTxt(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);

int main(){

	double H, L, nH, nL, K, t, c, tot, alfa, ro, time, step;

	readFromTxt(&H, &L, &nH, &nL, &K, &t, &c, &tot, &alfa, &ro, &time, &step);

	cout << H << endl;
	cout << L << endl;
	cout << nH << endl;
	cout << nL << endl;
	cout << K << endl;
	cout << t << endl;
	cout << c << endl;
	cout << tot << endl;
	cout << alfa << endl;
	cout << ro << endl;
	cout << time << endl;
	cout << step << endl;
	
	Grid *myGrid = generateGrid(H,L,nH,nL,t,K,c,tot,alfa,ro,time,step);
	myGrid->solve();
	

	system("pause");
	return 0;
}

Grid* generateGrid(double H, double L, double nH, double nL, double t, double K, double c, double tot, double alfa, double ro, double time, double step){

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
	return (new Grid(myNodes, myElements, K, c, t, tot, alfa, ro, time, step));
}

void readFromTxt(double* H, double* L, double* nH, double* nL, double* K, double* t, double* c, double* tot, double* alfa, double* ro, double* time, double* step) {
	ifstream inFile("Data.txt");
	string strOneLine;

	while (inFile)
	{
		getline(inFile, strOneLine);
		vector<string> list;
		size_t pos = 0;
		string token;
		while ((pos = strOneLine.find("=")) != string::npos) {
			token = strOneLine.substr(0, pos);
			list.push_back(token);
			strOneLine.erase(0, pos + 1);
		}
		list.push_back(strOneLine);
		
		if (list[0] == "H")
			*H = stod(list[1]);
		if (list[0] == "L")
			*L = stod(list[1]);
		if (list[0] == "nH")
			*nH = stod(list[1]);
		if (list[0] == "nL")
			*nL = stod(list[1]);
		if (list[0] == "K")
			*K = stod(list[1]);
		if (list[0] == "t")
			*t = stod(list[1]);
		if (list[0] == "c")
			*c = stod(list[1]);
		if (list[0] == "tot")
			*tot = stod(list[1]);
		if (list[0] == "alfa")
			*alfa = stod(list[1]);
		if (list[0] == "ro")
			*ro = stod(list[1]);
		if (list[0] == "time")
			*time = stod(list[1]);
		if (list[0] == "step")
			*step = stod(list[1]);
	}
	inFile.close();
}