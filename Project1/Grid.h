#pragma once
#include "Node.h"
#include "Element.h"
#include "UniversalElement.h"
class Grid {

public:

	vector <Node*> tabNodes;
	vector <Element*> tabElements;
	vector <UniversalElement*> tabUniversalElements;
	double K;
	double c;
	double t;
	double tot;
	double alfa;
	double ro;
	double time;
	double step;

	vector<vector<double>> globalH;
	vector<vector<double>> globalC;
	vector<double> globalP;

	vector<vector<double>> CdT;
	vector<vector<double>> HCdT;
	vector<double> t0;
	vector<double> t1;

	Grid();
	Grid(vector<Node*>, vector<Element*>, double, double, double, double, double, double, double, double);
	~Grid();

	void generateUniversalElements();
	void mapLocalToGlobal();
	void generateEquationElements();
	void updateTemperatures(vector<double>);
	void solve();

	void show(int);

};