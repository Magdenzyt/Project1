#pragma once
#include<vector>
#include<cmath>
#include"Element.h"

class UniversalElement {
public:

	Element* element;

	double K;
	double c;
	double ro;
	double alfa;
	double tot;
	double ksi[4];
	double eta[4];	
	
	double N[4][4];
	
	double xP[4];
	double yP[4];

	double dNdKsi[4][4];
	double dNdEta[4][4];

	double detJ[4];
	double Jacobian[4][4];

	double H[4][4];
	double C[4][4];
	
	double length[4];
	double HBC[4][4] = { {0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0} };

	double P[4] = { 0.0,0.0,0.0,0.0 };



	UniversalElement(Element*, double, double, double, double, double);
	~UniversalElement();

	void generateN();
	void generateUniversalXY();
	void generatedNdKsidNdEta();
	void generateJacobian();
	void generateH();
	void generateC();
	void generateHBC();
	void generateP();

};