#include "Grid.h"
#include<iostream>
#include<algorithm>
using namespace std;

vector < vector<double> > multiplyMatrices(vector < vector<double> > m1, vector < vector<double> > m2) {
	const int n = m1.size();     // a rows
	const int m = m1[0].size();  // a cols
	const int p = m2[0].size();  // b cols

	vector <vector<double>> result(n, vector<double>(p, 0));
	for (auto j = 0; j < p; ++j)
	{
		for (auto k = 0; k < m; ++k)
		{
			for (auto i = 0; i < n; ++i)
			{
				result[i][j] += m1[i][k] * m2[k][j];
			}
		}
	}
	return result;
}

vector <double> multiplyVectorMatrix(vector < vector<double> > m1, vector<double> v1) {
	int i, j;
	int m = m1.size();
	int n = v1.size();

	vector<double> result(m);

	for (i = 0; i < m; i++) {
		result[i] = 0.0;
		for (j = 0; j < n; j++)
			result[i] += m1[i][j] * v1[j];
	}
	return result;
}

void calculateInverse(vector< vector<double> > &A) {

	int n = A[0].size();
	
	vector<double> fragTmp;
	vector<vector<double>> tmp;
	for (int i = 0; i < n*2; i++) {
		fragTmp.push_back(0);
	}
	for (int i = 0; i < n; i++) {	
		tmp.push_back(fragTmp);
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			tmp[i][j] = A[i][j];
		}
	}

	for (int i = 0; i<n; i++)
	{
		for (int j = n; j<2 * n; j++)
		{
			if (i == j - n)
				tmp[i][j] = 1;
			else
				tmp[i][j] = 0;
		}
	}
	for (int i = 0; i<n; i++)
	{
		double t = tmp[i][i];
		for (int j = i; j<2 * n; j++)
			tmp[i][j] = tmp[i][j] / t;
		for (int j = 0; j<n; j++)
		{
			if (i != j)
			{
				t = tmp[j][i];
				for (int k = 0; k<2 * n; k++)
					tmp[j][k] = tmp[j][k] - t*tmp[i][k];
			}
		}
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			A[i][j] = tmp[i][j+n];
		}
	}
	/*for (int i = 0; i < n; i++) {
		for (int j = 0; j < n*2; j++) {
			cout << tmp[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << A[i][j] << " ";
		}
		cout << endl;
	}*/
}

vector<double> vectorMinusVector(vector<double> v1, vector<double> v2) {
	vector<double> result (v1.size(),0);
	for (int i = 0; i < v1.size(); i++) {
		result[i] = v1[i] - v2[i];
	}
	return result;
}

vector<double> vectorPlusVector(vector<double> v1, vector<double> v2) {
	vector<double> result(v1.size(), 0);
	for (int i = 0; i < v1.size(); i++) {
		result[i] = v1[i] + v2[i];
	}
	return result;
}

Grid::Grid() {}

Grid::Grid(vector<Node*> tabNodes, vector<Element*> tabElements, double K, double c, double t, double tot, double alfa, double ro, double time, double step) {
	
	this->tabNodes = tabNodes;
	this->tabElements = tabElements;
	this->K = K;
	this->c = c;
	this->t = t;
	this->tot = tot;
	this->alfa = alfa;
	this->ro = ro;
	this->time = time;
	this->step = step;

	for (int i = 0; i < tabNodes.size(); i++) {

		t0.push_back(0);
		t1.push_back(0);
	}
	vector<double> tmpCdT;
	vector<double> tmpHCdT;
	for (int i = 0; i < tabNodes.size(); i++) {
		tmpCdT.push_back(0);
		tmpHCdT.push_back(0);
	}
	for (int i = 0; i < tabNodes.size(); i++) {
		CdT.push_back(tmpCdT);
		HCdT.push_back(tmpHCdT);
	}
	

	generateUniversalElements();
	mapLocalToGlobal();
	generateEquationElements();


}

Grid::~Grid() {}

void Grid::generateUniversalElements() {
	for (int i = 0; i < tabElements.size(); i++) {
		UniversalElement *univElem = new UniversalElement(tabElements[i], ro, c, alfa, tot,K);
		tabUniversalElements.push_back(univElem);
	}
}

void Grid::mapLocalToGlobal() {
	for (int i = 0; i < tabNodes.size(); i++) {
		globalP.push_back(0);
		globalC.push_back(vector<double>(tabNodes.size(), 0));
		globalH.push_back(vector<double>(tabNodes.size(), 0));
	}
	for (int k = 0; k < tabElements.size(); k++) {
		for (int i = 0; i < 4; i++) {
			globalP[tabElements[k]->elemNodes[i]->id] += tabUniversalElements[k]->P[i];
			for (int j = 0; j < 4; j++) {
				globalH[tabElements[k]->elemNodes[i]->id][tabElements[k]->elemNodes[j]->id] += tabUniversalElements[k]->H[i][j];
				globalC[tabElements[k]->elemNodes[i]->id][tabElements[k]->elemNodes[j]->id] += tabUniversalElements[k]->C[i][j];
			}
		}
	}
}

void Grid::generateEquationElements() {
	for (int i = 0; i < globalC[1].size(); i++) {
		for (int j = 0; j < globalC[1].size(); j++) {
			CdT[i][j] = globalC[i][j] / step;
		}
	}
	for (int i = 0; i < globalH[1].size(); i++) {
		for (int j = 0; j < globalH[1].size(); j++) {
			HCdT[i][j] = globalH[i][j] + CdT[i][j];
		}
	}

	for (int i = 0; i < tabNodes.size(); i++) {
		t0[i] = (tabNodes[i]->t);
	}
	calculateInverse(HCdT);
}

void Grid::updateTemperatures(vector<double> t1) {
	for (int i = 0; i < tabNodes.size(); i++) {
		tabNodes[i]->t = t1[i];
	}
}

void Grid::solve() {
	for (int i = 0; i < tabNodes.size(); i++) {
		t0[i] = tabNodes[i]->t;
	}
	
	for (int i = 0; i <(time / step); i++) {
		
		vector<double> CdTT0P = vectorPlusVector(multiplyVectorMatrix(CdT,t0),globalP);
		t1 = multiplyVectorMatrix(HCdT,CdTT0P);

		cout << "-------------------Iteracja " << i + 1 << "-----------------" << endl;
		for (int j = 0; j < t1.size(); j++)
			cout << t1[j] << " ";
		cout << endl;
		auto min = min_element(t1.begin(), t1.end());
		auto max = max_element(t1.begin(), t1.end());
		cout << "Min: " << *min << "Max: " << *max << endl;
		updateTemperatures(t1);
		for (int j = 0; j < t1.size(); j++)
			t0[j] = t1[j];
	}	
}

void Grid::show(int it) {
	cout << "Iteracja: " << it << endl;
	cout << "H: " << endl;
	for (int i = 0; i < globalH[0].size(); i++) {
		for (int j = 0; j < globalH[0].size(); j++) {
			cout << globalH[i][j] <<" ";
		}
		cout << endl;
	}
	cout << "C: " << endl;
	for (int i = 0; i < globalC[0].size(); i++) {
		for (int j = 0; j < globalC[i].size(); j++) {
			cout << globalC[i][j] << " ";
		}
		cout << endl;
	}
}