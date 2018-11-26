#include "UniversalElement.h"
#include <math.h>

UniversalElement::UniversalElement(Element* myElement, double ro, double c, double alfa, double tot, double K) {

	this->element = myElement;
	this->ro = ro;
	this->c = c;
	this->alfa = alfa;
	this->tot = tot;
	this->K = K;

	this->ksi[0] = -1 / sqrt(3);
	this->ksi[1] = 1 / sqrt(3);
	this->ksi[2] = 1 / sqrt(3);
	this->ksi[3] = -1 / sqrt(3);

	this->eta[0] = -1 / sqrt(3);
	this->eta[1] = -1 / sqrt(3);
	this->eta[2] = 1 / sqrt(3);
	this->eta[3] = 1 / sqrt(3);

	generateN();
	generateUniversalXY();
	generatedNdKsidNdEta();
	generateJacobian();
	generateH();
	generateC();
	generateHBC();
	generateP();

}

UniversalElement::~UniversalElement() {}

void UniversalElement::generateN() {
	for (int i = 0; i < 4; i++) {
		this->N[0][i] = 0.25*(1 - ksi[i])*(1 - eta[i]);
		this->N[1][i] = 0.25*(1 + ksi[i])*(1 - eta[i]);
		this->N[2][i] = 0.25*(1 + ksi[i])*(1 + eta[i]);
		this->N[3][i] = 0.25*(1 - ksi[i])*(1 + eta[i]);
	}
}

void UniversalElement::generateUniversalXY() {
	for (int i = 0; i < 4; i++) {
		this->xP[i] = element->elemNodes[0]->x*N[0][i] + element->elemNodes[1]->x*N[1][i] + element->elemNodes[2]->x*N[2][i] + element->elemNodes[3]->x*N[3][i];
		this->yP[i] = element->elemNodes[0]->y*N[0][i] + element->elemNodes[1]->y*N[1][i] + element->elemNodes[2]->y*N[2][i] + element->elemNodes[3]->y*N[3][i];
	}
}

void UniversalElement::generatedNdKsidNdEta() {
	for (int i = 0; i < 4; i++) {
		this->dNdKsi[i][0] = -0.25*(1 - this->eta[i]);
		this->dNdKsi[i][1] = 0.25*(1 - this->eta[i]);
		this->dNdKsi[i][2] = 0.25*(1 + this->eta[i]);
		this->dNdKsi[i][3] = -0.25*(1 + this->eta[i]);
	}

	for (int i = 0; i < 4; i++) {
		this->dNdEta[i][0] = -0.25*(1 - this->ksi[i]);
		this->dNdEta[i][1] = -0.25*(1 + this->ksi[i]);
		this->dNdEta[i][2] = 0.25*(1 + this->ksi[i]);
		this->dNdEta[i][3] = 0.25*(1 - this->ksi[i]);
	}
}

void UniversalElement::generateJacobian() {
	double tabJ[4][4];

	for (int i = 0; i < 4; i++) {
		tabJ[i][0] = dNdKsi[i][0] * element->elemNodes[0]->x + dNdKsi[i][1] * element->elemNodes[1]->x + dNdKsi[i][2] * element->elemNodes[2]->x + dNdKsi[i][3] * element->elemNodes[3]->x;
		tabJ[i][1] = dNdKsi[i][0] * element->elemNodes[0]->y + dNdKsi[i][1] * element->elemNodes[1]->y + dNdKsi[i][2] * element->elemNodes[2]->y + dNdKsi[i][3] * element->elemNodes[3]->y;
		tabJ[i][2] = dNdEta[i][0] * element->elemNodes[0]->x + dNdEta[i][1] * element->elemNodes[1]->x + dNdEta[i][2] * element->elemNodes[2]->x + dNdEta[i][3] * element->elemNodes[3]->x;
		tabJ[i][3] = dNdEta[i][0] * element->elemNodes[0]->y + dNdEta[i][1] * element->elemNodes[1]->y + dNdEta[i][2] * element->elemNodes[2]->y + dNdEta[i][3] * element->elemNodes[3]->y;
	}

	for (int i = 0; i < 4; i++) {
		detJ[i] = tabJ[i][0] * tabJ[i][3] - tabJ[i][1] * tabJ[i][2];
	}

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			Jacobian[i][j] = tabJ[i][j] / detJ[i];
		}
	}
}

void UniversalElement::generateH() {
	double dNdX[4][4];
	double dNdY[4][4];

	double dNdXdNdXTdetJ[4][4][4];
	double dNdYdNdYTdetJ[4][4][4];
	double multiplyK[4][4][4];

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			dNdX[i][j] = Jacobian[i][0] * dNdKsi[i][j] + Jacobian[i][1] * dNdEta[i][j];
		}
	}

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			dNdY[i][j] = Jacobian[i][2] * dNdKsi[i][j] + Jacobian[i][3] * dNdEta[i][j];
		}
	}

	for (int k = 0; k < 4; k++) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				dNdXdNdXTdetJ[k][i][j] = dNdX[k][j] * dNdX[k][i] * detJ[k];
			}
		}
	}

	for (int k = 0; k < 4; k++) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				dNdYdNdYTdetJ[k][i][j] = dNdY[k][j] * dNdY[k][i] * detJ[k];
			}
		}
	}

	for (int k = 0; k < 4; k++) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				multiplyK[k][i][j] = K*(dNdXdNdXTdetJ[k][i][j] + dNdYdNdYTdetJ[k][i][j]);
			}
		}
	}

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			H[i][j] = multiplyK[0][i][j] + multiplyK[1][i][j] + multiplyK[2][i][j] + multiplyK[3][i][j];
		}
	}
}

void UniversalElement::generateC() {
	double cRoNNTdetJ[4][4][4];

	for (int k = 0; k < 4; k++) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				cRoNNTdetJ[k][i][j] = c*ro*detJ[k] * N[k][i] * N[k][j];
			}
		}
	}

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			C[i][j] = cRoNNTdetJ[0][i][j] + cRoNNTdetJ[1][i][j] + cRoNNTdetJ[2][i][j] + cRoNNTdetJ[3][i][j];
		}
	}
}

void UniversalElement::generateHBC() {

	int index0, index1;
	double pN[2][4];
	double pc0[4][4];
	double pc1[4][4];
	double sum[4][4];

	double points[8][2] = { { -1 / sqrt(3),-1 },
							{ 1 / sqrt(3),-1 },
							{ 1,-1 / sqrt(3) },
							{ 1, 1 / sqrt(3) },
							{ 1 / sqrt(3), 1 },
							{ -1 / sqrt(3),1 },
							{ -1,1 / sqrt(3) },
							{ -1,-1 / sqrt(3) } };

	for (int k = 0; k < 4; k++) {
		index0 = k;
		index1 = k + 1;
		if (index1 >= element->elemNodes.size())
			index1 = 0;
		length[k] = sqrt(pow(element->elemNodes[index0]->x - element->elemNodes[index1]->x, 2) + pow(element->elemNodes[index0]->y - element->elemNodes[index1]->y, 2));
		if (element->elemNodes[index0]->bc && element->elemNodes[index1]->bc) {
			pN[0][0] = 0.25*(1 - points[k * 2][0])*(1 - points[k * 2][1]);
			pN[0][1] = 0.25*(1 + points[k * 2][0])*(1 - points[k * 2][1]);
			pN[0][2] = 0.25*(1 + points[k * 2][0])*(1 + points[k * 2][1]);
			pN[0][3] = 0.25*(1 - points[k * 2][0])*(1 + points[k * 2][1]);

			pN[1][0] = 0.25*(1 - points[k * 2 + 1][0])*(1 - points[k * 2 + 1][1]);
			pN[1][1] = 0.25*(1 + points[k * 2 + 1][0])*(1 - points[k * 2 + 1][1]);
			pN[1][2] = 0.25*(1 + points[k * 2 + 1][0])*(1 + points[k * 2 + 1][1]);
			pN[1][3] = 0.25*(1 - points[k * 2 + 1][0])*(1 + points[k * 2 + 1][1]);

			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					pc0[i][j] = pN[0][j] * pN[0][i] * alfa;
				}
			}

			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					pc1[i][j] = pN[1][j] * pN[1][i] * alfa;
				}
			}

			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					sum[i][j] = (pc0[i][j] + pc1[i][j]) * (length[k] / 2);
				}
			}

			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					HBC[i][j] += sum[i][j];
				}
			}
		}
	}

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			H[i][j] = H[i][j] + HBC[i][j];
		}
	}

}

void UniversalElement::generateP() {
	int index0, index1;
	double pN[2][4];
	double points[8][2] = { { -1 / sqrt(3),-1 },
							{ 1 / sqrt(3),-1 },
							{ 1,-1 / sqrt(3) },
							{ 1, 1 / sqrt(3) },
							{ 1 / sqrt(3), 1 },
							{ -1 / sqrt(3),1 },
							{ -1,1 / sqrt(3) },
							{ -1,-1 / sqrt(3) } };
	for (int k = 0; k < 4; k++) {
		index0 = k;
		index1 = k + 1;
		if (index1 >= element->elemNodes.size())
			index1 = 0;
		if (element->elemNodes[index0]->bc && element->elemNodes[index1]->bc) {
			pN[0][0] = 0.25*(1 - points[k * 2][0])*(1 - points[k * 2][1]);
			pN[0][1] = 0.25*(1 + points[k * 2][0])*(1 - points[k * 2][1]);
			pN[0][2] = 0.25*(1 + points[k * 2][0])*(1 + points[k * 2][1]);
			pN[0][3] = 0.25*(1 - points[k * 2][0])*(1 + points[k * 2][1]);

			pN[1][0] = 0.25*(1 - points[k * 2 + 1][0])*(1 - points[k * 2 + 1][1]);
			pN[1][1] = 0.25*(1 + points[k * 2 + 1][0])*(1 - points[k * 2 + 1][1]);
			pN[1][2] = 0.25*(1 + points[k * 2 + 1][0])*(1 + points[k * 2 + 1][1]);
			pN[1][3] = 0.25*(1 - points[k * 2 + 1][0])*(1 + points[k * 2 + 1][1]);

			for (int i = 0; i< 2; i++) {
				for (int j = 0; j < 4; j++) {
					P[j] += (pN[i][j] * alfa * tot * (length[j] / 2.0));
				}
			}

		}

	}
}