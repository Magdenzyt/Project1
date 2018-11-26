#pragma once
#include<iostream>
using namespace std;
class Node {
public:

	int id;
	bool bc;
	double x, y, t;

	Node();
	Node(int id, double x, double y, double t);
	~Node();

	void show();
};