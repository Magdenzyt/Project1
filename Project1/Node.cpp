#include "Node.h"

Node::Node() {
	this->id = 0;
	this->x = 0;
	this->y = 0;
	this->t = 0;
}

Node::Node(int id, double x, double y, double t) {
	this->id = id;
	this->bc = false;
	this->x = x;
	this->y = y;
	this->t = t;
}

Node::~Node() {}

void Node::show() {
	cout << "Node id:" <<id<<" x:" << x << " y:" << y << " t:" << t << endl;
}