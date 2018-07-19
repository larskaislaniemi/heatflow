#include <iostream>

#include "Settings.hpp"
#include "Coord.hpp"
#include "Error.hpp"
#include "Grid.hpp"

using namespace std;

Coord::Coord(Settings *S) {
	this->S = S;
	this->coords = new int[S->dims];
	for (int dir = 0; dir < S->dims; dir++) {
		this->coords[dir] = 0;
	}
}

Coord::~Coord() {
	delete this->coords;
}

void Coord::copyfrom(Coord *C) {
	if (C->S->dims != this->S->dims) throw Error();
	for (int i = 0; i < C->S->dims; i++) {
		this->coords[i] = C->coords[i];
	}
}

void Coord::print() {
	for (int i = 0; i < this->S->dims; i++) {
		cout << this->coords[i] << "\t";
	}
}

