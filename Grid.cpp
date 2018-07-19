#include <iostream>
#include <sstream>
#include <fstream>
#include "Settings.hpp"
#include "Coord.hpp"
#include "Error.hpp"
#include "Grid.hpp"

using namespace std;

void Grid::print() {
	if (this->S->dims > 2) cerr << "Showing dimension 0 and 1." << endl;
	cout << endl;
	for (int i = 0; i < this->S->nx[0]; i++) {
		for (int j = 0; j < this->S->nx[1]; j++) {
			Coord c = Coord(this->S);
			c.coords[0] = i;
			c.coords[1] = j;
			cout << this->at(&c) << " ";
		}
		cout << endl;
	}
}

double Grid::diff(int *coord, int dir, int order) {
	double ret;
	
	if (dir >= this->S->dims) throw Error();
	if (coord[dir] <= 0) throw IndexError();
	if (coord[dir] >= this->S->nx[dir]-1) throw IndexError();
	
	if (order == 1) {
		ret = (this->val[coord[dir]+1] - this->val[coord[dir]-1]) / this->S->dx[dir];
	} else if (order == 2) {
		ret = (this->val[coord[dir]+1] - 2*this->val[coord[dir]] + this->val[coord[dir]-1]) / (this->S->dx[dir] * this->S->dx[dir]);
	} else {
		throw Error();
	}
	
	return ret;
}

int Grid::index(Coord *coord) {
	int index = 0;
	
	// index = x0 + x1*nx0 + x2*nx0*nx1 + ...
	for (int i = 0; i < this->S->dims; i++) {
		int thisterm;
		thisterm = coord->coords[i];
		for (int j = i-1; j >= 0; j--) {
			thisterm = thisterm * this->S->nx[j];
		}
		index = index + thisterm;
	}
	
	if (index >= size) throw IndexError();
	
	return index;
}

void Grid::set(Coord *cFrom, Coord *cTo, double val, int dim, Coord *C) {
	int delC = 0;
	if (C == 0) {
		delC = 1;
		C = new Coord(this->S);
		cout << "Setting all from (";
		cFrom->print();
		cout << ") to (";
		cTo->print();
		cout << ") to value " << val << endl;
	}
	for (int i = 0; i < this->S->nx[dim]; i++) {
		C->coords[dim] = i;
		if (dim == this->S->dims-1) {
			int inside = 1;
			for (int j = 0; j < this->S->dims; j++) {
				if (C->coords[j] > cTo->coords[j] || 
					C->coords[j] < cFrom->coords[j]) {
					inside = 0;
					break;
				}
			}
			if (inside > 0) {
				this->set(C, val);
			}
		} else {
			this->set(cFrom, cTo, val, dim+1, C);
		}
	}	
	if (delC) delete C;
}

void Grid::set(double val) {
	for (int i = 0; i < this->size; i++) {
		this->val[i] = val;
	}
}

void Grid::set(Coord *coord, double val) {
	this->val[this->index(coord)] = val;
}

void Grid::set(int pos, double val) {
	this->val[pos] = val;
}

double Grid::at(Coord *coord) {
	return this->val[this->index(coord)];
}

double Grid::at(int pos) {
	return this->val[pos];
}

double Grid::max() {
	double m = 0.0;
	for (int i = 0; i < this->size; i++) {
		if (this->val[i] > m) m = this->val[i];
	}
	return m;
}

double Grid::min() {
	double m = this->val[0];
	for (int i = 1; i < this->size; i++) {
		if (this->val[i] < m) m = this->val[i];
	}
	return m;
}

Grid::Grid(Settings *S) {
	this->S = S;
	this->size = 1;
	for (int i = 0; i < S->dims; i++) {
		this->size = this->size * S->nx[i];
	}
	this->val = new double[this->size];
}

Grid::~Grid() {
	delete val;
}


int Grid::evaleach(int (*f)(Grid *g, Coord *C), Coord *C, int dim) {
	int ret = 0;
	int delC = 0;
	if (C == 0) {
		C = new Coord(this->S);
		delC = 1;
	}
	for (int i = 0; i < this->S->nx[dim]; i++) {
		C->coords[dim] = i;
		if (dim == this->S->dims-1) {
			ret = ret | f(this, C);
		} else {
			ret = (ret | evaleach(f, C, dim+1));
		}
	}	
	if (delC) delete C;
	return ret;
}

void Grid::readfromcsv(char *filename) {
	Coord *C;
	C = new Coord(this->S);
	
	int linesest, linesread;
	linesest = 1;
	for (int i = 0; i < this->S->dims; i++) {
		linesest *= this->S->nx[i];
	}
	cerr << "Expecting " << linesest << " lines" << endl;
	linesread = 0;
	
	ifstream infile(filename); 
	string line = "";
	while (getline(infile, line)) {
		int dim;
		int done;
		double val;
		stringstream strstr(line);
		string word = "";
		linesread++;

		dim = -1;
		done = 0;
		while (getline(strstr,word, ',')) {
			dim++;
			if (dim < this->S->dims) {
				// coord field
				istringstream(word) >> C->coords[dim];
			} else if (dim == this->S->dims) {
				// value field
				istringstream(word) >> val;
				this->set(C, val);
				done = 1;
			} else {
				cerr << "Too many fields at the row." << endl;
				throw Error();
			}
		}
	}
	cerr << "Read " << linesread << " (expected " << linesest << ")" << endl;
	if (linesread < linesest) throw Error();
	infile.close();
	
	delete C;
}

void Grid::savetofile(char *fieldname, int (*lineformatfunction)(Grid *g, Coord *C)) {
	char *filename = 0;
	ofstream file;
	ostream *store_prev_ostream;
	
	filename = new char[1024];
	sprintf(filename, "%s.%d.csv", fieldname, this->S->cur_tstep);
	
	file.open(filename);
	if (file.is_open()) {
		store_prev_ostream = this->S->cur_output;
		this->S->cur_output = &file;

		Coord C(S);
		this->evaleach(lineformatfunction, &C);
	
		this->S->cur_output = store_prev_ostream;
	
		file.close();
	} else {
		throw Error();
	}
	delete filename;
}

void Grid::set(Grid *g) {
	if (g == 0) throw Error();
	if (this->size != g->size) throw Error();
	for (int i = 0; i < g->size; i++) {
		this->val[i] = g->val[i];
	}
}
