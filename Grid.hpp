#ifndef _GRID_HPP
#define _GRID_HPP
class Settings;

#include <iostream>
#include <fstream>
#include <string>
#include "Settings.hpp"
#include "Coord.hpp"


class Grid {
	public:
	double *val;	// actual values of the nodes
	Settings *S;
	int size;		// total size (nx1 * nx2 * nx3 * ...)
	
	Grid(Settings *S);
	~Grid();
	
	// at() : return value at given position
	double at(Coord *coord);
	double at(int pos);
	
	// set() : set value at given position
	void set(Coord *coord, double val);
	void set(int pos, double val);
	void set(double val);
	void set(Coord *cFrom, Coord *cTo, double val, int dim = 0, Coord *C = 0);
	void set(Grid *g);
	
	// convert Coord to a location in array val
	int index(Coord *coord);
	
	// return global maximum and minimum value
	double max();
	double min();
	
	double diff(int *coord, int dir, int order = 1);
	void print();
	
	// Save the contents of the grid into a file.
	// Format may be modified by lineformatfunction 
	// (see printcsv function for example). 
	// N.B. Reading initial values assumes  that 
	// lineformatfunction printcsv() has been used.
	void savetofile(char *fieldname, int (*lineformatfunction)(Grid *g, Coord *C));
	
	// apply function f to each node in the grid,
	// starting from coordinate C, dimension level dim
	// (default: apply to all)
	int evaleach(int (*f)(Grid *g, Coord *C), Coord *C = 0, int dim = 0);
	
	int printcsv(Coord *C);
	
	// read the contents of a file into the grid
	void readfromcsv(char *filename);
};
#endif
