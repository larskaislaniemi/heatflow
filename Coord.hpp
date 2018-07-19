#ifndef _COORD_HPP
#define _COORD_HPP
#include "Settings.hpp"

class Coord {
	public:
	int *coords;	// values of the coordinates
	Settings *S;
	
	void print();
	void copyfrom(Coord *C);
	
	Coord(Settings *S);
	~Coord();
};

#endif
