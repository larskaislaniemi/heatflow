#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include "Settings.hpp"
#include "Coord.hpp"
#include "Grid.hpp"
#include "Error.hpp"

using namespace std;

Settings::Settings() {
	this->readSettings();
	this->dimnames[0] = 'z';
	this->dimnames[1] = 'x';
	this->dimnames[2] = 'y';
	
	if(this->dims > MAXDIMS) throw Error();

	this->cur_tstep = -1;
	cur_output = &cout;
}

Settings::~Settings() {
	delete par_Tfield;
	delete par_denfield;
	delete par_cpfield;
	delete par_kfield;
	delete par_beta1field;
	delete par_beta2field;
	delete par_beta3field;
	delete par_hfield;

}

void Settings::readValue(int *to) {
	string str;
	char ch;
	while (cin.get(ch)) {
		if (ch == ' ' || ch < 28) {
			cin.putback(ch);
			getline(cin, str);
		} else {
			cin.putback(ch);
			getline(cin, str);
			istringstream iss(str);
			iss >> (*to) >> std::ws;
			return;
		}
	}
}

void Settings::readValue(double *to) {
	string str;
	char ch;
	while (cin.get(ch)) {
		if (ch == ' ' || ch < 28) {
			cin.putback(ch);
			getline(cin, str);
		} else {
			cin.putback(ch);
			getline(cin, str);
			istringstream iss(str);
			iss >> (*to) >> std::ws;
			return;
		}
	}
}

void Settings::readSettings() {
	int areas;
	int Tinifile;
	double *xmins;
	double *xmaxs;
	double *T0, *k0, *den0, *cp0, *h0, *beta1, *beta2, *beta3;
	
	Tinifile = 0;

	cerr << "Dims" << endl;
	readValue(&(this->dims));
	cerr << "Got " << this->dims << endl;
	if (dims < 1 || dims > MAXDIMS) throw ConfigError();

	cerr << "nxs" << endl;
	for (int i = 0; i < dims; i++) {
		readValue(&(this->nx[i]));
	}
	
	cerr << "Ls" << endl;
	for (int i = 0; i < dims; i++) {
		readValue(&(this->L[i]));
	}
	
	for (int i = 0; i < this->dims; i++) {
		this->dx[i] = this->L[i] / this->nx[i];
	}

	cerr << "dtcoeff" << endl;
	readValue(&(this->dtcoeff));
	
	cerr << "output interval" << endl;
	readValue(&(this->outputInterval));
	
	cerr << "runtime" << endl;
	readValue(&(this->runtime));
	this->runtime = this->runtime * (double)(52. * 7. * 24. * 60. * 60.);
	
	cerr << "bctypes" << endl;
	for (int i = 0; i < this->dims; i++) {
		cerr << "left" << endl;
		readValue(&(this->bctypes[2*i]));
		cerr << "right" << endl;
		readValue(&(this->bctypes[2*i+1]));
	}
	
	cerr << "bcvalues" << endl;
	for (int i = 0; i < this->dims; i++) {
		cerr << "left" << endl;
		readValue(&(this->bcvalues[2*i]));
		cerr << "right" << endl;
		readValue(&(this->bcvalues[2*i+1]));
	}

	cerr << "erosion speed" << endl;
	for (int i = 0; i < this->dims; i++) {
		cerr << "left" << endl;
		readValue(&(this->erosionspeed[2*i]));
		cerr << "right" << endl;
		readValue(&(this->erosionspeed[2*i+1]));
	}	
	
	cerr << "Tini from file?" << endl;
	readValue(&(Tinifile));
	
	cerr << "areas" << endl;
	readValue(&(areas));
	cerr << areas << endl;
	
	xmins = new double[dims*areas];
	xmaxs = new double[dims*areas];
	T0 = new double[areas];
	k0 = new double[areas];
	den0 = new double[areas];
	cp0 = new double[areas]; 
	h0 = new double[areas];
	beta1 = new double[areas];
	beta2 = new double[areas];
	beta3 = new double[areas];
	
	for (int i = 0; i < areas; i++) {
		for (int j = 0; j < dims; j++) {
			cerr << "xmin for area " << i << ", dimension " << j << endl;
			readValue(&(xmins[i*dims + j]));
			cerr << "xmax for area " << i << ", dimension " << j << endl;
			readValue(&(xmaxs[i*dims + j]));
		}
		cerr << "T0 for area " << i << endl;
		readValue(&(T0[i]));
		cerr << "k0 for area " << i << endl;
		readValue(&(k0[i]));
		cerr << "den0 for area " << i << endl;
		readValue(&(den0[i]));
		cerr << "cp0 for area " << i << endl;
		readValue(&(cp0[i]));
		cerr << "h0 for area " << i << endl;
		readValue(&(h0[i]));
		cerr << "beta1 for area " << i << endl;
		readValue(&(beta1[i]));
		cerr << "beta2 for area " << i << endl;
		readValue(&(beta2[i]));
		cerr << "beta3 for area " << i << endl;
		readValue(&(beta3[i]));
	}
	
	this->par_kfield = new Grid(this);
	this->par_Tfield = new Grid(this);
	this->par_denfield = new Grid(this);
	this->par_cpfield = new Grid(this);
	this->par_hfield = new Grid(this);
	this->par_beta1field = new Grid(this);
	this->par_beta2field = new Grid(this);
	this->par_beta3field = new Grid(this);
	this->cur_Tfield = this->par_Tfield;
	
	if (Tinifile > 0) {
		this->par_Tfield->readfromcsv((char *)TINIFILENAME);
	}

	for (int i = 0; i < areas; i++) {
		Coord cFrom(this);
		Coord cTo(this);
		for (int j = 0; j < dims; j++) {
			cFrom.coords[j] = floor(xmins[i*dims + j] / this->dx[j]);
			cTo.coords[j] = floor(xmaxs[i*dims + j] / this->dx[j]);
		}
		if (T0[i] >= 0) {
			this->par_Tfield->set(&cFrom, &cTo, T0[i]);
		}
		this->par_kfield->set(&cFrom, &cTo, k0[i]);
		this->par_denfield->set(&cFrom, &cTo, den0[i]);
		this->par_cpfield->set(&cFrom, &cTo, cp0[i]);
		this->par_hfield->set(&cFrom, &cTo, h0[i]);
		this->par_beta1field->set(&cFrom, &cTo, beta1[i]);
		this->par_beta2field->set(&cFrom, &cTo, beta2[i]);
		this->par_beta3field->set(&cFrom, &cTo, beta3[i]);
	}
	
	delete T0;
	delete k0;
	delete den0;
	delete cp0;
	delete h0;
	delete beta1;
	delete beta2;
	delete beta3;
}
