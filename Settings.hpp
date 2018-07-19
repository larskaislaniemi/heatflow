#ifndef _SETTINGS_HPP
#define _SETTINGS_HPP
#include "Grid.hpp"

#define MAXDIMS 3
#define TINIFILENAME "Tini.csv"

class Settings {
	public:
	int dims;		// dimensions of the problem
	char dimnames[MAXDIMS];	// not used currently
	int nx[MAXDIMS];		// number of nodes for each dim
	double L[MAXDIMS];		// physical size of each dim
	double dx[MAXDIMS];		// dx for each dim
	int bctypes[2*MAXDIMS];	// types of B.C. for each dim, both boundaries 
	double bcvalues[2*MAXDIMS];	// values of B.C. for each dim, both boundaries 
	double erosionspeed[2*MAXDIMS];	// erosion speed for each dim, both boundaries 
	double dtcoeff;			// dt coefficient (1.0 = stability criterion step size)
	double dt;				// current effective timestep size
	double runtime;			// total model run time (sec)
	double time;			// current time after model start (sec)
	int outputInterval;		// how often results are outputted to files
	
	int cur_tstep;			// current timestep (0,1,2,3,...)
	std::ostream *cur_output;	// pointer to current output stream (std::ostream), do not delete!
	Grid *cur_Tfield;	// pointer to current temperature field (do not delete!)
	Grid *par_kfield;	// parameter field for k0 values 
	Grid *par_Tfield;	// parameter field (initial values) for temperature
	Grid *par_denfield;	// parameter field for density0 values 
	Grid *par_cpfield;	// parameter field for cp0 values 
	Grid *par_hfield;	// parameter field for heat production
	Grid *par_beta1field;// parameter field for beta1 (in k(T)) values
	Grid *par_beta2field;// parameter field for beta2 (in k(T)) values
	Grid *par_beta3field;// parameter field for beta3 (in k(T)) values
	
	void readSettings();
	void readValue(int *to);
	void readValue(double *to);
	
	Settings();
	
	~Settings();
};

	
#endif
