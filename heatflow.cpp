/* Heatflow -- program to solve the heat equation in 1/2/3d
 * with variable and temperature dependent physical properties.
 * 
 * Heat equation is solved with finite difference method
 * (forward euler in time), discretized to account for 
 * temperature dependent material properties (see H.D.Baehr &
 * K. Stephan, Heat and Mass Transfer, Ch 2.4.4).
 * 
 * 2012-05-07 Lars Kaislaniemi <lars.kaislaniemi@iki.fi>
 * 
 * 
 * Version history:
 * 
 * 2012-05-11
 *
 * Added erosion.
 *
 * 2012-05-07
 * 
 * First version. Tested in 1D for consistency (maintains heat
 * flow; heat production sums up to total heat flow; etc.). Tested
 * in 2D qualitatively.
 * 
 * 
 * TODO: Benchmark 1D with k=k(T) against analytical solutions
 * 
 * 
 * Distributed under BSD license. See LICENSE.TXT for information.
 * 
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <stdlib.h>

#include "Settings.hpp"
#include "Coord.hpp"
#include "Grid.hpp"
#include "Error.hpp"


using namespace std;


double timestep(Grid *den, Grid *cp, Grid *k, Settings *S) {
	/* Return the largest possible stable timestep.
	 * Very conservative, uses global minimum and 
	 * maximum values of k to find the minimum/maximum
	 * ratios of k[i+1]/k[i]
	 */
	double dt;
	double kmin, kmax;
	
	for (int dir = 0; dir < S->dims; dir++) {
		double thisdt;
		kmin = k->min();
		kmax = k->max();

		//cerr << "kmin / kmax are " << kmin << " / " << kmax << endl;
		
		if (kmin == 0) {
			thisdt = S->dx[dir] * S->dx[dir] * cp->min() * den->min() / 
				(2 * kmax * (1/(1+kmin/kmax)));
		} else if (kmax == 0) {
			throw Error();
		} else {			
			thisdt = S->dx[dir] * S->dx[dir] * cp->min() * den->min() / 
				(2 * kmax * (1/(1+kmin/kmax) + 1/(1+kmax/kmin)));
		}
					
		if (dir == 0) dt = thisdt;
		else if (thisdt < dt) dt = thisdt;
	}
	
	return S->dtcoeff * dt;
}

double kT(double k0, double b1, double b2, double b3, double T) {
	double k;
	/* The function k(T) */
	k = b1*k0 + b2*T + b3*T*T;
	return k;
}

int kf(Grid *g, Coord *C) {
	/* Updates grid g with k(T) values
	 * at current temperature 
	 */
	if (g->S->cur_Tfield == 0) throw NoFieldError();
	if (g->S->par_kfield == 0) throw NoFieldError();
	if (g->S->par_beta1field == 0) throw NoFieldError();
	if (g->S->par_beta2field == 0) throw NoFieldError();
	if (g->S->par_beta3field == 0) throw NoFieldError();
	
	g->set(C, kT(g->S->par_kfield->at(C), g->S->par_beta1field->at(C), 
	g->S->par_beta2field->at(C), g->S->par_beta3field->at(C), g->S->cur_Tfield->at(C)));
	
	return 0;
}

void heatdiff(int dim, Grid *T, Grid *Tn, Grid *k, Grid *den, Grid *cp, Settings *S, Coord *C = 0) {
	/* Calculates next timestep of heat diffusion equation (without the 
	 * heat source component). This is recursive to account for each direction --
	 * calculations are made only on the "lowest" level.
	 */
	int delC = 0;
	if (C == 0) {
		C = new Coord(S);
		delC = 1;
	}
	for (int i = 0; i < S->nx[dim]; i++) {
		C->coords[dim] = i;
		if (dim == S->dims-1) {
			int isbound = 0;

			for (int dir = 0; dir < S->dims; dir++) {
				if (C->coords[dir] == 0 || C->coords[dir] == S->nx[dir]-1) {
					int indexadd;
					// boundary
					isbound = 9999;
					if (C->coords[dir] == 0) {
						isbound = 1;	// "left"
						indexadd = 0;
					} else {
						isbound = 2;	// "right"
						indexadd = 1;
					}
					if (S->bctypes[2*dir + indexadd] == 1) {
						// constant T
						Tn->set(C, S->bcvalues[2*dir + indexadd]);
					} else if (S->bctypes[2*dir + indexadd] == 2) {
						// constant flux
						
						/* to account for variable k=k(T)
						 * the B.C. has been extended:
						 * (i.e., on left boundary:)
						 * 
						 *  {1}  T[i+1] = - q * dz / k(T[i]) + T[i] ,
						 * where
						 *  {2}  k(T[i]) ~= k[i] + dk/dT|i * (T[i+1] - T[i])
						 *
						 * (Here T[i+1] refers to the temperature of any chosen 
						 * point close to point i, and the temperature field
						 * between these two points is assumed linear.)
						 *
						 * Let DT = T[i+1] - T[i]
						 * 
						 * so {1} and {2} become
						 * 
						 *  {3} DT * k[i] + DT^2 * dk/dT|i + q * dz = 0
						 * 
						 * which is a quadratic equation and can be solved 
						 * to find the DT (and thus T[i] when T[i+1] is
						 * known)
						 * 
						 * If k(T) is constant in respect to T, the 
						 * "traditional" formulation may be used:
						 * 
						 *  {4} T[i+1] = -q * dz / k + T[i]
						 */
						
						Coord *Cside = new Coord(S);
						Cside->copyfrom(C);
						if (isbound == 1) {
							// left
							Cside->coords[dir] = Cside->coords[dir] + 1;
						} else {
							// right
							Cside->coords[dir] = Cside->coords[dir] - 1;
						}

						double T0, T1, dT1, dT2, dTf, discr, T2;
						double ki, ki1, q, dz, dkdt;
						
						T0 = T->at(C);
						T1 = T->at(Cside);
						q = S->bcvalues[2*dir + indexadd];
						dz = S->dx[dir];
						
						if (T1-T0 > 0) {
							T2 = T0 + 1;
						} else {
							T2 = T0 - 1;
						}
						ki = kT(S->par_kfield->at(C), S->par_beta1field->at(C), 
						     S->par_beta2field->at(C), S->par_beta3field->at(C), T0);
						ki1 = kT(S->par_kfield->at(C), S->par_beta1field->at(C), 
						      S->par_beta2field->at(C), S->par_beta3field->at(C), T2);
						dkdt = (ki1-ki) / (T2 - T0);
						
						discr = ki*ki - 4 * dkdt * q * dz;
						
						if (dkdt != 0) {
							if (discr < 0) throw NoSolutionError();
						
							dT1 = (-ki + sqrt(discr)) / (2 * dkdt);
							dT2 = (-ki - sqrt(discr)) / (2 * dkdt);
							
							if (abs(dT1) < abs(dT2)) dTf = dT1;
							else dTf = dT2;
							
							//cout << "dT1 = " << dT1 << ", dT2 = " << dT2 << ", dTf = " << dTf << endl;
							if (isbound == 1) {
								// left
								Tn->set(C, T->at(Cside) + dTf);
							} else {
								// right
								Tn->set(C, T->at(Cside) - dTf);
							}
						} else {
							// k is (locally) constant (with respect to T)
							// use constant T solution for B.C.
							if (isbound == 1) {
								Tn->set(C, q * dz / ki + T->at(Cside));
							} else {
								Tn->set(C, -q * dz / ki + T->at(Cside));
							}
						}

						delete(Cside);
					} else if (S->bctypes[2*dir + indexadd] == 3) {
						// constant T with eroding surface
						double surf_level, surf_level_resid, eff_bcvalue;
						int surf_level_dx;
						double T1, T2;
						Coord *C0 = new Coord(S);
						double erosion_speed = S->erosionspeed[2*dir + indexadd];
						erosion_speed = erosion_speed / (1e6 * 52. * 7. * 24. * 60. * 60.);
						
						C0->copyfrom(C);
						
						surf_level = S->time * erosion_speed;
						
						surf_level_dx = (int)floor(surf_level / S->dx[dir]);
						surf_level_resid = surf_level - floor(surf_level / S->dx[dir]);
						
						if (isbound == 1) {
							for (int i = 0; i <= surf_level_dx; i++) {
								C0->coords[dir] = i;
								Tn->set(C0, S->bcvalues[2*dir + indexadd]);
								T->set(C0, S->bcvalues[2*dir + indexadd]);
							}
						} else {
							for (int i = S->nx[dir]-1; i >= S->nx[dir]-1 - surf_level_dx; i--) {
								C0->coords[dir] = i;
								Tn->set(C0, S->bcvalues[2*dir + indexadd]);
								T->set(C0, S->bcvalues[2*dir + indexadd]);
							}
						}
						
						/* Check this interpolation of BC value - now disabled *
						Coord *X0 = new Coord(S);*/
						Coord *X1 = new Coord(S);
						/*Coord *X2 = new Coord(S);
						X0->copyfrom(C0);
						X1->copyfrom(C0);
						X2->copyfrom(C0);
						if (isbound == 1) {
							// left
							X1->coords[dir] = X1->coords[dir] + 1;
							X2->coords[dir] = X2->coords[dir] + 2;
						} else {
							// right
							X1->coords[dir] = X1->coords[dir] - 1;
							X2->coords[dir] = X2->coords[dir] - 2;
						}						
						T1 = S->cur_Tfield->at(X1);
						T2 = S->cur_Tfield->at(X2);
						*/
						eff_bcvalue = S->bcvalues[2*dir + indexadd];
						//eff_bcvalue = -(T2 - S->bcvalues[2*dir + indexadd])*(S->dx[dir]) / (S->dx[dir] + surf_level_resid) + T2;
						//cerr << "effective bc value = " << eff_bcvalue << endl;
						
						
						Tn->set(X1, eff_bcvalue);
						
						/*delete C0;
						delete X0;*/
						delete X1;
						/*delete X2;*/
					}
				}
			}

			if (isbound == 0) {
				// not at the boundary at all,
				// calculate terms in all directions
				for (int dir = 0; dir < S->dims; dir++) {
					double thisterm = 0;
					Coord *Cleft = new Coord(S);
					Coord *Cright = new Coord(S);
					Cleft->copyfrom(C);
					Cright->copyfrom(C);
					Cleft->coords[dir] = C->coords[dir] - 1;
					Cright->coords[dir] = C->coords[dir] + 1;
					thisterm = 2 * k->at(C) * S->dt / (S->dx[dir]*S->dx[dir]);
					thisterm /= (den->at(C) * cp->at(C));
					thisterm *= ((T->at(Cright) - T->at(C)) / (1 + k->at(C) / k->at(Cright)) - 
								 (T->at(C) - T->at(Cleft)) / (1 + k->at(C) / k->at(Cleft)) );
					Tn->set(C, Tn->at(C) + thisterm); 
					delete(Cleft);
					delete(Cright);
				}
			}

			if (isbound == 0) Tn->set(C, Tn->at(C) + T->at(C));

		} else {
			heatdiff(dim+1, T, Tn, k, den, cp, S, C);
		}
		
	}
	if (delC) delete C;
}

void heatprod(int dim, Grid *T, Grid *den, Grid *cp, Settings *S, Coord *C = 0) {
	/* Add the internal heat production to the given node C */
	
	int delC = 0;
	if (C == 0) {
		C = new Coord(S);
		delC = 1;
	}
	for (int i = 0; i < S->nx[dim]; i++) {
		C->coords[dim] = i;
		if (dim == S->dims-1) {
			int isbound = 0;

			for (int dir = 0; dir < S->dims; dir++) {
				if (C->coords[dir] == 0) {
					// "left" boundary
					isbound = 1;
				}
				if (C->coords[dir] == S->nx[dir]-1) {
					// "right" boundary
					isbound = 1;
				}
			}

			if (isbound == 0) {
				double h = T->S->par_hfield->at(C);
				double c = T->S->par_cpfield->at(C);
				double d = T->S->par_denfield->at(C);
				T->set(C, T->at(C) + T->S->dt * h / (c*d));
			}

		} else {
			heatprod(dim+1, T, den, cp, S, C);
		}
		
	}
	if (delC) delete C;
}


int printcsv(Grid *g, Coord *C) {
	/* Print the location and value of the given node C
	 * in the grid g to the current output.
	 * Used by the file writing procedures.
	 */
	for (int i = 0; i < g->S->dims; i++) {
		*(g->S->cur_output) << C->coords[i] << ",";
	}
	*(g->S->cur_output) << g->at(C) << endl;
	
	return 0;
}


int main(int argc, char **argv) {
	Settings *S = new Settings();
	
	Grid k(S);	
	Grid den(S);
	Grid cp(S);
	Grid T(S);
	Grid Tn(S);
	
	k.set((double)0);
	T.set(S->par_Tfield);
	Tn.set((double)0);
	
	S->time = 0;
	
	cerr << "Running for " << S->runtime / (double)(52*7*24*60*60) << " years" << endl;
	
	while (S->time < S->runtime) {
		S->cur_tstep++;
		S->cur_Tfield = &T;

		// add heat production to each node
		heatprod(0, &T, &den, &cp, S);
		
		// copy the k0 values from settings
		// to the current k grid and then
		// calculate the values of k=k(T)
		k.set(S->par_kfield);
		k.evaleach(&kf);

		// same for density
		den.set(S->par_denfield);
		// NOT IMPLEMENTED k.evaleach(&denf);

		// same for heat capacity
		cp.set(S->par_cpfield);
		// NOT IMPLEMENTED k.evaleach(&cpf);
		

		// find the stable timestep using
		// updated values of (den, cp, ) k 
		// and T
		S->dt = timestep(&den, &cp, &k, S);

		cerr << "\r" << 100.*(double)S->time / (double)S->runtime << " %, timestep " << S->cur_tstep << " (dt = " << S->dt / (double)(60*60*24*52*7) << " years)\t\t";

		// calculate the diffusion
		heatdiff(0, &T, &Tn, &k, &den, &cp, S);
		
		S->time = S->time + S->dt;
		for (int i = 0; i < Tn.size; i++) {
			T.set(i, Tn.at(i));
		}
		Tn.set((double)0);
		if (S->cur_tstep % S->outputInterval == 0) {
			T.savetofile((char*)("T"), &printcsv);
			k.savetofile((char*)("k"), &printcsv);
		}
	}
	
	delete S;
	
	cout << endl;
	return 0;
}


