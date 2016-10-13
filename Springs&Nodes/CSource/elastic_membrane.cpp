/*****************************************************************************
 * 2D ELASTIC CELL MEMBRANE WITH RESTING TENSION AND SURFACE ADHESION				 *
 * TIME INTEGRATOR: FORWARD EULER (1st ORDER METHOD)                         *
 *                                                                           *
 * C.Copos 10/01/2016                                                        *
 ****************************************************************************/

#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <limits> 
#include <time.h> 
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <getopt.h>
#include <stdio.h>
using namespace std;

#include "structs_mem.h"

// TIME DISCRETIZATION
const double TMAX = 10.0;
const double tstep = 0.0001;
static int stop = floor(TMAX/tstep);
//static int stop = 1; /* RUNNING STATIC TESTS */

// CELL PARAMETERS
const double k_cell = 1000;
const double xi = 0.05;

// MEMBRANE/CORTEX PARAMETERS
double k_mem = 1.0; 					// membrane spring stiffness (units of force)
double gamma_mem = 1.0; 			// resting tension in membrane (units of force) 

// ADHESION BOND PARAMETERS
double k_adh = 0.5; 				
double crit_adh = 1.0;
const double k_contact = 100.0;

// SURFACE PARAMETERS
const double substrate_height = 0.1;
const double attach_dl0 = 0.4;

//---------------------------------------------------------------------------------------------------------//

// COMPUTE DIFFERENCE IN CLOCKTIME IN SECONDS
double diffclock(clock_t s, clock_t e) {
    double diffticks = s-e;
    double diffms = (diffticks)/CLOCKS_PER_SEC;
    return diffms;
}
// END COMPUTE DIFF IN CLOCKTIME IN SECONDS

// COMPUTE AREA OF CELL
double computeArea(int Npts, vertex Nodes[]) {
	double area = 0.;	
	
	int i, j=Npts-1;
	for(i=0; i<Npts; i++) {
		area += -(Nodes[j].def.x + Nodes[i].def.x)*(Nodes[j].def.y-Nodes[i].def.y);
		j = i; 
	}
  return area*0.5;
}
// END COMPUTE AREA

// COMPUTE ELASTIC FORCES
void computeElasticForces(int Npts, vertex Nodes[], vector ef[]) {
	int k, left, right;
	double ref_dist_L, ref_dist_R, def_dist_L, def_dist_R;
	double ds;
	vector t_L, t_R; // tangent vectors

	for(k=0; k<Npts; k++) {
		// setup neighbors and membrane length information
		if(k==0) { left = Npts-1; right = k+1; }
    else if(k==Npts-1) { left = k-1; right = 0; }
    else { left = k-1; right = k+1; }

    ref_dist_L = sqrt( pow(Nodes[k].ref.x - Nodes[left].ref.x,2.0) + pow(Nodes[k].ref.y - Nodes[left].ref.y,2.0) );
    ref_dist_R = sqrt( pow(Nodes[k].ref.x - Nodes[right].ref.x,2.0) + pow(Nodes[k].ref.y - Nodes[right].ref.y,2.0) );
    def_dist_L = sqrt( pow(Nodes[k].def.x - Nodes[left].def.x,2.0) + pow(Nodes[k].def.y - Nodes[left].def.y,2.0) );
    def_dist_R = sqrt( pow(Nodes[k].def.x - Nodes[right].def.x,2.0) + pow(Nodes[k].def.y - Nodes[right].def.y,2.0) );
	  Nodes[k].characteristicL = 0.5*(def_dist_L + def_dist_R);
	  ds = Nodes[k].characteristicL;

    t_L.x = (Nodes[left].def.x - Nodes[k].def.x)/def_dist_L;
    t_L.y = (Nodes[left].def.y - Nodes[k].def.y)/def_dist_L;
    t_R.x = (Nodes[right].def.x - Nodes[k].def.x)/def_dist_R;
    t_R.y = (Nodes[right].def.y - Nodes[k].def.y)/def_dist_R;

    // compute elastic forces
    ef[k].x = ((gamma_mem + k_mem*((def_dist_L - ref_dist_L)/ref_dist_L))*t_L.x + (gamma_mem + k_mem*((def_dist_R-ref_dist_R)/ref_dist_R))*t_R.x)/ds;
    ef[k].y = ((gamma_mem + k_mem*((def_dist_L - ref_dist_L)/ref_dist_L))*t_L.y + (gamma_mem + k_mem*((def_dist_R-ref_dist_R)/ref_dist_R))*t_R.y)/ds;
	}
}
// END COMPUTE ELASTIC FORCES

// COMPUTE CONTACT FORCES
void computeContactForces(int Npts, vertex Nodes[], vector f_contact[]) {
  int k;
  double v_dist;

  for(k=0; k<Npts; k++) {
    if(Nodes[k].bottom == 1) {
      /*v_dist = fabs(Nodes[k].def.y-2.6) - channel_height/2.0;
        if(v_dist>0) {
          f_contact[k].y = -k_contact*fabs(v_dist)*(Nodes[k].def.y-2.5)/fabs(Nodes[k].def.y-2.5);
        }
      */  // clean up on June 24
      v_dist = fabs(Nodes[k].def.y)-substrate_height;
      if(v_dist<0) {
        f_contact[k].y = k_contact*fabs(v_dist)*Nodes[k].def.y/fabs(Nodes[k].def.y);
      }
    }
  }
}
// END COMPUTE CONTACT FORCES

// COMPUTE P0 TERM
double computeP0Term(int Npts, vertex Nodes[], vector f[]) {
	int k, left, right;
	double p0_top, p0_bottom, p0;		
	double def_dist_L, def_dist_R;
	vector t_R, t_L, n;

	p0_top = 0;
	p0_bottom = 0;	
	for(k=0; k<Npts; k++) {
    // setup neighbors and membrane length information
    if(k==0) { left = Npts-1; right = k+1; }
    else if(k==Npts-1) { left = k-1; right = 0; }
    else { left = k-1; right = k+1; }

    def_dist_L = sqrt( pow(Nodes[k].def.x - Nodes[left].def.x,2.0) + pow(Nodes[k].def.y - Nodes[left].def.y,2.0) );
    def_dist_R = sqrt( pow(Nodes[k].def.x - Nodes[right].def.x,2.0) + pow(Nodes[k].def.y - Nodes[right].def.y,2.0) );

    t_L.x = (Nodes[left].def.x - Nodes[k].def.x)/def_dist_L;
    t_L.y = (Nodes[left].def.y - Nodes[k].def.y)/def_dist_L;
    t_R.x = (Nodes[right].def.x - Nodes[k].def.x)/def_dist_R;
    t_R.y = (Nodes[right].def.y - Nodes[k].def.y)/def_dist_R;

    n.x = 0.5*(-t_R.y*def_dist_R + t_L.y*def_dist_L)/Nodes[k].characteristicL;
    n.y = 0.5*(t_R.x*def_dist_R - t_L.x*def_dist_L)/Nodes[k].characteristicL;

    p0_top += (f[k].x/Nodes[k].characteristicL)*n.x*Nodes[k].characteristicL + (f[k].y/Nodes[k].characteristicL)*n.y*Nodes[k].characteristicL;
    p0_bottom += Nodes[k].characteristicL;
	}
 
  p0 = p0_top/p0_bottom;
	return p0;
}
// END COMPUTE P0 TERM

// COMPUTE PRESSURE CONTRIBUTION
void computePressureForces(int Npts, vertex Nodes[], vector pf[], double p) {
	int k, left, right;
	double def_dist_L, def_dist_R;
	vector t_L, t_R, n;

  // add pressure force
  for(k=0; k<Npts; k++) {
  	if(k==0) { left = Npts-1; right = k+1; }
    else if(k==Npts-1) { left = k-1; right = 0; }
    else { left = k-1; right = k+1; }

    def_dist_L = sqrt( pow(Nodes[k].def.x - Nodes[left].def.x,2.0) + pow(Nodes[k].def.y - Nodes[left].def.y,2.0) );
    def_dist_R = sqrt( pow(Nodes[k].def.x - Nodes[right].def.x,2.0) + pow(Nodes[k].def.y - Nodes[right].def.y,2.0) );

   	t_L.x = (Nodes[left].def.x-Nodes[k].def.x)/def_dist_L;
    t_L.y = (Nodes[left].def.y-Nodes[k].def.y)/def_dist_L;
    t_R.x = (Nodes[right].def.x - Nodes[k].def.x)/def_dist_R;
    t_R.y = (Nodes[right].def.y - Nodes[k].def.y)/def_dist_R;

    // compute pressure forces using normal = dT/ds / ||dT/ds||
    n.x = 0.5*(-t_R.y*def_dist_R + t_L.y*def_dist_L)/Nodes[k].characteristicL;
    n.y = 0.5*(t_R.x*def_dist_R - t_L.x*def_dist_L)/Nodes[k].characteristicL;

    // careful: these are force densities
   	pf[k].x = -p*n.x;
    pf[k].y = -p*n.y;
	}
}
// END COMPUTE PRESSURE CONTRIBUTION

// COMPUTE ADHESION FORCES
void computeAdhesionForces(int Npts, vertex Nodes[], vector f_adh[]) {
	int k;
	double dist;

  for(k=0; k<Npts; k++) {
		if(Nodes[k].sticky == 1) {
			dist = sqrt( pow(Nodes[k].def.x - Nodes[k].adh_link.x,2) + pow(Nodes[k].def.y - Nodes[k].adh_link.y,2) );

      f_adh[k].x = k_adh*((dist-attach_dl0)/attach_dl0)*(Nodes[k].adh_link.x - Nodes[k].def.x)/dist;
      f_adh[k].y = k_adh*((dist-attach_dl0)/attach_dl0)*(Nodes[k].adh_link.y - Nodes[k].def.y)/dist;
		}
	}
}
// END COMPUTE ADHESION FORCES

//-------------------------------------------------------------------------------------------------//

// PROGRESSION IN TIME
void progress(int Npts, vertex Nodes[]) {
	int j, k;

	vector ef[Npts];	
	vector pf[Npts];
	vector f[Npts];
	vector f_adh[Npts];
	vector f_contact[Npts];
	vector v[Npts];

  // file handling
  ofstream f1;
	f1.precision(16);
	f1.open("membrane.txt");

	// zero the force
  for(k=0; k<Npts; k++) {
		Nodes[k].characteristicL = 0;
		Nodes[k].force.x = 0;
  	Nodes[k].force.y = 0;
		f[k].x = 0;
		f[k].y = 0;
  } 

  // compute initial area
	j = Npts-1;
  double A0 = 0;
	for(k=0; k<Npts; k++) {
    A0 += -0.5*(Nodes[j].ref.x + Nodes[k].ref.x)*(Nodes[j].ref.y-Nodes[k].ref.y);
    j = k;
  } 

	// MARCHING IN TIME
	double p, currentA, p0;

	for(int t=0; t<stop; t++) {
		// ZERO OUT STUFF
		for(k=0; k<Npts; k++) {
			ef[k].x = 0;
			ef[k].y = 0;
      f_adh[k].x = 0;
      f_adh[k].y = 0;
      f_contact[k].x = 0;
			f_contact[k].y = 0;
			pf[k].x = 0;
			pf[k].y = 0;
		
			v[k].x = 0;
			v[k].y = 0;

			Nodes[k].bottom = 0; 
			if(Nodes[k].def.y < 0.5) { Nodes[k].bottom = 1; }
		}

		if(t==0) {
			for(k=0; k<Npts; k++) {
				f1 << t*tstep << " " << Nodes[k].def.x << " " << Nodes[k].def.y << " " << Nodes[k].ref.x << " " << Nodes[k].ref.y << " " << ef[k].x << " " << ef[k].y << " " << f_adh[k].x << " " << f_adh[k].y << " " << f_contact[k].x << " " << f_contact[k].y << " " << pf[k].x << " " << pf[k].y << " " << f[k].x << " " << f[k].y << " " << Nodes[k].characteristicL << " " << Nodes[k].sticky << " " << v[k].x << " " << v[k].y <<  endl;
			}
		}

		// COMPUTE ELASTIC FORCES
		computeElasticForces(Npts, Nodes, ef);

		// COMPUTE CONTACT FORCES
		computeContactForces(Npts, Nodes, f_contact);

		// COMPUTE ADHESION FORCE DENSITIES
		computeAdhesionForces(Npts, Nodes, f_adh);

		// RUPTURE ADHESIONS ABOVE FORCE THRESHOLD
		for(k=0; k<Npts; k++) {
			if( sqrt(pow(f_adh[k].x,2)+pow(f_adh[k].y,2))>crit_adh ) {
				f_adh[k].x = 0;
				f_adh[k].y = 0;
				Nodes[k].sticky = 0;
			}

			// compute total forces (f labels forces not force densities) 
			f[k].x = (ef[k].x + f_adh[k].x + f_contact[k].x)*Nodes[k].characteristicL;
      f[k].y = (ef[k].y + f_adh[k].y + f_contact[k].y)*Nodes[k].characteristicL;	
		}
	
		// COMPUTE PRESSURE CONTRIBUTION
		p0 = computeP0Term(Npts, Nodes, f);
		currentA = computeArea(Npts, Nodes);
  	p = p0 + k_cell*log(A0/currentA);
		computePressureForces(Npts, Nodes, pf, p);
		for(k=0; k<Npts; k++) {
      // careful: pf are force densities (need to convert them to forces)
      f[k].x += pf[k].x*Nodes[k].characteristicL;
      f[k].y += pf[k].y*Nodes[k].characteristicL;
		}

		// MOVE MEMBRANE POINTS
		for(k=0; k<Npts; k++) {
			Nodes[k].def.x = Nodes[k].def.x + tstep*f[k].x/(xi*Nodes[k].characteristicL);
			Nodes[k].def.y = Nodes[k].def.y + tstep*f[k].y/(xi*Nodes[k].characteristicL);
   		v[k].x = f[k].x/(xi*Nodes[k].characteristicL);
			v[k].y = f[k].y/(xi*Nodes[k].characteristicL);
		}

			if(t%10000==0 && t!=0) {
			printf("%.3f p0=%.5f, p=%.5f, A0=%.5f, A=%.5f\n",t*tstep,p0,p,A0,currentA);
			for(k=0; k<Npts; k++) {
				f1 << t*tstep << " " << Nodes[k].def.x << " " << Nodes[k].def.y << " " << Nodes[k].ref.x << " " << Nodes[k].ref.y << " " << ef[k].x << " " << ef[k].y << " " << f_adh[k].x << " " << f_adh[k].y << " " << f_contact[k].x << " " << f_contact[k].y << " " << pf[k].x << " " << pf[k].y << " " << f[k].x << " " << f[k].y << " " << Nodes[k].characteristicL << " " <<  Nodes[k].sticky << " " << v[k].x << " " << v[k].y << endl;  
			}
		}

	}	

	f1.close();
}

// DRAW START CONFIGURATION
void startCurve() {
  
  // file handling
  ifstream f1;
	f1.open("SideCellN162.txt");

  // determine length
  int Npoints = -1; 

  string b1;
  while( !f1.eof() ) { 
		getline(f1, b1); 
		Npoints++;
  }
  f1.close();

	f1.open("SideCellN162.txt");

  vertex Nodes[Npoints];
  int counter = 0;
	int M = 0;  // number of points on surface
  double c1, c2;
	while(f1 >> c1 >> c2) {
		Nodes[counter].init.x = c1; 
    Nodes[counter].init.y = c2; 
    Nodes[counter].ref.x = c1; 
    Nodes[counter].ref.y = c2; 
    Nodes[counter].def.x = c1; 
    Nodes[counter].def.y = c2; 

		// adhesion and surface information
		Nodes[counter].adh_link.x = std::numeric_limits<double>::quiet_NaN();
		Nodes[counter].adh_link.y = std::numeric_limits<double>::quiet_NaN();
		Nodes[counter].bottom = 0;
		Nodes[counter].sticky = 0;

		// check bottom and introduce bottom wall and bottom channel 
		if(Nodes[counter].init.y<3.0) {
			Nodes[counter].bottom = 1;
			Nodes[counter].adh_link.x = Nodes[counter].ref.x;
     	Nodes[counter].adh_link.y = -0.5;
			Nodes[counter].sticky= 1; 
			M++;
		}
		counter++;
  }
  f1.close();

	// progress points forward in time
  printf("%d points and %d points on the bottom / contact with surface.\n",Npoints,M);
	progress(Npoints, Nodes);
}

int main() {
	printf("Pill-shaped membrane simulation:\n");
  clock_t begin = clock();
  startCurve();
  clock_t end = clock();
  printf("Total computation time (s): %.10f\n", double(diffclock(end,begin)));

  return 0;
}
