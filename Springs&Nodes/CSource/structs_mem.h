#ifndef STRUCTS_H
#define STRUCTS_H

// VECTOR STRUCT
typedef struct vector{
  double x; // x-component
  double y; // y-component
} vector;

// MEMBRANE VERTEX STRUCT
typedef struct vertex{
  vector init; 						// initial reference coordinates
  vector ref; 						// reference coordinates
  vector def; 						// deformed/current coordinates 
  vector force; 					// force vector
	double characteristicL; // characteristic length unit
	int bottom; 						// 1 if this is a `bottom' point and 0 otherwise
	vector adh_link;				// ligand location on substrate
	int sticky;							// 1 for an `active' adhesion bond and 0 otherwise
} vertex;
#endif
