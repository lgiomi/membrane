#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define TOL 1e-10
#define RNMX (1.0-EPS)
#define LINESIZE 100
#define MAX_SIZE 100000
#define MAX_NEIGHBORS 500
#define PI 3.141592653589793238462643383279502884

#ifndef TYPES_DECLARATION

typedef struct {

	int label;

	long num_of_neighbors;
	long neighbor[MAX_NEIGHBORS];
	
  	double x;
	double y;
	double z;

	double gridx;
	double gridy;
	
	double hx;
	double hy;
	double hz;

	double h2;
	double kg;

	double h2_avg;
	double kg_avg;

	double phi;

	double nx;
	double ny;
	double nz;
	
	double area;
	double weight[MAX_NEIGHBORS];
	
} Vertex;

typedef struct {
	
	long v1;
	long v2;
	long v3;
	
} Triangle;

typedef struct {

  	long d;
  	long h;
  	long m;
  	long s;

} CPU_Time;

#endif
#define TYPES_DECLARATION

