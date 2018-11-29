// This library defines the data structures: Vertex, Triangle and CPU_Time

#define LINESIZE 100
#define MAX_SIZE 110000
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

