/******************************************************************
**                                                               **
**                           MEMBRANE                            **
**                                                               **
** Program: membrane.c                                           **
** Version: 2.1                                                  **
** Last Change: Nov 28 2016                                      **
**                                                               **
*******************************************************************/

#include "membrane.h"

#define  sign(x) ((x < 0) ? -1 : 1)
#define  min(x,y) ((x < y) ? x : y)
#define  max(x,y) ((x > y) ? x : y)
#define  mod(x,y) ((x % y + y) % y)

time_t t1, t2;

Triangle triangle[MAX_SIZE];
Vertex vertex[MAX_SIZE];

long num_of_meshpoint;
long num_of_triangles;
long num_of_iteration;
long method;

int num_of_domains;

double total_area;
double sigma;
double DT;

void rk2();
void rk4();
void run();
void end();
void init();
void euler();
void one_step();
void init_random();
void get_geometry();
void track_domains();
void find_neighbors();
void get_rhs(double *);
void progress_bar(long);
void import_mesh(char *);
void write_hi(FILE *, long);
void export_conf(long, long);
void get_time(time_t, CPU_Time *);
void export_graphic_complex(FILE *);

/*******************************************************************/

main() 
{	
	init();
	time(&t1);
	run();
	end();
	
  	return EXIT_SUCCESS;	
}

/*******************************************************************/

void init()
{
	double box_size, run_time;
	char f_name[32];
	long cflag;

	printf("Input mesh file ");
	scanf("%s",f_name);

	import_mesh(f_name);
	get_geometry();
	
	printf("Input run time ");  
  	scanf("%lg",&run_time);

  	printf("Input iterations ");  
  	scanf("%ld",&num_of_iteration);

	DT = run_time/num_of_iteration;

	printf("Choose integration method "); 
  	printf("\n");
  	printf("\n\t1. Euler");
  	printf("\n\t2. Runge-Kutta (2 steps)");
  	printf("\n\t3. Runge-Kutta (4 steps)");
  	printf("\n\n");
  	printf("Input flag ");
  	scanf("%ld",&method);
  
	printf("Input sigma "); 
  	scanf("%lg",&sigma);

	init_random();
}

/*******************************************************************/

void import_mesh(char *f_name)
{
	long i, v1, v2, v3, chi, dummy, num_of_edges=0;

	
	FILE *f_in = fopen(f_name,"r");

	/*
	fscanf(f_in,"%ld",&num_of_meshpoint);
	fscanf(f_in,"%ld",&num_of_triangles);
	
	vertex = malloc(num_of_meshpoint*sizeof(Vertex));
	triangle = malloc(num_of_triangles*sizeof(Triangle));
	
	for (i=0; i<num_of_meshpoint; i++){
		vertex[i].num_of_neighbors = 0;
		fscanf(f_in,"%lg%lg%lg",
		&vertex[i].x,
		&vertex[i].y,
		&vertex[i].z);
	}
	
	for (i=0; i<num_of_triangles; i++){
		fscanf(f_in,"%ld%ld%ld",&v1,&v2,&v3);
		triangle[i].v1 = v1-1;
		triangle[i].v2 = v2-1;
		triangle[i].v3 = v3-1;
	}
	*/
	
	char line[LINESIZE];
	
	fgets(line,LINESIZE-1,f_in); 
	fgets(line,LINESIZE-1,f_in); 
	fgets(line,LINESIZE-1,f_in); 
	fgets(line,LINESIZE-1,f_in); 
	fgets(line,LINESIZE-1,f_in); 
	fgets(line,LINESIZE-1,f_in); 
	fgets(line,LINESIZE-1,f_in); 
	fgets(line,LINESIZE-1,f_in);
	
	fscanf(f_in,"%ld",&num_of_meshpoint);

	if (num_of_meshpoint>MAX_SIZE){
		printf("Error: the number of vertices (%ld) exceeds MAX_SIZE (%d)\n",num_of_meshpoint,MAX_SIZE);
		exit(0);
	}
		
	for (i=0; i<num_of_meshpoint; i++){
		fscanf(f_in,"%ld%lg%lg%lg",
		&dummy,
		&vertex[i].x,
		&vertex[i].y,
		&vertex[i].z);
	}
	
	fgets(line,LINESIZE-1,f_in); 
	fgets(line,LINESIZE-1,f_in); 
	fgets(line,LINESIZE-1,f_in); 
	
	fscanf(f_in,"%ld",&num_of_triangles);
		
	if (num_of_triangles>MAX_SIZE){
		printf("Error: the number of triangles (%ld) exceeds MAX_SIZE (%d)\n",num_of_triangles,MAX_SIZE);
		exit(0);
	}	
		
	for (i=0; i<num_of_triangles; i++){
		fscanf(f_in,"%ld%ld%ld%ld%ld%ld%ld%ld",
		&dummy,
		&dummy,
		&dummy,
		&dummy,
		&dummy,
		&v1,
		&v2,
		&v3);
		triangle[i].v1 = v1-1;
		triangle[i].v2 = v2-1;
		triangle[i].v3 = v3-1;
	}
				
	fclose(f_in);		
				
	find_neighbors();			
				
	for (i=0; i<num_of_meshpoint; i++){
		num_of_edges += vertex[i].num_of_neighbors;
	}
	
	if (!num_of_edges%2){
		printf("Error: bad triangulation, 2E = %ld\n",num_of_edges);
		exit(0);
	}

	chi = num_of_meshpoint-num_of_edges/2+num_of_triangles;
	
	if (chi!=2){
		printf("Error: bad triangulation, chi = %ld\n",chi);
		exit(0);
	} 
}	

/*******************************************************************/

void find_neighbors()
{
	int is_neighbor(long,long);
	long i, max_neighbors=0; 
	
	for (i=0; i<num_of_meshpoint; i++){
		vertex[i].num_of_neighbors = 0;
	}
		
	for (i=0; i<num_of_triangles; i++){
		if (!is_neighbor(triangle[i].v1,triangle[i].v2)){
			vertex[triangle[i].v1].neighbor[vertex[triangle[i].v1].num_of_neighbors] = triangle[i].v2;
			vertex[triangle[i].v1].num_of_neighbors++;
			if (vertex[triangle[i].v1].num_of_neighbors==MAX_NEIGHBORS){
				max_neighbors = triangle[i].v1;
				break;
			}
		}
		if (!is_neighbor(triangle[i].v1,triangle[i].v3)){
			vertex[triangle[i].v1].neighbor[vertex[triangle[i].v1].num_of_neighbors] = triangle[i].v3;
			vertex[triangle[i].v1].num_of_neighbors++;
			if (vertex[triangle[i].v1].num_of_neighbors==MAX_NEIGHBORS){
				max_neighbors = triangle[i].v1;
				break;
			}
		}
		if (!is_neighbor(triangle[i].v2,triangle[i].v1)){
			vertex[triangle[i].v2].neighbor[vertex[triangle[i].v2].num_of_neighbors] = triangle[i].v1;
			vertex[triangle[i].v2].num_of_neighbors++;
			if (vertex[triangle[i].v2].num_of_neighbors==MAX_NEIGHBORS){
				max_neighbors = triangle[i].v2;
				break;
			}
		}
		if (!is_neighbor(triangle[i].v2,triangle[i].v3)){
			vertex[triangle[i].v2].neighbor[vertex[triangle[i].v2].num_of_neighbors] = triangle[i].v3;
			vertex[triangle[i].v2].num_of_neighbors++;
			if (vertex[triangle[i].v2].num_of_neighbors==MAX_NEIGHBORS){
				max_neighbors = triangle[i].v2;
				break;
			}
		}
		if (!is_neighbor(triangle[i].v3,triangle[i].v1)){
			vertex[triangle[i].v3].neighbor[vertex[triangle[i].v3].num_of_neighbors] = triangle[i].v1;
			vertex[triangle[i].v3].num_of_neighbors++;
			if (vertex[triangle[i].v3].num_of_neighbors==MAX_NEIGHBORS){
				max_neighbors = triangle[i].v3;
				break;
			}
		}
		if (!is_neighbor(triangle[i].v3,triangle[i].v2)){
			vertex[triangle[i].v3].neighbor[vertex[triangle[i].v3].num_of_neighbors] = triangle[i].v2;
			vertex[triangle[i].v3].num_of_neighbors++;
			if (vertex[triangle[i].v3].num_of_neighbors==MAX_NEIGHBORS){
				max_neighbors = triangle[i].v3;
				break;
			}
		}
	}
	
	if (max_neighbors){
		printf("Error: vertex %ld reached the maximum number of neighbors\n",max_neighbors);
		exit(0);
	}
}

/*******************************************************************/

int is_neighbor(long i, long j)
{
	long k;
	
	for (k=0; k<vertex[i].num_of_neighbors; k++){
		if (j==vertex[i].neighbor[k]) return 1;
	}
	
	return 0;
}

/*******************************************************************/

void get_geometry()
{
	double x[3], y[3], dx, dy, dz, norm, cota, cotb, base, height;
	double theta_this, theta_next, area_prev, area_next;
	
	long i, j, tmp, next, prev, this,  obtuse[MAX_SIZE], num_of_obtuse=0;
	long which_triangle(long,long,long);
	
	int swap(long *,long *), swapped;
	
	total_area = 0;	
			
	for (i=0; i<num_of_meshpoint; i++){
	
		// Construct an approximation for the tangent plane at each vertex
		
		this = vertex[i].neighbor[0];
		next = vertex[i].neighbor[1];
		
		x[0] = vertex[this].x-vertex[i].x;
		x[1] = vertex[this].y-vertex[i].y;
		x[2] = vertex[this].z-vertex[i].z;
		
		norm = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
		
		x[0] /= norm;
		x[1] /= norm;
		x[2] /= norm;
					
		y[0] = vertex[next].x-vertex[i].x;
		y[1] = vertex[next].y-vertex[i].y;
		y[2] = vertex[next].z-vertex[i].z;
			
		norm = x[0]*y[0]+x[1]*y[1]+x[2]*y[2];
			
		y[0] -= norm*x[0];
		y[1] -= norm*x[1];
		y[2] -= norm*x[2];
		
		norm = sqrt(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]);
		
		y[0] /= norm;
		y[1] /= norm;
		y[2] /= norm;
	
		// Sort neighbors count-clockwise on the tangent plane	
				
		do {
			
			swapped = 0;
			
			for (j=0; j<vertex[i].num_of_neighbors-1; j++){
				
				dx = (vertex[i].x-vertex[vertex[i].neighbor[j]].x)*x[0] 
					+(vertex[i].y-vertex[vertex[i].neighbor[j]].y)*x[1]
					+(vertex[i].z-vertex[vertex[i].neighbor[j]].z)*x[2];
					
				dy = (vertex[i].x-vertex[vertex[i].neighbor[j]].x)*y[0] 
					+(vertex[i].y-vertex[vertex[i].neighbor[j]].y)*y[1]
					+(vertex[i].z-vertex[vertex[i].neighbor[j]].z)*y[2];						
				
				theta_this = atan2(dy,dx);
				
				dx = (vertex[i].x-vertex[vertex[i].neighbor[j+1]].x)*x[0] 
					+(vertex[i].y-vertex[vertex[i].neighbor[j+1]].y)*x[1]
					+(vertex[i].z-vertex[vertex[i].neighbor[j+1]].z)*x[2];

				dy = (vertex[i].x-vertex[vertex[i].neighbor[j+1]].x)*y[0] 
					+(vertex[i].y-vertex[vertex[i].neighbor[j+1]].y)*y[1]
					+(vertex[i].z-vertex[vertex[i].neighbor[j+1]].z)*y[2];

				theta_next = atan2(dy,dx);
				
				if (theta_this > theta_next){
					swapped = swap(&vertex[i].neighbor[j],&vertex[i].neighbor[j+1]);
				}
			}
		} while(swapped);	
		
	// Calculate the Laplacian weights and the mean curvature
		
		vertex[i].hx = 0;
		vertex[i].hy = 0;
		vertex[i].hz = 0;
		vertex[i].area = 0;
	
		for (j=0; j<vertex[i].num_of_neighbors; j++){
			
			this = vertex[i].neighbor[j];
			prev = vertex[i].neighbor[mod(j-1,vertex[i].num_of_neighbors)];
			next = vertex[i].neighbor[mod(j+1,vertex[i].num_of_neighbors)];
		
			// Calculate cot(a)
			
			x[0] = vertex[this].x-vertex[prev].x;
			x[1] = vertex[this].y-vertex[prev].y;
			x[2] = vertex[this].z-vertex[prev].z;
			
			base = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
				
			x[0] /= base;
			x[1] /= base;
			x[2] /= base;
						
			y[0] = vertex[i].x-vertex[prev].x;
			y[1] = vertex[i].y-vertex[prev].y;
			y[2] = vertex[i].z-vertex[prev].z;
				
			norm = x[0]*y[0]+x[1]*y[1]+x[2]*y[2];
				
			y[0] -= norm*x[0];
			y[1] -= norm*x[1];
			y[2] -= norm*x[2];
			
			height = sqrt(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]);
			
			y[0] /= height;
			y[1] /= height;
			y[2] /= height;
					
			dx = (vertex[i].x-vertex[prev].x)*x[0]
				+(vertex[i].y-vertex[prev].y)*x[1]
				+(vertex[i].z-vertex[prev].z)*x[2];
				
			dy = (vertex[i].x-vertex[prev].x)*y[0]
				+(vertex[i].y-vertex[prev].y)*y[1]
				+(vertex[i].z-vertex[prev].z)*y[2];				
		    	
			cota = dx/dy;
		
			area_prev = base*height/2;
		
			if (cota<0){
				obtuse[num_of_obtuse] = which_triangle(i,this,prev); 
				num_of_obtuse++;
			}
			
			// Calculate cot(b)
			
			x[0] = vertex[this].x-vertex[next].x;
			x[1] = vertex[this].y-vertex[next].y;
			x[2] = vertex[this].z-vertex[next].z;
			
			base = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
				
			x[0] /= base;
			x[1] /= base;
			x[2] /= base;
			
			y[0] = vertex[i].x-vertex[next].x;
			y[1] = vertex[i].y-vertex[next].y;
			y[2] = vertex[i].z-vertex[next].z;
				
			norm = x[0]*y[0]+x[1]*y[1]+x[2]*y[2];
				
			y[0] -= norm*x[0];
			y[1] -= norm*x[1];
			y[2] -= norm*x[2];
			
			height = sqrt(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]);
			
			y[0] /= height;
			y[1] /= height;
			y[2] /= height;
			
			dx = (vertex[i].x-vertex[next].x)*x[0]
				+(vertex[i].y-vertex[next].y)*x[1]
				+(vertex[i].z-vertex[next].z)*x[2];
				
			dy = (vertex[i].x-vertex[next].x)*y[0]
				+(vertex[i].y-vertex[next].y)*y[1]
				+(vertex[i].z-vertex[next].z)*y[2];
			
			cotb = dx/dy;
			
			area_next = base*height/2;
			
			if (cotb<0){
				obtuse[num_of_obtuse] = which_triangle(i,this,next); 
				num_of_obtuse++;
			}
			
			// Calculate the un-normalized Laplacian weights
			
			vertex[i].weight[j] = (cota+cotb)/2;
			
			// Calculate the un-normalized mean curvature vector
			
			dx = vertex[i].x-vertex[this].x;
			dy = vertex[i].y-vertex[this].y;
			dz = vertex[i].z-vertex[this].z;
				
			vertex[i].hx += dx*(cota+cotb)/4;
			vertex[i].hy += dy*(cota+cotb)/4;
			vertex[i].hz += dz*(cota+cotb)/4;
			
			// Calculate the area of a vertex
			
			//vertex[i].area += (dx*dx+dy*dy+dz*dz)*(cota+cotb)/8;
			
			vertex[i].area += ( cota>0 ? (dx*dx+dy*dy+dz*dz)*cota/8 : area_prev/4);
			vertex[i].area += ( cotb>0 ? (dx*dx+dy*dy+dz*dz)*cotb/8 : area_next/4);
		}
		
		// Normalized the Laplacian weights 
		
		for (j=0; j<vertex[i].num_of_neighbors; j++){
			
			vertex[i].weight[j] /= vertex[i].area;	
		}
		
		// Normalize the mean curvature vector
				
		vertex[i].hx /= vertex[i].area;
		vertex[i].hy /= vertex[i].area;
		vertex[i].hz /= vertex[i].area;
		
		total_area += vertex[i].area;
	}
	
	FILE *f_ou;
	
	f_ou = fopen("ge.dat","w");
	
	for (i=0; i<num_of_meshpoint; i++){
		
		fprintf(f_ou,"%ld\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t",
			i,
			vertex[i].x,
			vertex[i].y,
			vertex[i].z,
			vertex[i].area,
			vertex[i].hx*vertex[i].hx+vertex[i].hy*vertex[i].hy+vertex[i].hz*vertex[i].hz);
			
		for (j=0; j<vertex[i].num_of_neighbors; j++){
			fprintf(f_ou,"%.10f\t",vertex[i].weight[j]);
		}	
		
		fprintf(f_ou,"\n");
	}

	fclose(f_ou);

	/*
	f_ou = fopen("debug.dat","w");

	for (i=229; i<=239; i++){
	
		double h = sqrt(vertex[i].hx*vertex[i].hx+vertex[i].hy*vertex[i].hy+vertex[i].hz*vertex[i].hz);
	
		//if (fabs(h-1)<0.1) continue;
			
		fprintf(f_ou,"vertex[%ld] = {%lg,%lg,%lg}\n\n",i,vertex[i].x,vertex[i].y,vertex[i].z);
		fprintf(f_ou,"\th %lg\n",h);
		fprintf(f_ou,"\ta %lg\n",vertex[i].area);
		for (j=0; j<vertex[i].num_of_neighbors; j++){
			fprintf(f_ou,"\tneighbor: %ld {%lg,%lg,%lg} w %lg\n",
			vertex[i].neighbor[j],
			vertex[vertex[i].neighbor[j]].x,
			vertex[vertex[i].neighbor[j]].y,
			vertex[vertex[i].neighbor[j]].z,
			vertex[i].weight[j]);	
		}
		fprintf(f_ou,"\n");
	}

	fclose(f_ou);
	
	f_ou = fopen("ob.dat","w");
	
	for (i=0; i<num_of_obtuse; i++){
		fprintf(f_ou,"%ld {%lg,%lg,%lg} {%lg,%lg,%lg} {%lg,%lg,%lg}\n",
		obtuse[i],
		vertex[triangle[obtuse[i]].v1].x,
		vertex[triangle[obtuse[i]].v1].y,
		vertex[triangle[obtuse[i]].v1].z,
		vertex[triangle[obtuse[i]].v2].x,
		vertex[triangle[obtuse[i]].v2].y,
		vertex[triangle[obtuse[i]].v2].z,
		vertex[triangle[obtuse[i]].v3].x,
		vertex[triangle[obtuse[i]].v3].y,
		vertex[triangle[obtuse[i]].v3].z);
	}
	
	fclose(f_ou);
	*/

	printf("\n");
	printf("\tVertices %ld\n",num_of_meshpoint);
	printf("\tTriangles %ld\n",num_of_triangles);
	printf("\tObtuse triangles %ld\n",num_of_obtuse);
	printf("\tTotal surface area %lg\n",total_area);
	printf("\n");

}

/*******************************************************************/

int swap(long *x, long *y)
{
	long z = (*x);
	
	(*x) = (*y); (*y) = z;
	
	return 1;
}

/*******************************************************************/

long which_triangle(long a, long b, long c)
{
	long i;
	
	for (i=0; i<num_of_triangles; i++){
		if (!(triangle[i].v1==a)||!(triangle[i].v1==b)||!(triangle[i].v1==c)) continue;
		if (!(triangle[i].v2==a)||!(triangle[i].v2==b)||!(triangle[i].v2==c)) continue;
		if (!(triangle[i].v3==a)||!(triangle[i].v3==b)||!(triangle[i].v3==c)) continue;
		return i;
	}
}

/*******************************************************************/

void init_random()
{
	double area, noise=0.1, ran2(long *);
	long i, seed;
		
	printf("Input area percentage ");
	scanf("%lg",&area);
	
	printf("Input random seed ");
    scanf("%ld",&seed);

	for (i=0; i<num_of_meshpoint; i++){
		vertex[i].phi = (2*area-1)+noise*(1-2*ran2(&seed));
	}
}

/*******************************************************************/

void get_rhs(double *rhs)
{
	double laplace(long), epsilon=0.1, lagrange=0;
	long i; 	
	
	// Calculate the Lagrange multiplier
	
	for (i=0; i<num_of_meshpoint; i++){
		lagrange += vertex[i].area*(sigma*laplace(i)-vertex[i].phi*(vertex[i].phi*vertex[i].phi-1)/(epsilon*epsilon));		
	}		
	
	lagrange /= total_area;
	
	// Calculate the right-hand side of the Allen-Cahn equation
	
	for (i=0; i<num_of_meshpoint; i++){
		
		rhs[i] = sigma*laplace(i)-vertex[i].phi*(vertex[i].phi*vertex[i].phi-1)/(epsilon*epsilon)-lagrange;
	}
}

/*******************************************************************/

double laplace(long i)
{
	double lb=0;
	long j;
	
	for (j=0; j<vertex[i].num_of_neighbors; j++){
		lb += vertex[i].weight[j]*(vertex[vertex[i].neighbor[j]].phi-vertex[i].phi); 
	}
	
	return lb;
}

/*******************************************************************/

void run()
{
	long t, export=50;
	char f_na[32];
	
	FILE *f_hi, *f_ou;
	
	f_hi = fopen("hi.dat","w");	
	
	for (t=0; t<num_of_iteration; t++){
		write_hi(f_hi,t);
		/*
		if (t%export==0){
			sprintf(f_na,"gc_%ld.m",t);
			f_ou = fopen(f_na,"w");
			export_graphic_complex(f_ou);
			fclose(f_ou);
		}
		*/
		//track_domains();
		progress_bar(t);
		one_step();
	}

	fclose(f_hi);
}

/*******************************************************************/

void one_step()
{
	switch(method){
  	case 1:
		euler();
    	break;
  	case 2:
		rk2();
    	break;
  	case 3:
		rk4();
    	break;
	}
}

/*******************************************************************/

void euler()
{
	double h=DT;
	double rhs[MAX_SIZE];
	long i;
		
	get_rhs(rhs);
	
	for (i=0; i<num_of_meshpoint; i++){
		vertex[i].phi += h*rhs[i];
	}
}

/*******************************************************************/

void rk2()
{
	double h=DT/2;
	double rhs0[MAX_SIZE]; 
	double rhs1[MAX_SIZE];
	double rhs2[MAX_SIZE];
	long i;
	
	for (i=0; i<num_of_meshpoint; i++){
		rhs0[i] = vertex[i].phi; 
	}
	
	get_rhs(rhs1);
	
	for (i=0; i<num_of_meshpoint; i++){
		vertex[i].phi = rhs0[i]+DT*rhs1[i];
	}
	
	get_rhs(rhs2);
	
	for (i=0; i<num_of_meshpoint; i++){
		vertex[i].phi = rhs0[i]+h*(rhs1[i]+rhs2[i]);
	}
}

/*******************************************************************/

void rk4()
{
	double h=DT/6; 
	double rhs0[MAX_SIZE];
	double rhs1[MAX_SIZE];
	double rhs2[MAX_SIZE];
	double rhs3[MAX_SIZE];
	double rhs4[MAX_SIZE];
	long i;
	
	for (i=0; i<num_of_meshpoint; i++){
		rhs0[i] = vertex[i].phi; 
	}
	
	get_rhs(rhs1);
	
	for (i=0; i<num_of_meshpoint; i++){
		vertex[i].phi = rhs0[i]+0.5*DT*rhs1[i];
	}
	
	get_rhs(rhs2);
	
	for (i=0; i<num_of_meshpoint; i++){
		vertex[i].phi = rhs0[i]+0.5*DT*rhs2[i];
	}
	
	get_rhs(rhs3);
	
	for (i=0; i<num_of_meshpoint; i++){
		vertex[i].phi = rhs0[i]+DT*rhs3[i];
	}
	
	get_rhs(rhs4);
	
	for (i=0; i<num_of_meshpoint; i++){
		vertex[i].phi = rhs0[i]+h*(rhs1[i]+2*rhs2[i]+2*rhs3[i]+rhs4[i]);
	}
}

/*******************************************************************/

void track_domains()
{
	int l, *labels, num_of_labels=1, position(int,int,int*);
	long i;
	
	void label_components(Vertex *,int *);	
	
	Vertex *pos, *neg;

	pos = malloc(num_of_meshpoint*sizeof(Vertex));
	neg = malloc(num_of_meshpoint*sizeof(Vertex));

	memcpy(pos,vertex,num_of_meshpoint*sizeof(Vertex));
	memcpy(neg,vertex,num_of_meshpoint*sizeof(Vertex));
	
	for (i=0; i<num_of_meshpoint; i++){
		if (vertex[i].phi>0){
			pos[i].label = 1;
			neg[i].label = 0;
		}
		if (vertex[i].phi<0){
			pos[i].label = 0;
			neg[i].label = 1;
		}
	}

	label_components(pos,&num_of_labels);	
	label_components(neg,&num_of_labels);
	
	for (i=0; i<num_of_meshpoint; i++){
		if (vertex[i].phi>0) vertex[i].label = pos[i].label;
		if (vertex[i].phi<0) vertex[i].label = neg[i].label;
	}
	
	labels = calloc(num_of_labels,sizeof(int));
	
	num_of_domains = 0;
	
	for (i=0; i<num_of_meshpoint; i++){
		
		l = position(vertex[i].label,num_of_domains,labels);
				
		if (l==0){
			num_of_domains++;
			labels[num_of_domains] = vertex[i].label;
		}
	}
	
	free(labels);
	free(pos);
	free(neg);
} 

/*******************************************************************/

int position(int x, int imax, int *list)
{
	int i;
	
	for (i=1; i<=imax; i++){
		if (x==list[i]) return i;
	}
	
	return 0;
}

/*******************************************************************/

void label_components(Vertex *v, int *num_of_labels)
{
	int i, j, l, lmin, lmax, *parent, new_label;
	int find(int,int *), unio(int,int,int *);
	
	parent = calloc(num_of_meshpoint,sizeof(int));
		
	// Assign each vertex the same label of its first labeled neighbor or create a new one	
		
	for (i=0; i<num_of_meshpoint; i++){
			
		if (v[i].label==0) continue;
				
		new_label=1;
		
		for (j=0; j<v[i].num_of_neighbors; j++){
			
			l = v[v[i].neighbor[j]].label;
						
			if( l==0 || l==1 ) continue;
			
			v[i].label = l;
			new_label = 0;
			break;
		}
		
		if (new_label){
	
			(*num_of_labels)++; 
			v[i].label = (*num_of_labels);
		}
	}

	// Check if two neighbors have different labels
	
	for (i=0; i<num_of_meshpoint; i++){
	
		if (v[i].label==0) continue;
	
		for (j=0; j<v[i].num_of_neighbors; j++){
			
			l = v[v[i].neighbor[j]].label;
			
			if (l==0) continue;
			
			if (l!=v[i].label){
				lmin = min(l,v[i].label);
				lmax = max(l,v[i].label);
				unio(lmin,lmax,parent); 
			}
		}
	}

	// Merge all labels

	for (i=0; i<num_of_meshpoint; i++){
		
		if (v[i].label==0) continue;
		
		v[i].label = find(v[i].label,parent);
	}
		
	free(parent);
}

/*******************************************************************/

int find(int x, int *parent)
{
	int j = x;
	
	while(parent[j]) j = parent[j];
	
	return j;
}

/*******************************************************************/

int unio(int x, int y, int *parent)
{
	int j=x, k=y;
	
	while (parent[j]) j = parent[j];
	while (parent[k]) k = parent[k];
	
	if (j!=k) parent[k] = j;
}

/*******************************************************************/

void progress_bar(long t)
{
	long i, ticks, percent, num_of_ticks=20;
	char progress[32];

	percent = 100*t/num_of_iteration;
	ticks = num_of_ticks*percent/100;
	
	for (i=0; i<20; i++){
		progress[i] = ((i<=ticks)?'#':' ');
	}
	progress[20] = '\0';
	
	printf("Progress: [%s] %ld%%\r",progress,percent);
	fflush(stdout);
}

/*******************************************************************/

void end()
{	
	CPU_Time cpu_time;
	
	FILE *f_ou;
	
	time(&t2);
  	get_time(t2-t1,&cpu_time);
	
	track_domains();
	
	export_conf(-1,-1);			
			
	printf("\n\n");
	printf("\tNumber of vertices %ld\n",num_of_meshpoint);
	printf("\tNumber of domains %d\n",num_of_domains);
	printf("\tTotal area %lg\n",total_area);
	printf("\tTime step %g\n",DT);
	printf("\tSigma %g\n",sigma);
	printf("\tCPU Time %ld:%ld:%ld:%ld\n",
	 cpu_time.d,cpu_time.h,
	 cpu_time.m,cpu_time.s);
  	printf("\n");
	
	f_ou = fopen("last.m","w");
	export_graphic_complex(f_ou);
	fclose(f_ou);
	
	f_ou = fopen("ou.dat","w");	
	fprintf(f_ou,"Number of vertices %ld\n",num_of_meshpoint);
	fprintf(f_ou,"Number of domains %d\n",num_of_domains);
	fprintf(f_ou,"Total area %lg\n",total_area);
	fprintf(f_ou,"Time step %g\n",DT);
	fprintf(f_ou,"Sigma %g\n",sigma);
	fprintf(f_ou,"CPU Time %ld:%ld:%ld:%ld\n",
	 cpu_time.d,cpu_time.h,
	 cpu_time.m,cpu_time.s);
	fclose(f_ou);

}

/*******************************************************************/

void write_hi(FILE *f_ou, long t)
{
	double phi_tot=0;
	long m, i, j;
 	
	for (i=0; i<num_of_meshpoint; i++){
		phi_tot += vertex[i].area*vertex[i].phi;
	}
	
	fprintf(f_ou,"%.10f\t%.10f\t%.10f\t%d\n",
	t*DT,
	vertex[0].phi,
	phi_tot,
	num_of_domains);
	
	// Legend
	
	// 1) Time
	// 2) Phi
	// 3) Integral of Phi
	// 4) Number of domains
}

/*******************************************************************/

void export_conf(long t, long period)
{
	char f_na[32];
	long i, j;

	FILE *f_ou;

	if (t%period!=0) return;
	sprintf(f_na,"s-t%ld.dat",t);
	if (t<0) sprintf(f_na,"last.dat");
	
	f_ou = fopen(f_na,"w");
	
	for (i=0; i<num_of_meshpoint; i++){
		fprintf(f_ou,"%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%d\n",
		vertex[i].x,
		vertex[i].y,
		vertex[i].z,
		vertex[i].hx,
		vertex[i].hy,
		vertex[i].hz,
		vertex[i].phi,
		vertex[i].label);
	}
	fclose(f_ou);
}

/*******************************************************************/

void export_graphic_complex(FILE *f_ou)
{
	double r, g, b;
	long i, j;
	
	fprintf(f_ou,"GraphicsComplex[{");

	for (i=0; i<num_of_meshpoint; i++){
		fprintf(f_ou,"{%.10f,%.10f,%.10f},",
		vertex[i].x,
		vertex[i].y,
		vertex[i].z);
	}

	fseek(f_ou,-1, SEEK_CUR);

	fprintf(f_ou,"},Polygon[{");
	
	for (i=0; i<num_of_triangles; i++){		
		
		fprintf(f_ou,"{%ld,%ld,%ld},",
		triangle[i].v1+1,
		triangle[i].v2+1,
		triangle[i].v3+1);
	}
	
	fseek(f_ou,-1, SEEK_CUR);
	
	fprintf(f_ou,"}],VertexColors->{");
	
	for (i=0; i<num_of_meshpoint; i++){
		r = 0.5*(1+atan(vertex[i].phi));
		g = 0.5*(1-atan(vertex[i].phi));
		b = 0;
		fprintf(f_ou,"RGBColor[{%lg,%lg,%lg}],",r,g,b);
		//if (vertex[i].phi<0){
		//	fprintf(f_ou,"Black,");
		//} 
		//else {
		//	fprintf(f_ou,"Hue[%lg],",1.0*vertex[i].label/num_of_domains);
		//}
	}

	fseek(f_ou,-1, SEEK_CUR);
	
	fprintf(f_ou,"}]");
}

/*******************************************************************/

void get_time(time_t t, CPU_Time *cpu_time)
{
  long sec_in_day, sec_in_ora, sec_in_min;

  sec_in_day = 86400;
  sec_in_ora = 3600;
  sec_in_min = 60;

  cpu_time->d = t/sec_in_day; t = t%sec_in_day;
  cpu_time->h = t/sec_in_ora; t = t%sec_in_ora;
  cpu_time->m = t/sec_in_min; t = t%sec_in_min;
  cpu_time->s = t;
}

/*******************************************************************/

double ran2(long *idum)
{
  	int j;
  	long k;
  	static long idum2=123456789;
  	static long iy=0;
  	static long iv[NTAB];
  	double temp;
 
  	if (*idum <= 0) { 
    	if (-(*idum) < 1) *idum=1; 
    	else *idum = -(*idum);
    	idum2 = (*idum);
    	for (j=NTAB+7; j>=0; j--) {
      		k = (*idum)/IQ1;
      		*idum = IA1*(*idum-k*IQ1)-k*IR1;
      		if (*idum < 0) *idum += IM1;
      		if (j < NTAB) iv[j] = *idum;
    	}
    	iy = iv[0];
  	}
  	k = (*idum)/IQ1; 
  	*idum = IA1*(*idum-k*IQ1)-k*IR1; 
  	if (*idum < 0) *idum += IM1; 
  	k = idum2/IQ2;
  	idum2 = IA2*(idum2-k*IQ2)-k*IR2; 
  	if (idum2 < 0) idum2 += IM2;
  	j = iy/NDIV; 
  	iy = iv[j]-idum2;
  	iv[j] = *idum; 
  	if (iy < 1) iy += IMM1;
  	if ((temp = AM*iy) > RNMX) return RNMX; 
  	else return temp;
}
