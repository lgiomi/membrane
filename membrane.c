// A solver for the non-local Allen-Cahn equation on curved surfaces
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

// Load program-specific libraries from subdirectory lib/
#include "lib/structures.h"
#include "lib/basic.h"
#include "lib/get_geometry.h"
#include "lib/solver.h"

// Ad-hoc min & max functions
#define  min(x,y) ((x < y) ? x : y)
#define  max(x,y) ((x > y) ? x : y)

//
// GLOBAL VARIABLES
//
// Time variable, to compute the start and end times of the simulation
time_t T1, T2;
// Various flags that determine whether specific settings were called from the line-command options
long N_ITERATIONS;
int C_FLAG=0,O_FLAG=1,I_FLAG=0,G_FLAG=0;
int V_FLAG=0;
int CUTOFF_FLAG=0;
double AREA_FLAG=0,VOLUME_FLAG=0;
double AVG_FLAG=0;
long INT_FLAG,N_EXPORT;
// Time-step related variables
double DT,DT_AUTO,RUN_TIME;
// Phi, \epsilon and tolerance of the evolution
double EPSILON,TOLERANCE;
// Extra parameters/flags
double HK_CUTOFF;
double G_SIGMA;
double CONSERVED=1;
double SCALE_FACTOR=1;
double LAGRANGE,PHI_AVG;
double CURRENT_TIME;
int N_DOMAINS;
long SEED;
int HALT_NOW=0;
//Couplings for the polynomial GL free energy
double GAMMA_H,GAMMA_H2,GAMMA_KG;
double K_BARRIER=1.;
//Couplings of the full thermodynamic model used in my work on antimixing
double T_V,J_V,Lk_V,Lkb_V,Mk_V,Mkb_V;
// Coordinates of center and direction of axis for the 2D cylindrical projection 
double C_X,C_Y,C_Z;
double A_X,A_Y,A_Z;
// Geometry-related global quantiteis
double TOTAL_AREA,TOTAL_VOLUME;
double WILLMORE_ENERGY,EULER_CHI;
double L_MIN,L_MAX,L_AVG;
double H2_MIN,H2_MAX,KG_MIN,KG_MAX;
double H2_AVG_MIN,H2_AVG_MAX,KG_AVG_MIN,KG_AVG_MAX;

//
// GLOBAL FUNCTIONS
//
void init();
void run();
void end();

void init_random_mixed();
void init_random_domains();
void init_import(char *);

void one_step();
void get_rhs(double *);

void euler();
void rk2();
void rk4();
void heun_euler();
void rkf45();

void get_geometry();

void track_domains();

void export_histogram(FILE *, long);
void export_graphic_complex(FILE *, long);
void export_dat(FILE *);

/*******************************************************************/

int main(int argc, char *argv[])
{	
	init(argc, argv);
	time(&T1);
	run();
	end();
	
  	return EXIT_SUCCESS;	
}

/*******************************************************************/

void init(int argc, char *argv[])
{
	double norm;
	char f_name[32];
	char import_name[32];
	int n;                           
   
	if(argc>1){
	RUN_TIME=0;
	EPSILON=0;
	N_ITERATIONS=1;
	N_EXPORT=0;
	GAMMA_H=0;
	GAMMA_H2=0;
	GAMMA_KG=0;
	HK_CUTOFF=0;
	INT_FLAG=0;
	SEED=0;
	TOLERANCE=0;
	A_X=0;
	A_Y=0;
	A_Z=1;

	
	for( n = 1; n < argc; n++ ){         /* Scan through args. */
		switch((int)argv[n][0]){     /* Check for option character "-" */
			case '-':	                  
					switch((int)argv[n][1]){
						case 'm':
								snprintf(f_name, sizeof(f_name), "%s", argv[n+1]);
								printf("Mesh file\t\t: %s\n",f_name);
								n++;
								break;
						case 'h':
								help();
								break;
						case 'L':
								O_FLAG=atol(argv[n+1]);
								printf("Output level\t\t: %d\n",O_FLAG);
								n++;
								break;
						case 't':
								RUN_TIME=atof(argv[n+1]);
								printf("Run time\t\t: %lg\n",RUN_TIME);
								n++;
								break;
						case 'v':
								VOLUME_FLAG=atof(argv[n+1]);
								printf("Rescale total volume to\t: %lg\n",VOLUME_FLAG);
								n++;
								break;
						case 'a':
								AREA_FLAG=atof(argv[n+1]);
								printf("Rescale total area to\t: %lg\n",AREA_FLAG);
								n++;
								break;
						case 'k':
								K_BARRIER=atof(argv[n+1]);
								printf("Potential barrier\t: %lg\n",K_BARRIER);
								n++;
								break;
						case 'e':
								EPSILON=atof(argv[n+1]);
								printf("Epsilon\t\t\t: %lg\n",EPSILON);
								n++;
								break;
						case 'i':
								N_ITERATIONS=atol(argv[n+1]);
								if(N_ITERATIONS==-1){printf("Will compute time step automatically\n");}
								else{printf("Total number of iterations: %ld\n",N_ITERATIONS);};
								n++;
								break;
						case 'x':
								N_EXPORT=atol(argv[n+1]);
								printf("Export configurations\t: every %ld steps\n",N_EXPORT);
								n++;
								break;
						case 'T':
								TOLERANCE=atof(argv[n+1]);
								printf("Tolerance set to\t: %g\n",TOLERANCE);
								n++;
								break;
						case 'l':
								CONSERVED=0;
								printf("Order parameter will not be CONSERVED\n");
								break;
						case 'M':
								AVG_FLAG=1;
								printf("Will use NN-averaged values for curvatures\n");
								break;
						case 'C':
								GAMMA_H=atof(argv[n+1]);
								GAMMA_H2=atof(argv[n+2]);
								GAMMA_KG=atof(argv[n+3]);
								printf("Couplings\t\t: (%lg H, %lg H^2, %lg KG )\n",GAMMA_H,GAMMA_H2,GAMMA_KG);
								n+=3;
								break;
						case 'P':
								C_X=atof(argv[n+1]);
								C_Y=atof(argv[n+2]);
								C_Z=atof(argv[n+3]);
								printf("Center of sphere\t: (%lg,%lg,%lg)\n",C_X,C_Y,C_Z);
								C_FLAG=1;
								n+=3;
								break;
						case 'A':
								A_X=atof(argv[n+1]);
								A_Y=atof(argv[n+2]);
								A_Z=atof(argv[n+3]);
								norm=sqrt(A_X*A_X+A_Y*A_Y+A_Z*A_Z);
								A_X/=norm;
								A_Y/=norm;
								A_Z/=norm;
								printf("North pole direction\t: (%lg,%lg,%lg)\n",A_X,A_Y,A_Z);
								C_FLAG=1;
								n+=3;
								break;
						case 'I':
								switch((int)argv[n+1][0]){
									case '1': INT_FLAG=1;break;
									case '2': INT_FLAG=2;break;
									case '3': INT_FLAG=3;break;
									case '4': INT_FLAG=4;break;
									case '5': INT_FLAG=5;break;
									default : printf("Illegal integration method  %c\n",(int)argv[n+1][0]);
										  printf("\nType ./membrane -h for help\n"); 
										  exit(0);
										  break;
								}
								printf("Integration method\t: %ld\n",INT_FLAG);		
								n++;
								break;
						case 'R':
								snprintf(import_name, sizeof(import_name), "%s", argv[n+1]);
								printf("Initial configuration\t: %s\n",import_name);
								I_FLAG=1;
								n++;
								break;
						case 'r':
								SEED=atol(argv[n+1]);
								PHI_AVG=atof(argv[n+2]);
								printf("Random initial data\t: (seed %ld, total concentration %lg)\n",SEED,PHI_AVG);
								I_FLAG=2;
								n+=2;
								break;
						case 'p':
								SEED=atol(argv[n+1]);
								PHI_AVG=atof(argv[n+2]);
								printf("Fill the mesh with ±1 randomly over NNN domains, with seed %ld and total concentration %lg \n",SEED,PHI_AVG);
								I_FLAG=3;
								n+=2;
								break;     
						case 'g':
								G_SIGMA=atof(argv[n+1]);
								G_FLAG=1;
								printf("Gaussian noise enabled (with dispersion %lg)\n",G_SIGMA);
								n++;
								break;
						case 'u':
								T_V=atof(argv[n+1]);
								J_V=atof(argv[n+2]);
								Lk_V=atof(argv[n+3]);
								Lkb_V=atof(argv[n+4]);
								Mk_V=atof(argv[n+5]);
								Mkb_V=atof(argv[n+6]);
								printf("Mean-field free energy, parameters T=%lg, J=%lg, Lk=%lg, Lkb=%lg, Mk=%lg, Mkb=%lg\n",T_V,J_V,Lk_V,Lkb_V,Mk_V,Mkb_V);
								V_FLAG=1;
								n+=6;
								break;
						case 'w':
								HK_CUTOFF=atof(argv[n+1]);
								printf("Cut-off for curvatures\t: %lg\n",HK_CUTOFF);
								CUTOFF_FLAG=1;
								n++;
								break;
						case 'F':
								printf("Use quintic interaction potential (unbounded, but does not affect Maxwell construction)\n");
								V_FLAG=2;
								break;
        
						default:  
							printf("Illegal option code \"-%c\"\n",(int)argv[n][1]);
							printf("Type ./membrane -h for help\n"); 
							exit(0);
							break;
					}
                 			break;
			default:  
				printf("\nError: give input in the following format:\n"); 
				print_cmd_line();
				printf("Type ./membrane -h for help\n"); 
				exit(0);
                 		break;
			}	
		}
		DT = RUN_TIME/N_ITERATIONS;
		if( RUN_TIME==0 ||  INT_FLAG == 0 || (N_ITERATIONS == 1 && INT_FLAG < 4) || I_FLAG == 0 || AREA_FLAG*VOLUME_FLAG>0 || AREA_FLAG < 0 || VOLUME_FLAG <0){
				printf("\nError: not enough arguments given or wrong parameter values. Write input in the following format:\n"); 
				print_cmd_line();
				printf("Type ./membrane -h for help\n"); 
				exit(0);
		}
		if(N_EXPORT == 0){
			printf("No -x option given, assuming no intermediate output.\n"); 
		}
	}
	else{
		printf("\nYou did not specify enough arguments to parse. Please give input in the following format:\n\n"); 
		print_cmd_line();
		printf("Type ./membrane -h for help\n"); 
		exit(0);
	}
	import_mesh(f_name);
	get_geometry();
	if(TOLERANCE==0 && INT_FLAG>=4){
		printf("You selected -I %ld but no tolerance was set: defaulting to 1E-6.\n",INT_FLAG);
		TOLERANCE=1E-6;
	}
	if(I_FLAG==1){
		init_import(import_name);
	}
	else if(I_FLAG==2){
		init_random_mixed();
	} 
	else if(I_FLAG==3){
		init_random_domains();
	} 
	;
}

/*******************************************************************/

void get_geometry()
{

	int a,b,c;
	double lab,lbc,lca;
	double cota,cotb,cotc;
	double theta_a,theta_b,theta_c;
	double sp,area_t;
	int this;
	double dx, dy, dz, base, x[3],y[3],z[3],rm[3];
	double gc[3];
	double xm[3],xM[3];

	L_MIN=1E10;
	L_MAX=0;
	L_AVG=0;

	H2_MIN=1E10;
	H2_MAX=-1E10;

	KG_MIN=1E10;
	KG_MAX=-1E10;

	H2_AVG_MIN=1E10;
	H2_AVG_MAX=-1E10;

	KG_AVG_MIN=1E10;
	KG_AVG_MAX=-1E10;

	gc[0]=0;
	gc[1]=0;
	gc[2]=0;

	xm[0]=0;
	xm[1]=0;
	xm[2]=0;

	xM[0]=0;
	xM[1]=0;
	xM[2]=0;

	long i, j, obtuse[MAX_SIZE], num_of_obtuse=0;
	
	TOTAL_AREA = 0;	
	TOTAL_VOLUME = 0;
	WILLMORE_ENERGY = 0;	
	EULER_CHI= 0;	

	FILE *f_ou;

	// Setting all vertex variables to zero before looping through triangles

	for (i=0; i<N_GRID_POINTS; i++){
		vertex[i].hx = 0;
		vertex[i].hy = 0;
		vertex[i].hz = 0;
		vertex[i].area = 0;
		vertex[i].kg=2*PI;
		for (j=0; j<vertex[i].num_of_neighbors; j++){
			vertex[i].weight[j] = 0;
		}
	}

	// Loop through triangles computing vertices weights, mean and Gaussian curvatures and areas

	for (i=0; i<N_GRID_TRIANGLES; i++){

		a=triangle[i].v1;
		b=triangle[i].v2;
		c=triangle[i].v3;

		lab=sqrt(pow((vertex[a].x-vertex[b].x),2.)+pow((vertex[a].y-vertex[b].y),2.)+pow((vertex[a].z-vertex[b].z),2.));
		lbc=sqrt(pow((vertex[b].x-vertex[c].x),2.)+pow((vertex[b].y-vertex[c].y),2.)+pow((vertex[b].z-vertex[c].z),2.));
		lca=sqrt(pow((vertex[c].x-vertex[a].x),2.)+pow((vertex[c].y-vertex[a].y),2.)+pow((vertex[c].z-vertex[a].z),2.));

		L_AVG+=.5*(lab+lbc+lca);

		sp=.5*(lab+lbc+lca);
		area_t=sqrt(sp*(sp-lab)*(sp-lbc)*(sp-lca));

		cota=(-lbc*lbc+lca*lca+lab*lab)/(4.*area_t);
		cotb=(-lca*lca+lbc*lbc+lab*lab)/(4.*area_t);
		cotc=(-lab*lab+lca*lca+lbc*lbc)/(4.*area_t);

		if(cota<0){
			obtuse[num_of_obtuse]=i;
			num_of_obtuse++;
			vertex[a].area+=.5*area_t;
			vertex[b].area+=.25*area_t;
			vertex[c].area+=.25*area_t;
		}
		else if(cotb<0){
			obtuse[num_of_obtuse]=i;
			num_of_obtuse++;
			vertex[a].area+=.25*area_t;
			vertex[b].area+=.5*area_t;
			vertex[c].area+=.25*area_t;
		}
		else if(cotc<0){
			obtuse[num_of_obtuse]=i;
			num_of_obtuse++;
			vertex[a].area+=.25*area_t;
			vertex[b].area+=.25*area_t;
			vertex[c].area+=.5*area_t;
		}
		else{
			vertex[a].area+=.125*(lca*lca*cotb+lab*lab*cotc);
			vertex[b].area+=.125*(lbc*lbc*cota+lab*lab*cotc);
			vertex[c].area+=.125*(lca*lca*cotb+lbc*lbc*cota);
		};

		theta_a=atan2(4.*area_t,-lbc*lbc+lca*lca+lab*lab);
		theta_b=atan2(4.*area_t,-lca*lca+lbc*lbc+lab*lab);
		theta_c=atan2(4.*area_t,-lab*lab+lca*lca+lbc*lbc);

		vertex[a].kg-=theta_a;
		vertex[b].kg-=theta_b;
		vertex[c].kg-=theta_c;

		for (j=0; j<vertex[a].num_of_neighbors; j++){
			if(vertex[a].neighbor[j]==b)vertex[a].weight[j] += .5*cotc;
			if(vertex[a].neighbor[j]==c)vertex[a].weight[j] += .5*cotb;
		}
		for (j=0; j<vertex[b].num_of_neighbors; j++){
			if(vertex[b].neighbor[j]==c)vertex[b].weight[j] += .5*cota;
			if(vertex[b].neighbor[j]==a)vertex[b].weight[j] += .5*cotc;
		}
		for (j=0; j<vertex[c].num_of_neighbors; j++){
			if(vertex[c].neighbor[j]==a)vertex[c].weight[j] += .5*cotb;
			if(vertex[c].neighbor[j]==b)vertex[c].weight[j] += .5*cota;
		}

		vertex[a].hx += .25*(cotb*(vertex[c].x-vertex[a].x)+cotc*(vertex[b].x-vertex[a].x));
		vertex[a].hy += .25*(cotb*(vertex[c].y-vertex[a].y)+cotc*(vertex[b].y-vertex[a].y));
		vertex[a].hz += .25*(cotb*(vertex[c].z-vertex[a].z)+cotc*(vertex[b].z-vertex[a].z));

		vertex[b].hx += .25*(cotc*(vertex[a].x-vertex[b].x)+cota*(vertex[c].x-vertex[b].x));
		vertex[b].hy += .25*(cotc*(vertex[a].y-vertex[b].y)+cota*(vertex[c].y-vertex[b].y));
		vertex[b].hz += .25*(cotc*(vertex[a].z-vertex[b].z)+cota*(vertex[c].z-vertex[b].z));

		vertex[c].hx += .25*(cota*(vertex[b].x-vertex[c].x)+cotb*(vertex[a].x-vertex[c].x));
		vertex[c].hy += .25*(cota*(vertex[b].y-vertex[c].y)+cotb*(vertex[a].y-vertex[c].y));
		vertex[c].hz += .25*(cota*(vertex[b].z-vertex[c].z)+cotb*(vertex[a].z-vertex[c].z));

	}

	L_AVG*=2./NUM_OF_EDGES;

	if(O_FLAG>=1){
		f_ou = fopen("triangles.dat","w");
		for (i=0; i<N_GRID_TRIANGLES; i++){
			fprintf(f_ou,"%ld\t%ld\t%ld\n",
				triangle[i].v1,
				triangle[i].v2,
				triangle[i].v3);
		}
		fclose(f_ou);
	}

	if(A_X*A_X<1.){

		x[0]=0;
		x[1]=-A_Z/sqrt(A_Z*A_Z+A_Y*A_Y);
		x[2]=A_Y/sqrt(A_Z*A_Z+A_Y*A_Y);

		z[0]=x[1]*A_Z-x[2]*A_Y;
		z[1]=x[2]*A_X-x[0]*A_Z;
		z[2]=x[0]*A_Y-x[1]*A_X;

		z[0]/=sqrt(z[0]*z[0]+z[1]*z[1]+z[2]*z[2]);
		z[1]/=sqrt(z[0]*z[0]+z[1]*z[1]+z[2]*z[2]);
		z[2]/=sqrt(z[0]*z[0]+z[1]*z[1]+z[2]*z[2]);
	}else{

		x[0]=0;
		x[1]=-1;
		x[2]=0;

		z[0]=0;
		z[1]=0;
		z[2]=1;

	};

	//printf("Orthonomal basis for spherical projection: n=(%lg,%lg,%lg), m=(%lg,%lg,%lg), o=(%lg,%lg,%lg)\n",A_X,A_Y,A_Z,x[0],x[1],x[2],z[0],z[1],z[2]);

	// One further loop through vertices to normalize things

	for (i=0; i<N_GRID_POINTS; i++){

		vertex[i].nx = -vertex[i].hx;
		vertex[i].ny = -vertex[i].hy;
		vertex[i].nz = -vertex[i].hz;

		base = sqrt(vertex[i].nx*vertex[i].nx+vertex[i].ny*vertex[i].ny+vertex[i].nz*vertex[i].nz);

		if(base>0){
		vertex[i].nx /= base;	
		vertex[i].ny /= base;
		vertex[i].nz /= base;
		};

		vertex[i].hx /= vertex[i].area;
		vertex[i].hy /= vertex[i].area;
		vertex[i].hz /= vertex[i].area;

		vertex[i].h2 = vertex[i].hx*vertex[i].hx+vertex[i].hy*vertex[i].hy+vertex[i].hz*vertex[i].hz;

		vertex[i].kg /= vertex[i].area;

		vertex[i].h2_avg = vertex[i].h2*vertex[i].area;
		vertex[i].kg_avg = vertex[i].kg*vertex[i].area;

		for (j=0; j<vertex[i].num_of_neighbors; j++){
			vertex[i].weight[j] /= vertex[i].area;	
		}

		// Compute projective angles

		if(C_FLAG==1){

			dx = vertex[i].x-C_X;
			dy = vertex[i].y-C_Y;
			dz = vertex[i].z-C_Z;

			dx/=sqrt(dx*dx+dy*dy+dz*dz);
			dy/=sqrt(dx*dx+dy*dy+dz*dz);
			dz/=sqrt(dx*dx+dy*dy+dz*dz);

			base=dx*A_X+dy*A_Y+dz*A_Z;

			rm[0]=dx-base*A_X;
			rm[1]=dy-base*A_Y;
			rm[2]=dz-base*A_Z;

			vertex[i].gridx=atan2(rm[0]*x[0]+rm[1]*x[1]+rm[2]*x[2],rm[0]*z[0]+rm[1]*z[1]+rm[2]*z[2]);
			vertex[i].gridy=atan2(base,sqrt(rm[0]*rm[0]+rm[1]*rm[1]+rm[2]*rm[2]));

		}

		TOTAL_AREA+=vertex[i].area;
		WILLMORE_ENERGY += vertex[i].area*vertex[i].h2;
		EULER_CHI += .5/PI*vertex[i].kg*vertex[i].area;

	}

	// Try to point all normal vectors outwards: this does only do a single loop through vertices.
	// It is not safe to assume that after just this iteration any mesh will be nicely oriented. 
	// However, it has been working so far, but remember, the general solution is more complicated
	// Moreover, use this loop to compute averaged values of H2 and KG

	for (i=0; i<N_GRID_POINTS; i++){

		area_t = 0;

		for (j=0; j<vertex[i].num_of_neighbors; j++){

			a = vertex[i].neighbor[j];

			base = vertex[i].nx*vertex[a].nx+vertex[i].ny*vertex[a].ny+vertex[i].nz*vertex[a].nz;

			if(base<0){vertex[a].nx*=-1;vertex[a].ny*=-1;vertex[a].nz*=-1;};

			vertex[i].h2_avg+=vertex[a].h2*vertex[a].area;
			vertex[i].kg_avg+=vertex[a].kg*vertex[a].area;
			area_t+=vertex[a].area;

		}

		vertex[i].h2_avg/=area_t+vertex[i].area;
		vertex[i].kg_avg/=area_t+vertex[i].area;
		
	}

	// Compute the total volume

	for (i=0; i<N_GRID_POINTS; i++){
		TOTAL_VOLUME+=1./3.*vertex[i].area*(vertex[i].nx*vertex[i].x+vertex[i].ny*vertex[i].y+vertex[i].nz*vertex[i].z);
	}

	// If rescaling was chosen, rescale all Vertex attributes accordingly, rescale mesh details and recompute area and volume

	if(AREA_FLAG>0 || VOLUME_FLAG>0){

		if(AREA_FLAG>0)SCALE_FACTOR=pow(AREA_FLAG/TOTAL_AREA,1./2.);
		if(VOLUME_FLAG>0)SCALE_FACTOR=pow(VOLUME_FLAG/TOTAL_VOLUME,1./3.);

		TOTAL_AREA=0;
		TOTAL_VOLUME=0;
		
		for (i=0; i<N_GRID_POINTS; i++){

			vertex[i].x*=SCALE_FACTOR;
			vertex[i].y*=SCALE_FACTOR;
			vertex[i].z*=SCALE_FACTOR;

			vertex[i].hx/=SCALE_FACTOR;
			vertex[i].hy/=SCALE_FACTOR;
			vertex[i].hz/=SCALE_FACTOR;

			vertex[i].h2/=pow(SCALE_FACTOR,2.);
			vertex[i].kg/=pow(SCALE_FACTOR,2.);

			vertex[i].h2_avg/=pow(SCALE_FACTOR,2.);
			vertex[i].kg_avg/=pow(SCALE_FACTOR,2.);

			vertex[i].area*=pow(SCALE_FACTOR,2.);

			for (j=0; j<vertex[i].num_of_neighbors; j++){
				vertex[i].weight[j]/=pow(SCALE_FACTOR,2.);
			}

			TOTAL_VOLUME+=1./3.*vertex[i].area*(vertex[i].nx*vertex[i].x+vertex[i].ny*vertex[i].y+vertex[i].nz*vertex[i].z);
			TOTAL_AREA+=vertex[i].area;
		}
	
		L_AVG*=SCALE_FACTOR;

	};

	if(CUTOFF_FLAG>0){

		for (i=0; i<N_GRID_POINTS; i++){

			vertex[i].h2=min(vertex[i].h2,HK_CUTOFF);
			vertex[i].kg=min(vertex[i].kg,HK_CUTOFF);
			vertex[i].kg=max(vertex[i].kg,-HK_CUTOFF);

			vertex[i].h2_avg=min(vertex[i].h2_avg,HK_CUTOFF);
			vertex[i].kg_avg=min(vertex[i].kg_avg,HK_CUTOFF);
			vertex[i].kg_avg=max(vertex[i].kg_avg,-HK_CUTOFF);

		}
	
	};

	// Compute the geometrical center and the boxing size 

	for (i=0; i<N_GRID_POINTS; i++){

		gc[0]+=1./N_GRID_POINTS*vertex[i].x;
		gc[1]+=1./N_GRID_POINTS*vertex[i].y;
		gc[2]+=1./N_GRID_POINTS*vertex[i].z;

		xm[0]=vertex[i].x<xm[0]?vertex[i].x:xm[0];
		xm[1]=vertex[i].y<xm[1]?vertex[i].y:xm[1];
		xm[2]=vertex[i].z<xm[2]?vertex[i].z:xm[2];

		xM[0]=vertex[i].x>xM[0]?vertex[i].x:xM[0];
		xM[1]=vertex[i].y>xM[1]?vertex[i].y:xM[1];
		xM[2]=vertex[i].z>xM[2]?vertex[i].z:xM[2];

	}

	// If -O 1 or higher, save geometry.dat
	
	if(O_FLAG>=1)f_ou = fopen("geometry.dat","w");

	// Moreover, compute min/max lengths and curvatures

	for (i=0; i<N_GRID_POINTS; i++){

		if(O_FLAG>=1){if(C_FLAG==0){
			fprintf(f_ou,"%ld\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10lg\t%.10lg\t%.10lg\t%.10lg\t",
				i,
				vertex[i].x,
				vertex[i].y,
				vertex[i].z,
				vertex[i].area,
				vertex[i].h2,
				vertex[i].kg,
				vertex[i].h2_avg,
				vertex[i].kg_avg,
				vertex[i].nx,
				vertex[i].ny,
				vertex[i].nz);
		}else if(C_FLAG==1){
			fprintf(f_ou,"%ld\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10lg\t%.10lg\t%.10lg\t%.10lg\t%.10lg\t%.10lg\t",
				i,
				vertex[i].x,
				vertex[i].y,
				vertex[i].z,
				vertex[i].area,
				vertex[i].h2,
				vertex[i].kg,
				vertex[i].h2_avg,
				vertex[i].kg_avg,
				vertex[i].gridx,
				vertex[i].gridy,
				vertex[i].nx,
				vertex[i].ny,
				vertex[i].nz);
		};}
			
		for (j=0; j<vertex[i].num_of_neighbors; j++){
			if(O_FLAG>=1)fprintf(f_ou,"%.10f\t",vertex[i].weight[j]);

			this = vertex[i].neighbor[j];
		
			x[0] = vertex[i].x-vertex[this].x;
			x[1] = vertex[i].y-vertex[this].y;
			x[2] = vertex[i].z-vertex[this].z;
			base=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);

			if(L_MAX<base){L_MAX=base;};
			if(L_MIN>base){L_MIN=base;};

			if(H2_MAX<vertex[i].h2){H2_MAX=vertex[i].h2;};
			if(H2_MIN>vertex[i].h2){H2_MIN=vertex[i].h2;};

			if(KG_MAX<vertex[i].kg){KG_MAX=vertex[i].kg;};
			if(KG_MIN>vertex[i].kg){KG_MIN=vertex[i].kg;};

			if(H2_AVG_MAX<vertex[i].h2_avg){H2_AVG_MAX=vertex[i].h2_avg;};
			if(H2_AVG_MIN>vertex[i].h2_avg){H2_AVG_MIN=vertex[i].h2_avg;};

			if(KG_AVG_MAX<vertex[i].kg_avg){KG_AVG_MAX=vertex[i].kg_avg;};
			if(KG_AVG_MIN>vertex[i].kg_avg){KG_AVG_MIN=vertex[i].kg_avg;};
				
		}	
			
		if(O_FLAG>=1)fprintf(f_ou,"\n");
	}	
	if(O_FLAG>=1)fclose(f_ou);

	if(O_FLAG==2){
		f_ou = fopen("mean_curvature.m","w");
		export_graphic_complex(f_ou,2);
		fclose(f_ou);

		f_ou = fopen("gaussian_curvature.m","w");
		export_graphic_complex(f_ou,3);
		fclose(f_ou);

	}

	// Print summary of mesh computations
	printf("\n");

	printf("\tVertices\t\t\t: %ld\n",N_GRID_POINTS);
	printf("\tTriangles\t\t\t: %ld\n",N_GRID_TRIANGLES);
	printf("\tObtuse triangles\t\t: %ld (%.1ld%% of total)\n",num_of_obtuse/2,100*num_of_obtuse/2/N_GRID_TRIANGLES);

	printf("\tGeometrical center\t\t: (%2.3f,%2.3f,%2.3f)\n",gc[0],gc[1],gc[2]);
	printf("\tCoordinate dimensions\t\t: (%2.3f,%2.3f,%2.3f)\n",xM[0]-xm[0],xM[1]-xm[1],xM[2]-xm[2]);

	printf("\tTotal surface area\t\t: %lg (typical length %lg)\n",TOTAL_AREA,sqrt(TOTAL_AREA));
	printf("\tTotal enclosed volume\t\t: %lg (typical length %lg)\n",TOTAL_VOLUME,pow(TOTAL_VOLUME,1./3.));

	printf("\tWillmore energy\t\t\t: %lg (asphericity %lg)\n",WILLMORE_ENERGY,WILLMORE_ENERGY/4/PI-1);
	printf("\tEuler characteristic\t\t: %lg\n",EULER_CHI);

	printf("\tMin/max/avg edge length \t: (%lg,%lg,%lg)\n",L_MIN,L_MAX,L_AVG);

	printf("\tMin/max/avg H2 \t\t\t: (%lg,%lg,%lg)\n",H2_MIN,H2_MAX,WILLMORE_ENERGY/TOTAL_AREA);
	printf("\tMin/max/avg KG \t\t\t: (%lg,%lg,%lg)\n",KG_MIN,KG_MAX,2*PI*EULER_CHI/TOTAL_AREA);

	printf("\tMin/max H2 (NN averaged)\t: (%lg,%lg)\n",H2_AVG_MIN,H2_AVG_MAX);
	printf("\tMin/max KG (NN averaged)\t: (%lg,%lg)\n",KG_AVG_MIN,KG_AVG_MAX);

	// The interface profile anywhere on the surface should countain at least six grid points:
	if(EPSILON==0){
		EPSILON=sqrt(K_BARRIER)*L_MAX;
		printf("\tSetting EPSILON to\t\t: %f\n",EPSILON);
	}else{
		printf("\tProposed EPSILON\t\t: %lg (minimal), %lg (averaged) or %lg (conservative)\n",L_MAX*sqrt(K_BARRIER/2),L_AVG*sqrt(K_BARRIER/2),1.1*L_MAX*sqrt(K_BARRIER/2));
	}

	if(V_FLAG==0 || V_FLAG ==2){
		// Interface thickness can be computed independently of curvature only in the plynomial GL model
		printf("\tInterface thickness\t\t: %lg (roughly %2.1f%% - %2.1f%% of membrane size)\n",3*sqrt(2/K_BARRIER)*EPSILON,3*sqrt(2/K_BARRIER/TOTAL_AREA)*EPSILON*100.,3*sqrt(2/K_BARRIER)*EPSILON/pow(TOTAL_VOLUME,1./3.)*100.);
		// Line tension can be computed independently of curvature only in the polynomial GL model
		printf("\tExpected line tension\t\t: %lg\n",2./3.*sqrt(2*K_BARRIER)*EPSILON);
		// What you really want to set in -C is "adimensional", not just \Delta k. 
		// For this reason we have to rescale the couplings that will enter in the flow equation
		// It used to be that the \eta's were assumed in -C, and the GAMMA's had to be rescaled as
		//GAMMA_H*=2./3.*sqrt(2*K_BARRIER*TOTAL_AREA);
		//GAMMA_H2*=2./3.*sqrt(2*K_BARRIER*TOTAL_AREA);
		//GAMMA_KG*=2./3.*sqrt(2*K_BARRIER*TOTAL_AREA);
		// However now it is more convenient to rescale only with \sigma, so that by keeping fixed the -C values 
		// and rescaling the geometry will produce an effect
		GAMMA_H*=2./3.*sqrt(2*K_BARRIER);
		GAMMA_H2*=2./3.*sqrt(2*K_BARRIER);
		GAMMA_KG*=2./3.*sqrt(2*K_BARRIER);

		printf("\tCouplings entering EOMs\t\t: (%lg H, %lg H^2, %lg KG)\n",GAMMA_H,GAMMA_H2,GAMMA_KG);
	}

	// The couplings in the \epsilon \to 0 are unbounded.s
	// However, since numerically EPSILON is finite, this sets a bound on the magnitude of the curvature couplings in order to preserve the double-well structure of the potential	

	if(V_FLAG==0){

		
		//printf("\tExtremal value of deltas \t: (%lg H, %lg H^2, %lg - %lg KG)\n",.75/K_BARRIER*GAMMA_H*EPSILON*sqrt(H2_MAX),.75/K_BARRIER*GAMMA_H2*EPSILON*H2_MAX,.75/K_BARRIER*GAMMA_KG*EPSILON*KG_MIN,.75/K_BARRIER*GAMMA_KG*EPSILON*KG_MAX);	
	
		printf("\tAllowed coupling ranges\t\t: |g_H|<%lg, |Delta k|<%lg, |Delta k_b|<%lg\n",4./3.*K_BARRIER/EPSILON/sqrt(H2_MAX),4./3.*K_BARRIER/EPSILON/H2_MAX,4./3.*K_BARRIER/EPSILON/max(KG_MAX,sqrt(KG_MIN*KG_MIN)));  

		printf("\tAllowed coupling ranges (NN avg): |g_H|<%lg, |Delta k|<%lg, |Delta k_b|<%lg\n",4./3.*K_BARRIER/EPSILON/sqrt(H2_AVG_MAX),4./3.*K_BARRIER/EPSILON/H2_AVG_MAX,4./3.*K_BARRIER/EPSILON/max(KG_AVG_MAX,sqrt(KG_AVG_MIN*KG_AVG_MIN))); 
	}

	if(V_FLAG==2){

		printf("\tAllowed coupling ranges\t\t: |g_H|<%lg, |Delta k|<%lg, |Delta k_b|<%lg\n",1./2.*K_BARRIER/EPSILON/sqrt(H2_MAX),1./2.*K_BARRIER/EPSILON/H2_MAX,1./2.*K_BARRIER/EPSILON/max(KG_MAX,sqrt(KG_MIN*KG_MIN)));  

		printf("\tAllowed coupling ranges (NN avg): |g_H|<%lg, |Delta k|<%lg, |Delta k_b|<%lg\n",1./2.*K_BARRIER/EPSILON/sqrt(H2_AVG_MAX),1./2.*K_BARRIER/EPSILON/H2_AVG_MAX,1./2.*K_BARRIER/EPSILON/max(KG_AVG_MAX,sqrt(KG_AVG_MIN*KG_AVG_MIN))); 
	}

	// For diffusion processes are stable only for DT/DX^2 \simeq 1/2:
	DT_AUTO=.1*L_MIN*L_MIN/2;

	if(INT_FLAG>=4){
		DT=DT_AUTO;
		printf("\tSetting initial time-step to\t: %lg\n",DT);
	}else{
		printf("\tProposed initial time-step\t: %lg\n",DT_AUTO);  
	}

	printf("\n");

}


/*******************************************************************/

void init_random_mixed()
{
	double ran2(long *);
	long i;
		
	for (i=0; i<N_GRID_POINTS; i++){
		vertex[i].phi = (2*PHI_AVG-1)+NOISE*(1-2*ran2(&SEED));
		if( V_FLAG==1 && vertex[i].phi>1. ){vertex[i].phi=.99;};
		if( V_FLAG==1 && vertex[i].phi<-1. ){vertex[i].phi=-.99;};	
	}
}

/*******************************************************************/

void init_import(char *import_name)
{
	int i,lines=0;
	char line[LINESIZE];

	if( access(import_name, F_OK ) == -1 ) {
		printf("\nError: file %s does not exist\n",import_name);
		exit(0);
	};

	FILE *f_in = fopen(import_name,"r");

	while(fgets(line,LINESIZE-1,f_in) != NULL)
	{
 		lines++;
	};
	fclose(f_in);

	//printf("Total number of lines : %d\n",lines);

	if(lines==N_GRID_POINTS){
		PHI_AVG=0;
		FILE *f_in = fopen(import_name,"r");
		for (i=0; i<N_GRID_POINTS; i++){
			if(fscanf(f_in,"%lg%*[^\n]",&vertex[i].phi)){
				PHI_AVG+=vertex[i].area*vertex[i].phi;
			}else{
				printf("Error: failed to read mesh point in %s.\n",import_name);
				exit(0);
			};
		}

		fclose(f_in);
	}else{
		printf("Error, number of lines of %s does not match mesh size\n",import_name);
		exit(0);
	};
	PHI_AVG/=TOTAL_AREA;
	PHI_AVG=.5*(1+PHI_AVG);
	printf("Imported initial configuration from %s with mean concentration %g\n",import_name,PHI_AVG);
}

/*******************************************************************/

void init_random_domains()
{
	long i,j,j_this,k,k_this;
	double phi_1=-.9,phi_2=.9;
	double ran2(long *),PHI_AVG_temp=0;

	if(PHI_AVG < .5*(1.+phi_1) || PHI_AVG > .5*(1.+phi_2)){
		printf("Warning: provided C0 does not lie in the range [%2.3g,%2.3g]",.5*(1.+phi_1),5*(1.+phi_2));
		phi_1 = min(phi_1,2*PHI_AVG-1);
		phi_2 = max(phi_2,2*PHI_AVG-1);
		printf(", reshifting to [%2.3g,%2.3g]\n",.5*(1.+phi_1),5*(1.+phi_2));
	};

	for (i=0; i<N_GRID_POINTS; i++){
		vertex[i].phi = phi_1;
		PHI_AVG_temp +=.5*(1.+phi_1)*vertex[i].area/TOTAL_AREA;
	}

	while(PHI_AVG_temp < PHI_AVG){
		i = floor(ran2(&SEED)*N_GRID_POINTS);
		//printf("Random vertex %d\n",i);
		if(vertex[i].phi < 0){
			vertex[i].phi= phi_2;
			PHI_AVG_temp +=.5*(phi_2-phi_1)*vertex[i].area/TOTAL_AREA;
		}
		for (j=0; j<vertex[i].num_of_neighbors; j++){
			j_this = vertex[i].neighbor[j];
			if(vertex[j_this].phi < 0){
				vertex[j_this].phi= phi_2;
			        PHI_AVG_temp +=.5*(phi_2-phi_1)*vertex[j_this].area/TOTAL_AREA;
			}
			for (k=0; k<vertex[j_this].num_of_neighbors; k++){
				k_this = vertex[j_this].neighbor[k];
				if(vertex[k_this].phi < 0){
					vertex[k_this].phi= phi_2;
				        PHI_AVG_temp +=.5*(phi_2-phi_1)*vertex[k_this].area/TOTAL_AREA;
				}
			}
		}
	};

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

// Define the part of the free energy density which does not depend on gradients (i.e. the potential)

double V(long i)
{
	double V_0,V_I;
	double phi,Phi;

	phi = vertex[i].phi;

	if(V_FLAG == 0){

		V_0 = K_BARRIER*.25*pow(pow(phi,2.)-1,2.);
		V_I = .25*pow((phi+1),2.)*(2-phi);

		if(AVG_FLAG==0)return V_0+V_I*EPSILON*(GAMMA_H2*vertex[i].h2+GAMMA_KG*vertex[i].kg+GAMMA_H*sqrt(vertex[i].h2));
		if(AVG_FLAG!=0)return V_0+V_I*EPSILON*(GAMMA_H2*vertex[i].h2_avg+GAMMA_KG*vertex[i].kg_avg+GAMMA_H*sqrt(vertex[i].h2_avg));

	} else if(V_FLAG == 1){

		Phi = .5*(phi+1.);

		V_0 = T_V*(Phi*log(Phi)+(1.-Phi)*log(1.-Phi))+J_V*Phi*(1.-Phi);
		V_I = Phi*(1.-Phi);

		if(AVG_FLAG==0)return V_0+V_I*(Lk_V*vertex[i].h2+Lkb_V*vertex[i].kg)+Phi*(Mk_V*vertex[i].h2+Mkb_V*vertex[i].kg);
		if(AVG_FLAG!=0)return V_0+V_I*(Lk_V*vertex[i].h2_avg+Lkb_V*vertex[i].kg_avg)+Phi*(Mk_V*vertex[i].h2_avg+Mkb_V*vertex[i].kg_avg);

	} else if(V_FLAG == 2){

		V_0 = K_BARRIER*.25*pow(pow(phi,2.)-1,2.);
		V_I = .5+pow(phi,3.)-pow(phi,5.);

		if(AVG_FLAG==0)return V_0+V_I*EPSILON*(GAMMA_H2*vertex[i].h2+GAMMA_KG*vertex[i].kg+GAMMA_H*sqrt(vertex[i].h2));
		if(AVG_FLAG!=0)return V_0+V_I*EPSILON*(GAMMA_H2*vertex[i].h2_avg+GAMMA_KG*vertex[i].kg_avg+GAMMA_H*sqrt(vertex[i].h2_avg));

	};;
	
}

// Homogeneous part of the potential

double Vp(long i)
{
	double V_0;
	double phi,Phi;

	phi = vertex[i].phi;

	if(V_FLAG == 0){

		V_0 = K_BARRIER*.25*pow(pow(phi,2.)-1,2.);

	} else if(V_FLAG == 1){

		Phi = .5*(phi+1.);

		V_0 = T_V*(Phi*log(Phi)+(1.-Phi)*log(1.-Phi))+J_V*Phi*(1.-Phi);

	};

	return V_0;
}

// The reaction term (i.e. the first derivative of the potential)

double dV(long i)
{
	double dV_0,dV_I;
	double phi,Phi;

	phi = vertex[i].phi;

	if(V_FLAG == 0){

		dV_0 = K_BARRIER*phi*(phi*phi-1);
		dV_I = .75*(1-phi*phi);

		if(AVG_FLAG==0)return dV_0+dV_I*EPSILON*(GAMMA_H2*vertex[i].h2+GAMMA_KG*vertex[i].kg+GAMMA_H*sqrt(vertex[i].h2));
		if(AVG_FLAG!=0)return dV_0+dV_I*EPSILON*(GAMMA_H2*vertex[i].h2_avg+GAMMA_KG*vertex[i].kg_avg+GAMMA_H*sqrt(vertex[i].h2_avg));

	} else if(V_FLAG == 1){

		Phi = .5*(phi+1.);

		dV_0 = 2.*T_V*atanh(phi)-J_V*phi;

		if(AVG_FLAG==0)return dV_0+(-phi*Lk_V+Mk_V)*vertex[i].h2+(-phi*Lkb_V+Mkb_V)*vertex[i].kg;
		if(AVG_FLAG!=0)return dV_0+(-phi*Lk_V+Mk_V)*vertex[i].h2_avg+(-phi*Lkb_V+Mkb_V)*vertex[i].kg_avg;

	} else if(V_FLAG == 2){

		dV_0 = K_BARRIER*phi*(phi*phi-1);
		dV_I = -.5*pow(phi,2.)*(5.-6.*pow(phi,2.));

		if(AVG_FLAG==0)return dV_0+dV_I*EPSILON*(GAMMA_H2*vertex[i].h2+GAMMA_KG*vertex[i].kg+GAMMA_H*sqrt(vertex[i].h2));
		if(AVG_FLAG!=0)return dV_0+dV_I*EPSILON*(GAMMA_H2*vertex[i].h2_avg+GAMMA_KG*vertex[i].kg_avg+GAMMA_H*sqrt(vertex[i].h2_avg));

	};
	

}

/*******************************************************************/

void get_rhs(double *rhs)
{
	double laplace(long);
	double dV(long);
	long i;
	
	// Calculate the Lagrange multiplier
	
	LAGRANGE = 0;

	for (i=0; i<N_GRID_POINTS; i++){

		LAGRANGE += vertex[i].area*(EPSILON*laplace(i)-dV(i)/EPSILON);		
	}		
	
	LAGRANGE /= TOTAL_AREA;
	
	// Calculate the right-hand side of the Allen-Cahn equation
	
	for (i=0; i<N_GRID_POINTS; i++){
		if(G_FLAG==0){
			rhs[i] = EPSILON*laplace(i)-dV(i)/EPSILON-CONSERVED*LAGRANGE;
		}else{
			rhs[i] = EPSILON*laplace(i)-dV(i)/EPSILON-CONSERVED*LAGRANGE+gauss_ran2(&SEED,0,G_SIGMA);
		}
	}
}

/*******************************************************************/

void run()
{
	long t=0;
	char f_math_na[32],f_dat_na[32];
	
	FILE *f_hi, *f_ou;
	
	if(O_FLAG>=1)f_hi = fopen("histo.dat","w");	

	CURRENT_TIME=0;

	if(N_EXPORT>0 && C_FLAG==0 && O_FLAG!=2){
		printf("You chose to export configurations every %ld time steps, but did not specify a center nor chose -O 2. Assuming -P 0 0 0.\n",N_EXPORT);
		C_FLAG=1;
		C_X=0;
		C_Y=0;
		C_Z=0;
	}

	if(N_ITERATIONS==-1 || INT_FLAG == 4 || INT_FLAG == 5){
		while(CURRENT_TIME<RUN_TIME && HALT_NOW==0){

			CURRENT_TIME+=DT;

			if(O_FLAG>=1)export_histogram(f_hi,t);

			if (N_EXPORT>0 && t%N_EXPORT==0){

				if(O_FLAG == 2 && C_FLAG==1){
					sprintf(f_math_na,"gc_%08ld.m",t);
					f_ou = fopen(f_math_na,"w");
					export_graphic_complex(f_ou,1);
					fclose(f_ou);
				}
				if(C_FLAG==1 && O_FLAG>=1){
					sprintf(f_dat_na,"gc_%08ld.dat",t);
					f_ou = fopen(f_dat_na,"w");
					export_dat(f_ou);
					fclose(f_ou);
				};
			}
			t++;
			if(t%100==0){printf("Step: %ld\t\t Current time:%.6g\r",t,CURRENT_TIME);}
			one_step();
		}
		printf("Total time steps: %ld\nSimulation time: %g\n",t,CURRENT_TIME);
	}
	else{
		for (t=0; t<N_ITERATIONS; t++){

			CURRENT_TIME+=DT;

			if(O_FLAG>=1)export_histogram(f_hi,t);

			if (N_EXPORT>0 && t%N_EXPORT==0){
				if(O_FLAG == 2 && C_FLAG==1){
					sprintf(f_math_na,"gc_%08ld.m",t);
					f_ou = fopen(f_math_na,"w");
					export_graphic_complex(f_ou,1);
					fclose(f_ou);
				}
				if(C_FLAG==1 && O_FLAG>=1){
					sprintf(f_dat_na,"gc_%08ld.dat",t);
					f_ou = fopen(f_dat_na,"w");
					export_dat(f_ou);
					fclose(f_ou);
				};
			}
		
			//track_domains();
			progress_bar(t,N_ITERATIONS);
			one_step();
		}
	};
	if(O_FLAG>=1)fclose(f_hi);
}

/*******************************************************************/

void one_step()
{
	switch(INT_FLAG){
  	case 1:
		euler();
    	break;
  	case 2:
		rk2();
    	break;
  	case 3:
		rk4();
    	break;
  	case 4:
		heun_euler();
    	break;
  	case 5:
		rkf45();
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
	
	for (i=0; i<N_GRID_POINTS; i++){
		vertex[i].phi += h*rhs[i];
	}
}

/*******************************************************************/

void rk2()
{
	double rhs0[MAX_SIZE]; 
	double rhs1[MAX_SIZE];
	double rhs2[MAX_SIZE];
	long i;
	
	for (i=0; i<N_GRID_POINTS; i++){
		rhs0[i] = vertex[i].phi; 
	}
	
	get_rhs(rhs1);
	
	for (i=0; i<N_GRID_POINTS; i++){
		vertex[i].phi = rhs0[i]+DT*rhs1[i];
	}
	
	get_rhs(rhs2);
	
	for (i=0; i<N_GRID_POINTS; i++){
		vertex[i].phi = rhs0[i]+0.5*DT*(rhs1[i]+rhs2[i]);
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
	
	for (i=0; i<N_GRID_POINTS; i++){
		rhs0[i] = vertex[i].phi; 
	}
	
	get_rhs(rhs1);
	
	for (i=0; i<N_GRID_POINTS; i++){
		vertex[i].phi = rhs0[i]+0.5*DT*rhs1[i];
	}
	
	get_rhs(rhs2);
	
	for (i=0; i<N_GRID_POINTS; i++){
		vertex[i].phi = rhs0[i]+0.5*DT*rhs2[i];
	}
	
	get_rhs(rhs3);
	
	for (i=0; i<N_GRID_POINTS; i++){
		vertex[i].phi = rhs0[i]+DT*rhs3[i];
	}
	
	get_rhs(rhs4);
	
	for (i=0; i<N_GRID_POINTS; i++){
		vertex[i].phi = rhs0[i]+h*(rhs1[i]+2*rhs2[i]+2*rhs3[i]+rhs4[i]);
	}
}

/*******************************************************************/

void heun_euler()
{
	double rhs0[MAX_SIZE];
	double rhs1[MAX_SIZE];
	double rhs2[MAX_SIZE];
	double Q[MAX_SIZE];
	double Q_average=0,delta,rhs_average=0;
	double DTmin=DT_AUTO*1E-5,DTmax=DT_AUTO*1E5;
	long i;

	FILE *f_ou;

	for (i=0; i<N_GRID_POINTS; i++){
		rhs0[i] = vertex[i].phi; 
	}

	get_rhs(rhs1);
	
	for (i=0; i<N_GRID_POINTS; i++){
		vertex[i].phi = rhs0[i]+DT*rhs1[i];
	}
	
	get_rhs(rhs2);

	for (i=0; i<N_GRID_POINTS; i++){
		vertex[i].phi = rhs0[i]+0.5*DT*(rhs1[i]+rhs2[i]);
		Q[i] = (rhs1[i]-rhs2[i])/2.;
		if(Q[i]<0)Q[i]*=-1;
		Q_average+=vertex[i].area*Q[i]/TOTAL_AREA;
		rhs_average+=vertex[i].area*sqrt(rhs2[i]*rhs2[i])/TOTAL_AREA;
	}

	if(rhs_average<TOLERANCE)HALT_NOW=1;

	delta=pow(20.*TOLERANCE/Q_average,.5);

	if(delta<.1){DT*=.1;}
	else if(delta>4.){DT*=4.;}
	else {DT*=delta;};

	if(DT<DTmin)DT=DTmin;
	if(DT>DTmax)DT=DTmax;

	if(O_FLAG>=3){
		f_ou = fopen("Q_avg_debug.dat","a");
		fprintf(f_ou,"%.8E\t %.8E\t %.8E\t %.8E\t %.8E\n",CURRENT_TIME,Q_average,delta,DT,rhs_average);
		fclose(f_ou);
	}
}

/*******************************************************************/

void rkf45()
{
	double rhs0[MAX_SIZE];
	double rhs1[MAX_SIZE];
	double rhs2[MAX_SIZE];
	double rhs3[MAX_SIZE];
	double rhs4[MAX_SIZE];
	double rhs5[MAX_SIZE];
	double rhs6[MAX_SIZE];
	double Q[MAX_SIZE];
	double Q_average=0,delta,rhs_average=0;
	long i;

	FILE *f_ou;

	for (i=0; i<N_GRID_POINTS; i++){
		rhs0[i] = vertex[i].phi; 
	}
	
	get_rhs(rhs1);
	
	for (i=0; i<N_GRID_POINTS; i++){
		vertex[i].phi = rhs0[i]+1./4.*DT*rhs1[i];
	}
	
	get_rhs(rhs2);
	
	for (i=0; i<N_GRID_POINTS; i++){
		vertex[i].phi = rhs0[i]+1./32.*DT*(3.*rhs1[i]+9.*rhs2[i]);
	}
	
	get_rhs(rhs3);
	
	for (i=0; i<N_GRID_POINTS; i++){
		vertex[i].phi = rhs0[i]+1./2197.*DT*(1932.*rhs1[i]-7200.*rhs2[i]+7296.*rhs3[i]);
	}
	
	get_rhs(rhs4);
	
	for (i=0; i<N_GRID_POINTS; i++){
		vertex[i].phi = rhs0[i]+DT*(439./216.*rhs1[i]-8.*rhs2[i]+3680./513.*rhs3[i]-845./4104.*rhs4[i]);
	}

	get_rhs(rhs5);

	for (i=0; i<N_GRID_POINTS; i++){
		vertex[i].phi = rhs0[i]+DT*(-8./27.*rhs1[i]+2.*rhs2[i]-3544./2565.*rhs3[i]+1859./4104.*rhs4[i]-11./40.*rhs5[i]);
	}

	get_rhs(rhs6);

	for (i=0; i<N_GRID_POINTS; i++){
		vertex[i].phi = rhs0[i]+DT*(25./216.*rhs1[i]+1408./2565.*rhs3[i]+2197./4104.*rhs4[i]-1./5.*rhs5[i]); 				// This is the correct integration step at order 4
		//vertex[i].phi = rhs0[i]+DT*(16./135.*rhs1[i]+6656./12825.*rhs3[i]+28561./56430.*rhs4[i]-9./50.*rhs5[i]+2./55.*rhs6[i]);	// This is the correct integration step at order 5
		Q[i] = (1/360.*rhs1[i]-128./4275.*rhs3[i]-2197./75240.*rhs4[i]+1/50.*rhs5[i]+2./55.*rhs6[i]); 					// This is their difference, which serves as an estimate to the (local) error
		if(Q[i]<0)Q[i]*=-1;
		Q_average+=vertex[i].area*Q[i]/TOTAL_AREA; 											// We use <abs(Q)>_\Sigma as a global error estimate 
		//rhs_average+=vertex[i].area*pow(rhs6[i]*rhs6[i],.5)/TOTAL_AREA;								// Estimate of the distance from equilibrium with arbitrary power. 
		rhs_average+=vertex[i].area*rhs6[i]*rhs6[i]/TOTAL_AREA;										// Estimate of the distance from equilibrium as <rhs^2>_\Sigma
	}

	if(Q_average>TOLERANCE){
		for (i=0; i<N_GRID_POINTS; i++){
			vertex[i].phi = rhs0[i];
		}
		CURRENT_TIME-=DT;
		
	}

	delta=pow(TOLERANCE*DT/Q_average/2.,.25);
	DT*=delta;

	if(rhs_average<TOLERANCE/10.)HALT_NOW=1;

	if(O_FLAG>=3){
		f_ou = fopen("Q_avg_debug.dat","a");
		fprintf(f_ou,"%.8E\t %.8E\t %.8E\t %.8E\t %.8E\n",CURRENT_TIME,Q_average,delta,DT,rhs_average);
		fclose(f_ou);
	}


}

/*******************************************************************/

void track_domains()
{
	int l, *labels, num_of_labels=1, position(int,int,int*);
	long i;
	
	void label_components(Vertex *,int *);	
	
	Vertex *pos, *neg;

	pos = malloc(N_GRID_POINTS*sizeof(Vertex));
	neg = malloc(N_GRID_POINTS*sizeof(Vertex));

	memcpy(pos,vertex,N_GRID_POINTS*sizeof(Vertex));
	memcpy(neg,vertex,N_GRID_POINTS*sizeof(Vertex));
	
	for (i=0; i<N_GRID_POINTS; i++){
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
	
	for (i=0; i<N_GRID_POINTS; i++){
		if (vertex[i].phi>0) vertex[i].label = pos[i].label;
		if (vertex[i].phi<0) vertex[i].label = neg[i].label;
	}
	
	labels = calloc(num_of_labels,sizeof(int));
	
	N_DOMAINS = 0;
	
	for (i=0; i<N_GRID_POINTS; i++){
		
		l = position(vertex[i].label,N_DOMAINS,labels);
				
		if (l==0){
			N_DOMAINS++;
			labels[N_DOMAINS] = vertex[i].label;
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
	
	parent = calloc(N_GRID_POINTS,sizeof(int));
		
	// Assign each vertex the same label of its first labeled neighbor or create a new one	
		
	for (i=0; i<N_GRID_POINTS; i++){
			
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
	
	for (i=0; i<N_GRID_POINTS; i++){
	
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

	for (i=0; i<N_GRID_POINTS; i++){
		
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

void end()
{	
	CPU_Time cpu_time;
	double c0=0,phisq=0,kin=0,pot=0,pot_p=0,tph,tpa,tpl,phiH2=0,phiKG=0;
	double kappa_avg=0,kappa_sq_avg=0;
	int i;
	
	FILE *f_ou;
	
	time(&T2);
  	get_time(T2-T1,&cpu_time);
	
	track_domains();

	if(O_FLAG==1 || O_FLAG==2){
		f_ou = fopen("last.m","w");
		export_graphic_complex(f_ou,1);
		fclose(f_ou);

		f_ou = fopen("interface.m","w");
		export_graphic_complex(f_ou,4);
		fclose(f_ou);
	}
	
	if(C_FLAG>=0){
		f_ou = fopen("last.dat","w");
		export_dat(f_ou);
		fclose(f_ou);
	};

 	
	for (i=0; i<N_GRID_POINTS; i++){

		tph=vertex[i].phi;
		tpa=vertex[i].area;
		tpl=laplace(i);

		phisq 	+=tpa*tph*tph;
		c0 	+=tpa*tph;
		kin	-=tpa*.5*tph*tpl*EPSILON;
		pot	+=tpa*V(i)/EPSILON;
		pot_p	+=tpa*Vp(i)/EPSILON;
		if(AVG_FLAG==0){
		phiH2	+=tpa*tph*vertex[i].h2;
		phiKG	+=tpa*tph*vertex[i].kg;
		}else
		{
		phiH2	+=tpa*tph*vertex[i].h2_avg;
		phiKG	+=tpa*tph*vertex[i].kg_avg;
		};
		kappa_avg+=tpa*(tpl-dV(i)/EPSILON/EPSILON)*sqrt(2.*Vp(i));
		kappa_sq_avg+=tpa*pow(tpl-dV(i)/EPSILON/EPSILON,2.)*EPSILON/2.;
	}

	kin*=3/sqrt(2.*K_BARRIER);
	pot*=3/sqrt(2.*K_BARRIER);
	pot_p*=3/sqrt(2.*K_BARRIER);
	phisq/=TOTAL_AREA;
	c0/=TOTAL_AREA;
	phiH2/=TOTAL_AREA;
	//phiKG/=TOTAL_AREA;
	phiKG/=2;

	printf("\n");
	printf("\tNumber of domains %d\n",N_DOMAINS);
	if(INT_FLAG<=3){printf("\tTime step %g\n",DT);};
	printf("\tCPU Time %ld:%ld:%ld:%ld\n",
	 cpu_time.d,cpu_time.h,
	 cpu_time.m,cpu_time.s);
  	printf("\n");
	
	f_ou = fopen("final.dat","w");	
	fprintf(f_ou,"%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%d\t%.10f\t%.10f\t%.10f\n",
	CURRENT_TIME,
	kin,
	pot_p,
	pot,
	phisq-c0*c0,
	phiH2,//-c0*WILLMORE_ENERGY/TOTAL_AREA,
	phiKG,//-c0*2*PI*EULER_CHI/TOTAL_AREA,
	LAGRANGE,
	N_DOMAINS,
	c0,
	kappa_avg,
	kappa_sq_avg);
	fclose(f_ou);

}

/*******************************************************************/

void export_histogram(FILE *f_ou, long t)
{
	long m, i, j;
	double V(long);

	double c0=0,phisq=0,kin=0,pot=0,pot_p=0,tph,tpa,tpl,phiH2=0,phiKG=0;
	double kappa_avg=0,kappa_sq_avg=0;
 	
	for (i=0; i<N_GRID_POINTS; i++){

		tph=vertex[i].phi;
		tpa=vertex[i].area;
		tpl=laplace(i);

		phisq 	+=tpa*tph*tph;
		c0 	+=tpa*tph;
		kin	-=tpa*.5*tph*tpl*EPSILON;
		pot	+=tpa*V(i)/EPSILON;
		pot_p	+=tpa*Vp(i)/EPSILON;
		if(AVG_FLAG==0){
		phiH2	+=tpa*tph*vertex[i].h2;
		phiKG	+=tpa*tph*vertex[i].kg;
		}else
		{
		phiH2	+=tpa*tph*vertex[i].h2_avg;
		phiKG	+=tpa*tph*vertex[i].kg_avg;
		};
		kappa_avg+=tpa*(tpl-dV(i)/EPSILON/EPSILON)*sqrt(2.*Vp(i));
		kappa_sq_avg+=tpa*pow(tpl-dV(i)/EPSILON/EPSILON,2.)*EPSILON/2.;
	}

	kin*=3/sqrt(2.*K_BARRIER);
	pot*=3/sqrt(2.*K_BARRIER);
	pot_p*=3/sqrt(2.*K_BARRIER);
	phisq/=TOTAL_AREA;
	c0/=TOTAL_AREA;
	phiH2/=TOTAL_AREA;
	//phiKG/=TOTAL_AREA;
	phiKG/=2;
	
	fprintf(f_ou,"%.10ld\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%d\t%.10f\t%.10f\t%.10f\n",
	t,
	DT,
	CURRENT_TIME,
	kin,
	pot_p,
	pot,
	phisq-c0*c0,
	phiH2,//-c0*WILLMORE_ENERGY/TOTAL_AREA,
	phiKG,//-c0*2*PI*EULER_CHI/TOTAL_AREA,
	LAGRANGE,
	N_DOMAINS,
	c0,
	kappa_avg,
	kappa_sq_avg
	);
	
}


/*******************************************************************/

void export_graphic_complex(FILE *f_ou, long l)
{
	double r, g, b;
	long i, j;
	
	fprintf(f_ou,"GraphicsComplex[{");

	for (i=0; i<N_GRID_POINTS; i++){
		fprintf(f_ou,"{%.10f,%.10f,%.10f},",
		vertex[i].x,
		vertex[i].y,
		vertex[i].z);
	}

	fseek(f_ou,-1, SEEK_CUR);

	fprintf(f_ou,"},Polygon[{");
	
	for (i=0; i<N_GRID_TRIANGLES; i++){		
		
		fprintf(f_ou,"{%ld,%ld,%ld},",
		triangle[i].v1+1,
		triangle[i].v2+1,
		triangle[i].v3+1);
	}
	
	fseek(f_ou,-1, SEEK_CUR);
	
	fprintf(f_ou,"}],VertexColors->{");
	
	for (i=0; i<N_GRID_POINTS; i++){
		switch(l){
			case 1: r = 0.5*(1+atan(vertex[i].phi));
				g = 0.5*(1-atan(vertex[i].phi));
				b = 0;
				break;
			case 2: r = (vertex[i].h2-H2_MIN)/(H2_MAX-H2_MIN);
				g = 0*(vertex[i].h2-H2_MIN)/(H2_MAX-H2_MIN);
				b = 1-r-g;
				break;
			case 3: r = 0*(vertex[i].kg-KG_MIN)/(KG_MAX-KG_MIN);
				b = (vertex[i].kg-KG_MIN)/(KG_MAX-KG_MIN);
				g = 1-r-b;
				break;
			case 4: r = (vertex[i].phi<0.5 && vertex[i].phi>-0.5)?0:1;
				b = (vertex[i].phi<0.5 && vertex[i].phi>-0.5)?0:1;
				g = (vertex[i].phi<0.5 && vertex[i].phi>-0.5)?0:1;
				break;
			 }

		fprintf(f_ou,"RGBColor[{%lg,%lg,%lg}],",r,g,b);
	}

	fseek(f_ou,-1, SEEK_CUR);
	
	fprintf(f_ou,"}]");
}

/*******************************************************************/

void export_dat(FILE *f_ou)
{
	long i,bphi;
	
	for (i=0; i<N_GRID_POINTS; i++){
		(vertex[i].phi>0)?(bphi=1):(bphi=0);
		fprintf(f_ou,"%.12e\t%ld\n",
		vertex[i].phi,
		bphi);
	}
}
