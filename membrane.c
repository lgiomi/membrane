/******************************************************************
**                                                               **
**                           MEMBRANE                            **
**                                                               **
**			      V 2.3                              **
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
long export;

int num_of_domains;
double l_min,l_max,h2_min,h2_max,kg_min,kg_max;
double c_x,c_y,c_z;
int c_flag=0,o_flag=1,i_flag=0;
double current_time;

double total_area;
double willmore_energy,euler_chi;
double sigma=1,c_0,epsilon,DTauto;
double DT,run_time;
double gamma_h,gamma_h2,gamma_kg;
double lagrange;
double current_energy,previous_energy;
double energy_error_max=.1,energy_error_min=.001;

long seed;

int halt_now=0;

void euler();
void rk2();
void rk4();
void rkf45();
void heun_euler();

void run();
void end();
void init();

void one_step();
void init_random();
void import_initial(char *);
void get_geometry();
void track_domains();
void find_neighbors();
void get_rhs(double *);
void progress_bar(long);
void import_mesh(char *);
void write_hi(FILE *, long);
void get_time(time_t, CPU_Time *);
void export_graphic_complex(FILE *, long);
void export_dat(FILE *);
void help();
void print_cmd_line();

/*******************************************************************/

main(int argc, char *argv[])
{	
	init(argc, argv);
	time(&t1);
	run();
	end();
	
  	return EXIT_SUCCESS;	
}

/*******************************************************************/

void init(int argc, char *argv[])
{
	double box_size;
	char f_name[32];
	char import_name[32];
	int n;                           
   
	if(argc>1){
	run_time=0;
	epsilon=0;
	num_of_iteration=1;
	export=0;
	gamma_h=0;
	gamma_h2=0;
	gamma_kg=0;
	method=0;
	seed=0;
	
	for( n = 1; n < argc; n++ ){         /* Scan through args. */
		switch((int)argv[n][0]){     /* Check for option character "-" */
			case '-':	                  
					switch((int)argv[n][1]){
						case 'm':
								snprintf(f_name, sizeof(f_name), "%s", argv[n+1]);
								printf( "Mesh file: %s\n",f_name);
								n++;
								break;
						case 'h':
								help();
								break;
						case 'O':
								o_flag=atol(argv[n+1]);
								printf( "Output level flag: %d\n",o_flag);
								n++;
								break;
						case 't':
								run_time=atof(argv[n+1]);
								printf( "Run time: %lg\n",run_time);
								n++;
								break;
						case 'e':
								epsilon=atof(argv[n+1]);
								printf( "Epsilon : %lg\n",epsilon);
								n++;
								break;
						case 'i':
								num_of_iteration=atol(argv[n+1]);
								if(num_of_iteration==-1){printf("Will compute time step automatically\n");}
								else{printf( "Total number of iterations: %ld\n",num_of_iteration);};
								n++;
								break;
						case 'x':
								export=atol(argv[n+1]);
								printf( "Export field configuration every %ld steps\n",export);
								n++;
								break;
						case 'C':
								gamma_h=atof(argv[n+1]);
								gamma_h2=atof(argv[n+2]);
								gamma_kg=atof(argv[n+3]);
								printf( "Couplings: (%lg H, %lg H^2, %lg K_G )\n",gamma_h,gamma_h2,gamma_kg);
								n+=3;
								break;
						case 'P':
								c_x=atof(argv[n+1]);
								c_y=atof(argv[n+2]);
								c_z=atof(argv[n+3]);
								printf( "Center of stereographic projection (%lg,%lg,%lg)\n",c_x,c_y,c_z);
								c_flag=1;
								n+=3;
								break;
						case 'I':
								switch((int)argv[n+1][0]){
									case '1': method=1;break;
									case '2': method=2;break;
									case '3': method=3;break;
									case '4': method=4;break;
									case '5': method=5;break;
									default:printf( "Illegal integration method  %c\n",(int)argv[n+1][0]);
										printf( "\nType ./membrane -h for help\n"); 
										exit(0);
										break;
								}
								printf( "Integration method: %ld\n",method);		
								n++;
								break;
						case 'r':
								seed=atol(argv[n+1]);
								c_0=atof(argv[n+2]);
								if(c_0>=0&&c_0<=1){
									printf( "Random initial data with seed %ld and total concentration %lg\n",seed,c_0);
									i_flag+=2;
									n+=2;
								}else{
									printf( "Error, total relative concentration has to be between 0 and 1 (you passed %lg) \n",c_0);
									printf( "\nType ./membrane -h for help\n"); 
									exit(0);
								};
								break;
						case 'R':
								snprintf(import_name, sizeof(import_name), "%s", argv[n+1]);
								printf( "Import initial configuration from file: %s\n",import_name);
								i_flag+=1;
								n++;
								break;             
						default:  
							printf( "Illegal option code \"-%c\"\n",(int)argv[n][1]);
							printf( "\nType ./membrane -h for help\n"); 
							exit(0);
							break;
					}
                 			break;
			default:  
				printf( "\nError: give input in the following format:\n"); 
				print_cmd_line();
				printf( "\nType ./membrane -h for help\n"); 
				exit(0);
                 		break;
			}	
		}
		DT = run_time/num_of_iteration;
		if( run_time==0 ||  method == 0 || (num_of_iteration == 1 && method < 4) || epsilon == 0 || i_flag <1){
				printf( "\nError: not enough arguments given or wrong parameter values. Write input in the following format:\n"); 
				print_cmd_line();
				printf( "Type ./membrane -h for help\n"); 
				exit(0);
		}
		if(export == 0){
			printf( "No -x option given, assuming no intermediate output.\n"); 
		}
	}
	else{
	printf("You did not specify arguments to parse, assume input is from stdin\n");
	printf("Input mesh file ");	if(scanf("%s",f_name)==1){};
	printf("Input run time ");	if(scanf("%lg",&run_time)==1){};
  	printf("Input iterations ");	if(scanf("%ld",&num_of_iteration)==1){};
	DT = run_time/num_of_iteration;
	printf("Choose integration method\n\n\t1. Euler\n\t2. Runge-Kutta (2 steps)\n\t3. Runge-Kutta (4 steps)\n\t4. Adaptive Heun-Euler\n\t5. RKF45\n\n"); if(scanf("%ld",&method)==1){};
	printf("Input sigma "); if(scanf("%lg",&sigma)==1){};
	printf("Input epsilon "); if(scanf("%lg",&epsilon)==1){};
	printf("Input area percentage "); if(scanf("%lg",&c_0)){};
	printf("Input random seed "); if(scanf("%ld",&seed)){};
	i_flag=2;
	}
	import_mesh(f_name);
	get_geometry();
	if(i_flag==1){
		import_initial(import_name);
	}
	else if(i_flag==2){
		init_random();
	};
}

/*******************************************************************/

void print_cmd_line()
{
	printf("./membrane -m MESH_FILE -t RUN_TIME -I METHOD -e EPSILON  (-r SEED MEAN_CONCENTRATION | -R START_FILE ) [-O LEVEL] [-x STEPS] [-i TOTAL_ITERATIONS] [-C GAMMA_H GAMMA_H^2 GAMMA_KG] [-P CX CY CZ]\n\n");
}

/*******************************************************************/

void help()
{
	printf("\nThis program takes input values as in-line commands or from stdin (deprecated). In the former case the syntax is\n\n");	
	print_cmd_line();
	printf("Where:\n");
	printf("\t -m MESH_FILE\t: specifies the .msh file (works only with standard gmsh mesh format)\n");
	printf("\t -t RUN_TIME\t: specifies the upper bound on the total time of the simulation\n");
	printf("\t -I METHOD\t: specifies the integration method\n\t\t\t\t1: Euler\n\t\t\t\t2: RK2\n\t\t\t\t3: RK4\n\t\t\t\t4: RK2-Euler with adaptive step-size [NOT FULLY WORKING]\n\t\t\t\t5: RKF45 with adaptive stepsize\n");
	printf("\t -e EPSILON\t: specifies the value of epsilon (i.e. diffusion constant and interface thickness)\n");
	printf("\t -r #1 #2\t: specifies the seed #1 for the random intial condition and the desired mean concentration value #2\n");
	printf("\t -R FILE\t: specifies where to import from the initial configuration of the phase field on the mesh (does not work with -r)\n");
	printf("\t -O LEVEL\t: choose which output files will be printed\n\t\t\t\t0: 'histo.dat', 'last.dat', 'final.dat'\n\t\t\t\t1: previous + 'geometry.dat','last.m' + 'gc_#.dat' if -x is set [DEFAULT]\n\t\t\t\t2: previous + 'mean_curvature.m','gaussian_curvature.m' + 'gc_#.m' if -x is set\n\t\t\t\t3: as in '1' + debug files\n");
	printf("\t -x STEPS\t: if specified and #!=0, decides the frequency with which to export field configurations \n");
	printf("\t -i ITERATIONS\t: specifies the total number of iterations (-I 4 and -I 5 do not use this parameter)\n");
	printf("\t -C #1 #2 #3\t: specifies the values of the three cubic couplings with H, H^2 and K_G\n");
	printf("\t -P #1 #2 #3\t: specifies the coordinates of the center for two-dimensional stereographic projection of the surface\n");

	printf("\n");
	exit(1);
}

/*******************************************************************/

void import_mesh(char *f_name)
{
	long i, v1, v2, v3, chi, triangle_test, dummy, num_of_edges=0;

	if( access( f_name, F_OK ) == -1 ) {
		printf("\nError: file %s does not exist\n",f_name);
		exit(0);
	};

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
	
	for (i=0; i<8; i++){
	if(fgets(line,LINESIZE-1,f_in)){};
	}
	
	if(fscanf(f_in,"%ld",&num_of_meshpoint)){};

	if (num_of_meshpoint>MAX_SIZE){
		printf("Error: the number of vertices (%ld) exceeds MAX_SIZE (%d)\n",num_of_meshpoint,MAX_SIZE);
		exit(0);
	}
		
	for (i=0; i<num_of_meshpoint; i++){
		if(fscanf(f_in,"%ld%lg%lg%lg",
		&dummy,
		&vertex[i].x,
		&vertex[i].y,
		&vertex[i].z)){}else{printf("Failed to read mesh point.");};
	}
	
	for (i=0; i<3; i++){
	if(fgets(line,LINESIZE-1,f_in)){};
	}
	
	if(fscanf(f_in,"%ld",&num_of_triangles)){};
		
	if (num_of_triangles>MAX_SIZE){
		printf("Error: the number of triangles (%ld) exceeds MAX_SIZE (%d)\n",num_of_triangles,MAX_SIZE);
		exit(0);
	}	
		
	for (i=0; i<num_of_triangles; i++){
		if(fscanf(f_in,"%ld%ld%ld%ld%ld%ld%ld%ld",
		&dummy,
		&triangle_test,
		&dummy,
		&dummy,
		&dummy,
		&v1,
		&v2,
		&v3)){};
		if(triangle_test!=2){break;}; //Only elements with the second column equal to 2 are triangles: therefore, ignore anything else.
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
		printf("Error: bad triangulation or not a g=0 surface, chi = %ld\n",chi);
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
	double theta_this, theta_next, area_prev, area_next, theta_a, theta_b;


	l_min=1E10;
	l_max=0;

	h2_min=1E10;
	h2_max=-1E10;

	kg_min=1E10;
	kg_max=-1E10;

	long i, j, tmp, next, prev, this,  obtuse[MAX_SIZE], num_of_obtuse=0;
	long which_triangle(long,long,long);
	
	int swap(long *,long *), swapped;
	
	total_area = 0;	
	willmore_energy = 0;	
	euler_chi= 0;	
			
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
	
		// Sort neighbors on the tangent plane (orientation is arbitrary)
				
		do {
			
			swapped = 0;
			
			for (j=0; j<vertex[i].num_of_neighbors-1; j++){
				
				dx = (vertex[i].x-vertex[vertex[i].neighbor[j]].x)*x[0] 
					+(vertex[i].y-vertex[vertex[i].neighbor[j]].y)*x[1]
					+(vertex[i].z-vertex[vertex[i].neighbor[j]].z)*x[2];
					
				dy = (vertex[i].x-vertex[vertex[i].neighbor[j]].x)*y[0] 
					+(vertex[i].y-vertex[vertex[i].neighbor[j]].y)*y[1]
					+(vertex[i].z-vertex[vertex[i].neighbor[j]].z)*y[2];						
				
				theta_this = atan2(-dy,-dx);

				
				dx = (vertex[i].x-vertex[vertex[i].neighbor[j+1]].x)*x[0] 
					+(vertex[i].y-vertex[vertex[i].neighbor[j+1]].y)*x[1]
					+(vertex[i].z-vertex[vertex[i].neighbor[j+1]].z)*x[2];

				dy = (vertex[i].x-vertex[vertex[i].neighbor[j+1]].x)*y[0] 
					+(vertex[i].y-vertex[vertex[i].neighbor[j+1]].y)*y[1]
					+(vertex[i].z-vertex[vertex[i].neighbor[j+1]].z)*y[2];

				theta_next = atan2(-dy,-dx);

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
		vertex[i].kg=2*PI;
	
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

			// Calculate angle at vertex i in triangle_prev

			x[0] = vertex[i].x-vertex[prev].x;
			x[1] = vertex[i].y-vertex[prev].y;
			x[2] = vertex[i].z-vertex[prev].z;

			norm = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
				
			x[0] /= norm;
			x[1] /= norm;
			x[2] /= norm;

			y[0] = vertex[i].x-vertex[this].x;
			y[1] = vertex[i].y-vertex[this].y;
			y[2] = vertex[i].z-vertex[this].z;

			norm = sqrt(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]);

			y[0] /= norm;
			y[1] /= norm;
			y[2] /= norm;

			theta_a = atan2(sqrt(1.-pow(x[0]*y[0]+x[1]*y[1]+x[2]*y[2],2.)),x[0]*y[0]+x[1]*y[1]+x[2]*y[2]);

			// Calculate angle at vertex i in triangle_next

			x[0] = vertex[i].x-vertex[next].x;
			x[1] = vertex[i].y-vertex[next].y;
			x[2] = vertex[i].z-vertex[next].z;

			norm = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
				
			x[0] /= norm;
			x[1] /= norm;
			x[2] /= norm;

			y[0] = vertex[i].x-vertex[this].x;
			y[1] = vertex[i].y-vertex[this].y;
			y[2] = vertex[i].z-vertex[this].z;

			norm = sqrt(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]);

			y[0] /= norm;
			y[1] /= norm;
			y[2] /= norm;

			theta_b = atan2(sqrt(1.-pow(x[0]*y[0]+x[1]*y[1]+x[2]*y[2],2.)),x[0]*y[0]+x[1]*y[1]+x[2]*y[2]);

			// Subtract the computed angles to 2Pi in K_G

			vertex[i].kg-=(theta_a+theta_b)/2;

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

			vertex[i].area += ( cota>=0 && theta_a <= PI/2 && PI-atan(1/cota)-theta_a <= PI/2 ? (dx*dx+dy*dy+dz*dz)*cota/8 : (theta_a >= PI/2 ? area_prev/4 : area_prev/8));
			vertex[i].area += ( cotb>=0 && theta_b <= PI/2 && PI-atan(1/cotb)-theta_b <= PI/2 ? (dx*dx+dy*dy+dz*dz)*cotb/8 : (theta_b >= PI/2 ? area_next/4 : area_next/8));
		}
		
		// Normalized the Laplacian weights 
		
		for (j=0; j<vertex[i].num_of_neighbors; j++){
			
			vertex[i].weight[j] /= vertex[i].area;	
		}
		
		// Normalize the mean curvature vector
				
		vertex[i].hx /= vertex[i].area;
		vertex[i].hy /= vertex[i].area;
		vertex[i].hz /= vertex[i].area;

		vertex[i].kg /= vertex[i].area;

		// Compute projective angles

		if(c_flag==1){
			dx = vertex[i].x-c_x;
			dy = vertex[i].y-c_y;
			dz = vertex[i].z-c_z;
			norm=sqrt(dx*dx+dy*dy+dz*dz);
			vertex[i].gridx=atan2(dy,dx);
			vertex[i].gridy=atan2(dz,sqrt(dx*dx+dy*dy));
		}

		// Calculate the total squared curvature 

		vertex[i].h2 = vertex[i].hx*vertex[i].hx+vertex[i].hy*vertex[i].hy+vertex[i].hz*vertex[i].hz;
		
		total_area += vertex[i].area;
		willmore_energy += vertex[i].area*vertex[i].h2;
		euler_chi += .5/PI*vertex[i].area*vertex[i].kg;
	}
	
	FILE *f_ou;
	
	if(o_flag>=1){
		f_ou = fopen("geometry.dat","w");

		for (i=0; i<num_of_meshpoint; i++){

			if(c_flag==0){
				fprintf(f_ou,"%ld\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10lg\t",
					i,
					vertex[i].x,
					vertex[i].y,
					vertex[i].z,
					vertex[i].area,
					vertex[i].h2,
					vertex[i].kg);
			}else if(c_flag==1){
				fprintf(f_ou,"%ld\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10lg\t%.10lg\t%.10lg\t",
					i,
					vertex[i].x,
					vertex[i].y,
					vertex[i].z,
					vertex[i].area,
					vertex[i].h2,
					vertex[i].kg,
					vertex[i].gridx,
					vertex[i].gridy);
			};
				
			for (j=0; j<vertex[i].num_of_neighbors; j++){
				fprintf(f_ou,"%.10f\t",vertex[i].weight[j]);
	
				this = vertex[i].neighbor[j];
			
				x[0] = vertex[i].x-vertex[this].x;
				x[1] = vertex[i].y-vertex[this].y;
				x[2] = vertex[i].z-vertex[this].z;
				base=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
	
				if(l_max<base){l_max=base;};
				if(l_min>base){l_min=base;};
	
				if(h2_max<vertex[i].h2){h2_max=vertex[i].h2;};
				if(h2_min>vertex[i].h2){h2_min=vertex[i].h2;};
	
				if(kg_max<vertex[i].kg){kg_max=vertex[i].kg;};
				if(kg_min>vertex[i].kg){kg_min=vertex[i].kg;};
				
			}	
			
			fprintf(f_ou,"\n");
	}	
	
	fclose(f_ou);
	}

	if(o_flag==2){
		f_ou = fopen("mean_curvature.m","w");
		export_graphic_complex(f_ou,2);
		fclose(f_ou);

		f_ou = fopen("gaussian_curvature.m","w");
		export_graphic_complex(f_ou,3);
		fclose(f_ou);
	}

	/*
	f_ou = fopen("obtuse_debug.dat","w");
	
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
	printf("\tObtuse triangles %ld (%.1ld%% of total)\n",num_of_obtuse/2,100*num_of_obtuse/2/num_of_triangles);
	printf("\tTotal surface area %lg\n",total_area);
	printf("\tWillmore energy %lg (asphericity %lg)\n",willmore_energy,willmore_energy/4/PI-1);
	printf("\tEuler characteristic %lg\n",euler_chi);
	printf("\tMin/max length (%lg,%lg)\n",l_min,l_max);
	printf("\tMin/max H2 (%lg,%lg)\n",h2_min,h2_max);
	printf("\tMin/max KG (%lg,%lg)\n",kg_min,kg_max);
	printf("\tProposed initial time-step %lg\n",.1/sigma*l_min*l_min/2);// For diffusion processes on flat space, explicit discretization schemes (as ours) are stable only for sigma*DT/DX^2 < 1/2.
	printf("\n");

	/*f_ou = fopen("negative_kg_debug.m","w");

	fprintf(f_ou,"{");
	for (i=0; i<num_of_meshpoint; i++){
		if(vertex[i].kg<-1){
			fprintf(f_ou,"{{");
			for (j=0; j<vertex[i].num_of_neighbors; j++){
				if(j >0)fprintf(f_ou,",");
				fprintf(f_ou,"Line[{{%2.3f10,%2.3f,%2.3f},{%2.3f,%2.3f,%2.3f}}]",
					vertex[i].x,
					vertex[i].y,
					vertex[i].z,
					vertex[vertex[i].neighbor[j]].x,
					vertex[vertex[i].neighbor[j]].y,
					vertex[vertex[i].neighbor[j]].z
				);
			}
			fprintf(f_ou,"},");
			fprintf(f_ou,"%2.5f,%2.5f,%2.5f},",vertex[i].area,vertex[i].h2*vertex[i].area,vertex[i].kg*vertex[i].area);
		};
	}
	fseek(f_ou,-1, SEEK_CUR);
	fprintf(f_ou,"}\n");

	fclose(f_ou);*/

	//exit(1);
	

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
	double  noise=0.1, ran2(long *);
	long i;
		
	for (i=0; i<num_of_meshpoint; i++){
		vertex[i].phi = (2*c_0-1)+noise*(1-2*ran2(&seed));
	}
}

/*******************************************************************/

void import_initial(char *import_name)
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

	if(lines==num_of_meshpoint){
		c_0=0;
		FILE *f_in = fopen(import_name,"r");
		for (i=0; i<num_of_meshpoint; i++){
		if(fscanf(f_in,"%lg%*[^\n]",&vertex[i].phi)){
			c_0+=vertex[i].area*vertex[i].phi;
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
	c_0=.5+c_0/2.;
	printf("Imported initial configuration from %s with mean concentration %g\n",import_name,c_0);
}

/*******************************************************************/

void get_rhs(double *rhs)
{
	double laplace(long);
	double dV(long);
	long i; 	
	
	// Calculate the Lagrange multiplier
	
	lagrange = 0;

	for (i=0; i<num_of_meshpoint; i++){

		lagrange += vertex[i].area*(sigma*laplace(i)-dV(i)/(epsilon*epsilon));		
	}		
	
	lagrange /= total_area;
	
	// Calculate the right-hand side of the Allen-Cahn equation
	
	for (i=0; i<num_of_meshpoint; i++){
		
		rhs[i] = sigma*laplace(i)-dV(i)/(epsilon*epsilon)-lagrange;
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

// Define the potential and its functional derivative [CHECK THESE ARE CORRECT]

double dV(long i)
{
	double dV_0,dV_int;

	dV_0 = vertex[i].phi*(vertex[i].phi*vertex[i].phi-1);

	dV_int = 3./4.*(1-vertex[i].phi*vertex[i].phi)*epsilon;
	
	return dV_0+dV_int*(gamma_h2*vertex[i].h2+gamma_kg*vertex[i].kg+gamma_h*sqrt(vertex[i].h2));
}

double V(long i)
{
	double V_0,V_int;

	V_0 = .25*pow(pow(vertex[i].phi,2.)-1,2.);

	V_int = epsilon*pow((vertex[i].phi+1),2.)*(2-vertex[i].phi)/4.;
	
	return V_0+V_int*(gamma_h2*vertex[i].h2+gamma_kg*vertex[i].kg+gamma_h*sqrt(vertex[i].h2));
}

/*******************************************************************/

void run()
{
	long t=0;
	char f_math_na[32],f_dat_na[32];
	
	FILE *f_hi, *f_ou;
	
	f_hi = fopen("histo.dat","w");	

	current_time=0;

	DTauto=.1/sigma*l_min*l_min/2;

	if(method>=4){
		DT=DTauto;
		printf("Setting initial time step automatically to %f\n",DT);
	}
	
	if(export>0 && c_flag==0 && o_flag!=2){
		printf("You chose to export configurations every %ld time steps, but did not specify a sterographic center nor chose -O 2. Assuming -C 0 0 0.\n",export);
		c_flag=1;
		c_x=0;
		c_y=0;
		c_z=0;
	}

	if(num_of_iteration==-1 || method == 4 || method == 5){
		while(current_time<run_time && halt_now==0){

			current_time+=DT;

			write_hi(f_hi,t);

			if (export>0 && t%export==0){

				if(o_flag == 2 && c_flag==1){
					sprintf(f_math_na,"gc_%08ld.m",t);
					f_ou = fopen(f_math_na,"w");
					export_graphic_complex(f_ou,1);
					fclose(f_ou);
				}
				if(c_flag==1 && o_flag>=1){
					sprintf(f_dat_na,"gc_%08ld.dat",t);
					f_ou = fopen(f_dat_na,"w");
					export_dat(f_ou);
					fclose(f_ou);
				};
			}
			t++;
			if(t%100==0){printf("Step: %ld\t\t\t Current time:%.6g\r",t,current_time);}
			one_step();
		}
		printf("Total time steps: %ld\nSimulation time: %g\n",t,current_time);
	}
	else{
		for (t=0; t<num_of_iteration; t++){

			current_time+=DT;

			write_hi(f_hi,t);

			if (export>0 && t%export==0){
				if(o_flag == 2 && c_flag==1){
					sprintf(f_math_na,"gc_%08ld.m",t);
					f_ou = fopen(f_math_na,"w");
					export_graphic_complex(f_ou,1);
					fclose(f_ou);
				}
				if(c_flag==1 && o_flag>=1){
					sprintf(f_dat_na,"gc_%08ld.dat",t);
					f_ou = fopen(f_dat_na,"w");
					export_dat(f_ou);
					fclose(f_ou);
				};
			}
		
			//track_domains();
			progress_bar(t);
			one_step();
		}
	};
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
	
	for (i=0; i<num_of_meshpoint; i++){
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
	
	for (i=0; i<num_of_meshpoint; i++){
		rhs0[i] = vertex[i].phi; 
	}
	
	get_rhs(rhs1);
	
	for (i=0; i<num_of_meshpoint; i++){
		vertex[i].phi = rhs0[i]+DT*rhs1[i];
	}
	
	get_rhs(rhs2);
	
	for (i=0; i<num_of_meshpoint; i++){
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

void heun_euler()
{
	double rhs0[MAX_SIZE];
	double rhs1[MAX_SIZE];
	double rhs2[MAX_SIZE];
	double Q[MAX_SIZE];
	double tol=1E-6,Q_average=0;
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
		vertex[i].phi = rhs0[i]+0.5*DT*(rhs1[i]+rhs2[i]);
		Q[i] = (rhs1[i]-rhs2[i])/2.;
		if(Q[i]<0)Q[i]*=-1;
		Q_average+=vertex[i].area*Q[i]/total_area;
	}

	if(Q_average<tol)halt_now=1;


	if(o_flag>=3){
		FILE *f_ou;
		f_ou = fopen("Q_avg_debug.dat","a");
		fprintf(f_ou,"%.8f\t %2.8f\t\n",current_time,Q_average);
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
	double tol=1E-4,Q_average=0,delta;
	double DTmin=DTauto,DTmax=DTauto*100.;
	long i;

	FILE *f_ou;

	for (i=0; i<num_of_meshpoint; i++){
		rhs0[i] = vertex[i].phi; 
	}
	
	get_rhs(rhs1);
	
	for (i=0; i<num_of_meshpoint; i++){
		vertex[i].phi = rhs0[i]+1./4.*DT*rhs1[i];
	}
	
	get_rhs(rhs2);
	
	for (i=0; i<num_of_meshpoint; i++){
		vertex[i].phi = rhs0[i]+1./32.*DT*(3.*rhs1[i]+9.*rhs2[i]);
	}
	
	get_rhs(rhs3);
	
	for (i=0; i<num_of_meshpoint; i++){
		vertex[i].phi = rhs0[i]+1./2197.*DT*(1932.*rhs1[i]-7200.*rhs2[i]+7296.*rhs3[i]);
	}
	
	get_rhs(rhs4);
	
	for (i=0; i<num_of_meshpoint; i++){
		vertex[i].phi = rhs0[i]+DT*(439./216.*rhs1[i]-8.*rhs2[i]+3680./513.*rhs3[i]-845./4104.*rhs4[i]);
	}

	get_rhs(rhs5);

	for (i=0; i<num_of_meshpoint; i++){
		vertex[i].phi = rhs0[i]+DT*(-8./27.*rhs1[i]+2.*rhs2[i]-3544./2565.*rhs3[i]+1859./4104.*rhs4[i]-11./40.*rhs5[i]);
	}

	get_rhs(rhs6);

	for (i=0; i<num_of_meshpoint; i++){
		vertex[i].phi = rhs0[i]+DT*(25./216.*rhs1[i]+1408./2565.*rhs3[i]+2197./4104.*rhs4[i]-1./5.*rhs5[i]);
		Q[i] =  (1/360.*rhs1[i]-128./4275.*rhs3[i]+2197./75240.*rhs4[i]+1/50.*rhs5[i]+2./55.*rhs6[i]);
		if(Q[i]<0)Q[i]*=-1;
		Q_average+=vertex[i].area*Q[i]/total_area;
	}

	if(Q_average<tol/100.)halt_now=1;

	delta=pow(tol/Q_average/2.,.25);

	if(current_time>epsilon){
		if(delta<.9){DT+=-.001*DTauto;}
		else if(delta>1.){DT+=.0001*DTauto;}
		else {DT*=delta;};
	}

	if(DT<DTmin)DT=DTmin;
	if(DT>DTmax)DT=DTmax;

	if(o_flag>=3){
		f_ou = fopen("Q_avg_debug.dat","a");
		fprintf(f_ou,"%.8f\t %2.8f\t %2.8f\t %2.8f\n",current_time,Q_average,delta,DT);
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

	printf("\n");
	printf("\tNumber of domains %d\n",num_of_domains);
	if(method<=3){printf("\tTime step %g\n",DT);};
	printf("\tCPU Time %ld:%ld:%ld:%ld\n",
	 cpu_time.d,cpu_time.h,
	 cpu_time.m,cpu_time.s);
  	printf("\n");
	
	if(o_flag>=1){
		f_ou = fopen("last.m","w");
		export_graphic_complex(f_ou,1);
		fclose(f_ou);
	}
	
	if(c_flag>=0){
		f_ou = fopen("last.dat","w");
		export_dat(f_ou);
		fclose(f_ou);
	};
	
	f_ou = fopen("final.dat","w");	
	fprintf(f_ou,"Number of domains %d\n",num_of_domains);
	if(method<=3){fprintf(f_ou,"Time step %g\n",DT);};
	fprintf(f_ou,"CPU Time %ld:%ld:%ld:%ld\n",
	 cpu_time.d,cpu_time.h,
	 cpu_time.m,cpu_time.s);
	fclose(f_ou);

}

/*******************************************************************/

void write_hi(FILE *f_ou, long t)
{
	long m, i, j;
	double V(long);

	double c0=0,phisq=0,kin=0,pot=0,tph,tpa,phiH2=0,phiKG=0;
 	
	for (i=0; i<num_of_meshpoint; i++){

		tph=vertex[i].phi;
		tpa=vertex[i].area;

		phisq 	+=tpa*tph*tph;
		c0 	+=tpa*tph;
		kin	-=tpa*sigma*.5*tph*laplace(i);
		pot	+=tpa*V(i)/(epsilon*epsilon);
		phiH2	+=tpa*tph*vertex[i].h2;
		phiKG	+=tpa*tph*vertex[i].kg;
	}

	current_energy=kin+pot;
	kin/=total_area;
	pot/=total_area;
	phisq/=total_area;
	c0/=total_area;
	phiH2/=total_area;
	phiKG/=total_area;
	
	fprintf(f_ou,"%.10ld\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%d\n",
	t,
	DT,
	current_time,
	kin,
	pot,
	kin+pot,
	phisq-c0*c0,
	phiH2-c0*willmore_energy/total_area,
	phiKG-c0*2*PI*euler_chi/total_area,
	lagrange,
	num_of_domains);
	
}


/*******************************************************************/

void export_graphic_complex(FILE *f_ou, long l)
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
		switch(l){
			case 1: r = 0.5*(1+atan(vertex[i].phi));
				g = 0.5*(1-atan(vertex[i].phi));
				b = 0;
				break;
			case 2: r = (vertex[i].h2-h2_min)/(h2_max-h2_min);
				g = 0*(vertex[i].h2-h2_min)/(h2_max-h2_min);
				b = 1-r-g;
				break;
			case 3: r = 0*(vertex[i].kg-kg_min)/(kg_max-kg_min);
				b = (vertex[i].kg-kg_min)/(kg_max-kg_min);
				g = 1-r-b;
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
	
	for (i=0; i<num_of_meshpoint; i++){
		(vertex[i].phi>0)?(bphi=1):(bphi=0);
		fprintf(f_ou,"%.10f\t%ld\n",
		vertex[i].phi,
		bphi);
	}
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
