/******************************************************************
**                                                               **
**                           MEMBRANE                            **
**                                                               **
**			      V 3.0                              **
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

long num_of_meshpoint,num_of_triangles,num_of_iteration,num_of_edges=0;

double c_x,c_y,c_z;
double a_x,a_y,a_z;
long method;
long export;
int c_flag=0,o_flag=1,i_flag=0;
double a_flag=0,v_flag=0;
double tol;
double DT,DTauto,run_time;
double c_0,epsilon;
double gamma_h,gamma_h2,gamma_kg;
double k_barrier=1.;
double conserved=1;
double scale=1;

double total_area,willmore_energy,euler_chi,total_volume;
double l_min,l_max,l_avg,h2_min,h2_max,kg_min,kg_max;

double lagrange;
double current_time;
int num_of_domains;

long seed;
int halt_now=0;

void euler();
void rk2();
void rk4();
void heun_euler();
void rkf45();

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
	double box_size,norm;
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
	tol=0;
	a_x=0;
	a_y=0;
	a_z=1;

	
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
								o_flag=atol(argv[n+1]);
								printf("Output level\t\t: %d\n",o_flag);
								n++;
								break;
						case 't':
								run_time=atof(argv[n+1]);
								printf("Run time\t\t: %lg\n",run_time);
								n++;
								break;
						case 'v':
								v_flag=atof(argv[n+1]);
								printf("Rescale total volume to\t: %lg\n",v_flag);
								n++;
								break;
						case 'a':
								a_flag=atof(argv[n+1]);
								printf("Rescale total area to\t: %lg\n",a_flag);
								n++;
								break;
						case 'k':
								k_barrier=atof(argv[n+1]);
								printf("Potential barrier\t: %lg\n",k_barrier);
								n++;
								break;
						case 'e':
								epsilon=atof(argv[n+1]);
								printf("Epsilon\t\t\t: %lg\n",epsilon);
								n++;
								break;
						case 'i':
								num_of_iteration=atol(argv[n+1]);
								if(num_of_iteration==-1){printf("Will compute time step automatically\n");}
								else{printf("Total number of iterations: %ld\n",num_of_iteration);};
								n++;
								break;
						case 'x':
								export=atol(argv[n+1]);
								printf("Export configurations\t: every %ld steps\n",export);
								n++;
								break;
						case 'T':
								tol=atof(argv[n+1]);
								printf("Tolerance set to\t: %g\n",tol);
								n++;
								break;
						case 'l':
								conserved=0;
								printf("Order parameter will not be conserved\n");
								break;
						case 'C':
								gamma_h=atof(argv[n+1]);
								gamma_h2=atof(argv[n+2]);
								gamma_kg=atof(argv[n+3]);
								printf("Couplings\t\t: (%lg H, %lg H^2, %lg KG )\n",gamma_h,gamma_h2,gamma_kg);
								n+=3;
								break;
						case 'P':
								c_x=atof(argv[n+1]);
								c_y=atof(argv[n+2]);
								c_z=atof(argv[n+3]);
								printf("Center of sphere\t: (%lg,%lg,%lg)\n",c_x,c_y,c_z);
								c_flag=1;
								n+=3;
								break;
						case 'A':
								a_x=atof(argv[n+1]);
								a_y=atof(argv[n+2]);
								a_z=atof(argv[n+3]);
								norm=sqrt(a_x*a_x+a_y*a_y+a_z*a_z);
								a_x/=norm;
								a_y/=norm;
								a_z/=norm;
								printf("North pole direction\t: (%lg,%lg,%lg)\n",a_x,a_y,a_z);
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
									default:printf("Illegal integration method  %c\n",(int)argv[n+1][0]);
										printf("\nType ./membrane -h for help\n"); 
										exit(0);
										break;
								}
								printf("Integration method\t: %ld\n",method);		
								n++;
								break;
						case 'r':
								seed=atol(argv[n+1]);
								c_0=atof(argv[n+2]);
								if(c_0>=0&&c_0<=1){
									printf("Random initial data\t: (seed %ld, total concentration %lg)\n",seed,c_0);
									i_flag+=2;
									n+=2;
								}else{
									printf("Error, total relative concentration has to be between 0 and 1 (you passed %lg) \n",c_0);
									printf("\nType ./membrane -h for help\n"); 
									exit(0);
								};
								break;
						case 'R':
								printf(import_name, sizeof(import_name), "%s", argv[n+1]);
								printf("Initial configuration\t: %s\n",import_name);
								i_flag+=1;
								n++;
								break;             
						default:  
							printf("Illegal option code \"-%c\"\n",(int)argv[n][1]);
							printf("\nType ./membrane -h for help\n"); 
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
		DT = run_time/num_of_iteration;
		if( run_time==0 ||  method == 0 || (num_of_iteration == 1 && method < 4) || i_flag <1 || a_flag*v_flag>0 || a_flag < 0 || v_flag <0){
				printf("\nError: not enough arguments given or wrong parameter values. Write input in the following format:\n"); 
				print_cmd_line();
				printf("Type ./membrane -h for help\n"); 
				exit(0);
		}
		if(export == 0){
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
	if(tol==0 && method>=4){
		printf("You selected -I %ld but no tolerance was set: defaulting to 1E-6.\n",method);
		tol=1E-6;
	}
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
	printf("\t./membrane -m MESH_FILE -t RUN_TIME -I METHOD (-r SEED MEAN_CONCENTRATION | -R START_FILE ) [-e EPSILON] [-T TOL] [-L LEVEL] [-x STEPS] [-i TOTAL_ITERATIONS] [-C GAMMA_H GAMMA_H^2 GAMMA_KG] [-P CX CY CZ] [-A NX NY NZ] [-k BARRIER] [-l] [-a AREA] [-v VOL]\n\n");
}

/*******************************************************************/

void help()
{
	printf("\nThis program takes input from in-line commands. The syntax is (commands delimited by [] are optional):\n\n");	
	print_cmd_line();
	printf("Where:\n");
	printf("\t -m MESH_FILE\t: import mesh from file (works only with standard gmsh mesh format .msh)\n");
	printf("\t -t RUN_TIME\t: total simulation time. For adaptive stepsize is an (obligatory) upper bound\n");
	printf("\t -I METHOD\t: choose integration method\n\t\t\t\t1: Euler\n\t\t\t\t2: RK2\n\t\t\t\t3: RK4\n\t\t\t\t4: RK2-Euler with adaptive step-size [NOT FULLY WORKING]\n\t\t\t\t5: RKF45 with adaptive stepsize\n");
	printf("\t -r $1 $2\t: random intial condition with seed $1 and mean concentration $2\n");
	printf("\t -R FILE\t: import initial configuration from file\n");
	printf("\t -e EPSILON\t: set the value of epsilon. If not set is computed automatically from average edge length\n");
	printf("\t -T TOL\t\t: set the tolerance for adaptive step-size integration methods\n");
	printf("\t -L LEVEL\t: choose which output files will be printed\n\t\t\t\t0: 'last.dat', 'final.dat'\n\t\t\t\t1: previous +  'histo.dat', 'geometry.dat', 'last.m', 'interface.m', 'triangles.dat' + 'gc_#.dat' if -x is set [DEFAULT]\n\t\t\t\t2: previous + 'mean_curvature.m', 'gaussian_curvature.m' + 'gc_#.m' if -x is set\n\t\t\t\t3: as in '1' + debug files\n");
	printf("\t -x STEPS\t: if specified and nonzero, decides the frequency with which to export field configurations \n");
	printf("\t -i ITERATIONS\t: total number of iterations (-I 4 and -I 5 do not use this parameter)\n");
	printf("\t -C $1 $2 $3\t: specifies the values of the couplings with H, H^2 and K_G\n");
	printf("\t -P $1 $2 $3\t: specifies the (x,y,z) coordinates of the center for projection of the surface onto a unit sphere\n");
	printf("\t -A $1 $2 $3\t: if -P has been given, specifies the north pole direction w.r.t coordinate axis [DEFAULT (0,0,1)]\n");
	printf("\t -k BARRIER\t: set the height of the potential barrier (default is 1)\n");
	printf("\t -l \t\t: switch off the conservation of order parameter\n");
	printf("\t -a AREA\t: rescale the mesh so that the total area is AREA\n");
	printf("\t -v VOL\t\t: rescale the mesh to that the total volume is VOL\n");

	printf("\n");
	exit(1);
}

/*******************************************************************/

void import_mesh(char *f_name)
{
	long i, v1, v2, v3, chi, triangle_test, dummy;

	if( access( f_name, F_OK ) == -1 ) {
		printf("\nError: file %s does not exist\n",f_name);
		exit(0);
	};

	FILE *f_in = fopen(f_name,"r");
	
	char line[LINESIZE];

	for (i=0; i<8; i++){
		if(fgets(line,LINESIZE-1,f_in)){};
	}
	
	if(fscanf(f_in,"%ld",&num_of_meshpoint)){};

	if (num_of_meshpoint>MAX_SIZE){
		printf("ERROR: number of vertices (%ld) exceeds MAX_SIZE (%d)\n",num_of_meshpoint,MAX_SIZE);
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
		printf("ERROR: number of triangles (%ld) exceeds MAX_SIZE (%d)\n",num_of_triangles,MAX_SIZE);
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

	// If each edge is not counted twice we might have a problem	
	if (!num_of_edges%2){
		printf("ERROR: bad triangulation, 2E = %ld\n",num_of_edges);
		exit(0);
	}

	chi = num_of_meshpoint-num_of_edges/2+num_of_triangles;
	
	if (chi!=2){
		printf("WARNING: apparently not a genus zero surface, Euler characteristic is %ld\n",chi);
		//exit(0);
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

	int a,b,c;
	double lab,lbc,lca;
	double cota,cotb,cotc;
	double theta_a,theta_b,theta_c;
	double sp,area_t;
	int this;
	double dx, dy, dz, base, x[3],y[3],z[3],rm[3];

	l_min=1E10;
	l_max=0;
	l_avg=0;

	h2_min=1E10;
	h2_max=-1E10;

	kg_min=1E10;
	kg_max=-1E10;

	long i, j, obtuse[MAX_SIZE], num_of_obtuse=0;
	
	total_area = 0;	
	total_volume = 0;
	willmore_energy = 0;	
	euler_chi= 0;	

	FILE *f_ou;

	// Setting all vertex variables to zero before looping through triangles

	for (i=0; i<num_of_meshpoint; i++){
		vertex[i].hx = 0;
		vertex[i].hy = 0;
		vertex[i].hz = 0;
		vertex[i].area = 0;
		vertex[i].kg=2*PI;
		for (j=0; j<vertex[i].num_of_neighbors; j++){
			vertex[i].weight[j] = 0;
		}
	}

	// Loop through triangles computing vertices weights, mean and gaussian curvatures and areas

	for (i=0; i<num_of_triangles; i++){

		a=triangle[i].v1;
		b=triangle[i].v2;
		c=triangle[i].v3;

		lab=sqrt(pow((vertex[a].x-vertex[b].x),2.)+pow((vertex[a].y-vertex[b].y),2.)+pow((vertex[a].z-vertex[b].z),2.));
		lbc=sqrt(pow((vertex[b].x-vertex[c].x),2.)+pow((vertex[b].y-vertex[c].y),2.)+pow((vertex[b].z-vertex[c].z),2.));
		lca=sqrt(pow((vertex[c].x-vertex[a].x),2.)+pow((vertex[c].y-vertex[a].y),2.)+pow((vertex[c].z-vertex[a].z),2.));

		l_avg+=.5*(lab+lbc+lca);

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

	l_avg*=2./num_of_edges;

	if(o_flag>=1){
		f_ou = fopen("triangles.dat","w");
		for (i=0; i<num_of_triangles; i++){
			fprintf(f_ou,"%ld\t%ld\t%ld\n",
				triangle[i].v1,
				triangle[i].v2,
				triangle[i].v3);
		}
		fclose(f_ou);
	}

	if(a_x*a_x<1.){

		x[0]=0;
		x[1]=-a_z/sqrt(a_z*a_z+a_y*a_y);
		x[2]=a_y/sqrt(a_z*a_z+a_y*a_y);

		z[0]=x[1]*a_z-x[2]*a_y;
		z[1]=x[2]*a_x-x[0]*a_z;
		z[2]=x[0]*a_y-x[1]*a_x;

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

	//printf("Orthonomal basis for spherical projection: n=(%lg,%lg,%lg), m=(%lg,%lg,%lg), o=(%lg,%lg,%lg)\n",a_x,a_y,a_z,x[0],x[1],x[2],z[0],z[1],z[2]);

	// One further loop through vertices to normalize things

	for (i=0; i<num_of_meshpoint; i++){

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

		vertex[i].h2_avg = vertex[i].h2/(vertex[i].num_of_neighbors+1);
		vertex[i].kg_avg = vertex[i].kg/(vertex[i].num_of_neighbors+1);

		for (j=0; j<vertex[i].num_of_neighbors; j++){
			vertex[i].weight[j] /= vertex[i].area;	
		}

		// Compute projective angles

		if(c_flag==1){

			dx = vertex[i].x-c_x;
			dy = vertex[i].y-c_y;
			dz = vertex[i].z-c_z;

			dx/=sqrt(dx*dx+dy*dy+dz*dz);
			dy/=sqrt(dx*dx+dy*dy+dz*dz);
			dz/=sqrt(dx*dx+dy*dy+dz*dz);

			base=dx*a_x+dy*a_y+dz*a_z;

			rm[0]=dx-base*a_x;
			rm[1]=dy-base*a_y;
			rm[2]=dz-base*a_z;

			vertex[i].gridx=atan2(rm[0]*x[0]+rm[1]*x[1]+rm[2]*x[2],rm[0]*z[0]+rm[1]*z[1]+rm[2]*z[2]);
			vertex[i].gridy=atan2(base,sqrt(rm[0]*rm[0]+rm[1]*rm[1]+rm[2]*rm[2]));

		}

		total_area+=vertex[i].area;
		willmore_energy += vertex[i].area*vertex[i].h2;
		euler_chi += .5/PI*vertex[i].kg*vertex[i].area;

	}

	// Try to point all normal vectors outwards: this does only one full loop through vertices
	// It is not safe to assume that after just this iteration any mesh will be nicely oriented. 
	// However, it has been working so far, but remember, the general solution is more complicated
	// Moreover, use this loop to compute averaged values of H2 and KG

	for (i=0; i<num_of_meshpoint; i++){
		for (j=0; j<vertex[i].num_of_neighbors; j++){

			a = vertex[i].neighbor[j];

			base = vertex[i].nx*vertex[a].nx+vertex[i].ny*vertex[a].ny+vertex[i].nz*vertex[a].nz;

			if(base<0){vertex[a].nx*=-1;vertex[a].ny*=-1;vertex[a].nz*=-1;};

			vertex[i].h2_avg+=vertex[a].h2/(vertex[i].num_of_neighbors+1);
			vertex[i].kg_avg+=vertex[a].kg/(vertex[i].num_of_neighbors+1);

		}
	}

	// Compute the total volume

	for (i=0; i<num_of_meshpoint; i++){
		total_volume+=1./3.*vertex[i].area*(vertex[i].nx*vertex[i].x+vertex[i].ny*vertex[i].y+vertex[i].nz*vertex[i].z);
	}

	// If rescaling was chosen, rescale all Vertex attributes accordingly, rescale mesh details and recompute area and volume
	if(a_flag>0 || v_flag>0){

		if(a_flag>0)scale=pow(a_flag/total_area,1./2.);
		if(v_flag>0)scale=pow(v_flag/total_volume,1./3.);

		total_area=0;
		total_volume=0;
		
		for (i=0; i<num_of_meshpoint; i++){

			vertex[i].x*=scale;
			vertex[i].y*=scale;
			vertex[i].z*=scale;

			vertex[i].hx/=scale;
			vertex[i].hy/=scale;
			vertex[i].hz/=scale;

			vertex[i].h2/=pow(scale,2.);
			vertex[i].kg/=pow(scale,2.);

			vertex[i].h2_avg/=pow(scale,2.);
			vertex[i].kg_avg/=pow(scale,2.);

			vertex[i].area*=pow(scale,2.);

			for (j=0; j<vertex[i].num_of_neighbors; j++){
				vertex[i].weight[j]/=pow(scale,2.);
			}

			total_volume+=1./3.*vertex[i].area*(vertex[i].nx*vertex[i].x+vertex[i].ny*vertex[i].y+vertex[i].nz*vertex[i].z);
			total_area+=vertex[i].area;
		}
	
		l_avg*=scale;

	};
	
	if(o_flag>=1)f_ou = fopen("geometry.dat","w");

	// If -O 1 or higher, save geometry.dat
	// Moreover, compute min/max lengths and curvatures

	for (i=0; i<num_of_meshpoint; i++){

		if(o_flag>=1){if(c_flag==0){
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
		}else if(c_flag==1){
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
			if(o_flag>=1)fprintf(f_ou,"%.10f\t",vertex[i].weight[j]);

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
			
		if(o_flag>=1)fprintf(f_ou,"\n");
	}	
	if(o_flag>=1)fclose(f_ou);

	if(o_flag>=2){
		f_ou = fopen("mean_curvature.m","w");
		export_graphic_complex(f_ou,2);
		fclose(f_ou);

		f_ou = fopen("gaussian_curvature.m","w");
		export_graphic_complex(f_ou,3);
		fclose(f_ou);

	}

	// Print summary of mesh computations
	printf("\n");
	printf("\tVertices\t\t\t: %ld\n",num_of_meshpoint);
	printf("\tTriangles\t\t\t: %ld\n",num_of_triangles);
	printf("\tObtuse triangles\t\t: %ld (%.1ld%% of total)\n",num_of_obtuse/2,100*num_of_obtuse/2/num_of_triangles);

	printf("\tTotal surface area\t\t: %lg (typical length %lg)\n",total_area,sqrt(total_area));
	printf("\tTotal enclosed volume\t\t: %lg (typical length %lg)\n",total_volume,pow(total_volume,1./3.));

	printf("\tWillmore energy\t\t\t: %lg (asphericity %lg)\n",willmore_energy,willmore_energy/4/PI-1);
	printf("\tEuler characteristic\t\t: %lg\n",euler_chi);

	printf("\tMin/max/avg edge length \t: (%lg,%lg,%lg)\n",l_min,l_max,l_avg);
	printf("\tMin/max/avg H2 \t\t\t: (%lg,%lg,%lg)\n",h2_min,h2_max,willmore_energy/total_area);
	printf("\tMin/max/avg KG \t\t\t: (%lg,%lg,%lg)\n",kg_min,kg_max,2*PI*euler_chi/total_area);

	// The interface profile anywhere on the surface should countain at least six grid points:
	if(epsilon==0){
		epsilon=1.1*l_max*sqrt(k_barrier/2);
		printf("\tSetting epsilon to\t\t: %f\n",epsilon);
	}else{
		printf("\tProposed epsilon\t\t: %lg (minimal), %lg (averaged) or %lg (conservative)\n",l_max*sqrt(k_barrier/2),l_avg*sqrt(k_barrier/2),1.1*l_max*sqrt(k_barrier/2));
	}

	// Interface thickness 
	printf("\tInterface thickness\t\t: %lg (roughly %2.1f%% - %2.1f%% of membrane size)\n",3*sqrt(2)*epsilon/sqrt(k_barrier),3*sqrt(2)*epsilon/sqrt(k_barrier)/sqrt(total_area)*100.,3*sqrt(2)*epsilon/sqrt(k_barrier)/pow(total_volume,1./3.)*100.);

	// Line tension
	printf("\tExpected line tension\t\t: %lg\n",2./3.*sqrt(2*k_barrier)*epsilon);

	// Real couplings are obtained from -C by multiplying back to real values
	gamma_h*=2./3.*sqrt(2*k_barrier*total_area);
	gamma_h2*=2./3.*sqrt(2*k_barrier*total_area);
	gamma_kg*=2./3.*sqrt(2*k_barrier*total_area);
	printf("\tCouplings entering EOMs\t\t: (%lg H, %lg H^2, %lg KG)\n",gamma_h,gamma_h2,gamma_kg);

	// The couplings in the \epsilon \to 0 limit could take any value. 
	// However, since numerically epsilon is finite, this sets a bound on the magnitude of the curvature couplings in order to preserve the double-well structure of the potential			
	printf("\tAllowed coupling ranges\t\t: |Leibler|<%lg, |Delta k|<%lg, |Delta k_b|<%lg\n",4./3.*k_barrier/epsilon/sqrt(h2_max),4./3.*k_barrier/epsilon/h2_max,4./3.*k_barrier/epsilon/max(kg_max,sqrt(kg_min*kg_min)));  

	// For diffusion processes are stable only for DT/DX^2 \simeq 1/2:
	DTauto=.1*l_min*l_min/2;

	if(method>=4){
		DT=DTauto;
		printf("\tSetting initial time-step to\t: %lg\n",DT);
	}else{
		printf("\tProposed initial time-step\t: %lg\n",DTauto);  
	}

	printf("\n");

	//exit(1);
	
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
	c_0/=total_area;
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

		lagrange += vertex[i].area*(laplace(i)-dV(i)/(epsilon*epsilon));		
	}		
	
	lagrange /= total_area;
	
	// Calculate the right-hand side of the Allen-Cahn equation
	
	for (i=0; i<num_of_meshpoint; i++){
		
		rhs[i] = laplace(i)-dV(i)/(epsilon*epsilon)-conserved*lagrange;
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

// Define the potential and its functional derivative 

double V(long i)
{
	double V_0,V_int;

	V_0 = k_barrier*.25*pow(pow(vertex[i].phi,2.)-1,2.);

	V_int = .25*pow((vertex[i].phi+1),2.)*(2-vertex[i].phi)*epsilon;
	
	return V_0+V_int*(gamma_h2*vertex[i].h2+gamma_kg*vertex[i].kg+gamma_h*sqrt(vertex[i].h2));
}

double dV(long i)
{
	double dV_0,dV_int;

	dV_0 = k_barrier*vertex[i].phi*(vertex[i].phi*vertex[i].phi-1);

	dV_int = .75*(1-vertex[i].phi*vertex[i].phi)*epsilon;
	
	return dV_0+dV_int*(gamma_h2*vertex[i].h2+gamma_kg*vertex[i].kg+gamma_h*sqrt(vertex[i].h2));
}


/*******************************************************************/

void run()
{
	long t=0;
	char f_math_na[32],f_dat_na[32];
	
	FILE *f_hi, *f_ou;
	
	if(o_flag>=1)f_hi = fopen("histo.dat","w");	

	current_time=0;

	if(export>0 && c_flag==0 && o_flag!=2){
		printf("You chose to export configurations every %ld time steps, but did not specify a center nor chose -O 2. Assuming -C 0 0 0.\n",export);
		c_flag=1;
		c_x=0;
		c_y=0;
		c_z=0;
	}

	if(num_of_iteration==-1 || method == 4 || method == 5){
		while(current_time<run_time && halt_now==0){

			current_time+=DT;

			if(o_flag>=1)write_hi(f_hi,t);

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
			if(t%100==0){printf("Step: %ld\t\t Current time:%.6g\r",t,current_time);}
			one_step();
		}
		printf("Total time steps: %ld\nSimulation time: %g\n",t,current_time);
	}
	else{
		for (t=0; t<num_of_iteration; t++){

			current_time+=DT;

			if(o_flag>=1)write_hi(f_hi,t);

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
	if(o_flag>=1)fclose(f_hi);
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
	double Q_average=0,delta,rhs_average=0;
	double DTmin=DTauto/10000,DTmax=DTauto*10000.;
	long i;

	FILE *f_ou;

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
		rhs_average+=vertex[i].area*sqrt(rhs2[i]*rhs2[i])/total_area;
	}

	if(rhs_average<tol)halt_now=1;

	delta=tol/Q_average/2;

	if(delta<.1){DT*=.1;}
	else if(delta>4.){DT*=4.;}
	else {DT*=delta;};

	if(DT<DTmin)DT=DTmin;
	if(DT>DTmax)DT=DTmax;

	if(o_flag>=3){
		f_ou = fopen("Q_avg_debug.dat","a");
		fprintf(f_ou,"%.8E\t %.8E\t %.8E\t %.8E\t %.8E\n",current_time,Q_average,delta,DT,rhs_average);
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
	double DTmin=DTauto/10000,DTmax=DTauto*10000.;
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
		//vertex[i].phi = rhs0[i]+DT*(25./216.*rhs1[i]+1408./2565.*rhs3[i]+2197./4104.*rhs4[i]-1./5.*rhs5[i]); 			// This is the correct integration step at order 4
		vertex[i].phi = rhs0[i]+DT*(16./135.*rhs1[i]+6656./12825.*rhs3[i]+28561./56430.*rhs4[i]-9./50.*rhs5[i]+2./55.*rhs6[i]); // This is the correct integration step at order 5
		Q[i] = (1/360.*rhs1[i]-128./4275.*rhs3[i]-2197./75240.*rhs4[i]+1/50.*rhs5[i]+2./55.*rhs6[i]); 				// This is their difference, which serves as an estimate to the (local) error
		if(Q[i]<0)Q[i]*=-1;
		Q_average+=vertex[i].area*Q[i]/total_area; 										// We use <abs(Q)>_\Sigma as a global error estimate 
		//rhs_average+=vertex[i].area*pow(rhs6[i]*rhs6[i],.5)/total_area;							// Estimate of the distance from equilibrium with arbitrary power. 
		rhs_average+=vertex[i].area*rhs6[i]*rhs6[i]/total_area;									// Estimate of the distance from equilibrium as <rhs^2>_\Sigma
	}

	if(rhs_average<tol)halt_now=1;

	delta=pow(tol/Q_average/2.,.25);

	if(delta<.1){DT*=.1;}
	else if(delta>1.5){DT*=1.5;}
	else {DT*=delta;};

	if(DT<DTmin)DT=DTmin;
	if(DT>DTmax)DT=DTmax;

	if(o_flag>=3){
		f_ou = fopen("Q_avg_debug.dat","a");
		fprintf(f_ou,"%.8E\t %.8E\t %.8E\t %.8E\t %.8E\n",current_time,Q_average,delta,DT,rhs_average);
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
	double c0=0,phisq=0,kin=0,pot=0,tph,tpa,phiH2=0,phiKG=0;
	int i;
	
	FILE *f_ou;
	
	time(&t2);
  	get_time(t2-t1,&cpu_time);
	
	track_domains();

	if(o_flag>=1){
		f_ou = fopen("last.m","w");
		export_graphic_complex(f_ou,1);
		fclose(f_ou);

		f_ou = fopen("interface.m","w");
		export_graphic_complex(f_ou,4);
		fclose(f_ou);
	}
	
	if(c_flag>=0){
		f_ou = fopen("last.dat","w");
		export_dat(f_ou);
		fclose(f_ou);
	};

 	
	for (i=0; i<num_of_meshpoint; i++){

		tph=vertex[i].phi;
		tpa=vertex[i].area;

		phisq 	+=tpa*tph*tph;
		c0 	+=tpa*tph;
		kin	-=tpa*.5*tph*laplace(i)*epsilon;
		pot	+=tpa*V(i)/epsilon;
		phiH2	+=tpa*tph*vertex[i].h2;
		phiKG	+=tpa*tph*vertex[i].kg;
	}

	kin/=total_area;
	pot/=total_area;
	phisq/=total_area;
	c0/=total_area;
	phiH2/=total_area;
	phiKG/=total_area;

	printf("\n");
	printf("\tNumber of domains %d\n",num_of_domains);
	if(method<=3){printf("\tTime step %g\n",DT);};
	printf("\tCPU Time %ld:%ld:%ld:%ld\n",
	 cpu_time.d,cpu_time.h,
	 cpu_time.m,cpu_time.s);
  	printf("\n");
	
	f_ou = fopen("final.dat","w");	
	fprintf(f_ou,"%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%d\n",
	current_time,
	kin,
	pot,
	kin+pot,
	phisq-c0*c0,
	phiH2-c0*willmore_energy/total_area,
	phiKG-c0*2*PI*euler_chi/total_area,
	lagrange,
	num_of_domains);
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
		kin	-=tpa*.5*tph*laplace(i)*epsilon;
		pot	+=tpa*V(i)/epsilon;
		phiH2	+=tpa*tph*vertex[i].h2;
		phiKG	+=tpa*tph*vertex[i].kg;
	}

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
	
	for (i=0; i<num_of_meshpoint; i++){
		(vertex[i].phi>0)?(bphi=1):(bphi=0);
		fprintf(f_ou,"%.12e\t%ld\n",
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
