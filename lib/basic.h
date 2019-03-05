// Implementations of help functions and of the pseudorandom number generators.

#define NOISE .1
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
#define RNMX (1.0-EPS)

#define LINESIZE 100

void print_cmd_line();
void help();

double ran2(long *);
double gauss_ran2(long *,double, double);

void get_time(time_t, CPU_Time *);
void progress_bar(long, long);

/*******************************************************************/

void print_cmd_line()
{
	printf("membrane -m MESH_FILE -t RUN_TIME -I METHOD (-r SEED C0 | -p SEED C0 | -R START_FILE ) [-e EPSILON] [-T TOLERANCE] [-L LEVEL] [-x STEPS] [-i TOTAL_ITERATIONS] [-C GAMMA_H GAMMA_H^2 GAMMA_KG] [-P CX CY CZ] [-A NX NY NZ] [-k BARRIER] [-g SIGMA] [-l] [-M] [-a AREA] [-v VOL] [-u T J Lk Lkb Mk Mkb] [-w CUTOFF] [-F]\n\n");
}

/*******************************************************************/

void help()
{
	printf("\nThis program takes input from in-line commands. The syntax is (commands delimited by [] are optional):\n\n");	
	print_cmd_line();
	printf("Where:\n");
	printf("\t -m MESH_FILE\t: import mesh from file (works only with standard gmsh mesh format .msh)\n");
	printf("\t -t RUN_TIME\t: total simulation time. For adaptive stepsize is an (obligatory) upper bound\n");
	printf("\t -I METHOD\t: choose integration method\n\t\t\t\t1: Euler\n\t\t\t\t2: RK2\n\t\t\t\t3: RK4\n\t\t\t\t4: RK2-Euler with adaptive stepsize\n\t\t\t\t5: RKF45 with adaptive stepsize\n");
	printf("\t -r $1 $2\t: quasi-constant random initial configuration (seed $1) centered around mean value $2 and variance %g\n",NOISE);
	printf("\t -p $1 $2\t: initial configuration is set for fields at ±1 at random, seed $1 and total concentration $2\n");
	printf("\t -R FILE\t: import initial configuration from file\n");
	printf("\t -e EPSILON\t: set the value of EPSILON (if it is not set, it will be computed automatically from average edge length)\n");
	printf("\t -T TOLERANCE\t: set the tolerance for adaptive stepsize integration methods\n");
	printf("\t -L LEVEL\t: choose which output files will be printed\n\t\t\t\t0: 'last.dat', 'final.dat'\n\t\t\t\t1: previous +  'histo.dat', 'geometry.dat', 'last.m', 'interface.m', 'triangles.dat' + 'gc_#.dat' if -x is set [DEFAULT]\n\t\t\t\t2: previous + 'mean_curvature.m', 'gaussian_curvature.m' + 'gc_#.m' if -x is set\n\t\t\t\t3: as in '1' without Mathematica files + debug dat \n");
	printf("\t -x STEPS\t: if specified and nonzero, decides the frequency with which to export field configurations \n");
	printf("\t -i ITERATIONS\t: total number of iterations (-I 4 and -I 5 do not use this parameter)\n");
	printf("\t -C $1 $2 $3\t: specifies the values of the couplings with H, H^2 and K_G\n");
	printf("\t -P $1 $2 $3\t: specifies the (x,y,z) coordinates of the center for projection of the surface onto a unit sphere\n");
	printf("\t -A $1 $2 $3\t: if -P has been given, specifies the north pole direction w.r.t coordinate axis [DEFAULT (0,0,1)]\n");
	printf("\t -k BARRIER\t: set the height of the potential barrier (default is 1)\n");
	printf("\t -g SIGMA\t: add white Gaussian noise to the EOM with standard deviation SIGMA.\n");
	printf("\t -l \t\t: switch off the conservation of order parameter\n");
	printf("\t -M \t\t: use values of curvatures averaged over nearest neighbours\n");
	printf("\t -a AREA\t: rescale the mesh so that the total area is AREA\n");
	printf("\t -v VOLUME\t: rescale the mesh so that the total volume is VOLUME (overrides -a)\n");
	printf("\t -u $1-$6\t: use lattice-gas mean-field free energy with\n\t\t\t\t -temperature T ($1)\n\t\t\t\t -quadratic homogeneous coupling J ($2)\n\t\t\t\t -quadratic mean curvature squared interaction Lk ($3)\n\t\t\t\t -quadratic Gaussian curvature interaction Lkb ($4)\n\t\t\t\t -linear mean curvature squared interaction Mk ($5)\n\t\t\t\t -linear Gaussian curvature interaction Mkb ($6)\n");
	printf("\t -w CUTOFF\t: add a cutoff for the H^2 and K absolute values (eventually NN averaged).\n");
	printf("\t -F \t\t: use quintic curvature coupling, which preserves the Maxwell values at ±1\n");

	printf("\n");
	exit(1);
}

/*******************************************************************/
//  Long period (>2 E18) random number generator of
//  L'Ecuyer with Bays-Durham suffle and added
//  safeguards. Returns a uniform random deviate between
//  0.0 and 1.0 (exclusive of the endpoint values).
//  Call with idum a negative integer to initialise;
//  thereafter, do not alter idum between successive
//  deviates in a sequence. RNMX should approximate the
//  largest floating value that is less than 1.

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

/*******************************************************************/
// Returns a normally distributed deviate with zero mean and unit variance,
// using ran2(idum) as the source of uniform random deviates.

double gauss_ran2(long *idum, double mean, double stddev)
{
  static int iset = 0;
  static double gset;
  double fac,rsq,v1,v2, deviate;

  if (*idum < 0) iset = 0; // Reinitialize
  if (iset == 0) { 
    // 
    // We don't have an extra deviate handy, so pick two uniform numbers
    // in the square extending from -1 to +1 in each direction, and
    // see if they are in the unit circle
    //
    do {
      v1 = 2.0 * ran2(idum) - 1.0; 
      v2 = 2.0 * ran2(idum) - 1.0;
      rsq = v1 * v1 + v2 * v2; 
    } while (rsq >= 1.0 || rsq == 0.0); // and if they are not, try again.
    fac = sqrt(-2.0 * log(rsq)/rsq);

    //
    // Now make the Box-Muller transformation to get two normal deviates.
    // Return one and save the other for next time.
    //
    gset = v1*fac;
    iset = 1;      // Set flag.
    deviate = v2*fac;
  }
  else {       // We have an extra deviate handy,
    iset = 0;      // so unset the flag,
    deviate = gset; // and return it.
  }

  return deviate * stddev + mean;
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

void progress_bar(long t, long tot)
{
	long i, ticks, percent, num_of_ticks=20;
	char progress[32];

	percent = 100*t/tot;
	ticks = num_of_ticks*percent/100;
	
	for (i=0; i<20; i++){
		progress[i] = ((i<=ticks)?'#':' ');
	}
	progress[20] = '\0';
	
	printf("Progress: [%s] %ld%%\r",progress,percent);
	fflush(stdout);
}
