// This library defines the functions used to import a mesh and compute all the geometric quantities associated to it

void import_mesh(char *);
void find_neighbors();
int is_neighbor(long,long);

long num_of_meshpoint,num_of_triangles,num_of_edges=0;

Triangle triangle[MAX_SIZE];
Vertex vertex[MAX_SIZE];

/*******************************************************************/

void find_neighbors()
{

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

