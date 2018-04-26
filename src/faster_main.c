#include "define.h"

#define MASTER 0
#define DIM 2
#define X 0
#define Y 1
#define GENERATE_INPUT_FILE 1

typedef double vector[DIM];
const double G = 6.673e-11;

int my_rank;
int num_processes;
int img_width, img_height;

vector *velocities = NULL;

void gen_init(int types[], double masses[], vector positions[], vector my_velocities[], int numParticlesLight, int numParticleMedium, int numParticleHeavy, int num_particles, int chunk);
void state_update(unsigned char* image, int types[], double masses[], vector positions[], vector my_velocities[], int num_particles, int chunk, int img_width, int img_height);
void compute_f(int my_particles, double masses[], vector my_forces[], vector positions[], int num_particles, int chunk);
void update_p(int my_particles, double masses[], vector my_forces[], vector my_positions[], vector my_velocities[], int num_particles, int chunk, double time_substep);
void test_gen_output(double masses[], vector positions[], vector my_velocities[], int num_particles, int chunk, double time);

int main(int argc, char* argv[]) {
	srand (1234);

	int numParticlesLight, numParticleMedium, numParticleHeavy, num_particles;

	int chunk, have_extra, padding_num, num_steps, steps, num_sub_steps, sub_steps, my_particles, output_freq, *types;

	double *masses;

	vector* my_positions, * positions, * my_velocities, * my_forces;

	double time_substep, time;
	double min_time_of_substeps = -1, max_time_of_substeps = -1, avg_time_of_substeps = -1;
	double start_time, stop_time, total_sub_time;

	// MPI bookkeeping
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_processes);

  // Get input
	if(my_rank == MASTER){
		numParticlesLight = strtol(argv[1], NULL, 10);
		numParticleMedium = strtol(argv[2], NULL, 10);
		numParticleHeavy = strtol(argv[3], NULL, 10);
		num_steps = strtol(argv[4], NULL, 10);
		num_sub_steps = atoi(argv[5]);
		time_substep = strtod(argv[6], NULL);

		num_particles = numParticlesLight + numParticleMedium + numParticleHeavy;
		output_freq = 1;

		img_width = atoi(argv[7]);
		img_height = atoi(argv[8]);
	}

	MPI_Bcast(num_particles, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	MPI_Bcast(num_steps, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	MPI_Bcast(time_substep, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	MPI_Bcast(output_freq, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

	have_extra = num_particles % num_processes;
	padding_num = (have_extra == 0) ? 0 : num_processes - have_extra;
	num_particles += padding_num;
	chunk = num_particles / num_processes;

	types = (int *) malloc(num_particles*sizeof(int));
	masses = (double *) malloc(num_particles*sizeof(double));
	positions = (vector *) malloc(num_particles*sizeof(vector));
	my_forces = (vector *) malloc(chunk*sizeof(vector));
	my_positions = positions + my_rank*chunk;
	my_velocities = (vector *) malloc(chunk*sizeof(vector));

	if (my_rank == MASTER) {
		velocities = (vector *) malloc(num_particles*sizeof(vector));
	}

	unsigned char* image = (unsigned char*) malloc(img_width * img_height * 3 * sizeof(unsigned char));

	// Generate particles
	gen_init(types, masses, positions, my_velocities, numParticlesLight, numParticleMedium, numParticleHeavy, num_particles, chunk);

	for (steps = 1; steps <= num_steps; steps++) {
		for (sub_steps = 1; sub_steps <= num_sub_steps; sub_steps++) {
			time += time_substep;

			start_time = MPI_Wtime();

			// Each processor calculates their particle's force
			for (my_particles = 0; my_particles < chunk; my_particles++) {
				if (types[my_particles] == DUMMY) {
					continue;
				}
				compute_f(my_particles, masses, my_forces, positions, num_particles, chunk);
			}

			// Each processor updates their particle's position
			for (my_particles = 0; my_particles < chunk; my_particles++) {
				if (types[my_particles] == DUMMY) {
					continue;
				}
				update_p(my_particles, masses, my_forces, my_positions, my_velocities, num_particles, chunk, time_substep);
			}

			stop_time = MPI_Wtime();

			// Calculate timing
			total_sub_time = stop_time - start_time;
			if ((min_time_of_substeps == -1) && (max_time_of_substeps == -1) && (avg_time_of_substeps == -1)) {
				min_time_of_substeps = total_sub_time;
				max_time_of_substeps = total_sub_time;
				avg_time_of_substeps = total_sub_time;
			} else {
				if (total_sub_time > max_time_of_substeps) {
					max_time_of_substeps = total_sub_time;
				}
				if (total_sub_time < min_time_of_substeps) {
					min_time_of_substeps = total_sub_time;
				}
				avg_time_of_substeps += total_sub_time;
			}
		}

		MPI_Allgather(MPI_IN_PLACE, DIM * chunk, MPI_DOUBLE, positions, DIM * chunk, MPI_DOUBLE, MPI_COMM_WORLD);

		// Save BMP
		initilize_img(image, img_width, img_height);
		state_update(image, types, masses, positions, my_velocities, num_particles, chunk, img_width, img_height);

		// Save the image
		char img_name[100];
		char img_full_path[100];
		snprintf(img_name, sizeof(char) * 32, "_%05i.bmp", steps);

		strcpy(img_full_path, argv[9]);
		strcat(img_full_path, img_name);
		saveBMP(img_full_path, image, img_width, img_height);
	}

	test_gen_output(masses, positions, my_velocities, num_particles, chunk, time);

	// Print timing
	if(my_rank == MASTER) {
		avg_time_of_substeps /= (double) (num_steps * num_sub_steps);
		printf("%.6lf %.6lf %.6lf\n", min_time_of_substeps, max_time_of_substeps, avg_time_of_substeps);
	}

	// Free memory
	free(masses);
	free(positions);
	free(my_forces);
	free(my_velocities);
	if (my_rank == MASTER) {
		free(velocities);
	}
	MPI_Finalize();
	return 0;
}


void gen_init(int types[], double masses[], vector positions[], vector my_velocities[], int numParticlesLight, int numParticleMedium, int numParticleHeavy, int num_particles, int chunk){
	int part;
	int line;
	double mass = 5;
	double gap = 10;

	if (my_rank == MASTER) {
		// light particles
		for (part = 0; part < numParticlesLight; ++part) {
			types[part] = LIGHT;
			masses[part] = massLightMin + (massLightMax - massLightMin) * drand48();
			positions[part][X] = POS_MIN_X + (POS_MAX_X - POS_MIN_X) * drand48();
			positions[part][Y] = POS_MIN_Y + (POS_MAX_Y - POS_MIN_Y) * drand48();

			// generate velocities
			double max_v, min_v, rand_v, rand_deg;
			max_v = velocityLightMin;
			min_v = velocityLightMax;
			rand_v = min_v + (max_v - min_v) * drand48();
			rand_deg = TWO_PI * drand48();
			velocities[part][X] = rand_v * cos(rand_deg);
			velocities[part][Y] = rand_v * sin(rand_deg);
		}

		// medium particles
		for (part = numParticlesLight; part < numParticlesLight+numParticleMedium; ++part) {
			types[part] = MEDIUM;
			masses[part] = massMediumMin + (massMediumMax - massMediumMin) * drand48();
			positions[part][X] = POS_MIN_X + (POS_MAX_X - POS_MIN_X) * drand48();
			positions[part][Y] = POS_MIN_Y + (POS_MAX_Y - POS_MIN_Y) * drand48();

			// generate velocities
			double max_v, min_v, rand_v, rand_deg;
			max_v = velocityMediumMin;
			min_v = velocityMediumMax;
			rand_v = min_v + (max_v - min_v) * drand48();
			rand_deg = TWO_PI * drand48();
			velocities[part][X] = rand_v * cos(rand_deg);
			velocities[part][Y] = rand_v * sin(rand_deg);
		}

		// heavy particles
		for (part = numParticlesLight+numParticleMedium; part < numParticlesLight+numParticleMedium+numParticleHeavy; ++part) {
			types[part] = HEAVY;
			masses[part] = massHeavyMin + (massLightMax - massHeavyMin) * drand48();
			positions[part][X] = POS_MIN_X + (POS_MAX_X - POS_MIN_X) * drand48();
			positions[part][Y] = POS_MIN_Y + (POS_MAX_Y - POS_MIN_Y) * drand48();

			// generate velocities
			double max_v, min_v, rand_v, rand_deg;
			max_v = velocityHeavyMin;
			min_v = velocityHeavyMax;
			rand_v = min_v + (max_v - min_v) * drand48();
			rand_deg = TWO_PI * drand48();
			velocities[part][X] = rand_v * cos(rand_deg);
			velocities[part][Y] = rand_v * sin(rand_deg);
		}

		// padding particles
		for (part = numParticlesLight+numParticleMedium+numParticleHeavy; part < num_particles; ++part) {
			types[part] = DUMMY;
			masses[part] = 0;
			positions[part][X] = 0;
			positions[part][Y] = 0;

			// generate velocities
			velocities[part][X] = 0;
			velocities[part][Y] = 0;
		}

		// generating initial state of simulation for validation testing
		#ifdef GENERATE_INPUT_FILE
		FILE *fp_generate = fopen("../data/input_file.csv","w+");
		if (!fp_generate) {
		  printf("Cannot create input file <%s>\n", "../data/input_file");
		  exit(1);
		}

		for(line = 0; line < num_particles; line++){
			fprintf(fp_generate, "%lf,%d,%lf,%lf,%lf,%lf,%lf\n",
				0, line, masses[line],
				positions[line][X], positions[line][Y],
				velocities[line][X], velocities[line][Y]
			);
		}

		fclose(fp_generate);
		#endif
	}

	// send particle data to all the processors
	MPI_Bcast(types, num_particles, MPI_INT, MASTER, MPI_COMM_WORLD);
	MPI_Bcast(masses, num_particles, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	MPI_Bcast(positions, DIM*num_particles, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	MPI_Scatter(velocities, DIM * chunk, MPI_DOUBLE, my_velocities, DIM * chunk, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
}

void state_update(unsigned char* image, int types[], double masses[], vector positions[], vector my_velocities[], int num_particles, int chunk, int img_width, int img_height) {
	int part;
	MPI_Gather(my_velocities, DIM * chunk, MPI_DOUBLE, velocities, DIM * chunk, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	if (my_rank == MASTER) {
		for (part = 0; part < num_particles; part++) {
			update_img(image, positions[part][X], positions[part][Y], types[part], img_width, img_height);
		}
	}
}

void compute_f(int my_particles, double masses[], vector my_forces[], vector positions[], int num_particles, int chunk){
	int k, part;
	double m_g;
	vector f_part_k;
	double len, len_3, fact;

	part = my_rank * chunk + my_particles;
	my_forces[my_particles][X] = my_forces[my_particles][Y] = 0.0;

	for (k = 0; k < num_particles; k++) {
		if (k != part) {
			if (masses[k] <= 0) {
				continue;
			}
			f_part_k[X] = positions[part][X] - positions[k][X];
			f_part_k[Y] = positions[part][Y] - positions[k][Y];

			len = sqrt(pow(f_part_k[X], 2) + pow(f_part_k[Y], 2));
			len_3 = pow(len, 3);

			m_g = G * masses[part] * masses[k];
			fact = m_g / len_3;

			f_part_k[X] *= fact;
			f_part_k[Y] *= fact;

			my_forces[my_particles][X] += f_part_k[X];
			my_forces[my_particles][Y] += f_part_k[Y];
		}
	}
}

void update_p(int my_particles, double masses[], vector my_forces[], vector my_positions[], vector my_velocities[], int num_particles, int chunk, double time_substep) {
	int part;
	double fact;

	part = my_rank * chunk + my_particles;
	fact = time_substep / masses[part];

	my_positions[my_particles][X] += time_substep * my_velocities[my_particles][X];
	my_positions[my_particles][Y] += time_substep * my_velocities[my_particles][Y];
	my_velocities[my_particles][X] += fact * my_forces[my_particles][X];
	my_velocities[my_particles][Y] += fact * my_forces[my_particles][Y];
}

void test_gen_output(double masses[], vector positions[], vector my_velocities[], int num_particles, int chunk, double time) {
	int line;
	MPI_Gather(my_velocities, DIM * chunk, MPI_DOUBLE, velocities, DIM * chunk, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	if (my_rank == MASTER) {
		FILE *fp_generate_output = fopen("../data/output_file.csv","w+");
		if (!fp_generate_output) {
			printf("Cannot create output file <%s>\n", "..timing/output_file");
			exit(1);
		}

		for(line = 0; line < num_particles; line++){
		 	fprintf(fp_generate_output, "%lf,%d,%lf,%lf,%lf,%lf,%lf\n",
		 		time, line, masses[line],
		 		positions[line][X], positions[line][Y],
		 		velocities[line][X], velocities[line][Y]
		 	);
		}

		fclose(fp_generate_output);
	}
}
