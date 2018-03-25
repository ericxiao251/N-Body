#include "define.h"

int main(int argc, char* argv[]) {
	srand (1234);
	int num_processes, my_rank, ave_particle, have_extra;

	if (argc != 10) {
		printf("Usage: %s numParticlesLight numParticleMedium numParticleHeavy numSteps subSteps timeSubStep imageWidth imageHeight imageFilenamePrefix\n", argv[0]);
	}

	// MPI bookkeeping
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_processes);

	// Input
	int numParticlesLight = atof(argv[1]);
	int numParticleMedium = atof(argv[2]);
	int numParticleHeavy = atof(argv[3]);
	int total_p_cnt = numParticlesLight + numParticleMedium + numParticleHeavy;

	int numSteps = atof(argv[4]);
	int subSteps = atof(argv[5]);
	double timeSubStep = atof(argv[6]);

	int img_width = atof(argv[7]);
	int img_height = atof(argv[8]);
	unsigned char* image = (unsigned char*)malloc(img_width * img_height * 3 * sizeof(unsigned char));

	// Particle per slave node
	ave_particle = total_p_cnt / (num_processes - 1);
	have_extra = total_p_cnt % (num_processes - 1);

	/*--------------------------------- Master Node ------------------------------------------*/
	if (my_rank == 0) {
		// Initilize Particles
		int i;
		double *P_data, **P, *P_force_data, **P_force;

		P_data = (double *)malloc(total_p_cnt * PARTICLE_PROPERTIES_COUNT * sizeof(double));
		P = (double **) malloc(total_p_cnt * sizeof(double *));
		P_force_data = (double *)malloc(total_p_cnt * FORCE_SUM_PROPERTIES_COUNT * sizeof(double));
		P_force = (double **) malloc(total_p_cnt * sizeof(double *));

		for (i = 0; i < total_p_cnt; ++i) {
			P[i] = &P_data[PARTICLE_PROPERTIES_COUNT * i];
			P_force[i] = &P_force_data[FORCE_SUM_PROPERTIES_COUNT * i];
		}
		particles_gen(P, numParticlesLight, numParticleMedium, numParticleHeavy);


		// Start iteration
		// Print particles
		print_properties_h();
		print_all_particles(P, total_p_cnt);

		// Cyclic distribute particles
		cyclic_master_send(P, total_p_cnt, num_processes, ave_particle, have_extra);

		// Receive and sum forces;
		for (i = 0; i < total_p_cnt * FORCE_SUM_PROPERTIES_COUNT; ++i) {
			P_force_data[i] = 0.0;
		}
		cyclic_master_receive(P_force, num_processes);

		// Update P and img based on forces
		update(image, P, P_force, total_p_cnt, img_width, img_height);

		// Save the image
		saveBMP(argv[9], image, img_width, img_height);

		// Release Memory
		free(P);
		free(P_data);
		free(P_force_data);
		free(P_force);
	}
	/*--------------------------------- Slave Node ------------------------------------------*/
	else {
		int i, init_p_id, total_p_send;
		MPI_Status status;
		double *P_data, **P;

		// Receive particles sent from master node

		MPI_Recv(&init_p_id, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD, &status);
		MPI_Recv(&total_p_send, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD, &status);

		P_data = (double *)malloc(total_p_send * PARTICLE_PROPERTIES_COUNT * sizeof(double));
		P = (double **) malloc(total_p_send * sizeof(double *));
		
		for (i = 0; i < total_p_send; ++i) {
			P[i] = &P_data[PARTICLE_PROPERTIES_COUNT * i];
			MPI_Recv(&P[i][0], PARTICLE_PROPERTIES_COUNT, MPI_DOUBLE, 0, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD, &status);
			LOG(("Slave Node %d: Receive particle %f from Master Node\n", my_rank, P[i][ID_COL]));
		}

		// Calculate forces
		cyclic_slave_cal_force(P, my_rank, num_processes, ave_particle, have_extra, total_p_cnt, total_p_send);


		// Release Memory
		free(P);
		free(P_data);

	}

	free(image);

	MPI_Finalize();
	return 0;
}