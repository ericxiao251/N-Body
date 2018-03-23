#include "define.h"

int main(int argc, char* argv[]) {
	srand (1234);
	int i, j, num_of_process, my_rank, ave_particle, have_extra;

	if (argc != 10) {
		printf("Usage: %s numParticlesLight numParticleMedium numParticleHeavy numSteps subSteps timeSubStep imageWidth imageHeight imageFilenamePrefix\n", argv[0]);
	}

	// MPI bookkeeping
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_of_process);

	// Input
	int numParticlesLight = atof(argv[1]);
	int numParticleMedium = atof(argv[2]);
	int numParticleHeavy = atof(argv[3]);
	int totalNumParticle = numParticlesLight + numParticleMedium + numParticleHeavy;

	int numSteps = atof(argv[4]);
	int subSteps = atof(argv[5]);
	double timeSubStep = atof(argv[6]);

	int width = atof(argv[7]);
	int height = atof(argv[8]);
	unsigned char* image = (unsigned char*)malloc(width * height * sizeof(unsigned char));

	// Particle per slave node
	ave_particle = totalNumParticle / (num_of_process - 1);
	have_extra = totalNumParticle % (num_of_process - 1);

	/*--------------------------------- Master Node ------------------------------------------*/
	if (my_rank == 0) {
		// Initilize Particles
		double *P_data = (double *)malloc(totalNumParticle * PROPERTIES_COUNT * sizeof(double));
		double **P = (double **) malloc(totalNumParticle * sizeof(double *));
		for (i = 0; i < totalNumParticle; ++i) {
			P[i] = &P_data[PROPERTIES_COUNT * i];
		}
		particles_gen(P, numParticlesLight, numParticleMedium, numParticleHeavy);

		// Print particles
		print_properties_h();
		print_all_particles(P, totalNumParticle);

		// Cyclic distribute particles
		cyclic_master_send(P, totalNumParticle, num_of_process, ave_particle, have_extra);


		//Save the image
		// saveBMP(argv[9], image, width, height);

		// Release Memory
		free(P);
		free(P_data);
	}
	/*--------------------------------- Slave Node ------------------------------------------*/
	else {
		double *P_data, **P, *P_force_data, **P_force;
		int i, init_p_id, total_p_send;
		MPI_Status status;

		// Receive particles sent from master node

		MPI_Recv(&init_p_id, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD, &status);
		MPI_Recv(&total_p_send, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD, &status);

		P_data = (double *)malloc(total_p_send * PROPERTIES_COUNT * sizeof(double));
		P_force_data = (double *)malloc(total_p_send * PROPERTIES_COUNT * sizeof(double));
		P = (double **) malloc(total_p_send * sizeof(double *));
		P_force = (double **) malloc(total_p_send * sizeof(double *));
		for (i = 0; i < total_p_send; ++i) {
			P[i] = &P_data[PROPERTIES_COUNT * i];
			P_force[i] = &P_force_data[PROPERTIES_COUNT * i];
			MPI_Recv(&P[i][0], PROPERTIES_COUNT, MPI_DOUBLE, 0, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD, &status);
			LOG(("Slave Node %d: Receive particle %f from Master Node\n", my_rank, P[i][ID_COL]));
		}

		// Calculate forces
		cyclic_slave_cal_force(P, P_data, P_force, P_force_data, my_rank, num_of_process, ave_particle, have_extra, totalNumParticle, total_p_send);


		// Release Memory
		free(P);
		free(P_data);
		free(P_force_data);
		free(P_force);
	}

	free(image);

	MPI_Finalize();
	return 0;
}