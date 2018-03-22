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
		double *P_data, **P;

		// Receive particles sent from master node
		block_cyclic_slave_receive(P, P_data, my_rank, num_of_process, ave_particle, have_extra);

	}

	free(image);

	MPI_Finalize();
	return 0;
}