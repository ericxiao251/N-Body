#include "define.h"

int main(int argc, char* argv[]) {
	srand (1234);
	int i, j;

	if (argc != 10) {
		printf("Usage: %s numParticlesLight numParticleMedium numParticleHeavy numSteps subSteps timeSubStep imageWidth imageHeight imageFilenamePrefix\n", argv[0]);
	}

	MPI_Init(&argc,&argv);

	int p, my_rank;

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	//variables
	int numParticlesLight = atof(argv[1]);
	int numParticleMedium = atof(argv[2]);
	int numParticleHeavy = atof(argv[3]);

	int numSteps = atof(argv[4]);
	int subSteps = atof(argv[5]);
	double timeSubStep = atof(argv[6]);

	int width = atof(argv[7]);
	int height = atof(argv[8]);
	unsigned char* image = (unsigned char*)malloc(width * height * sizeof(unsigned char));

	//Root node
	if (my_rank == 0) {
		int totalNumParticle = numParticlesLight + numParticleMedium + numParticleHeavy;
		double *P_data = (double *)malloc(totalNumParticle * PROPERTIES_COUNT * sizeof(double));
		double **P = (double **) malloc(totalNumParticle * sizeof(double *));
		for (i = 0; i < totalNumParticle; ++i) {
			P[i] = &P_data[PROPERTIES_COUNT * i];
		}
		particles_gen(P, numParticlesLight, numParticleMedium, numParticleHeavy);
		print_properties_h();
		print_all_particles(P, totalNumParticle);

		//Save the image
		// saveBMP(argv[9], image, width, height);
		free(P);
		free(P_data);
	}
	//Other nodes
	else {

	}

	free(image);

	MPI_Finalize();
	return 0;
}