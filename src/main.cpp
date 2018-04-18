#include "define.h"

int main(int argc, char* argv[]) {
	srand (1234);
	int num_processes, my_rank, ave_particle, padding_num, total_p_cnt, have_extra;

	double min_time_of_substeps = -1, max_time_of_substeps = -1, avg_time_of_substeps = -1;
	double start_time, stop_time, total_sub_time;

	if (argc != 10) {
		printf("Usage: %s numParticlesLight numParticleMedium numParticleHeavy numSteps subSteps timeSubStep imageWidth imageHeight imageFilenamePrefix\n", argv[0]);
	}

	// MPI bookkeeping
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_processes);

	// Input
	int numParticlesLight = atoi(argv[1]);
	int numParticleMedium = atoi(argv[2]);
	int numParticleHeavy = atoi(argv[3]);
	total_p_cnt = numParticlesLight + numParticleMedium + numParticleHeavy;

	int numSteps = atoi(argv[4]);
	int subSteps = atoi(argv[5]);
	int timeSubStep = atoi(argv[6]);

	int img_width = atoi(argv[7]);
	int img_height = atoi(argv[8]);
	unsigned char* image = (unsigned char*) malloc(img_width * img_height * 3 * sizeof(unsigned char));

	// Particle per slave node
	have_extra = total_p_cnt % (num_processes - 1);
	padding_num = (have_extra == 0) ? 0 : (num_processes - 1) - have_extra;
	total_p_cnt += padding_num;
	ave_particle = total_p_cnt / (num_processes - 1);

	/*--------------------------------- Master Node ------------------------------------------*/
	if (my_rank == 0) {
		// Initilize Particles
		int i, j, k;
		int time = 0;
		double *P_data, **P, *P_force_data, **P_force;
		char img_name[20];
		char img_full_path[100];

		P_data = (double *) malloc(total_p_cnt * PARTICLE_PROPERTIES_COUNT * sizeof(double));
		P = (double **) malloc(total_p_cnt * sizeof(double *));
		P_force_data = (double *) malloc(total_p_cnt * FORCE_SUM_PROPERTIES_COUNT * sizeof(double));
		P_force = (double **) malloc(total_p_cnt * sizeof(double *));

		for (i = 0; i < total_p_cnt; ++i) {
			P[i] = &P_data[PARTICLE_PROPERTIES_COUNT * i];
			P_force[i] = &P_force_data[FORCE_SUM_PROPERTIES_COUNT * i];
		}
		particles_gen(P, numParticlesLight, numParticleMedium, numParticleHeavy, padding_num);

		//print_properties_h();
		print_all_particles(P, total_p_cnt, time);
		for (j = 0; j < numSteps; ++j) {
			for (k = 0; k < subSteps; ++k) {
				int regenerate_img = (k == subSteps-1);
				time += timeSubStep;
				//printf("============== Current Time is %lf =============\n", time);

				start_time = MPI_Wtime();

				// Cyclic distribute particles
				cyclic_master_send(P, total_p_cnt, num_processes, ave_particle, 0);

				// Receive and sum forces;
				for (i = 0; i < total_p_cnt * FORCE_SUM_PROPERTIES_COUNT; ++i) {
					P_force_data[i] = 0.0;
				}
				cyclic_master_receive(P_force, num_processes);

				// Update position and velocity based on force
				// only update img when necessary.
				update(image, P, P_force, total_p_cnt, img_width, img_height, timeSubStep, regenerate_img);

				stop_time = MPI_Wtime();

				// Calculate timing
				total_sub_time = stop_time - start_time;
				if ((min_time_of_substeps == -1) && (max_time_of_substeps == -1) && (avg_time_of_substeps == -1)) {
					min_time_of_substeps = total_sub_time;
					max_time_of_substeps = total_sub_time;
					avg_time_of_substeps = total_sub_time;
				} else {
					if (total_sub_time > max_time_of_substeps) max_time_of_substeps = total_sub_time;
					if (total_sub_time < min_time_of_substeps) min_time_of_substeps = total_sub_time;
					avg_time_of_substeps += total_sub_time;
				}

				print_all_particles(P, total_p_cnt, time);
			}
			// Save the image
			snprintf(img_name, sizeof(char) * 32, "_%05i.bmp", j);

			strcpy(img_full_path, argv[9]);
			strcat(img_full_path, img_name);
			saveBMP(img_full_path, image, img_width, img_height);
		}
		// Release Memory
		free(P);
		free(P_data);
		free(P_force_data);
		free(P_force);
		free(image);

		// Print timing
		avg_time_of_substeps /= (double) (numSteps * subSteps);
		// printf("%.6lf %.6lf %.6lf\n", min_time_of_substeps, max_time_of_substeps, avg_time_of_substeps);
		printf("%d,%d,%d,%d,%d,%.6lf,%.6lf,%.6lf\n",
			numParticlesLight, numParticleMedium, numParticleHeavy, total_p_cnt - padding_num,
			num_processes,
			min_time_of_substeps, max_time_of_substeps, avg_time_of_substeps
		);
	}
	/*--------------------------------- Slave Node ------------------------------------------*/
	else {
		int i, j, k, init_p_id, total_p_send;
		MPI_Status status;
		double *P_data, **P;

		// Receive particles sent from master node
		for (j = 0; j < numSteps; ++j) {
			for (k = 0; k < subSteps; ++k) {
				MPI_Recv(&init_p_id, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD, &status);
				MPI_Recv(&total_p_send, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD, &status);

				P_data = (double *)malloc(total_p_send * PARTICLE_PROPERTIES_COUNT * sizeof(double));
				P = (double **) malloc(total_p_send * sizeof(double *));

				for (i = 0; i < total_p_send; ++i) {
					P[i] = &P_data[PARTICLE_PROPERTIES_COUNT * i];
					MPI_Recv(&P[i][0], PARTICLE_PROPERTIES_COUNT, MPI_DOUBLE, 0, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD, &status);
					//LOG(("Slave Node %d: Receive particle %f from Master Node\n", my_rank, P[i][ID_COL]));
				}
				// Calculate forces
				cyclic_slave_cal_force(P, my_rank, num_processes, ave_particle, 0, total_p_cnt, total_p_send);

				// Release Memory
				free(P);
				free(P_data);
			}
		}

	}
	MPI_Finalize();

	return 0;
}
