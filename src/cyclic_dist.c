#include "define.h"

void cyclic_master_send(double **P, int num_particles, int num_processes, int avg_particle, int have_extra) {
	int i, pId, total_p_send;

	for (i = 1; i < num_processes; ++i) {
		pId = i - 1;
		
		total_p_send = (int)(num_particles / (num_processes - 1)) + ((i <= have_extra) ? 1 : 0);

		MPI_Send(&pId, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD);
		MPI_Send(&total_p_send, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD);
		while (pId < num_particles) {
			LOG(("Master Node: Send particle %f to node %d\n", P[pId][ID_COL], i));
			MPI_Send(&P[pId][0], PROPERTIES_COUNT, MPI_DOUBLE, i, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD);
			pId += num_processes - 1;
		}
	}
}


void cyclic_slave_cal_force(double **P, double *P_data, double **P_force, double *P_force_data, 
								int local_rank, int num_processes, int avg_particle, int have_extra, int totalP, int total_p_send) {
	MPI_Status status;
	int i, j, k, ii, init_p_id, to_slave, from_slave;
	double **P_received, *P_received_data;
	force_list_node **force_lists, **force_list_pnters;

	/*------------------ Figure out where to send&receive particles ---------------------*/
	to_slave = (local_rank == num_processes - 1) ? 1 : local_rank + 1;
	from_slave = (local_rank == 1) ? num_processes - 1 : local_rank - 1;

	/*------------------------------ Allocate memory ------------------------------------*/
	P_received_data = (double *)malloc(total_p_send * PROPERTIES_COUNT * sizeof(double));
	P_received = (double **) malloc(total_p_send * sizeof(double *));

	force_lists = (force_list_node **)malloc(sizeof(force_list_node *) * total_p_send);
	force_list_pnters = (force_list_node **)malloc(sizeof(force_list_node *) * total_p_send);

	for (i = 0; i < total_p_send; ++i) {
		P_received[i] = &P_received_data[PROPERTIES_COUNT * i];
		force_lists[i] = (force_list_node *)malloc(sizeof(force_list_node));
		force_lists[i]->next = NULL;
		force_list_pnters[i] = force_lists[i];
	}

	/*------------------------------ Initial Calculation Phase ----------------------------------------*/
	for (i = 0; i < total_p_send-1; ++i) {
		for (j = i + 1; j < total_p_send; ++j) {
			LOG(("Slave Node %d: Calculate F%d, %d\n", local_rank, (int)P[i][ID_COL], (int)P[j][ID_COL]));
			grav_force_particles(force_list_pnters[i], P[i], P[j]);
			force_list_pnters[i]->next = (force_list_node *)malloc(sizeof(force_list_node));
			force_list_pnters[i] = force_list_pnters[i]->next;
		}
	}

	/*------------------------------ Ring Passing Phase ----------------------------------------*/
	for (k = 0; k < num_processes - 2; ++k) {
		// Ring Pass Cycle i
		if (k == 0) {
			MPI_Send(&P[0][0], total_p_send * PROPERTIES_COUNT, MPI_DOUBLE, to_slave, SLAVE_TO_SLAVE_TAG, MPI_COMM_WORLD);
			
			LOG(("Slave Node %d: Send Particle ", local_rank));
			for (ii=0;ii<total_p_send;++ii) {
				LOG(("%f ", P[ii][ID_COL]));
			}
			LOG(("to Slave Node %d\n", to_slave));
		} else {
			MPI_Send(&P_received[0][0], total_p_send * PROPERTIES_COUNT, MPI_DOUBLE, to_slave, SLAVE_TO_SLAVE_TAG, MPI_COMM_WORLD);
			LOG(("Slave Node %d: Send Particle ", local_rank));
			for (ii=0;ii<total_p_send;++ii) {
				LOG(("%f ", P_received[ii][ID_COL]));
			}
			LOG(("to Slave Node %d\n", to_slave));
		}
		MPI_Recv(&P_received[0][0], total_p_send * PROPERTIES_COUNT, MPI_DOUBLE, from_slave, SLAVE_TO_SLAVE_TAG, MPI_COMM_WORLD, &status);
		LOG(("Slave Node %d: Receive Particle ", local_rank));
		for (ii=0;ii<total_p_send;++ii) {
			LOG(("%f ", P_received[ii][ID_COL]));
		}
		LOG(("from Slave Node %d\n", from_slave));


		for (i = 0; i < total_p_send; ++i) {
			for (j = 0; j < total_p_send; ++j) {
				if (P_received[j][ID_COL] <= P[i][ID_COL]) {
					continue;
				}
				grav_force_particles(force_list_pnters[i], P[i], P_received[j]);
				LOG(("Slave Node %d: Calculate F%d, %d\n", local_rank, (int)P[i][ID_COL], (int)P_received[j][ID_COL]));
				force_list_pnters[i]->next = (force_list_node *)malloc(sizeof(force_list_node));
				force_list_pnters[i] = force_list_pnters[i]->next;
			}
		}
	}
	for (i = 0; i < total_p_send; ++i) {
		print_force_list(force_lists[i], force_list_pnters[i]);
	}


	/*------------------------------ Clean up ---------------------------------------------*/
	free(P_received_data);
	free(P_received);
	free(force_list_pnters);
	for (i = 0; i < total_p_send; ++i) {
		while (force_lists[i] != NULL) {
			force_list_node *nxt = force_lists[i]->next;
			free(force_lists[i]);
			force_lists[i] = nxt;
		}
	}
	free(force_lists);
}