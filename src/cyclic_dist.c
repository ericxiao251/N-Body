#include "define.h"

/* ---------------------------- Private Functions ------------------------------*/

void cyclic_slave_send_back(force_list_node **force_lists, force_list_node **force_list_pnters, int force_cnt, int total_p_send, int local_rank) {
	int i, j;
	double *force_buffer;
	force_list_node *f, *tmp, *end_list;

	MPI_Send(&force_cnt, 1, MPI_INT, 0, SLAVE_TO_MASTER_TAG, MPI_COMM_WORLD);
	force_buffer = (double *)malloc(sizeof(double) * FORCE_PROPERTIES_COUNT);
	//LOG(("Slave Node %d: Send %d force(s) to Master Node...\n", local_rank, force_cnt));
	for (i = 0; i < total_p_send; ++i) {
		f = force_lists[i];
		end_list = force_list_pnters[i];
		while (f != NULL && f != end_list) {
			force_buffer[FROM_COL] = f->pid2;
			force_buffer[TO_COL] = f->pid1;
			force_buffer[X_COL] = f->f12_x;
			force_buffer[Y_COL] = f->f12_y;
			MPI_Send(force_buffer, FORCE_PROPERTIES_COUNT, MPI_DOUBLE, 0, SLAVE_TO_MASTER_TAG, MPI_COMM_WORLD);
			tmp = f;
			f = f->next;
			free(tmp);
		}
	}
	free(force_buffer);
}

/* ---------------------------- Global Functions ------------------------------*/

// Master functions
void cyclic_master_send(double **P, int num_particles, int num_processes, int avg_particle, int have_extra) {
	int i, pId, total_p_send;

	for (i = 1; i < num_processes; ++i) {
		pId = i - 1;

		total_p_send = (int)(num_particles / (num_processes - 1)) + ((i <= have_extra) ? 1 : 0);

		MPI_Send(&pId, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD);
		MPI_Send(&total_p_send, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD);
		while (pId < num_particles) {
			LOG(("Master Node: Send particle %f to node %d\n", P[pId][ID_COL], i));
			MPI_Send(&P[pId][0], PARTICLE_PROPERTIES_COUNT, MPI_DOUBLE, i, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD);
			pId += num_processes - 1;
		}
	}
}

void cyclic_master_receive(double **P_force, int num_processes) {
	MPI_Status status;
	int i, j, force_cnt, to_force, from_force;
	double *force_buffer;
	force_buffer = (double *)malloc(sizeof(double) * FORCE_PROPERTIES_COUNT);

	for (i = 1; i < num_processes; ++i) {
		MPI_Recv(&force_cnt, 1, MPI_INT, i, SLAVE_TO_MASTER_TAG, MPI_COMM_WORLD, &status);
		// LOG(("Master Node: Receive %d force(s) to from Slave Node %d...\n", force_cnt, i));
		for (j = 0; j < force_cnt; ++j) {
			MPI_Recv(force_buffer, FORCE_PROPERTIES_COUNT, MPI_DOUBLE, i, SLAVE_TO_MASTER_TAG, MPI_COMM_WORLD, &status);
			to_force = (int)force_buffer[TO_COL];
			from_force = (int)force_buffer[FROM_COL];

			P_force[to_force][X_COL] += force_buffer[X_COL];
			P_force[to_force][Y_COL] += force_buffer[Y_COL];

			P_force[from_force][X_COL] -= force_buffer[X_COL];
			P_force[from_force][Y_COL] -= force_buffer[Y_COL];
		}
	}
	free(force_buffer);
}


// Slave functions
void cyclic_slave_cal_force(double **P, int local_rank, int num_processes, int avg_particle, int have_extra, int totalP, int total_p_send) {
	MPI_Status status;
	int i, j, k, ii, init_p_id, to_slave, from_slave, force_cnt = 0;
	double **P_received, *P_received_data;
	force_list_node **force_lists, **force_list_pnters;

	/*------------------ Figure out where to send&receive particles ---------------------*/
	to_slave = (local_rank == num_processes - 1) ? 1 : local_rank + 1;
	from_slave = (local_rank == 1) ? num_processes - 1 : local_rank - 1;

	/*------------------------------ Allocate memory ------------------------------------*/
	P_received_data = (double *)malloc(total_p_send * PARTICLE_PROPERTIES_COUNT * sizeof(double));
	P_received = (double **) malloc(total_p_send * sizeof(double *));

	force_lists = (force_list_node **)malloc(sizeof(force_list_node *) * total_p_send);
	force_list_pnters = (force_list_node **)malloc(sizeof(force_list_node *) * total_p_send);

	for (i = 0; i < total_p_send; ++i) {
		P_received[i] = &P_received_data[PARTICLE_PROPERTIES_COUNT * i];
		force_lists[i] = (force_list_node *)malloc(sizeof(force_list_node));
		force_lists[i]->next = NULL;
		force_list_pnters[i] = force_lists[i];
	}

	/*------------------------------ Initial Calculation Phase ----------------------------------------*/
	for (i = 0; i < total_p_send-1; ++i) {
		for (j = i + 1; j < total_p_send; ++j) {
			if (P[i][WEIGHT_COL] == DUMMY_WEIGHT || P[j][WEIGHT_COL] == DUMMY_WEIGHT) {
				// skip padding particles
				continue;
			}
			//LOG(("Slave Node %d: Calculate F%d, %d\n", local_rank, (int)P[i][ID_COL], (int)P[j][ID_COL]));
			grav_force_particles(force_list_pnters[i], P[i], P[j]);
			force_list_pnters[i]->next = (force_list_node *)malloc(sizeof(force_list_node));
			force_list_pnters[i] = force_list_pnters[i]->next;
			++force_cnt;
		}
	}
	/*------------------------------ Ring Passing Phase ----------------------------------------*/
	for (k = 0; k < num_processes - 2; ++k) {
		// Ring Pass Cycle i
		if (k == 0) {
			MPI_Send(&P[0][0], total_p_send * PARTICLE_PROPERTIES_COUNT, MPI_DOUBLE, to_slave, SLAVE_TO_SLAVE_TAG, MPI_COMM_WORLD);

			//for (ii=0;ii<total_p_send;++ii) {
			//	LOG(("Slave Node %d: Send Particle %f to Slave Node %d\n", local_rank, P[ii][ID_COL], to_slave));
			//}
		} else {
			MPI_Send(&P_received[0][0], total_p_send * PARTICLE_PROPERTIES_COUNT, MPI_DOUBLE, to_slave, SLAVE_TO_SLAVE_TAG, MPI_COMM_WORLD);
			//for (ii=0;ii<total_p_send;++ii) {
			//	LOG(("Slave Node %d: Send Particle %f to Slave Node %d\n", local_rank, P_received[ii][ID_COL], to_slave));
			//}
		}
		MPI_Recv(&P_received[0][0], total_p_send * PARTICLE_PROPERTIES_COUNT, MPI_DOUBLE, from_slave, SLAVE_TO_SLAVE_TAG, MPI_COMM_WORLD, &status);

		for (i = 0; i < total_p_send; ++i) {
			for (j = 0; j < total_p_send; ++j) {
				if (P[i][WEIGHT_COL] == DUMMY_WEIGHT || P_received[j][WEIGHT_COL] == DUMMY_WEIGHT
						|| P_received[j][ID_COL] <= P[i][ID_COL]) {
					continue;
				}
				grav_force_particles(force_list_pnters[i], P[i], P_received[j]);
				/////////////////LOG(("Slave Node %d: Calculate F%d, %d\n", local_rank, (int)P[i][ID_COL], (int)P_received[j][ID_COL]));
				force_list_pnters[i]->next = (force_list_node *)malloc(sizeof(force_list_node));
				force_list_pnters[i] = force_list_pnters[i]->next;
				++force_cnt;
			}
		}
	}
	//for (i = 0; i < total_p_send; ++i) {
	//	print_force_list(force_lists[i], force_list_pnters[i]);
	//}

	/*------------------------------ Send back --------------------------------------------*/
	cyclic_slave_send_back(force_lists, force_list_pnters, force_cnt, total_p_send, local_rank);
	/*------------------------------ Clean up ---------------------------------------------*/
	free(P_received_data);
	free(P_received);
	free(force_list_pnters);
	free(force_lists);
}
