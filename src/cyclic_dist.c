#include "define.h"

void cyclic_master_send(double **P, int numParticles, int numProcesses, int avgParticle, int hasExtra) {
	int i, pId, totalPSend;

	for (i = 1; i < numProcesses; ++i) {
		pId = i - 1;
		
		totalPSend = (int)(numParticles / (numProcesses - 1)) + (i <= hasExtra) ? 1 : 0;

		MPI_Send(&pId, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD);
		MPI_Send(&totalPSend, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD);
		while (pId < numParticles) {
			LOG(("Send particle %d to node %d\n\n", pId, i));
			MPI_Send(&P[pId][0], PROPERTIES_COUNT, MPI_DOUBLE, i, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD);
			pId += numProcesses - 1;
		}
	}
}

void block_cyclic_slave_receive(double **P, double *P_data, int local_rank, int numProcesses, int avgParticle, int hasExtra) {
	MPI_Status status;
	int i, initPId, totalPSend;
	
	MPI_Recv(&initPId, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD, &status);
	MPI_Recv(&totalPSend, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD, &status);

	P_data = (double *)malloc(totalPSend * PROPERTIES_COUNT * sizeof(double));
	P = (double **) malloc(totalPSend * sizeof(double *));
	for (i = 0; i < totalPSend; ++i) {
		P[i] = &P_data[PROPERTIES_COUNT * i];
		MPI_Recv(&P[i][0], PROPERTIES_COUNT, MPI_DOUBLE, 0, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD, &status);
	}

}