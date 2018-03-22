//cyclic_dist.h
#ifndef __CYCLIC_DIST__
#define __CYCLIC_DIST__

void cyclic_master_send(double **P, int numParticles, int numProcesses, int avgParticle, int hasExtra);
void block_cyclic_slave_receive(double **P, double *P_data, int local_rank, int numProcesses, int avgParticle, int hasExtra);
void cyclic_slave(/*Todo*/);

#endif