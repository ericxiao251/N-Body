//cyclic_dist.h
#ifndef __CYCLIC_DIST__
#define __CYCLIC_DIST__

void cyclic_master_send(double **P, int num_particles, int num_processes, int avg_particle, int have_extra);
void cyclic_slave_cal_force(double **P, double *P_data, double **P_force, double *P_force_data, 
								int local_rank, int num_processes, int avg_particle, int have_extra, int totalP, int total_p_send);
void cyclic_slave(/*Todo*/);

#endif