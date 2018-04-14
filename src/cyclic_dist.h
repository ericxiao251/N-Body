//cyclic_dist.h
#ifndef __CYCLIC_DIST_H__
#define __CYCLIC_DIST_H__

// Master functions
void cyclic_master_send(double **P, int num_particles, int num_processes, int avg_particle, int have_extra);
void cyclic_master_receive(double **P_force, int num_processes);


void cyclic_slave_cal_force(double **P, int local_rank, int num_processes, int avg_particle, int have_extra, int totalP, int total_p_send);

#endif
