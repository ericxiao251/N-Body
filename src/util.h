# ifndef __UTIL_H__
# define __UTIL_H__

// Print functions
void print_properties_h(void);
void print_particle(double *p);
void print_all_particles(double **P, int numParticle);

// Particle generation
void particles_gen(double **P, int light_cnt, int medium_cnt, int heavy_cnt);



#endif