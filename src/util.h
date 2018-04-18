#ifndef __UTIL_H__
#define __UTIL_H__

// Print functions
void print_properties_h(void);
void print_particle(double *p, double time);
void print_all_particles(double **P, int numParticle, double time);

// Particle generation
void particles_gen(double **P, int light_cnt, int medium_cnt, int heavy_cnt, int padding_cnt);

// Update particle
void update(unsigned char* image, double **P, double **P_force,
	int total_p_cnt, int img_width, int img_height,
	int step_size, int regenerate_img);



#endif
