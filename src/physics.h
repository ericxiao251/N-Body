#ifndef __PHYSICS_H__
#define __PHYSICS_H__

#define epsilon 0.000000000000000222
#define TWO_PI 6.28318530718
#define GRAVITY 6.673E-11

typedef struct force_list_node {
	int pid1;
	int pid2;
	double f12_x;
	double f12_y;
	struct force_list_node *next;
}force_list_node;

// Physics Function
void grav_force(double *f12_x, double *f12_y, double mass1, double mass2, double s1_x, double s1_y, double s2_x, double s2_y);
void grav_force_particles(force_list_node* force, double *p1, double *p2);
void rev_force(double *f21_x, double *f21_y, double f12_x, double f12_y);

void print_force_list(force_list_node *list, force_list_node *end_list);

#endif