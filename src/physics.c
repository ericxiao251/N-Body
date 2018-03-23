#include "define.h"

void grav_force(double *f12_x, double *f12_y, double mass1, double mass2, double s1_x, double s1_y, double s2_x, double s2_y) {
	double dist_sq, scalar;

	dist_sq = pow(s1_x - s2_x, 2) + pow(s1_y - s2_y, 2);
	scalar = (GRAVITY * mass1 * mass2) / dist_sq;
	*f12_x = scalar * (s2_x - s1_x);
	*f12_y = scalar * (s2_y - s1_y);
}

void grav_force_particles(force_list_node* force, double *p1, double *p2) {
	double x, y;
	
	grav_force(&x, &y, p1[WEIGHT_COL], p2[WEIGHT_COL], p1[POS_X_COL], p1[POS_Y_COL], p2[POS_X_COL], p2[POS_Y_COL]);

	force->pid1 = (int)p1[ID_COL];
	force->pid2 = (int)p2[ID_COL];
	force->f12_x = x;
	force->f12_y = y;
}

void rev_force(double *f21_x, double *f21_y, double f12_x, double f12_y) {
	*f21_x = -f12_x;
	*f21_y = -f12_y;
}

void print_force_list(force_list_node *list, force_list_node *end_list) {
	force_list_node *p = list;
	double x, y;
	while(p != NULL && p != end_list) {
		x = p->f12_x; y = p->f12_y;
		LOG(("{\n"));
		LOG(("\tF%d, %d = %e,\n", p->pid1, p->pid2, sqrt(x*x+y*y)));
		LOG(("\tF_x = %e, F_y = %e\n", x, y));
		LOG(("}\n"));
		p = p->next;
	}
}