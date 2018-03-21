#include "define.h"
#include "properties.h"

void location_gen(double *x, double *y) {
	*x = POS_MIN_X + (POS_MAX_X - POS_MIN_X) * drand48();
	*y = POS_MIN_Y + (POS_MAX_Y - POS_MIN_Y) * drand48();
}

void velocity_gen(double *x, double *y, int type) {
	double max_v, min_v, rand_v, rand_deg;

	max_v = (type == LIGHT) ? velocityLightMin : 
					(type == MEDIUM) ? velocityMediumMin :
					(type == HEAVY) ? velocityHeavyMin : 0;
	min_v = (type == LIGHT) ? velocityLightMax :
					(type == MEDIUM) ? velocityMediumMax :
					(type == HEAVY) ? velocityHeavyMax : 0;

	rand_v = min_v + (max_v - min_v) * drand48();
	rand_deg = TWO_PI * drand48();
	*x = rand_v * cos(rand_deg);
	*y = rand_v * sin(rand_deg);
}

void weight_gen(double *w, int type) {
	*w = (type == LIGHT) ? massLightMin + (massLightMax - massLightMin) * drand48() :
		 (type == MEDIUM) ? massMediumMin + (massMediumMax - massMediumMin) * drand48() :
		 (type == HEAVY) ? massHeavyMin + (massHeavyMax - massHeavyMin) * drand48(): 0;
}

void particles_gen_by_type(double **P, int type, int cnt) {
	int i;
	for (i = 0; i < cnt; ++i) {
		location_gen(&P[i][POS_X_COL], &P[i][POS_Y_COL]);
		velocity_gen(&P[i][VOL_X_COL], &P[i][VOL_Y_COL], type);
		weight_gen(&P[i][WEIGHT_COL], type);
	}
}

/*------------------------- Global Functions -----------------------------------*/
void print_properties_h(void) {
	printf("Light Particle:\n");
	printf("  %f <= velocity <= %f\n", velocityLightMin, velocityLightMax);
	printf("  %f <=   mass   <= %f\n", massLightMin, massLightMax);
	printf("Medium Particle:\n");
	printf("  %f <= velocity <= %f\n", velocityMediumMin, velocityMediumMax);
	printf("  %f <=   mass   <= %f\n", massMediumMin, massMediumMax);
	printf("Heavy Particle:\n");
	printf("  %f <= velocity <= %f\n", velocityHeavyMin, velocityHeavyMax);
	printf("  %f <=   mass   <= %f\n", massHeavyMin, massHeavyMax);
}

void print_particle(double *p, int id) {
	printf("{\n");
	printf("\tparticle_id: %d,\n", id);
	printf("\tmass = %f,\n", p[WEIGHT_COL]);
	printf("\tposition: x = %f, ", p[POS_X_COL]);
	printf("y = %f,\n", p[POS_Y_COL]);
	printf("\tvelocity: x = %f, ", p[VOL_X_COL]);
	printf("y = %f,\n", p[VOL_Y_COL]);
	printf("}\n");
}

void print_all_particles(double **P, int numParticle) {
	int i;
	for (i = 0; i < numParticle; ++i) {
		print_particle(P[i], i);
	}
}

void particles_gen(double **P, int light_cnt, int medium_cnt, int heavy_cnt) {
	particles_gen_by_type(&P[0],                      LIGHT,  light_cnt);
	particles_gen_by_type(&P[light_cnt],              MEDIUM, medium_cnt);
	particles_gen_by_type(&P[light_cnt + medium_cnt], HEAVY,  heavy_cnt);
}