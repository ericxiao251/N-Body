#include "define.h"

void location_gen(double *x, double *y, int type) {
	*x = (type != DUMMY) ? POS_MIN_X + (POS_MAX_X - POS_MIN_X) * drand48() : 0.0;
	*y = (type != DUMMY) ? POS_MIN_Y + (POS_MAX_Y - POS_MIN_Y) * drand48() : 0.0;
}

void velocity_gen(double *x, double *y, int type) {
	double max_v, min_v, rand_v, rand_deg;
	if (type == DUMMY) {
		*x = 0.0;
		*y = 0.0;
		return;
	}

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
		 (type == HEAVY) ? massHeavyMin + (massHeavyMax - massHeavyMin) * drand48(): DUMMY_WEIGHT;
}

void particles_gen_by_type(double **P, int type, int cnt) {
	int i;
	static double id = 0.0;
	for (i = 0; i < cnt; ++i) {
		P[i][ID_COL] = id++;
		location_gen(&P[i][POS_X_COL], &P[i][POS_Y_COL], type);
		velocity_gen(&P[i][VOL_X_COL], &P[i][VOL_Y_COL], type);
		weight_gen(&P[i][WEIGHT_COL], type);
	}
}

/*------------------------- Global Functions -----------------------------------*/
void print_properties_h(void) {
	LOG(("Light Particle:\n"));
	LOG(("  %f <= velocity <= %f\n", velocityLightMin, velocityLightMax));
	LOG(("  %f <=   mass   <= %f\n", massLightMin, massLightMax));
	LOG(("Medium Particle:\n"));
	LOG(("  %f <= velocity <= %f\n", velocityMediumMin, velocityMediumMax));
	LOG(("  %f <=   mass   <= %f\n", massMediumMin, massMediumMax));
	LOG(("Heavy Particle:\n"));
	LOG(("  %f <= velocity <= %f\n", velocityHeavyMin, velocityHeavyMax));
	LOG(("  %f <=   mass   <= %f\n", massHeavyMin, massHeavyMax));
}

void print_particle(double *p) {
	LOG(("{\n"));
	LOG(("\tparticle_id: %f,\n", p[ID_COL]));
	LOG(("\tmass = %f,\n", p[WEIGHT_COL]));
	LOG(("\tposition: x = %f, ", p[POS_X_COL]));
	LOG(("y = %f,\n", p[POS_Y_COL]));
	LOG(("\tvelocity: x = %f, ", p[VOL_X_COL]));
	LOG(("y = %f, v = %f,\n", p[VOL_Y_COL], sqrt(p[VOL_X_COL] * p[VOL_X_COL] + p[VOL_Y_COL] * p[VOL_Y_COL])));
	LOG(("}\n"));
}

void print_all_particles(double **P, int numParticle) {
	int i;
	for (i = 0; i < numParticle; ++i) {
		print_particle(P[i]);
	}
}

void particles_gen(double **P, int light_cnt, int medium_cnt, int heavy_cnt, int padding_cnt) {
	particles_gen_by_type(&P[0],                      LIGHT,  light_cnt);
	particles_gen_by_type(&P[light_cnt],              MEDIUM, medium_cnt);
	particles_gen_by_type(&P[light_cnt + medium_cnt], HEAVY,  heavy_cnt);
	particles_gen_by_type(&P[light_cnt + medium_cnt + heavy_cnt], DUMMY,  padding_cnt);
}

void update(unsigned char* image, double **P, double **P_force, int total_p_cnt, int img_width, int img_height, double step_size) {
	int i, cnt = 0;
	
	initilize_img(image, img_width, img_height);
	for (i = 0; i < total_p_cnt; ++i) {
		if (P[i][WEIGHT_COL] == DUMMY_WEIGHT) {
			// skip padding particles
			continue;
		}
		update_p(P[i], P_force[(int)(P[i][ID_COL])], step_size);
		if (update_img(image, P[i], img_width, img_height) > 0) {
			++cnt;
		}
	}
	LOG(("UPDATE IMG: %d particles in range of the img.\n", cnt));
}