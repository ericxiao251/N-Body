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
		P[i][TYPE_COL] = (double)type;
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

void print_particle(double *p, double time) {
	LOG(("{\n"));
	LOG(("\tparticle_id: %f,\n", p[ID_COL]));
	LOG(("\tparticle_type: %f,\n", p[TYPE_COL]));
	LOG(("\tmass = %f,\n", p[WEIGHT_COL]));
	LOG(("\tposition: x = %f, ", p[POS_X_COL]));
	LOG(("y = %f,\n", p[POS_Y_COL]));
	LOG(("\tvelocity: x = %f, ", p[VOL_X_COL]));
	LOG(("y = %f, v = %f,\n", p[VOL_Y_COL], sqrt(p[VOL_X_COL] * p[VOL_X_COL] + p[VOL_Y_COL] * p[VOL_Y_COL])));
	LOG(("}\n"));

	LOG(("t=%f, p_id=%d, mass=%f, p_x=%f, p_y=%f, v_x=%f, v_y=%f\n",
			time, (int)p[ID_COL], p[WEIGHT_COL], p[POS_X_COL], p[POS_Y_COL],
			p[VOL_X_COL], p[VOL_Y_COL]));

	if (p[WEIGHT_COL] <= 0.0) {
		return;
	}

	LOG(("%lf,%d,%lf,%lf,%lf,%lf,%lf\n",
			time, (int)p[ID_COL], p[WEIGHT_COL], p[POS_X_COL], p[POS_Y_COL],
			p[VOL_X_COL], p[VOL_Y_COL]));

	return;
}

void print_all_particles(double **P, int numParticle, double time) {
	int i;
	LOG(("{\n"));
	LOG(("t,p_id,mass,p_x,p_y,v_x,v_y\n"));
	for (i = 0; i < numParticle; ++i) {
		print_particle(P[i], time);
	}
	LOG(("}\n"));
}

void particles_gen(double **P, int light_cnt, int medium_cnt, int heavy_cnt, int padding_cnt) {
	particles_gen_by_type(&P[0],                      LIGHT,  light_cnt);
	particles_gen_by_type(&P[light_cnt],              MEDIUM, medium_cnt);
	particles_gen_by_type(&P[light_cnt + medium_cnt], HEAVY,  heavy_cnt);
	particles_gen_by_type(&P[light_cnt + medium_cnt + heavy_cnt], DUMMY,  padding_cnt);
}

void update(unsigned char* image, double **P, double **P_force, int total_p_cnt,
	int img_width, int img_height, int step_size, int regenerate_img) {
	int i, cnt = 0;

	initilize_img(image, img_width, img_height);
	for (i = 0; i < total_p_cnt; ++i) {
		if (P[i][WEIGHT_COL] == DUMMY_WEIGHT) {
			// skip padding particles
			continue;
		}
		update_p(P[i], P_force[(int)(P[i][ID_COL])], step_size);
		// Only regenerate the whole image when necessary
		if (regenerate_img) {
			// int in_img_rng = update_img(image, P[i], img_width, img_height);
			int in_img_rng = update_img(image, P[i][POS_X_COL], P[i][POS_Y_COL], (int)P[i][TYPE_COL], img_width, img_height);
			if (in_img_rng) {
				// if the particle is still in image range, update cnt
				++cnt;
			}
		}
	}
	// print cnt for debug only
	//if (regenerate_img) {
	//	LOG(("UPDATE IMG: %d particles in range of the img.\n", cnt));
	//}
}
