#ifndef __DEFINE_H__
#define __DEFINE_H__

/*--------------------- printf() wrapper ----------------------*/
// Usage: LOG(("%s %d\n", "abc", 123));
#define ENABLE_LOG 1
#include "log.h"

/*--------------------- Include ----------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
//#include "vector3d.h"
#include "physics.h"
#include "image.h"
#include "properties.h"
#include "util.h"
#include "cyclic_dist.h"

/*--------------------- Data structure ----------------------*/
// Particle
#define PARTICLE_PROPERTIES_COUNT 7 // p_id, p_weight, p_pos_x, p_pos_y, p_vol_x, p_vol_y, p_type
#define ID_COL 0
#define WEIGHT_COL 1
#define POS_X_COL 2
#define POS_Y_COL 3
#define VOL_X_COL 4
#define VOL_Y_COL 5
#define TYPE_COL 6

// Force
#define FORCE_PROPERTIES_COUNT 4 // force_x, force_y, force_to(p_id), force_from(p_id)
#define FORCE_SUM_PROPERTIES_COUNT 2 // force_x, force_y
#define X_COL 0
#define Y_COL 1
#define TO_COL 2
#define FROM_COL 3



/*--------------------- Convension ----------------------*/
#define LIGHT 0
#define MEDIUM 1
#define HEAVY 2
#define DUMMY 3
#define DUMMY_WEIGHT 0.0
#define MASTER_TO_SLAVE_TAG 883216 // arbitrary tag for messages sent from master to slaves
#define SLAVE_TO_MASTER_TAG 886123 // arbitrary tag for messages sent from slaves to master
#define SLAVE_TO_SLAVE_TAG 666123 // arbitrary tag for messages sent from slaves to master

/*--------------------- Constants ----------------------*/
#define POS_MAX_X 1024.0
#define POS_MIN_X 0.0
#define POS_MAX_Y 1024.0
#define POS_MIN_Y 0.0

#define X_RNG (POS_MAX_X-POS_MIN_X)
#define Y_RNG (POS_MAX_Y-POS_MIN_Y)

#define LIGHT_R 0
#define LIGHT_G 0
#define LIGHT_B 255
#define MEDIUM_R 0
#define MEDIUM_G 255
#define MEDIUM_B 0
#define HEAVY_R 255
#define HEAVY_G 0
#define HEAVY_B 0

#endif