#ifndef __DEFINE_H__
#define __DEFINE_H__

// printf() wrapper
// Usage: LOG(("%s %d\n", "abc", 123));
#define ENABLE_LOG 1
#include "log.h"

// Include
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include "vector3d.h"
#include "savebmp.h"
//#include "properties.h"
#include "util.h"
#include "cyclic_dist.h"

// Data structure
#define PROPERTIES_COUNT 5 // p_weight, p_pos_x, p_pos_y, p_vol_x, p_vol_y
#define WEIGHT_COL 0
#define POS_X_COL 1
#define POS_Y_COL 2
#define VOL_X_COL 3
#define VOL_Y_COL 4

// Convension
#define LIGHT 0
#define MEDIUM 1
#define HEAVY 2
#define MASTER_TO_SLAVE_TAG 883216 // arbitrary tag for messages sent from master to slaves
#define SLAVE_TO_MASTER_TAG 886123 // arbitrary tag for messages sent from slaves to master

// Constants
#define POS_MAX_X 2550.0
#define POS_MIN_X 0.0
#define POS_MAX_Y 1024.0
#define POS_MIN_Y 0.0
#define epsilon 0.000000000000000222
#define TWO_PI 6.28318530718

#endif