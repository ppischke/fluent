#ifndef __SECONDARY_RDIWAKAR_H
#define __SECONDARY_RDIWAKAR_H

#include "udf.h"
#include "dpm.h"
#include "math_ext.h"
#include "model.h"

#include "collision.h"



#define RDIWAKAR_C_B1      6.0    /*   6.0  */
#define RDIWAKAR_C_B2     M_PI    /*  M_PI  */
#define RDIWAKAR_C_S1      0.5    /*   0.5  */
#define RDIWAKAR_C_S2     20.0    /*  20.0  */

#define RDIWAKAR_DIST_TYPE           DIST_ROSIN_RAMMLER
#define RDIWAKAR_DIST_PARAMETER      DIST_PARAMETER



void rdiwakar_initialize(Tracked_Particle*);
void rdiwakar_initialize_p(Particle*);

void rdiwakar_update(Tracked_Particle*);
bool rdiwakar_breakup(Tracked_Particle*);

real rdiwakar_collision_outcome(Tracked_Particle*, Particle*, const real, real*, collision_outcome_t*);
void rdiwakar_collision_handler(Tracked_Particle*, Particle*, real*, const real, const real, const collision_outcome_t);



#endif
