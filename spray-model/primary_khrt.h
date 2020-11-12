#ifndef __PRIMARY_KHRT_H
#define __PRIMARY_KHRT_H

#include "udf.h"
#include "dpm.h"
#include "model.h"
#include "collision.h"



/*  weber number limits
//  stripping and bag regimes
//  [comparisons of diesel spray liquid penetration and vapor fuel distributions with in-cylinder optical measurements]
//  [baumgarten,  tab. 4.1.]
*/
#define KHRT_ABSOLUTE_VELOCITY

#define KHRT_CONE_ANGLE      0.09
#define KHRT_SCTR_ANGLE      0.09
#define KHRT_DRAG_COEFF      0.47

#define RT_CD                0.2        /*  diameter constant:   def.  0.3  */      
#define RT_CT                1.0        /*  time-scale const.:   0.2...1.0  */
#define RT_CL                7.0        /*  liquid core length:  max. 14.9  */
#define RT_CE                1.0
#define RT_L0                0.0

#define RT_CUTOFF			 0.05

#define RT_DIST_TYPE        DIST_ROSIN_R
#define RT_DIST_PARAMETER   DIST_PARAMETER





void rt_initialize(Tracked_Particle*);
void rt_initialize_p(Particle*);

void rt_update (Tracked_Particle*);
bool rt_breakup(Tracked_Particle*);

/*
real rt_collision_outcome(Tracked_Particle*, Particle*, real*, collision_outcome_t*);
void rt_collision_handler(Tracked_Particle*, Particle*, real*, const real, const collision_outcome_t);
*/



#endif
