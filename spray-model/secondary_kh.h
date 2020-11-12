#ifndef __SECONDARY_KH_H
#define __SECONDARY_KH_H

/*  struct prototypes  */
struct kh_constants_struct; typedef struct kh_constants_struct kh_constants_t;

#include "udf.h"
#include "dpm.h"
#include "math_ext.h"

#include "collision.h"



/*  weber number limits
//  stripping and bag regimes
//  [comparisons of diesel spray liquid penetration and vapor fuel distributions with in-cylinder optical measurements]
//  [baumgarten,  tab. 4.1.]
*/
#define KH_WE_MIN             80.0
#define KH_WE_MAX              1.0e7
#define KH_WE_FORCE            0.0

#define KH_CUTOFF              0.05

#define KH_B0                  0.61
#define KH_B1                 40.00

#define KH_DIST_TYPE        DIST_ROSIN_R
#define KH_DIST_PARAMETER   DIST_PARAMETER



struct kh_constants_struct
{   /*  weber number limits  */
    real We_min;
    real We_max;
    real We_force_release;
    /*  limitating mass fractions  */
    real cutoff_min;
    real cutoff_max;
    /*  model constants   */
    real B0;
    real B1;
    /*  distribution parameters  */
    distribution_type_t
         distribution_type;
    real distribution_q;
};



void kh_initialize(Tracked_Particle*);
void kh_initialize_p(Particle*);

void kh_update(Tracked_Particle*, kh_constants_t*);
bool kh_breakup(Tracked_Particle*, kh_constants_t*, bool force);

real kh_collision_outcome(Tracked_Particle*, Particle*, const real, real*, collision_outcome_t*);
void kh_collision_handler(Tracked_Particle*, Particle*, real*, const real, const real, const collision_outcome_t);



#endif
