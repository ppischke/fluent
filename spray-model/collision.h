#ifndef __COLLISION_H
#define __COLLISION_H



#include "udf.h"
#include "dpm.h"
#include "model.h"


/*  data logging  */
#define COLLISION_LOG "collision.log"
/*
#define COLLISION_EXT "collision-%8.6f.log"
*/


/*  parcel diameter limiter  */
#define PARCEL_DIAMETER_TOL                     1.0e-4
#define PARCEL_DIAMETER_ITER                    20
#define PARCEL_DIAMETER_INIT_VFR                1.0e-3


/*  collision partner handling  */
#define COLLISION_IGNORE_SLEEPING
#define COLLISION_IGNORE_PROBABILITY            1.0e-6

#define COLLISION_EFFICIENCY
#define COLLISION_VISCOSITY

#define COLLISION_FRAGMENT_UNIFORMITY           0.0
#define COLLISION_FRAGMENT_NUMBER               1.0


/*  disable sub-algorithms  */
/*
#define COLLISION_DISABLE_GRADIENT
#define COLLISION_DISABLE_BESSEL
#define COLLISION_DISABLE_VOI
*/


typedef enum collision_outcome_enum
{
    COLLISION_OUTCOME_COALESCENCE = 0,
    COLLISION_OUTCOME_STRETCHING,
    COLLISION_OUTCOME_REFLEXIVE,
    COLLISION_OUTCOME_BOUNCING,
    /*  regime for breakup coupling  */
    COLLISION_OUTCOME_COUPLED,
    COLLISION_REGIMES
}
collision_outcome_t;



typedef enum collision_shadow_enum
{
    ORIGINAL_POSITION = 0,
    PERIODIC_SHADOW_HIGH,
    PERIODIC_SHADOW_LOW
}
collision_shadow_t;





real collision_model
    (Particle*, Particle*, real);
real droplet_collide
    (Tracked_Particle* collector, Particle* target, real, const real, const real, const real, const real, const collision_outcome_t);
real droplet_coalesce
    (Tracked_Particle* collector, Particle* target, real, const real, const real, const real);
    
void parcel_diameter_estimator
    (Particle** ptrs_minor, Particle** ptrs_major, const int n_minor, const int n_major);
    
real parcel_diameter_initialize
    (Tracked_Particle* p);
real parcel_diameter_initialize_p
    (Particle*);

bool allow_collision
    (const Particle*);
real collision_probability
    (const Particle*, const Particle*, const real P_max);
real collision_efficiency
    (const Particle*, const Particle*);
void collision_algorithm
    (Particle** ptrs_minor, Particle** ptrs_major, const int n_minor, const int n_major);

int linked_list_to_ptrs(Particle* linked_list, Particle** ptrs, const int n_max);
void create_tp_from_p(Tracked_Particle* tp, Particle* p);
void revert_tp_to_p(Tracked_Particle* tp, Particle* p);





#endif
