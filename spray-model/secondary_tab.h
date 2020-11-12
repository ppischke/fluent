#ifndef __SECONDARY_TAB_H
#define __SECONDARY_TAB_H

/*  struct prototypes  */
struct tab_constants_struct; typedef struct tab_constants_struct tab_constants_t;

#include "udf.h"
#include "dpm.h"
#include "math_ext.h"

#include "collision.h"
#include "secondary_kh.h"



#define TAB_Ck                           8.0
#define TAB_Cb                           0.5
#define TAB_CF                          (1.0/3.0)
#define TAB_Cd                           5.0
#define TAB_K                          (10.0/3.0)
/*
#define TAB_Cv                           0.7
*/

#define TAB_N_PARCELS            2
#define TAB_WE_MIN              12.0
#define TAB_WE_MAX              40.0

#define TAB_DIST_TYPE           DIST_ROSIN_R
#define TAB_DIST_PARAMETER      DIST_PARAMETER

#define TAB_DRAG_SCALAR 		1.489

/*
#define TAB_COLLIDE
#define TAB_COALESCE
*/



struct tab_constants_struct
{   /*  limiting Weber numbers  */
    real We_min;
    real We_max;
    /*  number of child parcels  */
    int  number_of_parcels;
    /*  distribution parameters  */
    distribution_type_t
         distribution_type;
    real distribution_q;
};



void tab_initialize(Tracked_Particle*);
void tab_initialize_p(Particle*);

void tab_update(Tracked_Particle*, tab_constants_t*);
bool tab_breakup(Tracked_Particle*, tab_constants_t*, kh_constants_t*);

real tab_collision_outcome(Tracked_Particle*, Particle*, const real, real*, collision_outcome_t*);
void tab_collision_handler(Tracked_Particle*, Particle*, real*, const real, const real, const collision_outcome_t);



#endif
