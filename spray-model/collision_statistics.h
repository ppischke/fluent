#ifndef __COLLISION_STATISTICS_H
#define __COLLISION_STATISTICS_H



#include "udf.h"
#include "dpm.h"
#include "time.h"
#include "model.h"
#include "collision.h"



typedef struct collision_statistics_struct
{
    real    volume;
    int     n_parcels;
    int     n_pairs;
    int     n_tries;
    int     n_readjust;
    real    collisions;
    real    c_regimes[COLLISION_REGIMES];
    real    momentum;
    real    m_regimes[COLLISION_REGIMES];
    real    dissipation;
    real    d_regimes[COLLISION_REGIMES];
    clock_t execution_time;
}
collision_statistics_t;



void statistics_reset();
void statistics_add_target(const int, const int);
void statistics_add_event(const real, const collision_outcome_t, const real momentum, const real dissipation, const real We_coll, const real B_coll, const real d_coll);
void statistics_add_tries(const int, const int, const int);
void statistics_add_readjust();
void statistics_print();



#endif
