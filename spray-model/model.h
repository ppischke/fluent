#ifndef __MODEL_H
#define __MODEL_H

#include "udf.h"
#include "dpm.h"
#include "math_ext.h"



#define MODEL_CONF                      "model.conf"



/*  model user variables  */
#define LISA_BREAKUP_TIME                0
#define LISA_BREAKUP_LENGTH              1
#define RT_PROGRESS                     LISA_BREAKUP_TIME
#define RT_WAVELENGTH                   LISA_BREAKUP_LENGTH
#define KH_ACCUMULATED_MASS              2
#define KH_ACCUMULATED_D2                3
#define KH_ACCUMULATED_D3                4
#define KH_LIMIT                         5
#define TB_ACCUMULATED_MASS             KH_ACCUMULATED_MASS
#define TB_ACCUMULATED_D2               KH_ACCUMULATED_D2
#define TB_ACCUMULATED_D3               KH_ACCUMULATED_D3
#define TB_LIMIT						KH_LIMIT
#define TAB_DISPLACEMENT                 6
#define TAB_VELOCITY                     7
#define TAB_EQUILIBRIUM                  8
#define RDIWAKAR_DIAMETER                9
#define RDIWAKAR_TIMESCALE              10
#define N_MODEL_VARS                    12
/*  other user variables  */
#define PARCEL_GENERATION               12
#define PARCEL_DIAMETER                 13
#define PARCEL_DIAMETER_REDUCED         14
#define PARCEL_NEIGHBORS                15
#define PARCEL_WEIGHT                   16
#define PARCEL_BESSEL                   17
#define PARCEL_JACOBIAN                 18
#define PARCEL_EIGENVECTORS             27
#define PARCEL_EIGENVALUES              36
#define FLOW_VISCOSITY					39
#define FLOW_DENSITY					40
#define FLOW_PRESSURE					41
#define N_USER_VARS                     42

/*  size limitations  */
#define BREAKUP_MIN_DIAM                0.00
#define BREAKUP_MIN_DIAM_RATIO          0.00
#define BREAKUP_MAX_MASS_RATIO          0.95

#define DIST_PARAMETER                  2.5



/*  parcel state  */
#define P_ID(particle)                  ((particle)->part_id)
#define P_CELL_ID(particle)             ((particle)->cCell.ct.c)
#define P_TYPE(particle)                ((particle)->type)
#define P_USER(particle)                ((particle)->user)

#define P_GENERATION(particle)          ((particle)->user[PARCEL_GENERATION])
#define P_INC_GENERATION(particle)      ((int)(P_GENERATION(particle) += 1))

#define P_IS_NEW(particle)              ((P_TIME(particle) == P_INIT_TIME(particle)))
#define P_IS_SLEEPING(particle)         (((int)P_GENERATION(particle)) <  0)
#define P_IS_PRIMARY(particle)          (((int)P_GENERATION(particle)) == 0)
#define P_IS_SECONDARY(particle)        (((int)P_GENERATION(particle)) >  0)

#define P_COLLISION_PERIODIC(particle)  ((particle)->wallfilm_thread_id)

#define P_PDRED(particle)               (P_USER(particle)[PARCEL_DIAMETER_REDUCED])
#define P_PDIAM(particle)               (P_USER(particle)[PARCEL_DIAMETER])
#define P_NEIGHBORS(particle)           (P_USER(particle)[PARCEL_NEIGHBORS])
#define P_WMEAN(particle)               (P_USER(particle)[PARCEL_WEIGHT])
#define P_BESSEL(particle)              (P_USER(particle)[PARCEL_BESSEL])
#define P_JACOBIAN(particle)            (P_USER(particle)+PARCEL_JACOBIAN)
#define P_EIGENVECTORS(particle)        (P_USER(particle)+PARCEL_EIGENVECTORS)
#define P_EIGENVALUES(particle)         (P_USER(particle)+PARCEL_EIGENVALUES)

#define P_FLOW_VISCOSITY(particle)		(P_USER(particle)[FLOW_VISCOSITY])
#define P_FLOW_DENSITY(particle)		(P_USER(particle)[FLOW_DENSITY])
#define P_FLOW_PRESSURE(particle)		(P_USER(particle)[FLOW_PRESSURE])



/*  mass/diameter conversion  */
#define P_MASS_TO_DIAM(particle)        (P_DIAM(particle)  = pow(6.0/M_PI * P_MASS(particle) /P_RHO(particle),  1.0/3.0))
#define P_MASS_TO_DIAM0(particle)       (P_DIAM0(particle) = pow(6.0/M_PI * P_MASS0(particle)/P_RHO0(particle), 1.0/3.0))

#define P_DIAM_TO_MASS(particle)        (P_MASS(particle)  = M_PI/6.0 * P_RHO(particle) *pow3(P_DIAM(particle)))
#define P_DIAM_TO_MASS0(particle)       (P_MASS0(particle) = M_PI/6.0 * P_RHO0(particle)*pow3(P_DIAM0(particle)))

/*  copy states  */
#define P_STATE(p)						(&((p)->state))
#define P_STATE0(p)						(&((p)->state0))
#define P_INIT_STATE(p)					(&((p)->init_state))

#define P_STATE_TO_STATE(src,tgt)       (memcpy(&((tgt)->state),  &((src)->state),  sizeof(particle_state_t)))
#define P_STATE_TO_STATE0(src,tgt)      (memcpy(&((tgt)->state0), &((src)->state),  sizeof(particle_state_t)))

#define P_STATE0_TO_STATE(src,tgt)      (memcpy(&((tgt)->state),  &((src)->state0), sizeof(particle_state_t)))
#define P_STATE0_TO_STATE0(src,tgt)     (memcpy(&((tgt)->state0), &((src)->state0), sizeof(particle_state_t)))

#define P_COPY_CELL_ID(src,tgt)         (memcpy(P_CCELL(tgt), P_CCELL(src), sizeof(CX_Cell_Id)))



void reset_variables(Tracked_Particle*);
void reset_variables_p(Particle*);

void reset_models(Tracked_Particle*);
void reset_models_p(Particle*);

bool parcel_initialize(Tracked_Particle*);
bool parcel_initialize_p(Particle*);

void parcel_suspend(Tracked_Particle*, const real, const real);
void parcel_suspend_p(Particle*, const real, const real);

int inc_generation(Tracked_Particle*);
int inc_generation_p(Particle*);

Particle* create_product_parcel(Tracked_Particle*, const real, const real, const real, const real);



/*  QUERY STATUS  */
bool is_new(const Tracked_Particle*);
bool is_new_p(const Particle*);

bool is_sleeping(const Tracked_Particle*);
bool is_sleeping_p(const Particle*);

bool is_primary(const Tracked_Particle*);
bool is_primary_p(const Particle*);

bool is_secondary(const Tracked_Particle*);
bool is_secondary_p(const Particle*);

bool is_tertiary(const Tracked_Particle*);
bool is_tertiary_p(const Particle*);

/*
bool is_suspended(const Tracked_Particle*);
bool is_suspended_p(const Particle*);
*/
#define is_suspended    is_sleeping
#define is_suspended_p  is_sleeping_p

    
    
/*  MASS TO DIAMETER  */
real diam_to_mass(particle_state_t*);
real mass_to_diam(particle_state_t*);



#endif
