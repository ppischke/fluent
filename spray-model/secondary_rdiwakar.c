#define __SECONDARY_RDIWAKAR_C

#include "udf.h"
#include "dpm.h"
#include "math.h"
#include "math_ext.h"
#include "model.h"

#include "secondary_rdiwakar.h"





void rdiwakar_init(real* user)
{
    /*  initialize user scalars  */
    user[RDIWAKAR_DIAMETER]  = 0.0;
    user[RDIWAKAR_TIMESCALE] = 0.0;

    return;
}

void rdiwakar_initialize(Tracked_Particle* parent)
{
    rdiwakar_init(parent->user);
    return;
}

void rdiwakar_initialize_p(Particle* parent)
{
    rdiwakar_init(parent->user);
    return;
}





void rdiwakar_update(Tracked_Particle* parent)
{   /*  states  */
    cphase_state_t* gas_state = &(parent->cphase);

    /*  liquid material properties  */
    const real rho_liq
        = DPM_RHO
            (parent, P_MATERIAL(parent), P_T(parent));
    const real st_liq
        = DPM_SURFTEN
            (parent);

    /*  gas (cell) material properties  */
    const real rho_gas
        = gas_state->rho;
    const real mu_gas
        = gas_state->mu;


    /*  parent state  */
    const real d0
        = P_DIAM(parent);
    const real r0
        = P_DIAM(parent)/2.0;
    const real u_rel_sqr
        = v_rel_sqr(P_VEL(parent), gas_state->V);


    /*  dimensionless numbers  */
    const real We
        = rho_gas * u_rel_sqr * r0 / st_liq;
    const real Re
        = rho_gas * sqrt(u_rel_sqr) * 2.0*r0 / mu_gas;


    /*  bag regime boundary  */
    real rBag = We/RDIWAKAR_C_B1;
    /*  stripping regime boundary  */
    real rStripping = pow2(We/sqrt(Re)/RDIWAKAR_C_S1);


    if(P_USER(parent)[RDIWAKAR_DIAMETER] == 0.0)
    {
        if(rStripping > 1.0)
        {
            /*  stripping regime stable diameter  */
            P_USER(parent)[RDIWAKAR_DIAMETER] = rosin_r_diam(d0/rStripping, RDIWAKAR_DIST_PARAMETER);

            /*  stripping regime time-scale  */
            P_USER(parent)[RDIWAKAR_TIMESCALE] = RDIWAKAR_C_S2*r0*sqrt(rho_liq/rho_gas/u_rel_sqr);
        }
    else
        if(rBag > 1.0)
        {
            /*  bag regime stable diameter  */
            P_USER(parent)[RDIWAKAR_DIAMETER] = rosin_r_diam(d0/rBag, RDIWAKAR_DIST_PARAMETER);

            /*  bag regime time-scale  */
            P_USER(parent)[RDIWAKAR_TIMESCALE] = RDIWAKAR_C_B2*sqrt(rho_liq*POW3(r0)/2.0/st_liq);
        }
    }
    else
    {
        if((rStripping < 1.0) && (rBag < 1.0))
        {   /*  stable  */
            P_USER(parent)[RDIWAKAR_DIAMETER] = 0.0;
            P_USER(parent)[RDIWAKAR_TIMESCALE] = 0.0;
        }
    }
    return;
}





bool rdiwakar_breakup(Tracked_Particle* parent)
{
    real dst = P_USER(parent)[RDIWAKAR_DIAMETER];
    real tau = P_USER(parent)[RDIWAKAR_TIMESCALE];

    real d0 = P_DIAM(parent);
    real m0 = P_MASS(parent);

    if((dst > 0.0) && (dst < d0))
    {
        /*  diameter  */
        real d = dst + (d0-dst)*exp(-P_DT(parent)/tau);
        /*  mass conservation  */
        real m = m0 * POW3(d/d0);

        P_DIAM(parent) = P_DIAM0(parent) = d;
        P_MASS(parent) = P_MASS0(parent) = m;

        /*  number in parcel  */
        P_N(parent) *= m0/m;

        return true;
    }
    return false;
}



real rdiwakar_collision_outcome(Tracked_Particle* collector, Particle* target, const real P_red, real* impact_parameter, collision_outcome_t* collision_outcome)
{
    return 0.0;
}



void rdiwakar_collision_handler(Tracked_Particle* collector, Particle* target, real* user_coalesced, const real P_red, const real impact_parameter, const collision_outcome_t collision_outcome)
{
    switch(collision_outcome)
    {
    case COLLISION_OUTCOME_COALESCENCE:
#ifdef RDIWAKAR_COALESCE    
        user_coalesced[RDIWAKAR_DIAMETER]  = 0.0;
        user_coalesced[RDIWAKAR_TIMESCALE] = 0.0;
#endif      
        break;

    default:
        break;
    }
    return;
}
