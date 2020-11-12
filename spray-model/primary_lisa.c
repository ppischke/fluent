#include "math.h"
#include "math_ext.h"
#include "udf.h"
#include "dpm.h"
#include "model.h"

#include "primary_lisa.h"

/*  include secondary breakup models for combined models  */
#include "secondary_tab.h"
#include "secondary_kh.h"
#include "secondary_khtab.h"
#include "secondary_rdiwakar.h"





/*  LISA DISPERSION EQUATION, SHORT WAVE, VISCID  */
real omega_viscid(const real k, const real nu_liq, const real rho_gas, const real rho_liq, const real st_liq, const real u_sq)
{
    real C_nu = pow2(k) *  nu_liq;
    real C_v  = pow2(k) * rho_gas/rho_liq * u_sq;
    real C_st = pow3(k) *  st_liq/rho_liq;

    real C = 4.0*pow2(C_nu) + C_v - C_st;

    if(C_v <= C_st)
        return  0.0;
    else
        return -2.0*C_nu + sqrt(C);
}





void lisa_initialize(Tracked_Particle* parent)
{
    /*  already broken up  */
    if(parent->state.time != parent->time_of_birth)
        return;

    /*  states  */
    particle_state_t* part_state = &(parent->state);
    particle_state_t* prev_state = &(parent->state0);

    /*  sheet properties  */
    real m0 = part_state->mass;
    real hb = part_state->diam * 0.5;

    real u_sq = v_sqr(part_state->V);

    /*  liquid material properties  */
    real rho_liq
        = DPM_RHO
            (parent, P_MATERIAL(parent), P_T(parent));
    real st_liq
        = DPM_SURFTEN
            (parent);
    real mu_liq
        = DPM_MU
            (parent);
    real nu_liq
        = mu_liq/rho_liq;

    /*  gas (cell) material properties  */
    real rho_gas = C_R(P_CELL(parent),P_CELL_THREAD(parent))
                    * C_T(P_CELL(parent),P_CELL_THREAD(parent)) / parent->state.temp;

    /*  dimensionless numbers  */
    real We_gas = rho_gas * u_sq * hb / st_liq; 
    real Re_liq = rho_liq * sqrt(u_sq)/ mu_liq * 4.0*hb;


    /*  short wave assumption  */
    /*  inviscid               */
    /*  -- predictor --------  */
    real k;
    real k_s
        = 2.0/3.0 * We_gas / hb;
    real omega;
    real omega_s
        = omega_viscid(k_s, nu_liq, rho_gas, rho_liq, st_liq, u_sq);

    /*  short wave assumption  */
    /*  viscid                 */
    /*  -- corrector --------  */
    real delta_k = 0.1*k_s;
    
    while(fabs(delta_k/k_s) > M_ITER_ACCURACY)
    {   /*  iterative solution */
        omega = omega_viscid(k = k_s - delta_k, nu_liq, rho_gas, rho_liq, st_liq, u_sq);

        if(omega > omega_s)
        {   /*  update instable frequency and wave number  */
            k_s = k;
            omega_s = omega;
        }
        else
        {   /*  inverse and shorten step  */
            delta_k *= -0.5;
        }
    }

    
    /*  breakup: breakup time */
    real breakup_time   = 12.0/omega_s;
    real breakup_length = breakup_time*v_abs(part_state->V);
#ifdef LISA_FIXED_BREAKUP_LENGTH
    breakup_length = LISA_FIXED_BREAKUP_LENGTH;
#else
    if (breakup_length < LISA_MINIMUM_BREAKUP_LENGTH*2.0*hb)
        breakup_length = LISA_MINIMUM_BREAKUP_LENGTH*2.0*hb;
#endif  
    parent->user[LISA_BREAKUP_TIME]
        = breakup_length / v_abs(part_state->V);
    parent->user[LISA_BREAKUP_LENGTH]
        = breakup_length;

        
    /*  breakup: sheet -> ligament    */
    real d_32;
    real d_lig  = sqrt(16.0*hb/k_s);
    real Oh_lig = mu_liq / sqrt(rho_liq*st_liq*d_lig);

    
    /*  breakup: ligament -> droplet  */    
    d_32 = d_lig * 1.88*pow(1.0 + 3.0*Oh_lig, 1.0/6.0);
#ifdef LISA_FIXED_DIAM
    /*  simplified: fixed diameter    */
    d_32 = LISA_FIXED_DIAM;
#endif
#ifdef LISA_FIXED_DIAM_RATIO
    /*  simplified: fixed diam ratio  */
    d_32 = LISA_FIXED_DIAM_RATIO * 2.0*hb;
#endif


    /*  sample drop diam.  */
    real d = prev_state->diam = part_state->diam = rosin_r_diam(d_32, LISA_DIST_PARAMETER);
    /*  mass conservation  */
    real m = prev_state->mass = diam_to_mass(part_state);
    /*  par. conservation  */
    parent->number_in_parcel *= m0/m;


    const real I = 1.0 + 0.16/pow(Re_liq,0.125)*gauss_rand_lim(DISTRIB_SIGMA);
    /*  drop dispersion    */
    int i;
    for(i=0; i<ND_3; ++i)
    {
        prev_state->V[i] = part_state->V[i] *= I;
    }   
#ifdef LISA_FIXED_DISPERS_ANGLE 
    for(i=0; i<ND_3; ++i)
    {
        prev_state->V[i] = part_state->V[i] += LISA_FIXED_DISPERS_ANGLE*sqrt(u_sq)*gauss_rand_lim(SIGMA);
    }
#endif
    
    /*  update statistics  */
	parent->gvtp.n_shed++;	    
#ifdef VERBOSE_BREAKUP					
    Message("LISA breakup (hb, drr, d):   %e   %e   %e\n", hb, d_32, d);
#endif	
    return;
}





void lisa_breakup(Tracked_Particle* parent)
{
    /*  already broken up  */
    if(!is_sleeping(parent))
        return;

    /*  breakup time has not reached yet  */
    if(parent->state.time - parent->time_of_birth < parent->user[LISA_BREAKUP_TIME])
        return;

    /*  breakup  */
    inc_generation(parent);
    return;
}





/*  LISA-breakup UDF  */
DEFINE_DPM_SCALAR_UPDATE(lisa, cell, thread, initialize, parent)
{
    if(initialize)
    {
        if(parcel_initialize(parent))
            lisa_initialize(parent);
    }
    else
    {
        lisa_breakup(parent);

        if(!is_suspended(parent))
        {   /*  drop deformation  */
            tab_update(parent, NULL);
        }
    }
    return;
}





DEFINE_DPM_SCALAR_UPDATE(lisa_tab, cell, thread, initialize, parent)
{
    if(initialize)
    {
        if(parcel_initialize(parent))
        {   /*  primary breakup model  */
            lisa_initialize(parent);
            /*  secondary breakup models  */
            tab_initialize(parent);
        }
    }
    else
    {   /*  primary breakup model  */
        lisa_breakup(parent);

        if(!is_suspended(parent))
        {   /*  secondary breakup models  */
            tab_update(parent, NULL);
#ifndef LISA_TERTIARY
            if(is_primary(parent))
#endif      /*  no tertiary breakup   */                
            {   
                tab_breakup(parent, NULL, false);
            }
        }
    }
    return;
}





DEFINE_DPM_SCALAR_UPDATE(lisa_khtab, cell, thread, initialize, parent)
{
    if(initialize)
    {
        if(parcel_initialize(parent))
        {   /*  primary breakup model  */
            lisa_initialize(parent);
            /*  secondary breakup models  */
            khtab_initialize(parent);
        }
    }
    else
    {   /*  primary breakup model  */
        lisa_breakup(parent);

        if(!is_suspended(parent))
        {   /*  secondary breakup models  */
            tab_update(parent, NULL);
#ifndef LISA_TERTIARY
            if(is_primary(parent))
#endif      /*  no tertiary breakup   */
            {   
                kh_update (parent, NULL);
                kh_breakup(parent, NULL, false);
                tab_breakup(parent, NULL, false);
            }
        }
    }
    return;
}





DEFINE_DPM_SCALAR_UPDATE(lisa_rdiwakar, cell, thread, initialize, parent)
{
    if(initialize)
    {
        if(parcel_initialize(parent))
        {   /*  primary breakup model  */
            lisa_initialize(parent);
            /*  drop deformation  */
            tab_initialize(parent);
            /*  secondary breakup model  */
            rdiwakar_initialize(parent);
        }
    }
    else
    {   /*  primary breakup model  */
        lisa_breakup(parent);

        if(!is_suspended(parent))
        {   /*  drop deformation  */
            tab_update(parent, NULL);
            return;

            /*  secondary breakup model  */
            if(is_primary(parent))
            {   /*  no tertiary breakup  */
                rdiwakar_update(parent);
                rdiwakar_breakup(parent);
            }
        }
    }
    return;
}
