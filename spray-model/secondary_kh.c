#include "udf.h"
#include "dpm.h"
#include "math.h"
#include "math_ext.h"
#include "model.h"

#include "secondary_kh.h"
#include "primary_khrt.h"





/*  default kh constants  */
kh_constants_t kh_default_constants =
{
    KH_WE_MIN,
    KH_WE_MAX,
    KH_WE_FORCE,
    KH_CUTOFF,
    KH_CUTOFF,
    KH_B0,
    KH_B1,
    KH_DIST_TYPE,
    KH_DIST_PARAMETER
};





void kh_init(real* user)
{
    /*  initialize user scalars  */	
    user[KH_ACCUMULATED_MASS] = 0.0;
    user[KH_ACCUMULATED_D3]   = 0.0;
    user[KH_ACCUMULATED_D2]   = 0.0;  
    user[KH_LIMIT]            = 0.5+uniform_rand();
    return;
}

void kh_initialize(Tracked_Particle* parent)
{
    kh_init(parent->user);
    return;
}

void kh_initialize_p(Particle* parent)
{
    kh_init(parent->user);
    return;
}





void kh_update(Tracked_Particle* parent, kh_constants_t* KH)
{
	/*  states  */
    cphase_state_t*   gas_state  = &(parent->cphase);
    particle_state_t* part_state = &(parent->state);
	particle_state_t* init_state = &(parent->init_state); 
        
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

    /*  gas (cell) material properties  */
    real rho_gas;
    if(is_sleeping(parent))
        rho_gas = C_R(P_CELL(parent),P_CELL_THREAD(parent))
                    * C_T(P_CELL(parent),P_CELL_THREAD(parent)) / parent->state.temp;
    else
        rho_gas = gas_state->rho;
    
    real u_rel_sqr;
#ifdef KHRT_ABSOLUTE_VELOCITY   
    /*  primary breakup    */
    if(is_sleeping(parent))
    {   /*  absolute velocity for primary breakup   */
        u_rel_sqr = v_sqr(init_state->V);        
    }
    else
#endif
    {   /*  relative velocity for secondary breakup */
        u_rel_sqr = v_rel_sqr(part_state->V, gas_state->V);
    }

	/*  parent state  */	
    real dt = parent->time_step;

    real m0 = part_state->mass - parent->user[KH_ACCUMULATED_MASS];    
    real r0 = pow(3.0/4.0/M_PI * m0/rho_liq, 1.0/3.0);

    /*  dimensionless numbers  */
    real We_gas = rho_gas *    u_rel_sqr    * r0 / st_liq;
    real We_liq = rho_liq *    u_rel_sqr    * r0 / st_liq;
    real Re_liq = rho_liq * sqrt(u_rel_sqr) * r0 / mu_liq;

    real We_d = 2.0*We_gas;

    /*  outside kh regime?  */
    if(We_d < KH_WE_MIN
    || We_d > KH_WE_MAX)
    {
        return;
    }

    real Oh = sqrt(We_liq) / Re_liq;    /*  Ohnesorge  */
    real Ta = sqrt(We_gas) * Oh;        /*  Taylor     */

    /*  dimensionless frequency   [Baumgarten, eq. 4.121]  */
    real omega = (0.34 + 0.38*pow(We_gas,1.5)) / (1.0 + Oh) / (1.0 + 1.4*pow(Ta,0.6)) * sqrt(st_liq / rho_liq / pow3(r0));
    /*  most unstable wavelength  [Baumgarten, eq. 4.122]  */
    real lambda = 9.02 * (1.0 + 0.45*sqrt(Oh)) * (1.0 + 0.4*pow(Ta,0.7)) / pow((1.0 + 0.865*pow(We_gas,1.67)), 0.6);
    
    /*  breakup time  */
    real tau = 3.788 * KH_B1/lambda/omega / (1.0 - KH_B0*lambda);
    if(tau <= 0.0)
        /*  parent must not grow  */
        return;

    /*  update parent state  */
    real r = r0 * exp(-dt/tau);
    real m = 4.0/3.0 * M_PI * rho_liq * pow3(r);
    if(r > r0)
        /*  parent must not grow  */
        return;

    /*  calculate sauter mean diameter  */
    real r2 = (pow2(r0) - pow2(r)) * 1.50/KH_B0/lambda;
    real r3 = (pow3(r0) - pow3(r));
        
    /*  accumulate user scalars  */	
    parent->user[KH_ACCUMULATED_MASS] += (m0-m);
    parent->user[KH_ACCUMULATED_D2]   += 4.0*r2;
    parent->user[KH_ACCUMULATED_D3]   += 8.0*r3;
    return;
}





bool kh_breakup(Tracked_Particle* parent, kh_constants_t* KH, bool force)
{    
    /*  sample cut-off mass   */
	real m = parent->init_state.mass*KH_CUTOFF*parent->user[KH_LIMIT];

    /*  sample mean diameter  */
    real d;
	real d0  = parent->state.diam;
    real smd = parent->user[KH_ACCUMULATED_D3]/parent->user[KH_ACCUMULATED_D2];   

    /*  scattering velocity  */
    real u_s = 0.0;
    real x_s = 0.0;
#ifdef KHRT_CONE_ANGLE
    if(is_sleeping(parent))
    {   
        u_s = tan(KHRT_SCTR_ANGLE) * v_abs(parent->state.V);
        x_s = tan(KHRT_CONE_ANGLE) * v_rel_abs(parent->state.pos, parent->init_state.pos) + parent->state.diam/2.0;
    }
#endif	

    while(parent->state.mass > 0.0 && (parent->user[KH_ACCUMULATED_MASS] > m || force))
    {				
        /*  sample product diameter  */
        switch(KH_DIST_TYPE)
        {        
        case DIST_CHI_SQR:
            /*  sample from chi-square distribution  */
            d = chi_sqr_diam(smd);
            break;

        case DIST_ROSIN_R:
        default:
            /*  sample from rosin-rammler distribution  */
            d = rosin_r_diam(smd, KH_DIST_PARAMETER);
            break;
        }                       
        /*  create product parcel  */       
        Particle* child = create_product_parcel(parent, m, d, u_s, x_s);

        if(child != NULL)
        {
            reset_models_p(child);
        }        

        /*  parent state  */ 
		if(parent->user[KH_ACCUMULATED_MASS] > m)
		{
			parent->user[KH_ACCUMULATED_D2]   *= 1.0-m/parent->user[KH_ACCUMULATED_MASS];
			parent->user[KH_ACCUMULATED_D3]   *= 1.0-m/parent->user[KH_ACCUMULATED_MASS];        
			parent->user[KH_ACCUMULATED_MASS] *= 1.0-m/parent->user[KH_ACCUMULATED_MASS];
		}
		parent->user[KH_LIMIT] = 0.5+uniform_rand();

		/*  sample cut-off mass  */
		m = parent->init_state.mass*KH_CUTOFF*parent->user[KH_LIMIT];
		
        /*  update statistics  */
        parent->gvtp.n_shed++;
#ifdef VERBOSE_BREAKUP				
		Message("KH breakup (d0, d32, d):     %e   %e   %e\n", d0, smd, d);
#endif		
    }
    /*  no breakup  */
    return false;
}





real kh_collision_outcome(Tracked_Particle* collector, Particle* target, const real P_red, real* impact_parameter, collision_outcome_t* collision_outcome)
{
    return 0.0;
}



void kh_collision_handler(Tracked_Particle* collector, Particle* target, real* user_coalesced, const real P_red, const real impact_parameter, const collision_outcome_t collision_outcome)
{
    /*  multiple collisions  */
    real k1 = 1.0;
    real k2 = 1.0;
    
    if(P_N(collector) > P_N(target))
        k1 = P_red;
    else
        k2 = P_red;

    switch(collision_outcome)
    {
    case COLLISION_OUTCOME_COALESCENCE:
        if(user_coalesced == NULL)
            break;
        /*  handle cóalescence  */        
        user_coalesced[KH_ACCUMULATED_D2]   = k1*collector->user[KH_ACCUMULATED_D2]   + k2*target->user[KH_ACCUMULATED_D2];
        user_coalesced[KH_ACCUMULATED_D3]   = k1*collector->user[KH_ACCUMULATED_D3]   + k2*target->user[KH_ACCUMULATED_D3];
		user_coalesced[KH_ACCUMULATED_MASS] = k1*collector->user[KH_ACCUMULATED_MASS] + k2*target->user[KH_ACCUMULATED_MASS];
        break;
    default:
        break;
    }
    return;
}
