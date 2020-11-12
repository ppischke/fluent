#define __SECONDARY_TAB_C

#include "udf.h"
#include "dpm.h"
#include "math.h"
#include "math_ext.h"
#include "model.h"

#include "secondary_tab.h"
#include "secondary_kh.h"





/*  default tab constants  */
tab_constants_t tab_default_constants =
{
    TAB_WE_MIN,
    TAB_WE_MAX,
    TAB_N_PARCELS,
    TAB_DIST_TYPE,
    TAB_DIST_PARAMETER
};





void tab_init(real* user)
{
    /*  initialize user scalars  */
    user[TAB_DISPLACEMENT] = 0.0;
    user[TAB_VELOCITY]     = 0.0;
    return;
}

void tab_initialize(Tracked_Particle* parent)
{
    tab_init(parent->user);
    return;
}

void tab_initialize_p(Particle* parent)
{
    tab_init(parent->user);
    return;
}





void tab_update(Tracked_Particle* parent, tab_constants_t* TAB)
{   
    /*  states  */
    cphase_state_t*   gas_state  = &(parent->cphase);
    particle_state_t* part_state = &(parent->state);

    /*  liquid material properties  */
    const real rho_liq
        = DPM_RHO
            (parent, P_MATERIAL(parent), P_T(parent));
    const real st_liq
        = DPM_SURFTEN
            (parent);
    const real mu_liq
        = DPM_MU
            (parent);

    /*  gas (cell) material properties  */
    const real rho_gas
        = gas_state->rho;

    /*  parent state  */
    const real r
        = part_state->diam / 2.0;
    const real u_rel_sqr
        = v_rel_sqr(part_state->V, gas_state->V);
    const real dt
        = parent->time_step;

        
    /*  displacement and velocity  */
    real y0_eq = TAB_CF/TAB_Cb/TAB_Ck * rho_gas*u_rel_sqr*r/st_liq;
                
    real y0 = parent->user[TAB_DISPLACEMENT] - y0_eq;
    real y1 = parent->user[TAB_VELOCITY];   
    
    
    /*  frequencies and time constants  */
    const real k0
        = TAB_Cd * mu_liq/rho_liq / pow2(r) / 2.0;
    const real omega0_sq
        = TAB_Ck * st_liq/rho_liq / pow3(r);

    const real omega_sq = pow2(k0) - omega0_sq;

    if(omega_sq < 0.0)
    {   /*  oscillating  */        
        const real omega = sqrt(-omega_sq);

        const real y_cos =  y0;
        const real y_sin = (y0*k0 + y1)/omega;

        /*  update TAB displacement  */
        y0 = exp(-k0*dt) * (y_cos * cos(omega*dt) + y_sin * sin(omega*dt));
        /*  update TAB velocity  */
        y1 = exp(-k0*dt) * (y_cos *-sin(omega*dt) + y_sin * cos(omega*dt))*omega - y0*k0;
    }
    else
    {   /*  non-oscillating  */
        const real kA = k0 + sqrt(omega_sq);
        const real kB = k0 - sqrt(omega_sq);                       

        const real yA = (-y1 - y0*kB)/(kA-kB);
        const real yB = ( y1 + y0*kA)/(kA-kB);

        /*  update TAB displacement and velocity  */
        y0 =  (   yA * exp(-kA*dt) +    yB * exp(-kB*dt));
        y1 = -(kA*yA * exp(-kA*dt) + kB*yB * exp(-kB*dt));
    }
    /*  update displacement and velocity  */
	parent->user[TAB_EQUILIBRIUM]  = y0_eq;
    parent->user[TAB_DISPLACEMENT] = y0_eq + y0;
    parent->user[TAB_VELOCITY]     = y1;
	
    return;
}





bool tab_breakup(Tracked_Particle* parent, tab_constants_t* TAB, kh_constants_t* KH)
{
    /*  TAB displacement  */
    const real y0    = parent->user[TAB_DISPLACEMENT];    
    const real y0_eq = parent->user[TAB_EQUILIBRIUM];    
		
		
    /*  sample cut-off mass   */
    real m = parent->state.mass/TAB_N_PARCELS;

    /*  sample mean diameter  */
	real d;
	real d0  = parent->state.diam;
    real smd = parent->state.diam / (1.0 + TAB_Ck*TAB_K/20.0*POW2(y0_eq));
	
    /*  scattering velocity   */
    real x_s = 0.0;		
#ifdef TAB_Cv
    real u_s = TAB_Cb*TAB_Cv*d0*y1;
#else
    real u_s = 0.0;
#endif
		

    /*  TAB regime  */
    const real st_liq
        = DPM_SURFTEN(parent);    
    const real rho_gas
        = parent->cphase.rho;
		
    const real We
        = rho_gas * v_rel_sqr(parent->state.V, parent->cphase.V) * d0 / st_liq;
	
    if(y0 >= 1.0 && We > TAB_WE_MIN && We < TAB_WE_MAX)
    {		        
		/*  TAB breakup  */
        while(parent->state.mass > 0.0)
        {
            /*  sample child diameter  */
            switch(TAB_DIST_TYPE)
            {            
            case DIST_CHI_SQR:
                /*  sample from chi-square distribution  */
                d = chi_sqr_diam(smd);
                break;

            case DIST_ROSIN_R:
            default:
                /*  sample from rosin-rammler distribution  */
                d = rosin_r_diam(smd, TAB_DIST_PARAMETER);
                break;
            }       
            /*  create child parcel  */
            Particle* child = create_product_parcel(parent, m, d, u_s, x_s);

            if(child != NULL)
            {
                reset_models_p(child);
            }
			
			/*	update statistics  */
			parent->gvtp.n_shed++;
#ifdef VERBOSE_BREAKUP					
			Message("TAB breakup (d0, d32, d):      %e   %e   %e\n", d0, smd, d);
#endif
        }
        return true;
    }
    /*  no breakup  */
    return false;
}





real tab_collision_outcome(Tracked_Particle* collector, Particle* target, const real P_red, real* impact_parameter, collision_outcome_t* collision_outcome)
{
    return 0.0;
}





void tab_collision_handler(Tracked_Particle* collector, Particle* target, real* user_coalesced, const real P_red, const real impact_parameter, const collision_outcome_t collision_outcome)
{
    switch(collision_outcome)
    {
    case COLLISION_OUTCOME_COALESCENCE:
        if(user_coalesced == NULL)
            break;
#ifdef TAB_COALESCENCE
        /*  reset displacement  */
        user_coalesced[TAB_DISPLACEMENT] = TAB_COALESCENCE;
        user_coalesced[TAB_VELOCITY]     = 0.0;
#endif
        break;


    default:
#ifdef TAB_COLLIDE  
        if(P_N(collector) > P_N(target))
        {
            /*  average displacement in collector parcel  */
            P_USER(collector)[TAB_DISPLACEMENT]
                *= (1.0 - P_N(target)*P_red/P_N(collector));
            P_USER(collector)[TAB_VELOCITY]
                *= (1.0 - P_N(target)*P_red/P_N(collector));
            /*  reset target  */
            P_USER(target)[TAB_DISPLACEMENT]
                 =  0.0;
            P_USER(target)[TAB_VELOCITY]
                 =  0.0;
        }
        else
        {
            /*  average displacement in target parcel  */
            P_USER(target)[TAB_DISPLACEMENT]
                *= (1.0 - P_N(collector)*P_red/P_N(target));
            P_USER(target)[TAB_VELOCITY]
                *= (1.0 - P_N(collector)*P_red/P_N(target));
            /*  reset collector  */
            P_USER(collector)[TAB_DISPLACEMENT]
                 =  0.0;
            P_USER(collector)[TAB_VELOCITY]
                 =  0.0;
        }
#endif
        break;
    }
    return;
}





/*  dynamic drag  */
DEFINE_DPM_DRAG(dynamic_tab, Re, particle)
{
    if(is_sleeping(particle))
        return 0.0;

    real y0 = particle->user[TAB_DISPLACEMENT];    
    /*  displacement  */
    if (y0 > 1.0)
        y0 = 1.0;
    if (y0 < 0.0)
        y0 = 0.0;
    	
	/*	diameter ratio  */
	real dr = 1.0 + TAB_Cb*y0;
		
    /*  scale drag  */
    real cd = SphereDragCoeff(Re)*(1.0+y0*TAB_DRAG_SCALAR) * POW2(dr);

    /*  limit drag  */
    if(cd > 0.0)
        return cd;
    else
        return 18.0;
}
