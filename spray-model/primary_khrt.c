#include "math.h"
#include "math_ext.h"
#include "udf.h"
#include "dpm.h"
#include "model.h"

#include "primary_khrt.h"

#include "secondary_tab.h"
#include "secondary_rdiwakar.h"
#include "secondary_kh.h"





void rt_init(real* user)
{
    /*  initialize user scalars  */
    user[RT_PROGRESS]   = 0.0;
    user[RT_WAVELENGTH] = 0.0;
    return;
}

void rt_initialize(Tracked_Particle* parent)
{
    rt_init(parent->user);
    return;
}

void rt_initialize_p(Particle* parent)
{
    rt_init(parent->user);
    return;
}





void rt_update(Tracked_Particle* parent)
{
    /*  states  */
    cphase_state_t*   gas_state  = &(parent->cphase);
    particle_state_t* part_state = &(parent->state);
    particle_state_t* init_state = &(parent->init_state);       
    
    /*  gas material properties  */
    real rho_gas;
    if(is_sleeping(parent))
        rho_gas = C_R(P_CELL(parent),P_CELL_THREAD(parent))
                    * C_T(P_CELL(parent),P_CELL_THREAD(parent)) / parent->state.temp;
    else
        rho_gas = gas_state->rho;

    /*  liq material properties  */
    real rho_liq
        = DPM_RHO
            (parent, P_MATERIAL(parent), P_T(parent));
	real st_liq
        = DPM_SURFTEN
            (parent);

    /*  relative velocity and penetration  */   
    real u_rel;
    real x_rel = v_rel_abs(part_state->pos, init_state->pos);

    if(x_rel < fmax(RT_L0, RT_CL*sqrt(rho_liq/rho_gas)*part_state->diam) && parent->user[RT_PROGRESS] <= 0.0)
    {
        rt_initialize(parent);
        return;
    }

#ifdef KHRT_ABSOLUTE_VELOCITY
    /*  primary breakup  */     
    if(is_sleeping(parent))
    {   /*  absolute velocity for primary breakup   */
        u_rel = v_abs(init_state->V);       
    }
    else
#endif
    {   /*  relative velocity for secondary breakup */      
        u_rel = v_rel_abs(part_state->V, gas_state->V);
    }
    real dt = parent->time_step;
    
#ifdef KHRT_DRAG_COEFF
    real a = 3.0/4.0 * pow2(u_rel)/init_state->diam * rho_gas/rho_liq * KHRT_DRAG_COEFF;
#else
    real a = fabs(u_rel * gas_state->mu/rho_liq / pow2(part_state->diam) * DragCoeff(parent));
#endif

    /*  unstable frequency  */
    real omega = sqrt(2.0/3.0 / sqrt(3.0*st_liq) * pow3(sqrt(a * (rho_liq - rho_gas))) / (rho_liq + rho_gas));
    /*  unstable wave length  */
    real lambda = 2.0*M_PI*RT_CD * sqrt(3.0*st_liq / a / (rho_liq - rho_gas));  
    
    /*  breakup time  */
    real tau = RT_CT/omega;
    
    /*  breakup progress  */
    if(lambda < part_state->diam)
    {
        parent->user[RT_PROGRESS]  += dt/tau;
        parent->user[RT_WAVELENGTH] = lambda;
    }
    else
    {
        rt_initialize(parent);
        return;
    }
    return;
}





bool rt_breakup(Tracked_Particle* parent)
{
   /*  breakup  */
    if(parent->user[RT_PROGRESS] >= 1.0)
    {                                   
        /*  sample mean diameter  */        
        real d;
        real d0  = parent->state.diam;
        real smd = parent->state.diam * pow(parent->user[RT_WAVELENGTH]/parent->state.diam, RT_CE);

        /*  scattering velocity   */
        real u_s = 0.0;
        real x_s = 0.0;
#ifdef KHRT_CONE_ANGLE
        if(is_sleeping(parent))
        {   
            u_s = tan(KHRT_SCTR_ANGLE) * v_abs(parent->state.V);            
            x_s = tan(KHRT_CONE_ANGLE) * v_rel_abs(parent->state.pos, parent->init_state.pos) + parent->state.diam/2.0;
        }
#endif   
           
        while(parent->state.mass > 0.0)
        {
			/*  sample cut-off mass   */
			real m = parent->init_state.mass*RT_CUTOFF*(0.5+uniform_rand()); 
			
            /*  sample product diameter  */
            switch(RT_DIST_TYPE)
            {        
            case DIST_CHI_SQR:
                /*  sample from chi-square distribution  */
                d = chi_sqr_diam(smd);
                break;

            case DIST_ROSIN_R:
            default:
                /*  sample from rosin-rammler distribution  */
                d = rosin_r_diam(smd, RT_DIST_PARAMETER);
                break;
            }       
            /*  create product parcel  */
            Particle* child = create_product_parcel(parent, m, d, u_s, x_s);

            if(child != NULL)
            {
                reset_models_p(child);    
            }
			
            /*  update statistics  */
			parent->gvtp.n_shed++;
#ifdef VERBOSE_BREAKUP				
            Message("RT breakup (d0, d32, d):     %e   %e   %e\n", d0, smd, d);
#endif			
        }
        return true;
    }
    /*  no breakup  */
    return false;
}





/*  KHRT-TAB UDF  */
DEFINE_DPM_SCALAR_UPDATE(khrt_tab, cell, thread, initialize, parent)
{
    if(initialize)
    {
        if(parcel_initialize(parent))
        {
            /*  breakup models  */
            kh_initialize(parent);
            rt_initialize(parent);
            tab_initialize(parent);
        }
    }
    else
    {
        if(is_sleeping(parent))
        {   /*  KHRT breakup for primary droplets   */
            kh_update (parent, NULL);
            kh_breakup(parent, NULL, false);
            rt_update (parent);
            rt_breakup(parent);
        }
        else
        if(!is_suspended(parent))
        {   /*  TAB drop deformation   */
            tab_update(parent, NULL);

            /*  TAB secondary breakup  */
#ifndef KHRT_TERTIARY
            if(is_primary(parent))
#endif      /*  no tertiary breakup  */
            {   
                tab_breakup(parent, NULL, NULL);
            }
        }
    }
    return;
}





/*  KHRT-TAB UDF  */
DEFINE_DPM_SCALAR_UPDATE(khrt_khtab, cell, thread, initialize, parent)
{
    if(initialize)
    {
        if(parcel_initialize(parent))
        {
            /*  breakup models  */
            kh_initialize(parent);
            rt_initialize(parent);
            tab_initialize(parent);
        }
    }
    else
    {
        if(is_sleeping(parent))
        {   /*  KHRT breakup for primary droplets   */
            kh_update (parent, NULL);
            kh_breakup(parent, NULL, false);
            rt_update (parent);
            rt_breakup(parent);
        }
        else
        if(!is_suspended(parent))
        {   /*  TAB drop deformation   */
            tab_update(parent, NULL);
            
            /*  KHTAB secondary breakup  */
#ifndef KHRT_TERTIARY
            if(is_primary(parent))
#endif      /*  no tertiary breakup  */
            {   
                kh_update  (parent, NULL);
                kh_breakup (parent, NULL, false);                
                tab_breakup(parent, NULL, NULL);
            }
        }
    }
    return;
}





/*  KHRT UDF  */
DEFINE_DPM_SCALAR_UPDATE(khrt, cell, thread, initialize, parent)
{
    if(initialize)
    {
        if(parcel_initialize(parent))
        {
            /*  breakup models  */
            kh_initialize(parent);
            rt_initialize(parent);
            tab_initialize(parent);
        }
    }
    else
    {
        if(is_sleeping(parent))
        {   /*  KHRT breakup for primary droplets   */
            kh_update (parent, NULL);
            kh_breakup(parent, NULL, false);
            rt_update (parent);
            rt_breakup(parent);
        }
        else
        if(!is_suspended(parent))
        {   /*  TAB drop deformation   */
            tab_update(parent, NULL);           
        }
    }
    return;
}
