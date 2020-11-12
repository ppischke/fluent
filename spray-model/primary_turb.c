#include "udf.h"
#include "dpm.h"
#include "math.h"
#include "math_ext.h"
#include "model.h"
#include "primary_turb.h"
#include "secondary_tab.h"





void tb_init(real* user)
{
    /*  initialize user scalars  */
    user[TB_ACCUMULATED_MASS] = 0.0;
    user[TB_ACCUMULATED_D3]   = 0.0;
    user[TB_ACCUMULATED_D2]   = 0.0; 
    user[TB_LIMIT]            = 0.5+uniform_rand(); 	
    return;
}

void tb_initialize(Tracked_Particle* parent)
{
    tb_init(parent->user);
    return;
}

void tb_initialize_p(Particle* parent)
{
    tb_init(parent->user);
    return;
}





void tb_update(Tracked_Particle* parent)
{           
    /*  liquid material properties  */
    const real rho_liq
        = DPM_RHO
            (parent, P_MATERIAL(parent), P_T(parent));
    const real mu_liq
        = DPM_MU
            (parent);

    /*  nozzle statee  */
    const real U0 = v_abs(parent->init_state.V);
    const real D0 = parent->init_state.diam;

    /*  parent state  */
    const real m0 = parent->state0.mass - parent->user[TB_ACCUMULATED_MASS];    
    const real d0 = pow(6.0/M_PI * m0/rho_liq, 1.0/3.0);
    if(m0 <= 0.0)
    {   /*  no further breakup  */
        return;
    }

    /*  time  */
    const real t  = parent->state.time - parent->init_state.time;
    const real dt = parent->time_step;
    
    /*  dimensionless numbers  */
    const real Re = rho_liq*U0*D0 / mu_liq;

    /*  initial turbulent kinetic energy    */
    const real k0 = pow2(U0 * 0.16*pow(Re,-1.0/8.0));    
    /*  initial turbulent dissipation rate  */ 
    const real e0 = TB_CMU/TB_CL * pow(k0,3.0/2.0)/D0;

    /*  turbulent time-scale  */
    const real tau_t = k0/e0;    
    /*  dimensionless time    */
    const real theta = 1.0 + (TB_C2E-1.0)*t/tau_t;

    /*  turbulent kinetic energy  */
    const real k = k0/pow(theta,1.0/(TB_C2E-1.0));  
    /*  turbulent dissipation rate  */
    const real e = e0/pow(theta,1.0/(TB_C2E-1.0)+1.0);
    
    /*  breakup time and length scale  */
    const real tau_b = TB_CTB * TB_CMU * k/e;
    const real len_b = TB_CLB * TB_CMU * k/e * sqrt(k)/M_PI_2;

    /*  update parent state  */
    const real d = fmax(d0 - len_b/tau_b*dt, 0.0);
    const real m = 1.0/6.0*M_PI * rho_liq * pow3(d);
        
    /*  calculate sauter mean diameter  */    
    const real D2 = (pow3(d0) - pow3(d)) / len_b;
    const real D3 = (pow3(d0) - pow3(d));

    /*  accumulate user scalars  */
    parent->user[TB_ACCUMULATED_MASS] += (m0-m);
    parent->user[TB_ACCUMULATED_D2]   +=  D2;
    parent->user[TB_ACCUMULATED_D3]   +=  D3;
    	
    if(d < len_b)    
    {	
		tb_breakup(parent, true);
    }
    return;
}





bool tb_breakup(Tracked_Particle* parent, bool force)
{   
  	/*  sample cut-off mass   */	
	real m = parent->init_state.mass*TB_CUTOFF*parent->user[TB_LIMIT];		
	
	/*  sample mean diameter  */
	real d;
	real d0  = parent->state.diam;
    real smd = parent->user[TB_ACCUMULATED_D3]/parent->user[TB_ACCUMULATED_D2];   
	
    /*  scattering velocity   */
#ifdef TB_CONE_ANGLE	
    real u_s = tan(TB_SCTR_ANGLE) * v_abs(parent->state.V);
    real x_s = tan(TB_CONE_ANGLE) * v_rel_abs(parent->state.pos, parent->init_state.pos) + parent->state.diam/2.0;         
#else
	real u_s = 0.0;
	real x_s = 0.0;
#endif	
	
    while(parent->state.mass > 0.0 && (parent->user[TB_ACCUMULATED_MASS] > m || force))
    {
        /*  sample product diameter  */
        switch(TB_DIST_TYPE)
        {        
        case DIST_CHI_SQR:
            /*  sample from chi-square distribution  */
            d = chi_sqr_diam(smd);
            break;

        case DIST_ROSIN_R:
        default:
            /*  sample from rosin-rammler distribution  */
            d = rosin_r_diam(smd, TB_DIST_PARAMETER);
            break;
        }       	    
        /*  create product parcel  */       
        Particle* child = create_product_parcel(parent, m, d, u_s, x_s);
		
        if(child != NULL)
        {
            reset_models_p(child);
        }        
		
        /*  parent state  */                
		if(parent->user[TB_ACCUMULATED_MASS] > m)
		{
			parent->user[TB_ACCUMULATED_D2]   *= 1.0-m/parent->user[TB_ACCUMULATED_MASS];
			parent->user[TB_ACCUMULATED_D3]   *= 1.0-m/parent->user[TB_ACCUMULATED_MASS];
			parent->user[TB_ACCUMULATED_MASS] *= 1.0-m/parent->user[TB_ACCUMULATED_MASS];
		}
		parent->user[TB_LIMIT] = 0.5+uniform_rand();
		
		/*  sample cut-off mass  */
		m = parent->init_state.mass*TB_CUTOFF*parent->user[TB_LIMIT];
		
        /*  update statistics  */
        parent->gvtp.n_shed++;		
#ifdef VERBOSE_BREAKUP		
        Message("Turb. breakup (d0, d32, d):     %e   %e   %e\n", d0, smd, d);         	        
#endif
    }
    /*  no breakup  */
    return false;
}






/*  KHRT UDF  */
DEFINE_DPM_SCALAR_UPDATE(turb_breakup, cell, thread, initialize, parent)
{
    if(initialize)
    {
        if(parcel_initialize(parent))
        {
            /*  breakup models  */
            tb_initialize(parent);
            tab_initialize(parent);
        }
    }
    else
    {
        if(is_sleeping(parent))
        {   /*  KHRT breakup for primary droplets   */
            tb_update (parent);
            tb_breakup(parent, false);
        }
        else
        if(!is_suspended(parent))
        {   /*  TAB drop deformation   */
            tab_update(parent, NULL);           
        }
    }
    return;
}
