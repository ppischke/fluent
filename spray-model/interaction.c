#include "udf.h"
#include "dpm.h"
#include "sg.h"
#include "metric.h"
#include "math.h"
#include "math_ext.h"
#include "model.h"
#include "turb.h"
#include "vaporization.h"
#include "interaction.h"





bool source_update = false;





/*  PARTICLE SOURCE IN CELL (PSIC)  */
DEFINE_DPM_SOURCE(interaction, cell, thread, source, strength, particle)
{   	
#if !RPHOST
	source_update = true;


	
    /*  FLOW PROPERTIES  */		
	P_FLOW_VISCOSITY(particle)    = particle->cphase.mu;
	P_FLOW_DENSITY(particle)      = particle->cphase.rho;
	P_FLOW_PRESSURE(particle)     = particle->cphase.pressure;
		
	
	
    if(is_sleeping(particle))
    {   /*  clear source terms  */        
        source->mass = 0.0;
        source->energy = 0.0;
        memset(source->momentum_s, 0,sizeof(real)*ND_3);
        memset(source->momentum_ap,0,sizeof(real)*ND_3);
        return;
    }	
	
	
	
    /*  materials and species  */
    Material* gas_material = THREAD_MATERIAL(P_CELL_THREAD(particle));
    Material* liq_material = P_MATERIAL(particle);

    /*  component indices  */
    int i;          /*  liquid phase  */
    int j;          /*  gas phase     */

    /*  states  */
    particle_state_t* part_state = &(particle->state);
    particle_state_t* prev_state = &(particle->state0);

    /*  mass and temperature  */
    real m0 = prev_state->mass;
    real m  = part_state->mass;
    real m_vap;
    real m_comp0;   /*  multicomponent only  */
    real m_comp;    /*  multicomponent only  */

    real T0 = prev_state->temp;
    real T  = part_state->temp;    	
	real T_gas = particle->cphase.temp;
    real T_vap;     /*  vaporization temperature  */
    real T_ref;		/*  reference temperature     */
	real h_vap;
	
    real q_s;
    real q_l;
	
    

    switch(particle->type)
    {
    /*  SINGLE-COMPONENT  */
    case DPM_TYPE_DROPLET:
        /*  component index */
        j = particle->injection->evap_species_index;
        /*  component evaporated  */
        m_vap = m0 - m;
        
        /*  vaporization temperature  */
        T_vap = MATERIAL_PROP(liq_material, PROP_boil_temp);
		/*  reference temperature  */
		T_ref = MATERIAL_PROP(gas_material->component[j], PROP_reference_temp);
		
        /*  latent heat  */
        h_vap = MATERIAL_PROP_INTEGRAL(gas_material->component[j], PROP_Cp, T_ref, T_vap) - MATERIAL_PROP(liq_material, PROP_latent_heat) + MATERIAL_PROP_INTEGRAL(liq_material, PROP_Cp, T_vap, T0);

        /*  update explicit mass source  */
        source->mass = source->species[j] = strength * m_vap;

        /*  update explicit energy source  */
        q_s = strength * m     * MATERIAL_PROP_INTEGRAL(liq_material, PROP_Cp, T, T0);    
        q_l = strength * m_vap * h_vap;

        source->energy = q_l + q_s;
        break;



    /*  MULTI-COMPONENT  */
    case DPM_TYPE_MULTICOMPONENT:

        /*  clear source terms  */
        source->mass = 0.0;    	
        source->energy = 0.0; 
	
        for(i=0; i<liq_material->n_components; ++i)
        {
            /*  component index */
            j = particle->component.fluid_index[i];
            /*  component evaporated  */
            m_comp0 = m0 * particle->component.state0[i];
            m_comp  = m  * particle->component.state[i];
            m_vap   = m_comp0 - m_comp;
            
            /*  vaporization temperature  */
            T_vap = MATERIAL_PROP(liq_material->component[i], PROP_boil_temp);
			/*  reference temperature  */
			T_ref = MATERIAL_PROP(gas_material->component[j], PROP_reference_temp);
			
            /*  latent heat  */
            h_vap = MATERIAL_PROP_INTEGRAL(gas_material->component[j], PROP_Cp, T_ref, T_vap) - MATERIAL_PROP(liq_material->component[i], PROP_latent_heat)
                  + MATERIAL_PROP_INTEGRAL(liq_material->component[i], PROP_Cp, T_vap, T0);

            /*  update explicit mass source  */            
			source->mass += source->species[j] = strength * m_vap;
			
            /*  update explicit energy source  */
            q_s = strength * m_comp * MATERIAL_PROP_INTEGRAL(liq_material->component[i], PROP_Cp, T, T0);            
            q_l = strength * m_vap  * h_vap;					
            
            source->energy += q_l + q_s;
        }
		/*  update implicit energy source  */
		source->energy_ap = -fmax(source->energy/(T-T_gas),0.0);
        break;
						


    /*  UNSUPPORTED PARTICLE TYPE  */
    default:
        break;
    }	
	


    
    /*  TIME AND LENGTH SCALES */    
    
    /*  particle time scale    */	
	const real tau_p = particle->state.rho * pow2(particle->state.diam) / particle->cphase.mu / DragCoeff(particle);        
    /*  integration time step  */
    const real dt = particle->state.time - particle->state0.time;
    
    /*  particle length scale  */
    const real l_p = particle->state.diam;
    /*  turbulent length scale */
    const real l_t = particle->cphase.cell_eqv_length;
    
    
    
    /*  MOMENTUM SOURCE TERMS  */
	for(i=0; i<ND_3; ++i)
	{
		source->momentum_s[i] = strength * (m0*particle->state0.V[i] - m*particle->state.V[i]);
                       
        if(tau_p > 0.0)
        {
            source->momentum_ap[i] = -strength*dt * m/tau_p;
        }
	}		
	

			
    /*  TURBULENT SOURCE TERMS */	
#ifdef KSG_PRODUCTION
	real V[ND_3];	
	v_rel(particle->cphase.V, particle->V_prime, V);
    
    if(tau_p > 0.0)
    {
		C_DPMS_KSG_S(cell,thread) += strength*dt * m/tau_p * v_rel_sqr(particle->state.V, V) * pow(l_p/l_t, 4.0/3.0);
	}
#endif

#ifdef KSG_DAMPING
	if(tau_p > 0.0)
	{
		C_DPMS_TAU_S(cell,thread) += strength*dt * m/tau_p * 2.0 * pow(l_p/l_t, 4.0/3.0);
	}
#endif
    return;		
#endif
}





DEFINE_SOURCE(source_ksg,c,t,dS,eqn)
{   
	/*  explicit  */
	real Sk = -C_DPMS_TAU(c,t)/C_VOLUME(c,t)*C_K(c,t) + C_DPMS_KSG(c,t)/C_VOLUME(c,t);
    /*  implicit  */
	dS[eqn] = -C_DPMS_TAU(c,t)/C_VOLUME(c,t); 	
	/*  unresolved kinetic energy  */
    return Sk;
}





DEFINE_ADJUST(source_relax, d)
{
	if(source_update)
	{	
		Thread *t;      
		thread_loop_c(t,d)
		{   
			cell_t c;               
			begin_c_loop(c,t)
			{   
				/*	store previous time-step volume fraction  */				
				C_DPMS_KSG(c,t)      = C_DPMS_KSG_S(c,t);
				C_DPMS_KSG_S(c,t)    = 0.0;								
                
				C_DPMS_TAU(c,t)      = C_DPMS_TAU_S(c,t);
				C_DPMS_TAU_S(c,t)    = 0.0;								
			}
			end_c_loop(c,t)
		}
	}
	source_update = false;		
}





DEFINE_EXECUTE_AT_END(source_reset)
{   
    Domain *d = Get_Domain(1);
    Thread *t;      
    thread_loop_c(t,d)
    {   
        cell_t c;               
        begin_c_loop(c,t)
        {   			
			/*	reset previous time-step sources  */			
			C_DPMS_KSG_S(c,t) = 0.0;
            C_DPMS_TAU_S(c,t) = 0.0;
            
            /*  store effective pressure gradient */
            int i;
            for(i=0; i<ND_3; ++i)
            {
                C_P_G_EFF(c,t,i) = C_P_G(c,t)[i] - (C_DPMS_MOM_S(c,t)[i] + (C_VEL(c,t,i)-C_VEL_M1(c,t,i))*C_DPMS_MOM_AP(c,t)[i])/C_VOLUME(c,t);
            }
        }
        end_c_loop(c,t)
    }
	source_update = false;
}
