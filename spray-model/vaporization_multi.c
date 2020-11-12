#include "udf.h"
#include "dpm.h"
#include "math.h"
#include "math_ext.h"
#include "model.h"
#include "vaporization.h"





/*  MULTICOMPONENT VAPORIZATION LAW
//      inert heating
//      inert cooling
//      vaporization
//      boiling
*/
void vaporization_multi(Tracked_Particle* particle)
{
#if !RPHOST 
    /*  time step  */
    real dt = particle->time_step;    

    /*  properties of states  */
    cphase_state_t*   gas_state  = &(particle->cphase);
    particle_state_t* part_state = &(particle->state);    

	
	
    /*  materials  */
    Material* gas_material = THREAD_MATERIAL(P_CELL_THREAD(particle));
    Material* liq_material = P_MATERIAL(particle);

    /*  species indices  */
    int n_gas = gas_material->n_components;
    int n_liq = liq_material->n_components;

    int i;  /*  species index liquid phase  */
    int j;  /*  species index gas phase     */    
    int i_of_j[VAPOR_MAX_SPECIES];
	
	

	/*  mark non-volatile components with -1  */
	memset(i_of_j, -1, sizeof(int)*VAPOR_MAX_SPECIES);
	/*  mark volatile components with liquid phase index  */
    for(i=0; i<n_liq; ++i)
    {	
        j = particle->component.fluid_index[i];        
        i_of_j[j] = i;
	}
            

            
    /*  properties of state  */        
    real p_gas = gas_state->pressure;    
    real T_gas = gas_state->temp;

    /*  droplet temperatures  */
	real T0 = part_state->temp;
    real T  = T0; 	
       
    /*  droplet mass  */    
	real m0 = part_state->mass;
    real m  = m0;
    
    /*  droplet diameter */
	real d0 = mass_to_diam(part_state);	
    real d  = d0;

    /*  molar fractions  */        
	real x_liq[VAPOR_MAX_SPECIES];
    real x_gas[VAPOR_MAX_SPECIES];
    real x_vap[VAPOR_MAX_SPECIES];
            
    /*  mass fractions  */		
    real*y_liq = particle->component.state;	
    real*y_gas = particle->cphase.yi;
    real y_vap[VAPOR_MAX_SPECIES];
	
	/*	partial masses  */
	real m_liq0[VAPOR_MAX_SPECIES];
	for(i=0; i<n_liq; ++i)
	{
		m_liq0[i] = m0*y_liq[i];
	}

                    
    
	/*  molar fractions of liquid phase  */
	real sum_y_over_mwi_liq = 0.0;
	/*  convert from mass fractions   */
	for(i=0; i<n_liq; ++i)
	{	
		j = particle->component.fluid_index[i];
		
		sum_y_over_mwi_liq += x_liq[i] = y_liq[i]/MATERIAL_PROP(gas_material->component[j], PROP_mwi);
	}
	for(i=0; i<n_liq; ++i)
	{
		x_liq[i] /= sum_y_over_mwi_liq;        
	}      

	
	
	/*  molar fractions of gas phase  */
	real sum_y_over_mwi_gas = 0.0;
	/*  convert from mass fractions   */
	for(j=0; j<n_gas; ++j)
	{        
		sum_y_over_mwi_gas += x_gas[j] = y_gas[j]/MATERIAL_PROP(gas_material->component[j], PROP_mwi);		
	}
	for(j=0; j<n_gas; ++j)
	{
		x_gas[j] /= sum_y_over_mwi_gas;
	}      
		
        
        
    /*  vapor pressure  */
    real p_crit = 0.0;
	/*  calculate from raoult's law   */
	for(i=0; i<n_liq; ++i)
	{
        p_crit += x_liq[i]*DPM_VAPOR_PRESSURE(particle, liq_material->component[i], INFINITY);
    }
    
	/*  molar fractions of volatiles  */
	real volatile_x_vap = 0.0;  /*  volatiles in vapor phase  */
	real volatile_x_gas = 0.0;  /*  volatiles in gas phase    */
    /*  calculate from raoult's law   */
	for(i=0; i<n_liq; ++i)
	{
		j = particle->component.fluid_index[i];
		
		volatile_x_vap += x_vap[j] = x_liq[i]*DPM_VAPOR_PRESSURE(particle, liq_material->component[i], T0) / fmin(p_gas, p_crit);                
		volatile_x_gas += x_gas[j];
	}     		
#ifdef VAPOR_MAX_FRACTION
	/*  limit maximum vapor fraction  */
	if(volatile_x_vap >= VAPOR_MAX_FRACTION)
	{
		for(i=0; i<n_liq; ++i)
		{	
			j = particle->component.fluid_index[i];
			
			x_vap[j]  *= VAPOR_MAX_FRACTION/volatile_x_vap;
		}
		volatile_x_vap = VAPOR_MAX_FRACTION;
	}
#endif
	

	
	/*  molar fractions of non-volatiles  */
	for(j=0; j<n_gas; ++j)
	{           
		i = i_of_j[j];
		/*  non-volatile?  */
		if(i==-1)                
		{
			x_vap[j] = x_gas[j] * (1.0-volatile_x_vap) / (1.0-volatile_x_gas);
		}                      
	}
	
	

	/*  mass fractions of vapor  */
	real sum_x_times_mwi_vap = 0.0;
	/*  convert from molar fractions  */
	for(j=0; j<n_gas; ++j)
	{		
		sum_x_times_mwi_vap += y_vap[j] = x_vap[j]*MATERIAL_PROP(gas_material->component[j], PROP_mwi);
	}
	for(j=0; j<n_gas; ++j)
	{
		y_vap[j] /= sum_x_times_mwi_vap;      
	}
	
	
	
	/*  mass fractions of volatiles  */
	real volatile_y_vap = 0.0;  /*  volatiles in vapor phase  */
	real volatile_y_gas = 0.0;  /*  volatiles in gas phase  */    
	/*  sum over volatile fractions  */
	for(i=0; i<n_liq; ++i)
	{    
		j = particle->component.fluid_index[i];
		
		volatile_y_vap += y_vap[j];
		volatile_y_gas += y_gas[j];		
	}      
			
	

	/*  gas (cell) material properties  */    
	real k_gas
		= gas_state->tCond;
	real cp_gas
		= gas_state->sHeat;
	real mu_gas
		= gas_state->mu;
	real rho_gas
		= gas_state->rho;    

	/*  diffusivity  */ 
	real D_gas[VAPOR_MAX_SPECIES];

	for(i=0; i<n_liq; ++i)
	{   
		D_gas[i] = MATERIAL_PROP_POLYNOMIAL(liq_material->component[i], PROP_binary_diffusivity, D_T_REF) * pow(T_gas/D_T_REF,1.5) / (p_gas/D_P_REF);
	}
	


	/*  dimensionless numbers  */
	/*  apply 1/3 rule here    */
	real Re = rho_gas * v_rel_abs(part_state->V, gas_state->V) * d0 / mu_gas;
		
	/*  mass transfer number   */
	real Bm = (volatile_y_vap - volatile_y_gas)/(1.0 - volatile_y_vap);    
			
	/*  mass transfer  */
	real Sc[VAPOR_MAX_SPECIES];
	real Sh[VAPOR_MAX_SPECIES];
	for(i=0; i<n_liq; ++i)
	{
		Sc[i] = mu_gas / rho_gas / D_gas[i];
		Sh[i] = correlation_abramzon_sirignano(Re, Sc[i], Bm);        
	}

	

	/*  component diffusive fluxes  */
	real j_flux[VAPOR_MAX_SPECIES];
	real sum_j_flux = 0.0;

	for(i=0; i<n_liq; ++i)
	{   
		j = particle->component.fluid_index[i];

		/*  mass-based diffusive flux  */
		j_flux[i] = Sh[i]/d0 * rho_gas*D_gas[i] * (y_vap[j]-y_gas[j]);

		/*  diffusive flux  */
		sum_j_flux += j_flux[i];
	}


	
	/*  component mass fluxes  */
	real m_flux[VAPOR_MAX_SPECIES];        
	real sum_m_flux = 0.0;

	for(i=0; i<n_liq; ++i)
	{   
		j = particle->component.fluid_index[i];        

		/*  mass-based multi-component stefan correction  */
		m_flux[i] = (j_flux[i] + sum_j_flux*y_vap[j]/(1.0-volatile_y_vap));
		
#ifdef VAPOR_ALLOW_CONDENSATION
		/*  no condensation on components evaporated  */
		if (m_flux[i]<0.0 && m_liq0[i]<=0.0)            
			m_flux[i]=0.0;
#else
		/*  no condensation  */
		if (m_flux[i]<0.0)            
			m_flux[i]=0.0;
#endif
		/*  mass flux  */
		sum_m_flux += m_flux[i];			
	}    
	


	/*  mass balances  */		
	real m_liq[VAPOR_MAX_SPECIES];
	real m_vap[VAPOR_MAX_SPECIES];    
	
	m = 0.0;
	for(i=0; i<n_liq; ++i)
	{			
		m_vap[i] = m_flux[i]*M_PI*POW2(d0)*dt;
		/*  component masses  */			
		if (m_vap[i] > m_liq0[i])
			m_vap[i] = m_liq0[i];
									
		/*  vaporized masses  */
		m += m_liq[i] = m_liq0[i] - m_vap[i];
	}

	
		
	/*  heat balances  */
	if(m > 0.0)
	{               
		/*  update mass fractions */
		for(i=0; i<n_liq; ++i)
		{               
			y_liq[i] = m_liq[i]/m;
		}
	}	
	
	

	/*  heat transfer  */        			
	real q_vap = 0.0;   
	real c_liq = 0.0;
		
	for(i=0; i<n_liq; ++i)
	{   
		j = particle->component.fluid_index[i];

		/*  reference temperature  */
		real T_vap = MATERIAL_PROP(liq_material->component[i], PROP_boil_temp);     

		/*  latent heat  */                
		real h_vap = MATERIAL_PROP_INTEGRAL(liq_material->component[i], PROP_Cp, T0, T_vap) 
				   + MATERIAL_PROP         (liq_material->component[i], PROP_latent_heat) 
				   + MATERIAL_PROP_INTEGRAL(gas_material->component[j], PROP_Cp, T_vap, T0);
				
		/*  latent heat  */        
		q_vap += m_vap[i] * h_vap;       
		
		/*  droplet sensible heat  */
		c_liq += m_liq[i] * MATERIAL_PROP_POLYNOMIAL(liq_material->component[i], PROP_Cp, T0);								
	}
	if((q_vap > 0.0) ^ (m0 > m))
		q_vap = 0.0;


	/*  mass transfer coefficient  */
	real g = sum_j_flux/(volatile_y_vap-volatile_y_gas);

	/*  heat transfer coefficient  */
	real h = g*cp_gas;			
	   
	int k;
	for(k=0; k<VAPOR_HTC_ITER; ++k)
	{
		/*  heat transfer number  */
		real Bt = Bm * fabs(g*cp_gas/h);
		/*  heat transfer  */
		real Pr = mu_gas * cp_gas / k_gas;    
		real Nu = correlation_abramzon_sirignano(Re, Pr, Bt);
		
		/*  heat transfer coefficient  */
		h = Nu * k_gas/d0;        			
	}
	h = h * M_PI*POW2(d0)*dt;	
				
	
	
	if(c_liq + h > 0.0)
	{
		/*  particle temperature  */    
		T = (c_liq*T0 + h*T_gas - q_vap) / (c_liq + h);
	}
	
	
	
	/*	droplet state  */    
	part_state->mass = m;        
	part_state->temp = T;    
	part_state->rho  = DPM_RHO(particle, P_MATERIAL(particle), T);  
	d = mass_to_diam(part_state);
	

	
    /*  droplet supercritical  */
#ifndef VAPOR_KEEP_SUPERCRITICAL    
    if(DPM_SURFTEN(particle) <= 0.0)
    {		
        part_state->mass = 0.0;
        part_state->diam = 0.0;  		
    }
#endif    
#endif
    return;
}
