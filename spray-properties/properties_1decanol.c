#include "udf.h"
#include "dpm.h"
#include "math.h"

#include "properties_mixture.h"



#define M_1_DECANOL  158.284;



real property_1_decanol_liq_viscosity(real T)
{   
	const real T_min = 285.00;
	const real T_max = 504.07;

	real viscosity;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return viscosity = exp(-80.656 + 6325.5/T + 9.646*log(T));	    
}



real property_1_decanol_liq_surface_tension(real T)
{	
/*	TO DO:
//	limit surface tension
*/

	const real T_min = 280.05;
	const real T_max = 687.30;
	
	real surface_tension;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return surface_tension = 0.051263*pow((1.0 - T/T_max), 1.0395);	
}



real property_1_decanol_liq_specific_heat(real T)
{
	const real T_min = 280.05;
	const real T_max = 504.07;

	real specific_heat;

	if (T < T_min)	
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

   	return specific_heat = (4.9885e6 - 5.2898e4*T + 216.35*pow2(T) - 0.37538*pow3(T) + 2.3674e-4*pow4(T)) / M_1_DECANOL;    
}



real property_1_decanol_liq_density(real T)
{
	const real T_min = 280.05;
	const real T_max = 687.30;

	real density;

	if (T < T_min)	
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return density = 0.37384 / (pow(0.224241, 1.0+pow(1.0-(T/T_max), 0.26646))) * M_1_DECANOL;
}



real property_1_decanol_heat_vaporization(real T)
{
	const real T_min = 280.05;
	const real T_max = 687.30;

	real heat_vaporization;

	if (T < T_min)
		T = T_min;
  	else
	if (T > T_max)
		T = T_max;

	return heat_vaporization = 1.1750e08 * pow(1.0-T/T_max, 0.65112) / M_1_DECANOL;
}



real property_1_decanol_vapor_pressure(real T)
{	
	const real T_crit = 687.30;
	const real T_min  = 280.05;

	real vapor_pressure;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_crit)
		T = T_crit;    

	return vapor_pressure = exp(250.59 - 19169/T - 32.903*log(T) + 1.4627e-5*pow2(T));
}



real property_1_decanol_vapor_diffusivity(real T, real p)
{


	real vapor_diffusivity;
	real temp = (158.284+29.0) / (29.0*158.284);

	return vapor_diffusivity = 1.013e-2 * pow(T , 1.75) * pow(temp , 0.5) / pow2(pow(19.7,1.0/3.0)+pow(215.93,1.0/3.0)) / p;
}



real property_1_decanol_vapor_thermal_cond(real T)
{
	const real T_min =  504.07;
	const real T_max =  1000.00;

	real vapor_thermal_cond;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return vapor_thermal_cond = -3.0720e-1 * pow(T, 0.489) / (1.0 - 6.75e4/T - 2.94e7/pow2(T));
}



real property_1_decanol_vapor_viscosity(real T)
{
	const real T_min =  280.05;
	const real T_max = 1000.00;

	real vapor_viscosity;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return vapor_viscosity = 5.5065e-8 * pow(T, 0.8341) / (1.0 + 79.56/T);
}



DEFINE_DPM_PROPERTY(one_decanol_liq_viscosity, c, t, p)
{
	return property_1_decanol_liq_viscosity(P_T(p));
}



DEFINE_DPM_PROPERTY(one_decanol_liq_surface_tension, c, t, p)
{	
	return property_1_decanol_liq_surface_tension(P_T(p));
}



DEFINE_DPM_PROPERTY(one_decanol_liq_specific_heat, c, t, p)
{
	return property_1_decanol_liq_specific_heat(P_T(p));
}



DEFINE_DPM_PROPERTY(one_decanol_liq_density, c, t, p)
{
	return property_1_decanol_liq_density(P_T(p));
}



DEFINE_DPM_PROPERTY(one_decanol_heat_vaporization, c, t, p)
{
	return property_1_decanol_heat_vaporization(P_T(p));
}



DEFINE_DPM_PROPERTY(one_decanol_vapor_pressure, c, t, p)
{	
	return property_1_decanol_vapor_pressure(P_T(p));
}




DEFINE_PROPERTY(one_decanol_vapor_diffusivity, c, t)
{
	return property_1_decanol_vapor_diffusivity(C_T(c,t), C_P(c,t));
}




DEFINE_PROPERTY(one_decanol_vapor_viscosity, c, t)
{
	return property_1_decanol_vapor_viscosity(C_T(c,t));
}



DEFINE_PROPERTY(one_decanol_vapor_thermal_cond, c, t)
{
	return property_1_decanol_vapor_thermal_cond(C_T(c,t));
}

