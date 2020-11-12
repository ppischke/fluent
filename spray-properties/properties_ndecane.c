#include "udf.h"
#include "dpm.h"
#include "math.h"
#include "misc.h"

#include "properties_mixture.h"



#define M_N_DECANE  142.285



real property_n_decane_liq_viscosity(real T)
{
	const real T_min = 243.51;
	const real T_max = 448.15;

	real viscosity;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return viscosity = exp(-16.468 + 1533.5/T + 0.7511*log(T));	
}



real property_n_decane_liq_surface_tension(real T)
{	
/*	TO DO:
//	limit surface tension
*/

	const real T_min = 243.51;
	const real T_max = 617.70;
	
	real surface_tension;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return surface_tension = 0.055435*pow((1.0-T/T_max), 1.3095);	
}



real property_n_decane_liq_specific_heat(real T)
{
	const real T_min = 243.51;
	const real T_max = 460.0;

	real specific_heat;

	if (T < T_min)	
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

   	return specific_heat = (278626.0 - 197.91*T + 1.0737*pow2(T)) / M_N_DECANE;
}



real property_n_decane_liq_density(real T)
{
	const real T_min = 243.51;
	const real T_max = 617.7;

	real density;

	if (T < T_min)	
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return density = 0.42831 / (pow(0.25745, 1.0+pow(1.0-(T/T_max), 0.28912))) * M_N_DECANE;
}



real property_n_decane_heat_vaporization(real T)
{
	const real T_min = 243.51;
	const real T_max = 617.7;

	real heat_vaporization;

	if (T < T_min)
		T = T_min;
  	else
	if (T > T_max)
		T = T_max;

	return heat_vaporization = 6.6126e7 * pow(1.0-T/T_max, 0.39797) / M_N_DECANE;
}



real property_n_decane_vapor_pressure(real T)
{	
	const real T_crit = 617.7;
	const real T_min  = 243.51;

	real vapor_pressure;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_crit)
		T = T_crit;

	return vapor_pressure = exp(112.73 - 9749.6/T - 13.245*log(T) + 7.1266e-6*pow2(T));
}



real property_n_decane_vapor_diffusivity(real T, real p)
{


	real vapor_diffusivity;

	return vapor_diffusivity = 1.013e-2*pow(T , 1.75)*pow((M_N_DECANE+29.0)/(29.0*M_N_DECANE) , 0.5)/(p*pow2(pow(19.7,1/3)+pow(209.82,1/3)));
}



real property_n_decane_vapor_thermal_cond(real T)
{
	const real T_min =  447.3;
	const real T_max = 1000.0;

	real vapor_thermal_cond;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return vapor_thermal_cond = -668.40 * pow(T, 0.9323) / (1.0 - 4.0710e9/T);
}



real property_n_decane_vapor_viscosity(real T)
{
	const real T_min =  243.51;
	const real T_max = 1000.00;

	real vapor_viscosity;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return vapor_viscosity = 2.6400e-8 * pow(T, 0.9487) / (1.0 + 71.000/T);
}



DEFINE_DPM_PROPERTY(n_decane_liq_viscosity, c, t, p)
{
	return property_n_decane_liq_viscosity(P_T(p));
}



DEFINE_DPM_PROPERTY(n_decane_liq_surface_tension, c, t, p)
{	
	return property_n_decane_liq_surface_tension(P_T(p));
}



DEFINE_DPM_PROPERTY(n_decane_liq_specific_heat, c, t, p)
{
	return property_n_decane_liq_specific_heat(P_T(p));
}



DEFINE_DPM_PROPERTY(n_decane_liq_density, c, t, p)
{
	return property_n_decane_liq_density(P_T(p));
}



DEFINE_DPM_PROPERTY(n_decane_heat_vaporization, c, t, p)
{
	return property_n_decane_heat_vaporization(P_T(p));
}



DEFINE_DPM_PROPERTY(n_decane_vapor_pressure, c, t, p)
{	
	return property_n_decane_vapor_pressure(P_T(p));
}




DEFINE_PROPERTY(n_decane_vapor_diffusivity, c, t)
{
	return property_n_decane_vapor_diffusivity(C_T(c,t));
}




DEFINE_PROPERTY(n_decane_vapor_viscosity, c, t)
{
	return property_n_decane_vapor_viscosity(C_T(c,t));
}



DEFINE_PROPERTY(n_decane_vapor_thermal_cond, c, t)
{
	return property_n_decane_vapor_thermal_cond(C_T(c,t));
}

