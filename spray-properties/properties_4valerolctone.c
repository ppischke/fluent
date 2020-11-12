#include "udf.h"
#include "dpm.h"
#include "math.h"

#include "properties_mixture.h"



#define M_4_VALEROLACTONE    100.117



real property_4_valerolactone_liq_viscosity(real T)
{
	const real T_min = 242.15;
	const real T_max = 572.15;

	real viscosity;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return viscosity = exp(-2.8539 + 563.34/T - 1.1445*log(T));	
}



real property_4_valerolactone_liq_surface_tension(real T)
{	
/*	TO DO:
//	limit surface tension
*/

	const real T_min = 242.15;
	const real T_max = 727.00;
	
	real surface_tension;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return surface_tension = 6.6815e-2*pow(1-T/T_max, 1.3153) ;	
}



real property_4_valerolactone_liq_specific_heat(real T)
{
	const real T_min = 363.50;
	const real T_max = 543.50;

	real specific_heat;

	if (T < T_min)	
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

   	return specific_heat = (9.5800e4 + 277.50*T) / M_4_VALEROLACTONE;
}



real property_4_valerolactone_liq_density(real T)
{
	const real T_min = 242.15;
	const real T_max = 727.00;

	real density;

	if (T < T_min)	
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return density = 0.98680 / (pow(0.27532, 1.0+pow(1.0-(T/T_max), 0.34900))) * M_4_VALEROLACTONE;
}



real property_4_valerolactone_heat_vaporization(real T)
{
	const real T_min = 242.15;
	const real T_max = 747.00;

	real heat_vaporization;

	if (T < T_min)
		T = T_min;
  	else
	if (T > T_max)
		T = T_max;

	return heat_vaporization = 6.1320e7 * pow(1.0-T/T_max, 0.29680) / M_4_VALEROLACTONE;
}



real property_4_valerolactone_vapor_pressure(real T)
{	
	const real T_crit = 242.15;
	const real T_min  = 727.00;

	real vapor_pressure;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_crit)
		T = T_crit;

	return vapor_pressure = exp(83.640 - 8.7843e3/T - 8.8878*log(T) + 4.5338e-6*pow(T , 2));
}



real property_4_valerolactone_vapor_diffusivity(real T, real p)
{


	real vapor_diffusivity;

	return vapor_diffusivity = 1.013e-2*pow(T , 1.75)*pow((M_4_VALEROLACTONE+29.0)/(29.0*M_4_VALEROLACTONE) , 0.5)/(p*pow2(pow(19.7,1/3)+pow(110.2,1/3)));
}



real property_4_valerolactone_vapor_thermal_cond(real T)
{
	const real T_min =  480.65;
	const real T_max = 1000.00;

	real vapor_thermal_cond;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return vapor_thermal_cond = 3.8300e-4 * pow(T, 0.87240) / (1.0 - 705.50/T);
}



real property_4_valerolactone_vapor_viscosity(real T)
{
	const real T_min =  242.15;
	const real T_max = 1000.00;

	real vapor_viscosity;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return vapor_viscosity = 8.0120e-8 * pow(T, 0.85930) / (1.0 + 61.350/T);
}



DEFINE_DPM_PROPERTY(four_valerolactone_liq_viscosity, c, t, p)
{
	return property_4_valerolactone_liq_viscosity(P_T(p));
}



DEFINE_DPM_PROPERTY(four_valerolactone_liq_surface_tension, c, t, p)
{	
	return property_4_valerolactone_liq_surface_tension(P_T(p));
}



DEFINE_DPM_PROPERTY(four_valerolactone_liq_specific_heat, c, t, p)
{
	return property_4_valerolactone_liq_specific_heat(P_T(p));
}



DEFINE_DPM_PROPERTY(four_valerolactone_liq_density, c, t, p)
{
	return property_4_valerolactone_liq_density(P_T(p));
}



DEFINE_DPM_PROPERTY(four_valerolactone_heat_vaporization, c, t, p)
{
	return property_4_valerolactone_heat_vaporization(P_T(p));
}



DEFINE_DPM_PROPERTY(four_valerolactone_vapor_pressure, c, t, p)
{	
	return property_4_valerolactone_vapor_pressure(P_T(p));
}



DEFINE_PROPERTY(four_valerolactone_vapor_diffusivity, c, t)
{
	return property_4_valerolactone_vapor_diffusivity(C_T(c,t), C_P(c,t));
}



DEFINE_PROPERTY(four_valerolactone_vapor_viscosity, c, t)
{
	return property_4_valerolactone_vapor_viscosity(C_T(c,t));
}



DEFINE_PROPERTY(four_valerolactone_vapor_thermal_cond, c, t)
{
	return property_4_valerolactone_vapor_thermal_cond(C_T(c,t));
}



