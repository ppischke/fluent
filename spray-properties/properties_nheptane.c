#include "udf.h"
#include "dpm.h"
#include "math.h"
#include "misc.h"

#include "properties_mixture.h"



#define M_N_HEPTANE    100.204



real property_n_heptane_liq_viscosity(real T)
{
	const real T_min = 182.57;
	const real T_max = 373.15;

	real viscosity;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return viscosity = exp(-24.451 + 1533.1/T + 2.0087*log(T));	
}



real property_n_heptane_liq_surface_tension(real T)
{	
/*	TO DO:
//	limit surface tension
*/

	const real T_min = 182.57;
	const real T_max = 540.20;
	
	real surface_tension;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return surface_tension = 5.4143e-2*pow(1-T/T_max, 1.2512) ;	
}



real property_n_heptane_liq_specific_heat(real T)
{
	const real T_min = 182.57;
	const real T_max = 520.00;

	real specific_heat;

	if (T < T_min)	
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

   	return specific_heat = (pow2(61.260)/(1-T/T_max) + 3.1441e5 - 2*61.260*1824.6*(1-T/T_max) + 61.260*2547.9*pow2(1 - T/T_max) - pow2(1824.6)*pow3(1 - T/T_max)/3 + 1824.6*2547.9*pow(1 - T/T_max,4)/2) / M_N_HEPTANE;
}



real property_n_heptane_liq_density(real T)
{
	const real T_min = 182.57;
	const real T_max = 540.20;

	real density;

	if (T < T_min)	
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return density = 0.61259 / (pow(0.26211, 1.0+pow(1.0-(T/T_max), 0.28141))) * M_N_HEPTANE;
}



real property_n_heptane_heat_vaporization(real T)
{
	const real T_min = 182.57;
	const real T_max = 540.20;

	real heat_vaporization;

	if (T < T_min)
		T = T_min;
  	else
	if (T > T_max)
		T = T_max;

	return heat_vaporization = 5.0014e7 * pow(1.0-T/T_max, 0.38795) / M_N_HEPTANE;
}



real property_n_heptane_vapor_pressure(real T)
{	
	const real T_crit = 182.57;
	const real T_min  = 540.20;

	real vapor_pressure;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_crit)
		T = T_crit;

	return vapor_pressure = exp(87.829 - 6996.4/T - 9.8802*log(T) + 7.2099e-6*pow(T , 2));
}



real property_n_heptane_vapor_diffusivity(real T, real p)
{


	real vapor_diffusivity;

	return vapor_diffusivity = 1.013e-2*pow(T , 1.75)*pow((M_N_HEPTANE+29.0)/(29.0*M_N_HEPTANE) , 0.5)/(p*pow2(pow(19.7,1/3)+pow(148.26,1/3)));
}



real property_n_heptane_vapor_thermal_cond(real T)
{
	const real T_min =  339.15;
	const real T_max = 1000.00;

	real vapor_thermal_cond;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return vapor_thermal_cond = -7.0028e-2 * pow(T, 0.38068) / (1.0 - 7049.9/T - 2.4005e6/pow2(T));
}



real property_n_heptane_vapor_viscosity(real T)
{
	const real T_min =  182.57;
	const real T_max = 1000.00;

	real vapor_viscosity;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return vapor_viscosity = 6.6720e-8 * pow(T, 0.82837) / (1.0 + 85.752/T);
}



DEFINE_DPM_PROPERTY(n_heptane_liq_viscosity, c, t, p)
{
	return property_n_heptane_liq_viscosity(P_T(p));
}



DEFINE_DPM_PROPERTY(n_heptane_liq_surface_tension, c, t, p)
{	
	return property_n_heptane_liq_surface_tension(P_T(p));
}



DEFINE_DPM_PROPERTY(n_heptane_liq_specific_heat, c, t, p)
{
	return property_n_heptane_liq_specific_heat(P_T(p));
}



DEFINE_DPM_PROPERTY(n_heptane_liq_density, c, t, p)
{
	return property_n_heptane_liq_density(P_T(p));
}



DEFINE_DPM_PROPERTY(n_heptane_heat_vaporization, c, t, p)
{
	return property_n_heptane_heat_vaporization(P_T(p));
}



DEFINE_DPM_PROPERTY(n_heptane_vapor_pressure, c, t, p)
{	
	return property_n_heptane_vapor_pressure(P_T(p));
}



DEFINE_PROPERTY(n_heptane_vapor_diffusivity, c, t)
{
	return property_n_heptane_vapor_diffusivity(C_T(c,t));
}



DEFINE_PROPERTY(n_heptane_vapor_viscosity, c, t)
{
	return property_n_heptane_vapor_viscosity(C_T(c,t));
}



DEFINE_PROPERTY(n_heptane_vapor_thermal_cond, c, t)
{
	return property_n_heptane_vapor_thermal_cond(C_T(c,t));
}



