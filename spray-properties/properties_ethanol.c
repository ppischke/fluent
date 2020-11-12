#include "udf.h"
#include "dpm.h"
#include "math.h"
#include "misc.h"

#include "properties_mixture.h"



#define M_ETHANOL    46.069



real property_ethanol_liq_viscosity(real T)
{
	const real T_min = 200.00;
	const real T_max = 440.00;

	real viscosity;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return viscosity = exp(7.8750 + 781.98/T - 3.0418*log(T));	
}



real property_ethanol_liq_surface_tension(real T)
{	
/*	TO DO:
//	limit surface tension
*/

	const real T_min = 273.15;
	const real T_max = 503.15;
	
	real surface_tension;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return surface_tension = 3.7640e-2 - 2.1570e-5*T - 1.0250e-7*pow2(T) ;	
}



real property_ethanol_liq_specific_heat(real T)
{
	const real T_min = 159.05;
	const real T_max = 390.00;

	real specific_heat;

	if (T < T_min)	
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

   	return specific_heat = (1.0264e5 - 139.63*T - 3.0341e-2*pow2(T) + 2.0386e-3*pow3(T)) / M_ETHANOL;
}



real property_ethanol_liq_density(real T)
{
	const real T_min = 159.05;
	const real T_max = 513.92;

	real density;

	if (T < T_min)	
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return density = 1.6480 / (pow(0.27627, 1.0+pow(1.0-(T/T_max), 0.23310))) * M_ETHANOL;
}



real property_ethanol_heat_vaporization(real T)
{
	const real T_min = 159.05;
	const real T_max = 513.92;

	real heat_vaporization;

	if (T < T_min)
		T = T_min;
  	else
	if (T > T_max)
		T = T_max;

	return heat_vaporization = 5.6900e7 * pow(1.0-T/T_max, 0.33590) / M_ETHANOL;
}



real property_ethanol_vapor_pressure(real T)
{	
	const real T_crit = 513.92;
	const real T_min  = 159.05;

	real vapor_pressure;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_crit)
		T = T_crit;

	return vapor_pressure = exp(74.475 - 7164.3/T - 7.3270*log(T) + 3.1340e-6*pow(T , 2));
}



real property_ethanol_vapor_diffusivity(real T, real p)
{


	real vapor_diffusivity;

	return vapor_diffusivity = 1.013e-2*pow(T , 1.75)*pow((M_ETHANOL+29.0)/(29.0*M_ETHANOL) , 0.5)/(p*pow2(pow(19.7,1/3)+pow(51.77,1/3)));
}



real property_ethanol_vapor_thermal_cond(real T)
{
	const real T_min =  293.15;
	const real T_max = 1000.00;

	real vapor_thermal_cond;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return vapor_thermal_cond = -1.0109e-2 * pow(T, 0.64750) / (1.0 - 7.3320e3/T - 2.6800e5/pow2(T));
}



real property_ethanol_vapor_viscosity(real T)
{
	const real T_min =  200.00;
	const real T_max = 1000.00;

	real vapor_viscosity;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return vapor_viscosity = 1.0613e-7 * pow(T, 0.80660) / (1.0 + 52.700/T);
}



DEFINE_DPM_PROPERTY(ethanol_liq_viscosity, c, t, p)
{
	return property_ethanol_liq_viscosity(P_T(p));
}



DEFINE_DPM_PROPERTY(ethanol_liq_surface_tension, c, t, p)
{	
	return property_ethanol_liq_surface_tension(P_T(p));
}



DEFINE_DPM_PROPERTY(ethanol_liq_specific_heat, c, t, p)
{
	return property_ethanol_liq_specific_heat(P_T(p));
}



DEFINE_DPM_PROPERTY(ethanol_liq_density, c, t, p)
{
	return property_ethanol_liq_density(P_T(p));
}



DEFINE_DPM_PROPERTY(ethanol_heat_vaporization, c, t, p)
{
	return property_ethanol_heat_vaporization(P_T(p));
}



DEFINE_DPM_PROPERTY(ethanol_vapor_pressure, c, t, p)
{	
	return property_ethanol_vapor_pressure(P_T(p));
}



DEFINE_PROPERTY(ethanol_vapor_diffusivity, c, t)
{
	return property_ethanol_vapor_diffusivity(C_T(c,t));
}



DEFINE_PROPERTY(ethanol_vapor_viscosity, c, t)
{
	return property_ethanol_vapor_viscosity(C_T(c,t));
}



DEFINE_PROPERTY(ethanol_vapor_thermal_cond, c, t)
{
	return property_ethanol_vapor_thermal_cond(C_T(c,t));
}



