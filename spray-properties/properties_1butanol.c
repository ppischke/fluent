#include "udf.h"
#include "dpm.h"
#include "math.h"

#include "properties_mixture.h"



#define M_1_BUTANOL    74.123



real property_1_butanol_liq_viscosity(real T)
{
	const real T_min = 190.00;
	const real T_max = 390.81;

	real viscosity;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return viscosity = exp(-35.426 + 3184.5/T - 3.2965*log(T) - 3.00e-27*pow(T,10));	
}



real property_1_butanol_liq_surface_tension(real T)
{	
/*	TO DO:
//	limit surface tension
*/

	const real T_min = 273.15;
	const real T_max = 413.15;
	
	real surface_tension;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return surface_tension = 4.9830e-2 - 8.5400e-5*T ;	
}



real property_1_butanol_liq_specific_heat(real T)
{
	const real T_min = 184.51;
	const real T_max = 390.81;

	real specific_heat;

	if (T < T_min)	
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

   	return specific_heat = (1.9120e5 - 730.40*T + 2.2998*pow2(T)) / M_1_BUTANOL;
}



real property_1_butanol_liq_density(real T)
{
	const real T_min = 184.51;
	const real T_max = 563.05;

	real density;

	if (T < T_min)	
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return density = 0.96500 / (pow(0.26660, 1.0+pow(1.0-(T/T_max), 0.24419))) * M_1_BUTANOL;
}



real property_1_butanol_heat_vaporization(real T)
{
	const real T_min = 184.51;
	const real T_max = 563.05;

	real heat_vaporization;

	if (T < T_min)
		T = T_min;
  	else
	if (T > T_max)
		T = T_max;

	return heat_vaporization = 6.7390e7 * pow(1.0-T/T_max, 0.17300 + 0.2915*T/T_max) / M_1_BUTANOL;
}



real property_1_butanol_vapor_pressure(real T)
{	
	const real T_crit = 563.05;
	const real T_min  = 184.51;

	real vapor_pressure;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_crit)
		T = T_crit;

	return vapor_pressure = exp(93.173 - 9185.9/T - 9.7464*log(T) + 4.7796e-18*pow(T , 6));
}



real property_1_butanol_vapor_diffusivity(real T, real p)
{


	real vapor_diffusivity;

	return vapor_diffusivity = 1.013e-2*pow(T , 1.75)*pow((M_1_BUTANOL+29.0)/(29.0*M_1_BUTANOL) , 0.5)/(p*pow2(pow(19.7,1/3)+pow(92.81,1/3)));
}



real property_1_butanol_vapor_thermal_cond(real T)
{
	const real T_min =  370.00;
	const real T_max = 800.00;

	real vapor_thermal_cond;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return vapor_thermal_cond = -4.4940e-2 * pow(T, 4.4600e-2) / (1.0 - 1355.2/T);
}



real property_1_butanol_vapor_viscosity(real T)
{
	const real T_min =  184.51;
	const real T_max = 1000.00;

	real vapor_viscosity;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return vapor_viscosity = 1.4031e-6 * pow(T, 0.4611) / (1.0 + 537.00/T);
}



DEFINE_DPM_PROPERTY(one_butanol_liq_viscosity, c, t, p)
{
	return property_1_butanol_liq_viscosity(P_T(p));
}



DEFINE_DPM_PROPERTY(one_butanol_liq_surface_tension, c, t, p)
{	
	return property_1_butanol_liq_surface_tension(P_T(p));
}



DEFINE_DPM_PROPERTY(one_butanol_liq_specific_heat, c, t, p)
{
	return property_1_butanol_liq_specific_heat(P_T(p));
}



DEFINE_DPM_PROPERTY(one_butanol_liq_density, c, t, p)
{
	return property_1_butanol_liq_density(P_T(p));
}



DEFINE_DPM_PROPERTY(one_butanol_heat_vaporization, c, t, p)
{
	return property_1_butanol_heat_vaporization(P_T(p));
}



DEFINE_DPM_PROPERTY(one_butanol_vapor_pressure, c, t, p)
{	
	return property_1_butanol_vapor_pressure(P_T(p));
}



DEFINE_PROPERTY(one_butanol_vapor_diffusivity, c, t)
{
	return property_1_butanol_vapor_diffusivity(C_T(c,t), C_P(c,t));
}



DEFINE_PROPERTY(one_butanol_vapor_viscosity, c, t)
{
	return property_1_butanol_vapor_viscosity(C_T(c,t));
}



DEFINE_PROPERTY(one_butanol_vapor_thermal_cond, c, t)
{
	return property_1_butanol_vapor_thermal_cond(C_T(c,t));
}


