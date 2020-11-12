#include "udf.h"
#include "dpm.h"
#include "math.h"

#include "properties_mixture.h"



#define M_TETRAHYDROFURFURYL_ALCOHOL    102.133



real property_tetrahydrofurfuryl_alcohol_liq_viscosity(real T)
{
	const real T_min = 220.00;
	const real T_max = 503.00;

	real viscosity;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return viscosity = exp(-7.9742 + 2745.4/T - 1.1468*log(T));	
}



real property_tetrahydrofurfuryl_alcohol_liq_surface_tension(real T)
{	
/*	TO DO:
//	limit surface tension
*/

	const real T_min = 294.75;
	const real T_max = 450.80;
	
	real surface_tension;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return surface_tension = 6.5073e-2 - 8.8258e-5*T - 1.6138e-8*pow(T,2.0);	
}



real property_tetrahydrofurfuryl_alcohol_liq_specific_heat(real T)
{
	const real T_min = 296.65;
	const real T_max = 469.50;

	real specific_heat;

	if (T < T_min)	
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

   	return specific_heat = (52700.0 + 435.80*T) / M_TETRAHYDROFURFURYL_ALCOHOL;
}



real property_tetrahydrofurfuryl_alcohol_liq_density(real T)
{
	const real T_min = 193.00;
	const real T_max = 639.00;

	real density;

	if (T < T_min)	
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return density = 0.97001 / (pow(0.2813, 1.0+pow(1.0-(T/T_max), 0.23837))) * M_TETRAHYDROFURFURYL_ALCOHOL;
}



real property_tetrahydrofurfuryl_alcohol_heat_vaporization(real T)
{
	const real T_min = 193.00;
	const real T_max = 639.00;

	real heat_vaporization;

	if (T < T_min)
		T = T_min;
  	else
	if (T > T_max)
		T = T_max;

	return heat_vaporization = 6.4109e7 * pow(1.0-T/T_max, 0.28538) / M_TETRAHYDROFURFURYL_ALCOHOL;
}



real property_tetrahydrofurfuryl_alcohol_vapor_pressure(real T)
{	
	const real T_crit = 639.00;
	const real T_min  = 193.00;

	real vapor_pressure;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_crit)
		T = T_crit;

	return vapor_pressure = exp(150.59 - 11574.0/T - 19.025*log(T) + 1.4141e-5*pow(T,2.0));
}



real property_tetrahydrofurfuryl_alcohol_vapor_diffusivity(real T, real p)
{


	real vapor_diffusivity;

	return vapor_diffusivity = 1.013e-2*pow(T , 1.75)*pow((M_TETRAHYDROFURFURYL_ALCOHOL+29.0)/(29.0*M_TETRAHYDROFURFURYL_ALCOHOL) , 0.5)/(p*pow(pow(19.7,1/3)+pow(114.82,1/3),2.0));
}



real property_tetrahydrofurfuryl_alcohol_vapor_thermal_cond(real T)
{
	const real T_min =  450.80;
	const real T_max = 1000.00;

	real vapor_thermal_cond;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return vapor_thermal_cond = 2.8570e-4 * pow(T, 0.90229) / (1.0 + 702.94/T);
}



real property_tetrahydrofurfuryl_alcohol_vapor_viscosity(real T)
{
	const real T_min =  193.00;
	const real T_max = 1000.00;

	real vapor_viscosity;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return vapor_viscosity = 1.6196e-7 * pow(T, 0.74453) / (1.0 + 145.79/T);
}



DEFINE_DPM_PROPERTY(tetrahydrofurfuryl_alcohol_liq_viscosity, c, t, p)
{
	return property_tetrahydrofurfuryl_alcohol_liq_viscosity(P_T(p));
}



DEFINE_DPM_PROPERTY(tetrahydrofurfuryl_alcohol_liq_surface_tension, c, t, p)
{	
	return property_tetrahydrofurfuryl_alcohol_liq_surface_tension(P_T(p));
}



DEFINE_DPM_PROPERTY(tetrahydrofurfuryl_alcohol_liq_specific_heat, c, t, p)
{
	return property_tetrahydrofurfuryl_alcohol_liq_specific_heat(P_T(p));
}



DEFINE_DPM_PROPERTY(tetrahydrofurfuryl_alcohol_liq_density, c, t, p)
{
	return property_tetrahydrofurfuryl_alcohol_liq_density(P_T(p));
}



DEFINE_DPM_PROPERTY(tetrahydrofurfuryl_alcohol_heat_vaporization, c, t, p)
{
	return property_tetrahydrofurfuryl_alcohol_heat_vaporization(P_T(p));
}



DEFINE_DPM_PROPERTY(tetrahydrofurfuryl_alcohol_vapor_pressure, c, t, p)
{	
	return property_tetrahydrofurfuryl_alcohol_vapor_pressure(P_T(p));
}



DEFINE_PROPERTY(tetrahydrofurfuryl_alcohol_vapor_diffusivity, c, t)
{
	return property_tetrahydrofurfuryl_alcohol_vapor_diffusivity(C_T(c,t), C_P(c,t));
}



DEFINE_PROPERTY(tetrahydrofurfuryl_alcohol_vapor_viscosity, c, t)
{
	return property_tetrahydrofurfuryl_alcohol_vapor_viscosity(C_T(c,t));
}



DEFINE_PROPERTY(tetrahydrofurfuryl_alcohol_vapor_thermal_cond, c, t)
{
	return property_tetrahydrofurfuryl_alcohol_vapor_thermal_cond(C_T(c,t));
}

