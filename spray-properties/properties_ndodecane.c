#include "udf.h"
#include "dpm.h"
#include "math.h"
#include "misc.h"

#include "properties_mixture.h"



#define M_N_DODECANE    170.338



real property_n_dodecane_liq_viscosity(real T)
{
	const real T_min = 263.57;
	const real T_max = 489.47;

	real viscosity;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return viscosity = exp(-20.607 + 1943.0/T + 1.3205*log(T));	
}



real property_n_dodecane_liq_surface_tension(real T)
{	
/*	TO DO:
//	limit surface tension
*/

	const real T_min = 263.57;
	const real T_max = 658.00;
	
	real surface_tension;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return surface_tension = 5.5493e-2*pow(1-T/T_max, 1.3262) ;	
}



real property_n_dodecane_liq_specific_heat(real T)
{
	const real T_min = 263.57;
	const real T_max = 330.00;

	real specific_heat;

	if (T < T_min)	
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

   	return specific_heat = (5.0821e5 - 1368.7*T + 3.1015*pow2(T)) / M_N_DODECANE;
}



real property_n_dodecane_liq_density(real T)
{
	const real T_min = 263.57;
	const real T_max = 658.00;

	real density;

	if (T < T_min)	
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return density = 0.35541 / (pow(0.25511, 1.0+pow(1.0-(T/T_max), 0.29368))) * M_N_DODECANE;
}



real property_n_dodecane_heat_vaporization(real T)
{
	const real T_min = 263.57;
	const real T_max = 658.00;

	real heat_vaporization;

	if (T < T_min)
		T = T_min;
  	else
	if (T > T_max)
		T = T_max;

	return heat_vaporization = 7.7337e7 * pow(1.0-T/T_max, 0.40681) / M_N_DODECANE;
}



real property_n_dodecane_vapor_pressure(real T)
{	
	const real T_crit = 658.00;
	const real T_min  = 263.57;

	real vapor_pressure;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_crit)
		T = T_crit;

	return vapor_pressure = exp(137.47 - 1.1976e4/T - 16.698*log(T) + 8.0906e-6*pow(T , 2));
}



real property_n_dodecane_vapor_diffusivity(real T, real p)
{


	real vapor_diffusivity;

	return vapor_diffusivity = 1.013e-2*pow(T , 1.75)*pow((M_N_DODECANE+29.0)/(29.0*M_N_DODECANE) , 0.5)/(p*pow2(pow(19.7,1/3)+pow(250.86,1/3)));
}



real property_n_dodecane_vapor_thermal_cond(real T)
{
	const real T_min =  489.47;
	const real T_max = 1000.00;

	real vapor_thermal_cond;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return vapor_thermal_cond = 5.7190e-6 * pow(T, 1.4699) / (1.0 - 579.40/T);
}



real property_n_dodecane_vapor_viscosity(real T)
{
	const real T_min =  263.57;
	const real T_max = 1000.00;

	real vapor_viscosity;

	if (T < T_min)
		T = T_min;
	else
	if (T > T_max)
		T = T_max;

	return vapor_viscosity = 6.3440e-8 * pow(T, 0.82870) / (1.0 + 219.50/T);
}



DEFINE_DPM_PROPERTY(n_dodecane_liq_viscosity, c, t, p)
{
	return property_n_dodecane_liq_viscosity(P_T(p));
}



DEFINE_DPM_PROPERTY(n_dodecane_liq_surface_tension, c, t, p)
{	
	return property_n_dodecane_liq_surface_tension(P_T(p));
}



DEFINE_DPM_PROPERTY(n_dodecane_liq_specific_heat, c, t, p)
{
	return property_n_dodecane_liq_specific_heat(P_T(p));
}



DEFINE_DPM_PROPERTY(n_dodecane_liq_density, c, t, p)
{
	return property_n_dodecane_liq_density(P_T(p));
}



DEFINE_DPM_PROPERTY(n_dodecane_heat_vaporization, c, t, p)
{
	return property_n_dodecane_heat_vaporization(P_T(p));
}



DEFINE_DPM_PROPERTY(n_dodecane_vapor_pressure, c, t, p)
{	
	return property_n_dodecane_vapor_pressure(P_T(p));
}



DEFINE_PROPERTY(n_dodecane_vapor_diffusivity, c, t)
{
	return property_n_dodecane_vapor_diffusivity(C_T(c,t));
}



DEFINE_PROPERTY(n_dodecane_vapor_viscosity, c, t)
{
	return property_n_dodecane_vapor_viscosity(C_T(c,t));
}



DEFINE_PROPERTY(n_dodecane_vapor_thermal_cond, c, t)
{
	return property_n_dodecane_vapor_thermal_cond(C_T(c,t));
}



