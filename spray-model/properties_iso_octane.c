#include "udf.h"
#include "dpm.h"
#include "math.h"
#include "math_ext.h"



#define M_ISO_OCTANE    114.231



real property_iso_octane_liq_viscosity(real T)
{
    const real T_min = 165.78;
    const real T_max = 543.8;

    real viscosity;

    if (T < T_min)
        T = T_min;
    else
    if (T > T_max)
        T = T_max;

    return viscosity = exp(-12.928 + 1137.5/T + 0.25725*log(T)); /*... - 3.6929e-28*pow(T,10.0)); */
}



real property_iso_octane_liq_surface_tension(real T)
{   
/*  TO DO:
//  limit surface tension
*/

    const real T_min = 165.78;
    const real T_max = 543.8;
    
    real surface_tension;

    if (T < T_min)
        T = T_min;
    else
    if (T > T_max)
        T = T_max;

    return surface_tension = 0.047675*pow((1.0 - T/T_max), 1.2018); 
}



real property_iso_octane_liq_specific_heat(real T)
{
    const real T_min = 165.78;
    const real T_max = 543.8;

    real specific_heat;

    if (T < T_min)  
        T = T_min;
    else
    if (T > T_max)
        T = T_max;

    return specific_heat = (95275.0 + 696.7*T - 1.3765*pow2(T) + 2.1734e-3*pow3(T)) / M_ISO_OCTANE;
}



real property_iso_octane_liq_density(real T)
{
    const real T_min = 165.78;
    const real T_max = 543.8;

    real density;

    if (T < T_min)  
        T = T_min;
    else
    if (T > T_max)
        T = T_max;

    return density = 0.59059 / (pow(0.27424, 1.0+pow(1.0-(T/T_max), 0.2847))) * M_ISO_OCTANE;
}



real property_iso_octane_heat_vaporization(real T)
{
    const real T_min = 165.78;
    const real T_max = 543.8;

    real heat_vaporization;

    if (T < T_min)
        T = T_min;
    else
    if (T > T_max)
        T = T_max;

    return heat_vaporization = 47711000.0 * pow(1.0-T/T_max, 0.37949) / M_ISO_OCTANE;
}



real property_iso_octane_vapor_pressure(real T)
{   
    const real T_crit = 543.8;
    const real T_min  = 165.78;

    real vapor_pressure;

    if (T < T_min)
        T = T_min;
    else
    if (T > T_crit)
        T = T_crit;

    return vapor_pressure = exp(84.912 - 6722.2/T - 9.5157*log(T) + 7.2244e-6*pow2(T));
}



real property_iso_octane_vapor_thermal_cond(real T)
{
    const real T_min =  355.15;
    const real T_max = 1000.00;

    real vapor_thermal_cond;

    if (T < T_min)
        T = T_min;
    else
    if (T > T_max)
        T = T_max;

    return vapor_thermal_cond = 1.7580e-5 * pow(T, 1.3114) / (1.0 + 392.9/T);
}



real property_iso_octane_vapor_viscosity(real T)
{
    const real T_min =  165.78;
    const real T_max = 1000.00;

    real vapor_viscosity;

    if (T < T_min)
        T = T_min;
    else
    if (T > T_max)
        T = T_max;

    return vapor_viscosity = 1.1070e-7 * pow(T, 0.7460) / (1.0 + 72.400/T);
}



DEFINE_DPM_PROPERTY(iso_octane_liq_viscosity, c, t, p, T)
{
    return property_iso_octane_liq_viscosity(T);
}



DEFINE_DPM_PROPERTY(iso_octane_liq_surface_tension, c, t, p, T)
{   
    return property_iso_octane_liq_surface_tension(T);
}



DEFINE_DPM_PROPERTY(iso_octane_liq_density, c, t, p, T)
{
    return property_iso_octane_liq_density(T);
}



DEFINE_DPM_PROPERTY(iso_octane_vapor_pressure, c, t, p, T)
{   
    return property_iso_octane_vapor_pressure(T);
}



DEFINE_PROPERTY(iso_octane_vapor_viscosity, c, t)
{
    return property_iso_octane_vapor_viscosity(C_T(c,t));
}



DEFINE_PROPERTY(iso_octane_vapor_thermal_cond, c, t)
{
    return property_iso_octane_vapor_thermal_cond(C_T(c,t));
}

