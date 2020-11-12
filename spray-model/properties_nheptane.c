#include "udf.h"
#include "dpm.h"
#include "math.h"
#include "math_ext.h"



#define M_N_HEPTANE  100.204



real property_n_heptane_liq_viscosity(real T)
{   
    const real T_min = 182.57;
    const real T_max = 540.20;

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
/*  TO DO:
//  limit surface tension
*/

    const real T_min = 182.57;
    const real T_max = 540.20;
    
    real surface_tension;

    if (T < T_min)
        T = T_min;
    else
    if (T > T_max)
        T = T_max;

    return surface_tension = 0.054143*pow((1.0 - T/T_max), 1.2512); 
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

    return density = 0.61259 / (pow(0.28211, 1.0+pow(1.0-(T/T_max), 0.28141))) * M_N_HEPTANE;
}



real property_n_heptane_vapor_pressure(real T)
{   
    const real T_crit = 540.20;
    const real T_min  = 182.57;

    real vapor_pressure;

    if (T < T_min)
        T = T_min;
    else
    if (T > T_crit)
        T = T_crit;    

    return vapor_pressure = exp(87.829 - 6996.4/T - 9.8802*log(T) + 7.2099e-6*pow2(T));
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

    return vapor_thermal_cond = 7.0028e-2 * pow(T, 0.38068) / (1.0 - 7049.9/T);
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



DEFINE_DPM_PROPERTY(n_heptane_liq_viscosity, c, t, p, T)
{
    return property_n_heptane_liq_viscosity(T);
}



DEFINE_DPM_PROPERTY(n_heptane_liq_surface_tension, c, t, p, T)
{   
    return property_n_heptane_liq_surface_tension(T);
}



DEFINE_DPM_PROPERTY(n_heptane_liq_density, c, t, p, T)
{
    return property_n_heptane_liq_density(T);
}



DEFINE_DPM_PROPERTY(n_heptane_vapor_pressure, c, t, p, T)
{   
    return property_n_heptane_vapor_pressure(T);
}



DEFINE_PROPERTY(n_heptane_vapor_viscosity, c, t)
{
    return property_n_heptane_vapor_viscosity(C_T(c,t));
}



DEFINE_PROPERTY(n_heptane_vapor_thermal_cond, c, t)
{
    return property_n_heptane_vapor_thermal_cond(C_T(c,t));
}

