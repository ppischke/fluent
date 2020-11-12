#include "udf.h"
#include "dpm.h"
#include "math.h"
#include "math_ext.h"



#define M_N_HEXANE  86.177;



real property_n_hexane_liq_viscosity(real T)
{   
    const real T_min = 177.84;
    const real T_max = 507.43;

    real viscosity;

    if (T < T_min)
        T = T_min;
    else
    if (T > T_max)
        T = T_max;

    return viscosity = exp(-19.116 + 1154.3/T + 1.2482*log(T));     
}



real property_n_hexane_liq_surface_tension(real T)
{   
/*  TO DO:
//  limit surface tension
*/

    const real T_min = 177.84;
    const real T_max = 507.43;
    
    real surface_tension;

    if (T < T_min)
        T = T_min;
    else
    if (T > T_max)
        T = T_max;

    return surface_tension = 0.056081*pow((1.0 - T/T_max), 1.2843); 
}



real property_n_hexane_liq_density(real T)
{
    const real T_min = 177.84;
    const real T_max = 507.43;

    real density;

    if (T < T_min)  
        T = T_min;
    else
    if (T > T_max)
        T = T_max;

    return density = 0.71470 / (pow(0.26500, 1.0+pow(1.0-(T/T_max), 0.27810))) * M_N_HEXANE;
}



real property_n_hexane_vapor_pressure(real T)
{   
    const real T_crit = 507.43;
    const real T_min  = 177.84;

    real vapor_pressure;

    if (T < T_min)
        T = T_min;
    else
    if (T > T_crit)
        T = T_crit;    

    return vapor_pressure = exp(165.47 - 8353.3/T - 23.927*log(T) + 2.9496e-2*T);
}



real property_n_hexane_vapor_thermal_cond(real T)
{
    const real T_min =  290.00;
    const real T_max =  480.00;

    real vapor_thermal_cond;

    if (T < T_min)
        T = T_min;
    else
    if (T > T_max)
        T = T_max;

    return vapor_thermal_cond = 3.6780e-6 * pow(T, 1.5677) / (1.0 + 360.00/T);
}



real property_n_hexane_vapor_viscosity(real T)
{
    const real T_min =  300.00;
    const real T_max = 1000.00;

    real vapor_viscosity;

    if (T < T_min)
        T = T_min;
    else
    if (T > T_max)
        T = T_max;

    return vapor_viscosity = 2.9400e-6 * pow(T, 0.35370) / (1.0 + 718.30/T);
}



DEFINE_DPM_PROPERTY(n_hexane_liq_viscosity, c, t, p, T)
{
    return property_n_hexane_liq_viscosity(T);
}



DEFINE_DPM_PROPERTY(n_hexane_liq_surface_tension, c, t, p, T)
{   
    return property_n_hexane_liq_surface_tension(T);
}



DEFINE_DPM_PROPERTY(n_hexane_liq_density, c, t, p, T)
{
    return property_n_hexane_liq_density(T);
}



DEFINE_DPM_PROPERTY(n_hexane_vapor_pressure, c, t, p, T)
{   
    return property_n_hexane_vapor_pressure(T);
}



DEFINE_PROPERTY(n_hexane_vapor_viscosity, c, t)
{
    return property_n_hexane_vapor_viscosity(C_T(c,t));
}



DEFINE_PROPERTY(n_hexane_vapor_thermal_cond, c, t)
{
    return property_n_hexane_vapor_thermal_cond(C_T(c,t));
}

