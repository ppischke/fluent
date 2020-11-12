#include "udf.h"
#include "dpm.h"
#include "math.h"
#include "math_ext.h"



#define M_N_DECANE  142.285                             /* Molmasse in g/mol */



real property_n_decane_liq_viscosity(real T)                        /* liq-Viskosität in Pa*s */
{
    const real T_min = 243.51;
    const real T_max = 617.70;

    real viscosity;

    if (T < T_min)
        T = T_min;
    else
    if (T > T_max)
        T = T_max;

    return viscosity = exp(-16.468 + 1533.5/T + 0.7511*log(T));
}



real property_n_decane_liq_surface_tension(real T)                  /* Oberflächenspannung in N/m */
{
/*  TO DO:
//  limit surface tension
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



real property_n_decane_liq_density(real T)                      /* liq-Dichte in kg/m^3 */
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



real property_n_decane_vapor_pressure(real T)                       /* Dampfdruck in Pa */
{
    const real T_crit = 617.7;
    const real T_min  = 243.51;

    real vapor_pressure;

    if (T < T_min)
        T = T_min;
    else
    if (T > T_crit)
        T = T_crit;

    return vapor_pressure = exp(112.73 - 9749.6/T - 13.245*log(T) + (7.1266e-6)*(pow2(T)));
}



real property_n_decane_vapor_viscosity(real T)                      /* gas-Dampfviskosität in Pa*s */
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



real property_n_decane_vapor_thermal_cond(real T)                   /* gas-Wärmeleitfähigkeit in W/mK */
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





DEFINE_DPM_PROPERTY(n_decane_liq_viscosity, c, t, p, T)
{
    return property_n_decane_liq_viscosity(T);
}



DEFINE_DPM_PROPERTY(n_decane_liq_surface_tension, c, t, p, T)
{
    return property_n_decane_liq_surface_tension(T);
}



DEFINE_DPM_PROPERTY(n_decane_liq_density, c, t, p, T)
{
    return property_n_decane_liq_density(T);
}



DEFINE_DPM_PROPERTY(n_decane_vapor_pressure, c, t, p, T)
{   
    return property_n_decane_vapor_pressure(T);
}




DEFINE_PROPERTY(n_decane_vapor_viscosity, c, t)
{
    return property_n_decane_vapor_viscosity(C_T(c,t));
}



DEFINE_PROPERTY(n_decane_vapor_thermal_cond, c, t)
{
    return property_n_decane_vapor_thermal_cond(C_T(c,t));
}
