#include "udf.h"
#include "dpm.h"
#include "math.h"
#include "math_ext.h"



#define M_N_1METHYL  142.200                                /* Molmasse in g/mol */



real property_n_1methyl_liq_viscosity(real T)                       /* liq-Viskosität in Pa*s */
{
    const real T_min = 242.67;
    const real T_max = 772.00;

    real viscosity;

    if (T < T_min)
        T = T_min;
    else
    if (T > T_max)
        T = T_max;

    return viscosity = exp(-91.835 + 5709.3/T + 11.735*log(T));
}



real property_n_1methyl_liq_surface_tension(real T)                 /* Oberflächenspannung in N/m */
{
/*  TO DO:
//  limit surface tension
*/

    const real T_min = 242.67;
    const real T_max = 772.00;

    real surface_tension;

    if (T < T_min)
        T = T_min;
    else
    if (T > T_max)
        T = T_max;

    return surface_tension = 0.07539*pow((1.0-T/T_max), 1.2925);
}



real property_n_1methyl_liq_density(real T)                     /* liq-Dichte in kg/m^3 */
{
    const real T_min = 242.67;
    const real T_max = 772.00;

    real density;

    if (T < T_min)
        T = T_min;
    else
    if (T > T_max)
        T = T_max;

    return density = 0.54778 / (pow(0.25357, 1.0+pow(1.0-(T/T_max), 0.28041))) * M_N_1METHYL;
}



real property_n_1methyl_vapor_pressure(real T)                      /* Dampfdruck in Pa */
{
    const real T_crit = 772.00;
    const real T_min  = 242.67;

    real vapor_pressure;

    if (T < T_min)
        T = T_min;
    else
    if (T > T_crit)
        T = T_crit;

    return vapor_pressure = exp(67.566 - 8737.0/T - 6.3362*log(T) + (1.6377e-6)*(pow2(T)));
}



real property_n_1methyl_vapor_viscosity(real T)                     /* gas-Dampfviskosität in Pa*s */
{
    const real T_min =  242.67;
    const real T_max = 1000.00;

    real vapor_viscosity;

    if (T < T_min)
        T = T_min;
    else
    if (T > T_max)
        T = T_max;

    return vapor_viscosity = 2.6213e-7 * pow(T, 0.64259) / (1.0 + 235.21/T);
}



real property_n_1methyl_vapor_thermal_cond(real T)                  /* gas-Wärmeleitfähigkeit in W/mK */
{
    const real T_min =  517.83;
    const real T_max = 1000.0;

    real vapor_thermal_cond;

    if (T < T_min)
        T = T_min;
    else
    if (T > T_max)
        T = T_max;

    return vapor_thermal_cond = 1.4947 * pow(T, -0.32034) / (1.0 - 1.1200e3/T + 2.6315e6/pow2(T));
}





DEFINE_DPM_PROPERTY(n_1methyl_liq_viscosity, c, t, p, T)
{
    return property_n_1methyl_liq_viscosity(T);
}



DEFINE_DPM_PROPERTY(n_1methyl_liq_surface_tension, c, t, p, T)
{
    return property_n_1methyl_liq_surface_tension(T);
}



DEFINE_DPM_PROPERTY(n_1methyl_liq_density, c, t, p, T)
{
    return property_n_1methyl_liq_density(T);
}



DEFINE_DPM_PROPERTY(n_1methyl_vapor_pressure, c, t, p, T)
{   
    return property_n_1methyl_vapor_pressure(T);
}



DEFINE_PROPERTY(n_1methyl_vapor_viscosity, c, t)
{
    return property_n_1methyl_vapor_viscosity(C_T(c,t));
}



DEFINE_PROPERTY(n_1methyl_vapor_thermal_cond, c, t)
{
    return property_n_1methyl_vapor_thermal_cond(C_T(c,t));
}
