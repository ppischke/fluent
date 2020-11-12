#include "udf.h"
#include "dpm.h"



const real M_air = 28.966;



real property_air_thermal_cond(real T)
{
    const real T_min =   70.00;
    const real T_max = 1000.00;

    real vapor_thermal_cond;

    if (T < T_min)
        T = T_min;
    else
    if (T > T_max)
        T = T_max;

    return vapor_thermal_cond = 3.1417e-4 * pow(T, 0.7786) / (1.0 - 0.7116/T + 2.1217e3/pow(T,2));
}



real property_air_viscosity(real T)
{
    const real T_min =   80.00;
    const real T_max = 1000.00;

    real vapor_thermal_cond;

    if (T < T_min)
        T = T_min;
    else
    if (T > T_max)
        T = T_max;

    return vapor_thermal_cond = 1.425e-6 * pow(T, 0.5039) / (1.0 + 108.30/T);
}



DEFINE_PROPERTY(air_thermal_cond, c, t)
{
    return property_air_thermal_cond(C_T(c,t));
}



DEFINE_PROPERTY(air_viscosity, c, t)
{
    return property_air_viscosity(C_T(c,t));
}
