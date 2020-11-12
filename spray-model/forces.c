#include "udf.h"
#include "dpm.h"
#include "math.h"
#include "math_ext.h"
#include "model.h"
#include "forces.h"
#include "interaction.h"





/*  drag coefficient of a circular disc  */
real DiscDragCoeff(const real Re, const real dr)
{
    real Re_dr = Re * dr;
    /*  Clift, Grace, Weber:  Bubbles, Drops, Particles  */
    if (Re_dr > 133.0)
        Re_dr = 133.0;

    const real cd = 64.0/M_PI/Re_dr * (1.0 + 0.138*pow(Re_dr,0.792));
    /*  fluent convention  */
    return 18.0/24.0*Re * cd;
}





/*  drag coefficient of a spheroid  */
real SpheroidDragCoeff(const real Re, const real dr)
{
    real Re_dr = Re * dr;
    /*  Clift, Grace, Weber:  Bubbles, Drops, Particles  */
    if (Re_dr > 1.0e4)
        Re_dr = 1.0e4;

    const real xd = -1.660 + 0.3958*log10(Re_dr) - 0.03*pow2(log10(Re_dr));
        
    const real cd =  108.0 * pow(Re_dr,xd);
    /*  fluent convention  */
    return 18.0/24.0*Re * cd;
}





DEFINE_DPM_BODY_FORCE(forces, particle, i)
{
	if(is_sleeping(particle))
		return  0.0;
	else
		return -C_P_G_EFF(P_CELL(particle),P_CELL_THREAD(particle),i)/particle->state.rho;
}
