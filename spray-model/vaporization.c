#include "udf.h"
#include "dpm.h"
#include "math.h"
#include "math_ext.h"
#include "model.h"
#include "vaporization.h"



/*  HEAT/MASS-TRANSFER CORRELATIONS
//      abramzon-sirignano
//      ackermann
*/
real correlation_abramzon_sirignano(const real Re_d, const real Pr, const real B)
{
    /*  heat transfer correlation by Abramzon and Sirignano  */

    /*  Sirignano: P. 37ff, correlation by Ranz-Marshall  */
    const real k = 0.848;
    /*  Sirignano: Eqns. 2.53ff  */
    real f;
    if(B == 0.0)
        f = 1.0;
    else
        f = log(1.0+B)/B;

    real Nu;
	if(Re_d < 5.0)
		Nu = 2.0*f + k*(pow(1.0+2.0*Re_d*Pr, 1.0/3.0) * fmax(1.0, pow(2.0*Re_d, 0.077)) - 1.0) / pow(1.0+B, 0.7);
	else
		Nu = 2.0*f + k*(pow(Re_d, 0.5) * pow(Pr, 1.0/3.0)) / pow(1.0+B, 0.7);

    /*  return Nusselt number  */
    return Nu;
}



real correlation_ackermann(const real Re_d, const real Pr, const real B)
{
    /*  heat transfer correlation with Ackermann correction  */
    real f;
    if(B == 0.0)
        f = 1.0;
    else
        f = log(1.0+B)/B;

    real Nu = f * (2.0 + (0.4*pow(Re_d, 1.0/2.0) + 0.06*pow(Re_d, 2.0/3.0)) * pow(Pr, 0.4));

    /*  return Nusselt number  */
    return Nu;
}





real T_prime(Tracked_Particle* particle)
{
    /*  velocity gradient  */
    real u_grad = v_sqr(particle->cphase.DelV[0])
				+ v_sqr(particle->cphase.DelV[1])
                + v_sqr(particle->cphase.DelV[2]);
                
    if(u_grad > 0.0)
    {
        /*  scalar gradient  */             
        real T_grad = 0.0; /* v_sqr(C_T_G(P_CELL(particle),P_CELL_THREAD(particle),i)); */
            
        /*  scale fluctuation  */
        real T_fluc = gauss_rand_lim(DISTRIB_SIGMA) * sqrt(particle->cphase.tke*T_grad/u_grad);               
        /*  limit fluctuation  */				
        return fmax(T_fluc, -particle->cphase.temp);                
    }	
    else
        return 0.0;
}



real y_prime(Tracked_Particle* particle, const int i)
{
    /*  velocity gradient  */
    real u_grad = v_sqr(particle->cphase.DelV[0])
				+ v_sqr(particle->cphase.DelV[1])
                + v_sqr(particle->cphase.DelV[2]);
                
    if(u_grad > 0.0)
    {
        /*  scalar gradient  */
        real y_grad = 0.0; /* v_sqr(C_YI_G(P_CELL(particle),P_CELL_THREAD(particle),i)); */
    
        /*  scale fluctuation  */
        real y_fluc = gauss_rand_lim(DISTRIB_SIGMA) * sqrt(particle->cphase.tke*y_grad/u_grad);
        /*  limit fluctuation  */
        return fmax(y_fluc, -particle->cphase.yi[i]);               
    }
    else
        return 0.0;
}





/*  VAPORIZATION LAWS  */
DEFINE_DPM_LAW(vaporization, particle, coupled)
{
    if(is_sleeping(particle))
        return;

    switch(particle->type)
    {
    /*  SINGLE-COMPONENT  */
    case DPM_TYPE_DROPLET:
        vaporization_single(particle);
        break;

    /*  MULTI-COMPONENT  */
    case DPM_TYPE_MULTICOMPONENT:
        vaporization_multi(particle);
        break;

    /*  UNSUPPORTED PARTICLE TYPE  */
    default:
        break;
    }
}
