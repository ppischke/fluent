#ifndef __VAPORIZATION_H
#define __VAPORIZATION_H

#include "udf.h"
#include "dpm.h"


/*  reference state for diffusion coefficient  */
#define D_P_REF                     100.0e3
#define D_T_REF                     300.0


/*  limiters  */
#define VAPOR_MAX_FRACTION          0.9999
#define VAPOR_MAX_SPECIES           MAX_SPE_EQNS

#define VAPOR_HTC_ITER              5


real correlation_abramzon_sirignano
        (const real Re_d, const real Pr, const real B);
real correlation_ackermann
        (const real Re_d, const real Pr, const real B);

real T_prime(Tracked_Particle* particle);
real y_prime(Tracked_Particle* particle, const int);
        
void vaporization_single(Tracked_Particle*);
void vaporization_multi(Tracked_Particle*);





#endif
