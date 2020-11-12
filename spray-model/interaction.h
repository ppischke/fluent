#ifndef __INTERACTION_H
#define __INTERACTION_H

#include "udf.h"
#include "dpm.h"





#define KSG_PRODUCTION
#define KSG_DAMPING





#define C_VEL(c,t,i)                (i==0?C_U(c,t): \
                                    (i==1?C_V(c,t): \
                                    (i==2?C_W(c,t):0.0)))
#define C_VEL_M1(c,t,i)             (i==0?C_U_M1(c,t): \
                                    (i==1?C_V_M1(c,t): \
                                    (i==2?C_W_M1(c,t):0.0)))
#define C_VEL_M2(c,t,i)             (i==0?C_U_M2(c,t): \
                                    (i==1?C_V_M2(c,t): \
                                    (i==2?C_W_M2(c,t):0.0)))

#define C_K_M1(c,t)                 C_STORAGE_R(c,t,SV_K_M1)
#define C_D_M1(c,t)                 C_STORAGE_R(c,t,SV_D_M1)

#define C_DPMS_MASS(c,t)            C_SOURCE_MASS_DPM(c,t)
#define C_DPMS_KSG(c,t)             C_UDMI(c,t,0)
#define C_DPMS_KSG_S(c,t)	        C_UDMI(c,t,1)
#define C_DPMS_TAU(c,t)             C_UDMI(c,t,2)
#define C_DPMS_TAU_S(c,t)	        C_UDMI(c,t,3)

#define C_P_G_EFF(c,t,i)	        C_UDMI(c,t,4+i)



#endif
