#ifndef __PRIMARY_LISA_H
#define __PRIMARY_LISA_H

#include "udf.h"
#include "dpm.h"
#include "model.h"



#define LISA_DIST_PARAMETER         DIST_PARAMETER
/*
#define LISA_FIXED_DISPERS_ANGLE      8.7e-3
*/
#define LISA_FIXED_BREAKUP_LENGTH   350.0e-6



void lisa_initialize(Tracked_Particle*);



#endif
