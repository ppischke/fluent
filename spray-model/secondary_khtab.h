#ifndef __SECONDARY_KHTAB_H
#define __SECONDARY_KHTAB_H

#include "udf.h"
#include "dpm.h"
#include "math_ext.h"

#include "collision.h"
#include "secondary_kh.h"
#include "secondary_tab.h"



void khtab_initialize(Tracked_Particle*);
void khtab_initialize_p(Particle*);

void khtab_update(Tracked_Particle*, kh_constants_t*, tab_constants_t*);
bool khtab_breakup(Tracked_Particle*, kh_constants_t*, tab_constants_t*, bool force);

real khtab_collision_outcome(Tracked_Particle*, Particle*, real const, real*, collision_outcome_t*);
void khtab_collision_handler(Tracked_Particle*, Particle*, real*, real const, const real, const collision_outcome_t);



#endif
