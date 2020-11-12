#include "udf.h"
#include "dpm.h"
#include "math_ext.h"

#include "secondary_khtab.h"

#include "secondary_kh.h"
#include "secondary_tab.h"



void khtab_initialize(Tracked_Particle* droplet)
{
    kh_initialize (droplet);
    tab_initialize(droplet);
    return;
}

void khtab_initialize_p(Particle* droplet)
{
    kh_initialize_p (droplet);
    tab_initialize_p(droplet);
    return;
}



void khtab_update(Tracked_Particle* droplet, kh_constants_t* KH, tab_constants_t* TAB)
{
    kh_update (droplet, KH);
    tab_update(droplet, TAB);
    return;
}



bool khtab_breakup(Tracked_Particle* droplet, kh_constants_t* KH, tab_constants_t* TAB, bool force)
{
    bool breakup = false;

    breakup |= kh_breakup (droplet, KH, force);
    breakup |= tab_breakup(droplet, TAB, KH);
    return breakup;
}



real khtab_collision_outcome(Tracked_Particle* collector, Particle* target, const real P_red, real* impact_parameter, collision_outcome_t* collision_outcome)
{
    return 0.0;
}



void khtab_collision_handler(Tracked_Particle* collector, Particle* target, real* user_coalesced, const real P_red, const real impact_parameter, const collision_outcome_t collision_outcome)
{
    kh_collision_handler
        (collector, target, user_coalesced, P_red, impact_parameter, collision_outcome);
    tab_collision_handler
        (collector, target, user_coalesced, P_red, impact_parameter, collision_outcome);
    return;
}
