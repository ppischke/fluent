#ifndef __PRIMARY_TURB_H
#define __PRIMARY_TURB_H

#include "udf.h"
#include "dpm.h"
#include "model.h"
#include "collision.h"





#define TB_CONE_ANGLE       0.09
#define TB_SCTR_ANGLE       0.09

#define TB_C2E              1.92  /*  standard k-epsilon model   */
#define TB_CMU              0.09  /*  standard k-epsilon model   */
#define TB_CL               0.04  /*  fully developed pipe flow  */
#define TB_CLB              2.0   /*  Huh-Gosman                 */
#define TB_CTB              1.0   /*  Huh-Gosman                 */

#define TB_DIST_TYPE        DIST_ROSIN_R
#define TB_DIST_PARAMETER   DIST_PARAMETER

#define TB_CUTOFF           0.10





void tb_initialize(Tracked_Particle*);
void tb_initialize_p(Particle*);

void tb_update (Tracked_Particle*);
bool tb_breakup(Tracked_Particle*, bool);

/*
real tb_collision_outcome(Tracked_Particle*, Particle*, real*, collision_outcome_t*);
void tb_collision_handler(Tracked_Particle*, Particle*, real*, const real, const collision_outcome_t);
*/



#endif
