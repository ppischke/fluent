#define __COLLISION_STATISTICS_C

#include "udf.h"
#include "dpm.h"
#include "dpm_parallel.h"
#include "math.h"
#include "math_ext.h"
#include "time.h"
#include "model.h"

#include "collision.h"
#include "collision_statistics.h"



collision_statistics_t collision_statistics;
/*  extended collision log  */
#ifdef COLLISION_EXT    
    FILE* ext;
#endif



void statistics_reset()
{
    memset(&collision_statistics, 0, sizeof(collision_statistics_t));

    collision_statistics.execution_time = clock();
    
    /*  extended collision log  */
#ifdef COLLISION_EXT
#if PARALLEL
    if(myid == 0)
#endif
    {   /*  name extended log  */
        char ext_log[80];
        sprintf(ext_log, COLLISION_EXT, solver_par.flow_time);
        
        /*  open extended log  */
        ext = fopen(ext_log,"a");
    }
#endif  
    return;
}



void statistics_add_tries(const int n_parcels, const int n_targets, const int n_tries)
{
    collision_statistics.n_parcels += n_parcels;
    collision_statistics.n_pairs   += n_parcels*n_targets;
    collision_statistics.n_tries   += n_tries;
    return;
}



void statistics_add_event(const real collisions, const collision_outcome_t collision_outcome, const real momentum, const real dissipation, const real We_coll, const real B_coll, const real d_rel)
{
    collision_statistics.momentum += momentum;
    collision_statistics.m_regimes[collision_outcome] += momentum;
    
    collision_statistics.dissipation += dissipation;    
    collision_statistics.d_regimes[collision_outcome] += dissipation;
    
    collision_statistics.collisions += collisions;
    collision_statistics.c_regimes[collision_outcome] += collisions;
    
#ifdef COLLISION_EXT
    fprintf(ext, "%6.4e\t%.1d\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\n",collisions,collision_outcome,momentum,dissipation,We_coll,B_coll,d_rel);
#endif
    return;
}



void statistics_add_readjust()
{
    collision_statistics.n_readjust++;
    return;
}



void statistics_print()
{
#if RP_NODE
    real storage[COLLISION_REGIMES];  /*  storage for PRF_GRSUM  */    
#endif
    
    collision_statistics.execution_time = clock()-collision_statistics.execution_time;

    real exectime     = (double)collision_statistics.execution_time*1.0e-3;

    int  n_parcels    = PRF_GISUM1(collision_statistics.n_parcels);
    int  n_pairs      = PRF_GISUM1(collision_statistics.n_pairs);
    int  n_tries      = PRF_GISUM1(collision_statistics.n_tries);
    int  n_readjust   = PRF_GISUM1(collision_statistics.n_readjust);    

    real collisions   = PRF_GRSUM1(collision_statistics.collisions);
    real momentum     = PRF_GRSUM1(collision_statistics.momentum);
    real dissipation  = PRF_GRSUM1(collision_statistics.dissipation);
    
    real c_regimes[COLLISION_REGIMES];
    memcpy(c_regimes, collision_statistics.c_regimes, sizeof(real)*COLLISION_REGIMES);
    PRF_GRSUM(c_regimes, COLLISION_REGIMES, storage);
    
    real m_regimes[COLLISION_REGIMES];
    memcpy(m_regimes, collision_statistics.m_regimes, sizeof(real)*COLLISION_REGIMES);  
    PRF_GRSUM(m_regimes, COLLISION_REGIMES, storage);
    
    real d_regimes[COLLISION_REGIMES];
    memcpy(d_regimes, collision_statistics.d_regimes, sizeof(real)*COLLISION_REGIMES);  
    PRF_GRSUM(d_regimes, COLLISION_REGIMES, storage);
    
    if(n_pairs<=0)
        return;

#if PARALLEL
    if(myid == 0)
#endif
    {
        /*  log to console  */
        Message("Collision statistics:\n");
        Message("  execution time (ms):   %6.4e\n",exectime);
        Message("  number of parcels:     %6d\n",n_parcels);
        Message("  number of pairs:       %6d\n",n_pairs);
        Message("  number of tries:       %6d\n",n_tries);
        Message("  number of adjustments: %6d\n",n_readjust);
        
        Message("  number of collisions:  %6.4e\n",collisions);
        if(collisions > 0.0)
        {
            Message("    coalescence:           %4.2f\n",c_regimes[COLLISION_OUTCOME_COALESCENCE]/collisions);
            Message("    stretching sep.:       %4.2f\n",c_regimes[COLLISION_OUTCOME_STRETCHING] /collisions);
            Message("    reflexive sep.:        %4.2f\n",c_regimes[COLLISION_OUTCOME_REFLEXIVE]  /collisions);
            Message("    bouncing:              %4.2f\n",c_regimes[COLLISION_OUTCOME_BOUNCING]   /collisions);
            Message("    breakup coupled:       %4.2f\n",c_regimes[COLLISION_OUTCOME_COUPLED]    /collisions);
        }
        Message("  momentum exchanged:    %6.4e\n",momentum);
        if(momentum > 0.0)
        {
            Message("    coalescence:           %4.2f\n",m_regimes[COLLISION_OUTCOME_COALESCENCE]/momentum);
            Message("    stretching sep.:       %4.2f\n",m_regimes[COLLISION_OUTCOME_STRETCHING] /momentum);
            Message("    reflexive sep.:        %4.2f\n",m_regimes[COLLISION_OUTCOME_REFLEXIVE]  /momentum);
            Message("    bouncing:              %4.2f\n",m_regimes[COLLISION_OUTCOME_BOUNCING]   /momentum);
            Message("    breakup coupled:       %4.2f\n",m_regimes[COLLISION_OUTCOME_COUPLED]    /momentum);
        }
        Message("  energy dissipated:     %6.4e\n",dissipation);
        if(dissipation > 0.0)
        {
            Message("    coalescence:           %4.2f\n",d_regimes[COLLISION_OUTCOME_COALESCENCE]/dissipation);
            Message("    stretching sep.:       %4.2f\n",d_regimes[COLLISION_OUTCOME_STRETCHING] /dissipation);
            Message("    reflexive sep.:        %4.2f\n",d_regimes[COLLISION_OUTCOME_REFLEXIVE]  /dissipation);
            Message("    bouncing:              %4.2f\n",d_regimes[COLLISION_OUTCOME_BOUNCING]   /dissipation);
            Message("    breakup coupled:       %4.2f\n",d_regimes[COLLISION_OUTCOME_COUPLED]    /dissipation);
        }
        
        /*  log to file  */
        FILE* log = fopen(COLLISION_LOG,"a");

        if(!ftell(log))
            fprintf(log, "time\tparcels\tpairs\ttries\tadj\tcoll\tmomentum\tdissipation\tclsc\tstr\trflx\tbnce\tcoupled\tclsc\tstr\trflx\tbnce\tcoupled\tclsc\tstr\trflx\tbnce\tcoupled\texec. time\n");

        fprintf(log, "%6.4f\t%d\t%d\t%d\t%d\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\n",
            solver_par.flow_time*1.0e3,
            n_parcels,            
            n_pairs,
            n_tries,
            n_readjust,
            collisions,
            momentum,
            dissipation,
            c_regimes[COLLISION_OUTCOME_COALESCENCE],
            c_regimes[COLLISION_OUTCOME_STRETCHING],
            c_regimes[COLLISION_OUTCOME_REFLEXIVE],
            c_regimes[COLLISION_OUTCOME_BOUNCING],
            c_regimes[COLLISION_OUTCOME_COUPLED],
            m_regimes[COLLISION_OUTCOME_COALESCENCE],
            m_regimes[COLLISION_OUTCOME_STRETCHING],
            m_regimes[COLLISION_OUTCOME_REFLEXIVE],
            m_regimes[COLLISION_OUTCOME_BOUNCING],
            m_regimes[COLLISION_OUTCOME_COUPLED],
            d_regimes[COLLISION_OUTCOME_COALESCENCE],
            d_regimes[COLLISION_OUTCOME_STRETCHING],
            d_regimes[COLLISION_OUTCOME_REFLEXIVE],
            d_regimes[COLLISION_OUTCOME_BOUNCING],
            d_regimes[COLLISION_OUTCOME_COUPLED],           
            exectime);

        fclose(log);
    }
    
    /*  extended collision log  */  
#ifdef COLLISION_EXT
#if PARALLEL
    if(myid == 0)
#endif
    {   /*  close extended log  */
        fclose(ext);
    }
#endif  
    return;
}
