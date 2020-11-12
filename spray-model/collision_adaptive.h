#ifndef __COLLISION_ADAPTIVE_H
#define __COLLISION_ADAPTIVE_H



#include "udf.h"
#include "dpm.h"



/*  
//  target number of parcels according to Pischke, 2014
//  431*27/8 to ensure proper mesh size even at the boundary
//  of the minor control volume
*/
#define COLLISION_ADAPT_LIMIT                   1450
/*  arbitrary values  */
#define COLLISION_ADAPT_MINOR                   100
#define COLLISION_ADAPT_MAJOR                   2700




struct adaptive_cell_struct;



typedef struct adaptive_cell_struct
{
    /*  pointers  */
    struct adaptive_cell_struct* parent;
    struct adaptive_cell_struct* children;
    Particle** ptrs_minor;
    Particle** ptrs_major;
    
    int n_minor;
    int n_major;
    
    /*  geometry  */
    real centroid[ND_3];    
    real resolution;
}
adaptive_cell_t;



#ifdef __COLLISION_ADAPTIVE_C
/*  for internal use only  */
bool sort_particle_to_list(Particle*, adaptive_cell_t*);
bool sort_particles_to_cell(Particle**, const int, adaptive_cell_t*);
void free_children(adaptive_cell_t*);
#endif



int  adaptive_build(Particle** ptrs, const int);
real adaptive_minor_volume(const int);
real adaptive_major_volume(const int);
int  adaptive_search(const int, Particle*** ptrs_minor, Particle*** ptrs_major, int* n_minor, int* n_major);
void adaptive_destroy();



#endif
