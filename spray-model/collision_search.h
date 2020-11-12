#ifndef __COLLISION_SEARCH_H
#define __COLLISION_SEARCH_H



#include "udf.h"
#include "dpm.h"



#define COLLISION_MESH_RESOLUTION            1.0e-3
#define COLLISION_MESH_CELLS                 250



struct dualcell_list_struct;
struct dualcell_cell_struct;



typedef struct dualcell_item_struct
{
    Particle* particle;
    struct dualcell_item_struct* next[3];
    int search_id;
    int flag;
}
dualcell_item_t;



typedef struct dualcell_cell_struct
{
    /*  pointers     */
    struct dualcell_item_struct* list[3];
    int search_id;
}
dualcell_cell_t;



#ifdef __COLLISION_SEARCH_C
/*  for internal use only  */
void alloc_dualcell_mesh(const int n);
void free_dualcell_mesh();
void sort_particle(Particle*);
int  cell_id_from_index(int *);
int  cell_id_from_position(const real*, int*);
int  search_cell(const int, Particle**, const int, int, const int);
#endif



void dualcell_build(Particle** particles, const int n_particles);
real dualcell_minor_volume();
real dualcell_major_volume();
int  dualcell_search(const Particle*, Particle** ptrs_minor, Particle** ptrs_major, int* n_minor, int* n_major, const int n_max);
void dualcell_revert_periodicity(Particle** ptrs);
void dualcell_reset_search();
void dualcell_destroy();



#endif
