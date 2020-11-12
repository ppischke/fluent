#define __COLLISION_ADAPTIVE_C

#include "udf.h"
#include "dpm.h"
#include "math_ext.h"
#include "time.h"

#include "model.h"
#include "collision.h"
#include "collision_adaptive.h"





/*  cell stack  */
adaptive_cell_t*  p_root;
adaptive_cell_t** p_cells;
int n_cells;
int n_sorted;

/*  orientation */
real random_orientation[ND_3][ND_3];
real random_displacement[ND_3];





bool sort_particle_to_list(Particle* particle, adaptive_cell_t* cell)
{
    real position[ND_3];
    /*  mesh orientation   */
    multiply_row_vector_matrix(P_POS(particle), random_orientation, position);
    /*  mesh displacement  */       
    v_rel(position, random_displacement, position);     
    
    int i;
    for(i=0; i<ND_3; ++i)
    {   /*  particle outside cell?   */
        if(fabs(position[i] - cell->centroid[i]) > cell->resolution*1.5)
            return false;
    }
    /*  add particle to major cell  */
    cell->ptrs_major[cell->n_major++] = particle;

    for(i=0; i<ND_3; ++i)
    {   /*  particle outside cell?   */
        if(fabs(position[i] - cell->centroid[i]) > cell->resolution*0.5)
            return true;
    }
    /*  add particle to minor cell  */
    cell->ptrs_minor[cell->n_minor++] = particle;

    return true;
}



bool sort_particles_to_cell(Particle** ptrs, const int n_ptrs, adaptive_cell_t* cell)
{   
    /*  clear children  */
    cell->children = NULL;

    /*  alloc particle lists  */    
    cell->ptrs_minor = malloc((n_ptrs+1)*sizeof(Particle*));
    cell->n_minor = 0;
    
    cell->ptrs_major = malloc((n_ptrs+1)*sizeof(Particle*));    
    cell->n_major = 0;
    
    /*  build particle lists  */
    int i;
    for(i=0; i<n_ptrs; ++i)
    {
#ifdef COLLISION_IGNORE_SLEEPING
        if(!allow_collision(ptrs[i]))                
            continue;
#endif        

        sort_particle_to_list(ptrs[i], cell);
    }
    
    /*  terminate with NULL   */
    cell->ptrs_major[cell->n_major] = NULL;
    cell->ptrs_minor[cell->n_minor] = NULL;
            
    /*  minors?  */
    if(cell->n_minor <= 0)
    {
        /*  free particle lists  */
        free(cell->ptrs_minor);
        free(cell->ptrs_major);

        cell->n_minor = 0;
        cell->n_major = 0;

        return false;
    }   

    /*  refine?  */
    if(cell->n_major > COLLISION_ADAPT_MAJOR && cell->n_minor > COLLISION_ADAPT_MINOR)
    {
        /*  alloc children  */      
        cell->children = malloc(8*sizeof(adaptive_cell_t));     
    
        int j;
        for(j=0; j<8; ++j)
        {
            cell->children[j].parent = cell;
        
            /*  set child centroid  */          
            for(i=0; i<ND_3; ++i)
            {
                if(j&(1<<i))
                    cell->children[j].centroid[i] = cell->centroid[i]+cell->resolution*0.25;
                else
                    cell->children[j].centroid[i] = cell->centroid[i]-cell->resolution*0.25;
            }
            /*  set child length  */
            cell->children[j].resolution = cell->resolution*0.5;                            

            /*  sort particles to child  */                     
            sort_particles_to_cell(cell->ptrs_major, cell->n_major, cell->children+j);                      
        }
        /*  free particle lists  */
        free(cell->ptrs_minor);
        free(cell->ptrs_major);

        cell->n_minor = 0;
        cell->n_major = 0;
        
        return true;
    }   
    /*  coarsen  */
    if(cell->n_major < COLLISION_ADAPT_LIMIT && cell->parent)
    {
        /*  increase number of particles  */
        cell->n_major = cell->parent->n_major;
        
        /*  copy parent particle list     */
        memcpy(cell->ptrs_major, cell->parent->ptrs_major, cell->parent->n_major*sizeof(Particle*));        
    }
                
    n_sorted += cell->n_minor;

    /*  randomize cell list  */
    int i_cell = (int)floor(uniform_rand()*(real)n_cells);  
    
    p_cells[n_cells++]
        = p_cells[i_cell];  
    p_cells[i_cell]
        = cell;
    
    return true;
}



void free_children(adaptive_cell_t* cell)
{
    if(cell->n_minor > 0)
        /*  free minor particle list  */
        free(cell->ptrs_minor);

    if(cell->n_major > 0)
        /*  free major particle list  */
        free(cell->ptrs_major);
        
    if(cell->children)
    {   
        int j;
        for(j=0; j<8; ++j)
        {
            /*  free child cells  */
            free_children(cell->children+j);
        }
        /*  free child pointers   */
        free(cell->children);
    }   
}



int adaptive_build(Particle** ptrs, const int n)
{
    /*  alloc cell list  */
    n_sorted = 0;
    n_cells = 0;    
    p_cells = malloc(sizeof(adaptive_cell_t*)*n);   
    memset(p_cells, 0, sizeof(adaptive_cell_t*)*n);
    
    /*  alloc root cell  */
    p_root = malloc(sizeof(adaptive_cell_t));           
    memset(p_root, 0, sizeof(adaptive_cell_t));
    
                        
    /*  randomize mesh orientation  */
    create_random_base(random_orientation);

    Message("Random mesh orientation:\n");
    Message("  matrix:  %7.4f  %7.4f  %7.4f\n",random_orientation[0][0],random_orientation[1][0],random_orientation[2][0]);
    Message("           %7.4f  %7.4f  %7.4f\n",random_orientation[0][1],random_orientation[1][1],random_orientation[2][1]);
    Message("           %7.4f  %7.4f  %7.4f\n",random_orientation[0][2],random_orientation[1][2],random_orientation[2][2]);

    
    /*  determine mesh resolution  */
    p_root->resolution = 0.0;        
    
    real x_max[ND_3];
    real x_min[ND_3];

    int i;    
    for(i=0; i<ND_3; ++i)
    {
        x_max[i] = -INFINITY;
        x_min[i] =  INFINITY;
    }    
    
    int j;
    for(j=0; j<n; ++j)
    {
        real position[ND_3];
        /*  mesh orientation  */
        multiply_row_vector_matrix(P_POS(ptrs[j]), random_orientation, position);    
        /*  mesh dimensions   */
        for(i=0; i<ND_3; ++i) 
        {        
            x_max[i] = fmax(x_max[i], position[i]);
            x_min[i] = fmin(x_min[i], position[i]);
        }                           
    }        
        
    for(i=0; i<ND_3; ++i)
    {   
        p_root->resolution = fmax(p_root->resolution, x_max[i]-x_min[i])*2.0;
    }
    
    
    /*  randomize mesh origin  */  
    for(i=0; i<ND_3; ++i)
    {   
        random_displacement[i] = 0.5*(x_max[i]+x_min[i]) + 0.5*(uniform_rand()-0.5)*(p_root->resolution);
    }
    
    
    /*  sort particles  */
    sort_particles_to_cell(ptrs, n, p_root);
    
    
    /*  return  */
    Message("Search mesh built\n");
    Message("  cells:      %d\n",n_cells);
    Message("  particles:  %d\n",n);
    Message("  sorted:     %d\n",n_sorted);
    
    return n_cells;
}



int adaptive_search(const int cell_id, Particle*** ptrs_minor, Particle*** ptrs_major, int* n_minor, int* n_major)
{   /*  cell index out of bounds?  */
    if(cell_id >= n_cells)
    {
        *ptrs_minor = NULL;
        *ptrs_major = NULL;

        *n_minor = 0;
        *n_major = 0;
    
        return n_cells;
    }       
    else
    {   
        *ptrs_minor = p_cells[cell_id]->ptrs_minor;
        *ptrs_major = p_cells[cell_id]->ptrs_major;

        *n_minor = p_cells[cell_id]->n_minor;
        *n_major = p_cells[cell_id]->n_major;
    
        return n_cells;
    }
}



real adaptive_minor_volume(const int cell_id)
{   /*  cell index out of bounds?  */
    if(cell_id >= n_cells)
        return 0;       

    /*  volume of minor cell  */
    return POW3(p_cells[cell_id]->resolution);
}

real adaptive_major_volume(const int cell_id)
{   /*  cell index out of bounds?  */
    if(cell_id >= n_cells)
        return 0;       

    /*  volume of minor cell  */
    return POW3(p_cells[cell_id]->resolution)*27.0;
}



void adaptive_destroy()
{
    /*  free children  */   
    free_children(p_root);  
    
    /*  free cells  */
    free(p_cells);
    free(p_root);
    return;
}
