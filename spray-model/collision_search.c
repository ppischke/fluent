#define __COLLISION_SEARCH_C

#include "udf.h"
#include "dpm.h"
#include "math_ext.h"
#include "time.h"

#include "model.h"
#include "collision.h"
#include "collision_search.h"





/*  node stack  */
dualcell_cell_t*  p_cells;
dualcell_item_t*  p_items;
dualcell_item_t*  p_free;



/*  mesh size  */
real  mesh_orientation[ND_3][ND_3];
real  mesh_displacement[ND_3];
real  mesh_resolution;
int   mesh_divisions;

real(*periodic_matrix)[ND_3];  /*  rotation matrix  */
real  periodic_angle;          /*  rotation angle   */

/*  search id  */
int search_id;
int n_items;
int n_cells;



void alloc_dualcell_mesh(const int n)
{
    /*  number of elements and number of divisions   */
    n_cells = POW3(mesh_divisions);
    n_items = n;

    /*  allocate items and cells  */
    p_cells = malloc(sizeof(dualcell_cell_t)*n_cells);
    p_items = malloc(sizeof(dualcell_item_t)*n);
    p_free = p_items;

    /*  clear search mesh  */
    memset(p_cells, sizeof(dualcell_cell_t)*n_cells, 0);
    /*  clear linked list  */
    memset(p_items, sizeof(dualcell_item_t)*n, 0);


    /*  randomize mesh orientation  */
    int i;
    for(i=0; i<ND_3; ++i)
    {
        mesh_displacement[i] = uniform_rand()*mesh_resolution;
    }
    create_random_base(mesh_orientation);

    Message("Random mesh orientation:\n");
    Message("  matrix:  %7.4f  %7.4f  %7.4f\n",mesh_orientation[0][0],mesh_orientation[1][0],mesh_orientation[2][0]);
    Message("           %7.4f  %7.4f  %7.4f\n",mesh_orientation[0][1],mesh_orientation[1][1],mesh_orientation[2][1]);
    Message("           %7.4f  %7.4f  %7.4f\n",mesh_orientation[0][2],mesh_orientation[1][2],mesh_orientation[2][2]);
    Message("  vector:  %7.4f  %7.4f  %7.4f\n",mesh_displacement[0]/1.0e-3,mesh_displacement[1]/1.0e-3,mesh_displacement[2]/1.0e-3);

    return;
}



void free_dualcell_mesh()
{   /*  free stacks  */
    free(p_cells);
    free(p_items);
    return;
}




int cell_id_from_position(const real position[ND_3], int* dim_id)
{
    real transformed[ND_3];

    /*  mesh orientation   */
    multiply_row_vector_matrix(position, mesh_orientation, transformed);
    /*  mesh displacement  */
    v_sum(transformed, mesh_displacement, transformed);

    int i;
    for(i=0; i<ND_3; ++i)
    {
        dim_id[i] = (int)round(transformed[i]/mesh_resolution) + mesh_divisions/2;
    }
    /*  return cell id  */
    return cell_id_from_index(dim_id);
}



int cell_id_from_index(int* dim_id)
{
    int cell_id=0;
    int m=1;

    int i;
    for(i=0; i<ND_3; ++i)
    {
        /*  inside mesh?  */
        if(dim_id[i] < 0)
            return -1;
            /*
            dim_id[i] = 0;
            */
        if(dim_id[i] >= mesh_divisions)
            return -1;
            /*
            dim_id[i] = mesh_divisions-1;
            */

        /*  inside mesh? */
        cell_id += dim_id[i]*m;

        /*  dim id multiplicator  */
        m *= mesh_divisions;
    }
    return cell_id;
}



void sort_particle_to_cell(Particle* particle)
{
#ifdef COLLISION_IGNORE_SLEEPING
    if(!allow_collision(particle))                
        return;
#endif        
        
    int index[ND_3];
    real position[ND_3];

    collision_shadow_t shadow;
    collision_shadow_t shadows;
    /*  periodic shadows  */
    if(periodic_matrix)
        shadows = 3;
    else
        shadows = 1;


    p_free->particle = particle;
    p_free->flag = 0;
    p_free->search_id = 0;


    /*  reset shadow flag  */
    P_COLLISION_PERIODIC(particle) = ORIGINAL_POSITION;

    for(shadow=0; shadow<shadows; ++shadow)
    {
        /*  shadow positions  */
        switch(shadow)
        {
        case PERIODIC_SHADOW_HIGH:
            /*  transform shadow targets  */
            multiply_row_vector_matrix
                (P_POS(particle), periodic_matrix, position);
            break;

        case PERIODIC_SHADOW_LOW:
            /*  transform shadow targets  */
            multiply_matrix_col_vector
                (periodic_matrix, P_POS(particle), position);
            break;

        case ORIGINAL_POSITION:
            /*  untransformed target  */
            memcpy(position, P_POS(particle), sizeof(real)*ND_3);
        default:
            break;
        }


        /*  cell id  */
        int cell_id = cell_id_from_position(position,index);
        if (cell_id < 0)
            return;

        /*  prepend particle  */
        p_free->next[shadow] = p_cells[cell_id].list[shadow];

        /*  assign to cell  */
        p_cells[cell_id].list[shadow] = p_free;
    }
    /*  return  */
    p_free++;
    return;
}



int search_cell(const int cell_id, Particle** ptrs, const int n_max, int n, const int flag)
{
    /*  cell id out of bounds  */
    if(cell_id < 0)
        return n;

    collision_shadow_t shadow;
    collision_shadow_t shadows;
    /*  periodic shadows  */
    if(periodic_matrix)
        shadows = 3;
    else
        shadows = 1;


    for(shadow=0; shadow<shadows; ++shadow)
    {
        dualcell_item_t* item;

        for(item = p_cells[cell_id].list[shadow]; item != NULL; item = item->next[shadow])
        {
            /*  mark with flag  */
            if(flag)
            {
                if(item->flag == flag)
                    /*  item flagged  */
                    continue;
                else
                    /*  flag item     */
                    item->flag = flag;
            }


            /*  mark with search id  */
            if(item->search_id == search_id)
                /*  already in list  */
                continue;
            else
                /*  mark as found    */
                item->search_id = search_id;


            if(n<n_max)
            {
                /*  add particle to ptr list  */
                ptrs[n++] = item->particle;

                /*  shadow velocity  */
                real position[ND_3];
                real velocity[ND_3];

                /*  shadow  */
                switch(shadow)
                {
                case PERIODIC_SHADOW_HIGH:
                    /*  transform shadow targets  */
                    multiply_row_vector_matrix
                        (P_VEL(item->particle), periodic_matrix, velocity);
                    memcpy
                        (P_VEL(item->particle), velocity, sizeof(real)*ND_3);
                    multiply_row_vector_matrix
                        (P_POS(item->particle), periodic_matrix, position);
                    memcpy
                        (P_POS(item->particle), position, sizeof(real)*ND_3);
                    break;

                case PERIODIC_SHADOW_LOW:
                    /*  transform shadow targets  */
                    multiply_matrix_col_vector
                        (periodic_matrix, P_VEL(item->particle), velocity);
                    memcpy
                        (P_VEL(item->particle), velocity, sizeof(real)*ND_3);
                    multiply_matrix_col_vector
                        (periodic_matrix, P_POS(item->particle), position);
                    memcpy
                        (P_POS(item->particle), position, sizeof(real)*ND_3);
                    break;

                case ORIGINAL_POSITION:
                default:
                    break;
                }
                /*  set shadow id    */
                P_COLLISION_PERIODIC(item->particle) = shadow;
            }
        }
    }
    return n;
}





void dualcell_build(Particle** particles, const int n_particles)
{
    Domain* domain = Get_Domain(1);
    Thread* thread;

    search_id = 0;

    /*  periodic boundary conditions  */
    periodic_matrix = NULL;

    thread_loop_f(thread, domain)
    {   /*  find periodic thread  */
        if(PERIODIC_FACE_THREAD_P(thread))
        {   /*  assign rotation matrix  */
            periodic_angle  = THREAD_VAR(thread).per.angle;
            periodic_matrix = malloc(sizeof(real)*ND_3*ND_3);

            /*  randomize periodic direction  */
            if(uniform_rand()<0.5)
                /*  copy as is  */
                transpose_matrix(THREAD_VAR(thread).per.rm, periodic_matrix);
            else
                /*  transpose   */
                memcpy(periodic_matrix,THREAD_VAR(thread).per.rm,sizeof(real)*ND_3*ND_3);

            Message("Periodic boundary condition found\n");
            Message("  angle:  %8.4f\n", periodic_angle * 180.0/M_PI);
            Message("  matrix:  %7.4f  %7.4f  %7.4f\n", periodic_matrix[0][0], periodic_matrix[0][1], periodic_matrix[0][2]);
            Message("           %7.4f  %7.4f  %7.4f\n", periodic_matrix[1][0], periodic_matrix[1][1], periodic_matrix[1][2]);
            Message("           %7.4f  %7.4f  %7.4f\n", periodic_matrix[2][0], periodic_matrix[2][1], periodic_matrix[2][2]);
            break;
        }
    }

    /*  configuration  */
    mesh_divisions  = COLLISION_MESH_CELLS;
    mesh_resolution = COLLISION_MESH_RESOLUTION;

    /*  alloc search tree  */
    alloc_dualcell_mesh(n_particles);

    /*  loop over particles  */
    int i;
    for(i=0; i<n_particles; ++i)
    {   /*  sort to tree     */
        if(particles[i]->number_in_parcel <= 0.0)
            continue;   

        sort_particle_to_cell(particles[i]);
    }
    /*  return  */
    Message("Search mesh built\n");
    Message("  particles:  %d\n",n_particles);
    Message("  sorted:     %d\n",p_free - p_items);
    return;
}



int dualcell_search(const Particle* centroid, Particle** ptrs_minor, Particle** ptrs_major, int* n_minor, int* n_major, const int n_max)
{
    int minor_index[ND_3];
    int major_index[ND_3];
    int offset[ND_3];
    int cell_id;


    search_id++;
    /*  minor cell id  */
    cell_id = cell_id_from_position(P_POS(centroid),minor_index);

    /*  minor cell searched?  */
    if(cell_id < 0 || p_cells[cell_id].search_id != 0)
    {
        /*  already searched  */
        /*  terminate list    */
        ptrs_minor[0] = NULL;
        ptrs_major[0] = NULL;
        return *n_major = *n_minor = 0;
    }
    else
        /*  mark as searched  */
        p_cells[cell_id].search_id = search_id;

    /*  search minor cell  */
    *n_major = *n_minor = search_cell(cell_id, ptrs_minor, n_max, 0, true);

    /*  minor cell empty?  */
    if(*n_minor<=0)
    {
        /*  empty cell      */
        /*  terminate list  */
        ptrs_minor[0] = NULL;
        ptrs_major[0] = NULL;
        return *n_major = *n_minor = 0;
    }


    /*  major cells  */
    memcpy(ptrs_major, ptrs_minor, sizeof(Particle*)*(*n_minor));

    for(offset[0]=0; offset[0]<3; ++offset[0])
    for(offset[1]=0; offset[1]<3; ++offset[1])
    for(offset[2]=0; offset[2]<3; ++offset[2])
    {
        int i;
        for(i=0; i<ND_3; ++i)
        {
            if(offset[i] >= 2)
                major_index[i] = minor_index[i]-1;
            else
                major_index[i] = minor_index[i]+offset[i];
        }

        /*  major cell id  */
        cell_id = cell_id_from_index(major_index);
        /*  search major cell  */
        *n_major = search_cell(cell_id, ptrs_major, n_max, *n_major, false);
    }


    /*  terminate list  */
    ptrs_minor[*n_minor] = NULL;
    ptrs_major[*n_major] = NULL;
    return *n_major;
}



void dualcell_reset_search(const int n)
{
    int i;
    for(i=0; i<n_cells; ++i)
    {
        p_cells[i].search_id = 0;       
    }
    
    for(i=0; i<n_items; ++i)
    {
        p_items[i].search_id = 0;
        p_items[i].flag = false;
    }
}



void dualcell_revert_periodicity(Particle** ptrs)
{
    /*  periodic shadows  */
    if(!periodic_matrix)
        return;


    while((*ptrs) != NULL)
    {
        real v[ND_3];

        /*  revert shadow targets in periodic systems  */
        switch(P_COLLISION_PERIODIC(*ptrs))
        {
        case PERIODIC_SHADOW_HIGH:
            /*  position  */
            multiply_matrix_col_vector
                (periodic_matrix, P_POS(*ptrs), v);
            memcpy
                (P_POS(*ptrs), v, sizeof(real)*ND_3);
            /*  velocity  */
            multiply_matrix_col_vector
                (periodic_matrix, P_VEL(*ptrs), v);
            memcpy
                (P_VEL(*ptrs), v, sizeof(real)*ND_3);
            break;

        case PERIODIC_SHADOW_LOW:
            /*  position  */
            multiply_row_vector_matrix
                (P_POS(*ptrs), periodic_matrix, v);
            memcpy
                (P_POS(*ptrs), v, sizeof(real)*ND_3);
            /*  velocity  */
            multiply_row_vector_matrix
                (P_VEL(*ptrs), periodic_matrix, v);
            memcpy
                (P_VEL(*ptrs), v, sizeof(real)*ND_3);
            break;

        case ORIGINAL_POSITION:
        default:
            break;
        }
        P_COLLISION_PERIODIC(*(ptrs++)) = ORIGINAL_POSITION;
    }
}



real dualcell_minor_volume()
{   /*  volume of minor cell  */
    return POW3(mesh_resolution);
}

real dualcell_major_volume()
{   /*  volume of major cell  */
    return POW3(mesh_resolution)*27.0;
}



void dualcell_destroy()
{
    /*  free per. matrix  */
    if(periodic_matrix)
        free(periodic_matrix);

    /*  free search tree  */
    free_dualcell_mesh();
    return;
}
