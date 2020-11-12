#include "udf.h"
#include "dpm.h"
#include "math.h"
#include "math_ext.h"
#include "model.h"





void reset_variables_ptr(real* user)
{
    int i;
    for(i=0; i<N_USER_VARS; ++i)
        user[i] = 0.0;
    return;
}

void reset_variables(Tracked_Particle* particle)
{
    reset_variables_ptr(particle->user);
    return;
}

void reset_variables_p(Particle* particle)
{
    reset_variables_ptr(particle->user);
    return;
}


void reset_models_ptr(real* user)
{
    int i;
    for(i=0; i<N_MODEL_VARS; ++i)
        user[i] = 0.0;
    return;
}

void reset_models(Tracked_Particle* particle)
{
    reset_models_ptr(particle->user);
    return;
}

void reset_models_p(Particle* particle)
{
    reset_models_ptr(particle->user);
    return;
}





bool parcel_initialize(Tracked_Particle* particle)
{           
    if((int)particle->user[PARCEL_GENERATION] == 0 && particle->state.time == particle->time_of_birth)
    {
        particle->user[PARCEL_GENERATION] = -1;        				       
        return true;
    }
    return false;
}

bool parcel_initialize_p(Particle* particle)
{   
    if((int)particle->user[PARCEL_GENERATION] == 0 && particle->state.time == particle->time_of_birth)
    {
        particle->user[PARCEL_GENERATION] = -1;
        return true;
    }
    return false;
}





int inc_generation(Tracked_Particle* particle)
{
    return (int)(particle->user[PARCEL_GENERATION] += 1.0);
}

int inc_generation_p(Particle* particle)
{
    return (int)(particle->user[PARCEL_GENERATION] += 1.0);
}





Particle* create_product_parcel(Tracked_Particle* parent, const real cutoff_mass, const real d, const real u_s, const real x_s)
{
    particle_state_t* part_state = &(parent->state);
    particle_state_t* prev_state = &(parent->state0);
    particle_state_t* prod_state;
    
    
    /*  coordinate system    */
    real n[ND_3];
    if(is_sleeping(parent))
    {   /*  align with direction of motion for primary breakup  */  
        v_rel(parent->state.pos,parent->init_state.pos,n);      
    }
    else
    {   /*  align with relative velocity for secondary breakup  */
        v_rel(parent->state.V,parent->cphase.V,n);
    }   
    real E[ND_3][ND_3];                                             
    create_orthonormal(n, E);
    
    /*  scattering  */
    real R[ND_3];
    
    /*  scattering position  */
	gauss_rand_vector(R,DISTRIB_SIGMA);

    real X[ND_3];   
    X[0] = x_s*R[0];
    X[1] = x_s*R[1];
    X[2] = x_s*R[2];
    
    real X_s[ND_3];
    multiply_matrix_col_vector(E, X, X_s);
    
    /*  scattering velocity  */ 	
    scale(R, exp(-v_sqr(R)/2.0), R);  
	
    real U[ND_3];		
    U[0] = 0.0;
    U[1] = u_s*R[1];
    U[2] = u_s*R[2];    
    
    real U_s[ND_3];
    multiply_matrix_col_vector(E, U, U_s);                      
                    

    /*  cutoff mass  */
    if(cutoff_mass <= 0.0)
        return NULL;

    /*  create product parcel  */
    Particle* product = new_particle(parent->injection, true);
    /*  copy parent state  */
    copy_tp_to_p(parent, product);

    /*  assign product state   */
    prod_state = &(product->state);


	/*	append product to particle list  */
	splice_particle_to_list(product, &(parent->injection->p), parent->injection->p_tail);
	/*	track product in next time-step  */
	MARK_PARTICLE(product, P_FL_SPAWNED);

    
    if(cutoff_mass >= part_state->mass * BREAKUP_MAX_MASS_RATIO)
    {
        /*  update product parcel  */
        prod_state->diam = d;
        diam_to_mass(prod_state);
        product->number_in_parcel = parent->number_in_parcel*part_state->mass/prod_state->mass;

        /*  update parent parcel   */
        prev_state->mass -= part_state->mass;
        mass_to_diam(prev_state);

        part_state->mass = 0.0;
        part_state->diam = 0.0; 
    }
    else
    {
        /*  update product parcel  */
        prod_state->diam = d;
        diam_to_mass(prod_state);
        product->number_in_parcel = parent->number_in_parcel*cutoff_mass/prod_state->mass;

        /*  update parent parcel   */
        prev_state->mass -= cutoff_mass;
        mass_to_diam(prev_state);
        
        part_state->mass -= cutoff_mass;
        mass_to_diam(part_state);
    }
    
    /*  update generation  */
    inc_generation_p(product);    

    /*  copy initial state */
    memcpy(P_INIT_STATE(product),
			P_STATE(product), sizeof(particle_state_t));
    memcpy(P_COMPONENT_INIT(product),
            P_COMPONENT(product), sizeof(real)*P_N_COMPONENTS(product));

    /*  update velocities  */
    int i;
    for(i=0; i<ND_3; ++i)
    {   /*  scattering     */
        prod_state->V[i]   += U_s[i];
        prod_state->pos[i] += X_s[i];
    }   
						
    /*  return pointer to product parcel  */
    return product;
}





bool is_new(const Tracked_Particle* particle)
{
    return particle->state.time == particle->init_state.time;
}

bool is_new_p(const Particle* particle)
{
    return particle->state.time == particle->init_state.time;
}


bool is_sleeping(const Tracked_Particle* particle)
{
    return ((int)particle->user[PARCEL_GENERATION] < 0);
}

bool is_sleeping_p(const Particle* particle)
{
    return ((int)particle->user[PARCEL_GENERATION] < 0);
}


bool is_primary(const Tracked_Particle* particle)
{
    return ((int)particle->user[PARCEL_GENERATION] == 0);
}

bool is_primary_p(const Particle* particle)
{
    return ((int)particle->user[PARCEL_GENERATION] == 0);
}


bool is_secondary(const Tracked_Particle* particle)
{
    return ((int)particle->user[PARCEL_GENERATION] > 0);
}

bool is_secondary_p(const Particle* particle)
{
    return ((int)particle->user[PARCEL_GENERATION] > 0);
}





/*  convert mass to diameter and diameter to mass  */
real mass_to_diam(particle_state_t* state)
{
    return state->diam = pow((6.0/M_PI * state->mass/state->rho), 1.0/3.0);
}

real diam_to_mass(particle_state_t* state)
{
    return state->mass = M_PI/6.0 * state->rho*pow3(state->diam);
}
