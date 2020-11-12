#define __COLLISION_C

#include "udf.h"
#include "dpm.h"
#include "dpm_parallel.h"
#include "time.h"
#include "math.h"
#include "math_ext.h"
#include "math_sse.h"
#include "linalg.h"
#include "model.h"
#include "collision.h"
#include "collision_adaptive.h"
#include "collision_search.h"
#include "collision_statistics.h"
#include "secondary_tab.h"
#include "secondary_kh.h"
#include "secondary_khtab.h"
#include "secondary_rdiwakar.h"
#include "vaporization.h"



#define p1      collector
#define p2      target

#define n1      P_N(collector)
#define n2      P_N(target)
#define g1      P_GENERATION(collector)
#define g2      P_GENERATION(target)
#define d1      P_DIAM(collector)
#define d2      P_DIAM(target)
#define D1      P_PDIAM(collector)
#define D2      P_PDIAM(target)
#define R1      P_PDRED(collector)
#define R2      P_PDRED(target)
#define m1      P_MASS(collector)
#define m2      P_MASS(target)
#define s1     TP_COMPONENT(collector)
#define s2      P_COMPONENT(target)
#define T1      P_T(collector)
#define T2      P_T(target)
#define r1      P_RHO(collector)
#define r2      P_RHO(target)
#define u1      P_VEL(collector)
#define u2      P_VEL(target)
#define x1      P_POS(collector)
#define x2      P_POS(target)
#define user1   P_USER(collector)
#define user2   P_USER(target)
#define J       ((void*)P_JACOBIAN(collector))
#define EV      ((void*)P_EIGENVECTORS(collector))
#define ev      ((real*)P_EIGENVALUES(collector))

#define st      DPM_SURFTEN(collector)
#define mu_liq  DPM_MU(collector)
#define rho_liq P_RHO(collector)
#define mu_gas  P_FLOW_VISCOSITY(collector)
#define rho_gas P_FLOW_DENSITY(collector)
#define p_gas   P_FLOW_PRESSURE(collector)




real boundary_We_bouncing(const real d_rel, const real B, const real Oh_red, const real p_red)
{       
    /*  parameters  */
    const real phi = (3.351);    
    const real tau = (1.0-B) * (1.0+d_rel);
    
    real xi;
    if(tau > 1.0)
        xi = -0.25*(1.0+tau)*POW2(tau-2.0) + 1.0;
    else
        xi =  0.25*(3.0-tau)*POW2(tau);

    /*  weber number  */
    const real We = d_rel * (1.0+POW2(d_rel)) * (4.0*phi-12.0) / (xi * (1.0 - POW2(B)));

    return We + fmin(Oh_red*sqrt(p_red)/12.0,20.0);
}



real boundary_We_reflexive(const real d_rel, const real B, const real Oh_liq)
{
    /*  parameters  */
    const real xi = 0.5*B * (1.0+d_rel);
    
    if(xi<d_rel)
    {
        const real eta1 = 2.0*POW2( 1.0 -xi) * sqrt(      1.0  -POW2(xi)) - 1.0;
        const real eta2 = 2.0*POW2(d_rel-xi) * sqrt(POW2(d_rel)-POW2(xi)) - POW3(d_rel);

        /*  weber number  */
        const real We = 3.0 * (7.0*powf(1.0+pow3(d_rel),2.0/3.0) - 4.0*(1.0+pow2(d_rel)))  *  (d_rel * pow2(1.0+pow3(d_rel))) / (pow6(d_rel)*eta1 + eta2);

        if(We > 0.0)
            return We + 680.0*Oh_liq;
        else
            return INFINITY;        
    }
    else        
    {
        return INFINITY;
    }
}



real boundary_B_stretching(const real d_rel, const real We_coll, const real Oh_liq)
{
    if(d_rel > 0.0 && We_coll - 630*Oh_liq > 0.0)
        return sqrt(2.4/(We_coll - 630*Oh_liq) * (POW3(1.0/d_rel) - 2.4*POW2(1.0/d_rel) + 2.7/d_rel));
    else
        return INFINITY;
}





real collision_model(Particle* collector_p, Particle* target, real P_red)
{
#ifdef COLLISION_DISABLE_MODEL		
	statistics_add_event(P_red*fmin(P_N(collector_p),P_N(target)),COLLISION_OUTCOME_BOUNCING,0.0,0.0,0.0,0.0,0.0);
	return 0.0;
#else
    /*  tracked particle struct  */
    Tracked_Particle  collector_tp;
    Tracked_Particle* collector = &collector_tp;

    create_tp_from_p(collector, collector_p);


    /*  collision parameters   */
    real B_coll = sqrt(uniform_rand());
    real z;

    real Re_coll;
    real We_coll;
    real d_coll;
    real u_coll;

    real u_rel[ND_3];
    real x_rel[ND_3];    

    /*  effective collision velocity  */
    v_rel(u2, u1, u_rel);
    v_rel(x2, x1, x_rel);

    if(P_BESSEL(collector) > 0.0)
    {   
        /*  velocity fluctuations  */       
        real u_bar[ND_3];
        multiply_matrix_col_vector(J, x_rel, u_bar);        
        v_rel(u_rel, u_bar, u_rel);     
    
        /*  unbiased fluctuations  */   
        scale(u_rel, 1.0/P_BESSEL(collector), u_rel);       
    }   
    u_coll = v_abs(u_rel);

    if(d1>d2)
    {   /*  collector diameter is larger than target diameter (O'Rourke)  */        
        d_coll  = d2/d1;
        We_coll = rho_liq * d2 * POW2(u_coll) / st;
        Re_coll = sqrt(rho_liq*rho_gas) * d2 * u_coll / mu_liq;
    }
    else
    {   /*  target diameter is larger than collector diameter  */
        d_coll  = d1/d2;
        We_coll = rho_liq * d1 * POW2(u_coll) / st;
        Re_coll = sqrt(rho_liq*rho_gas) * d1 * u_coll / mu_liq;
	}

#ifdef COLLISION_VISCOSITY
    /*	liquid ohnesorge number   */
    real Oh_liq = mu_liq/sqrt(rho_liq*fmin(d1,d2)*st);
    /*	gas ohnesorge number  */
	real Oh_gas = mu_gas/sqrt(rho_gas*fmin(d1,d2)*st);
    /*  reduced pressure  */		
	real p_red  = p_gas*fmin(d1,d2)/st;
#else    
    real Oh_liq = 0.0;
    real Oh_gas = 0.0;
    real p_red  = 0.0;
#endif
	
    /*  collision regime boundaries  */    
    real We_bouncing  = boundary_We_bouncing (d_coll, B_coll,  Oh_liq/Oh_gas, p_red);
    real B_stretching = boundary_B_stretching(d_coll, We_coll, Oh_liq);
    real We_reflexive = boundary_We_reflexive(d_coll, B_coll,  Oh_liq);

    collision_outcome_t collision_outcome;

    		
    if(We_coll < We_bouncing)
    /*  bouncing  */
        collision_outcome = COLLISION_OUTCOME_BOUNCING;
    else
    if(B_coll  > B_stretching)
    /*  stretching  */
        collision_outcome = COLLISION_OUTCOME_STRETCHING;        
    else
    if(We_coll > We_reflexive)
    /*  reflexive  */
        collision_outcome = COLLISION_OUTCOME_REFLEXIVE;    
    else
    /*  coalescence  */        
        collision_outcome = COLLISION_OUTCOME_COALESCENCE;        


    /*  collision breakup coupling  */
    z = khtab_collision_outcome
        (collector, target, P_red, &B_coll, &collision_outcome);
    z = rdiwakar_collision_outcome
        (collector, target, P_red, &B_coll, &collision_outcome);


    real n_coll;
    /*  droplet interaction  */
    switch(collision_outcome)
    {
    case COLLISION_OUTCOME_COALESCENCE:
        /*  coalescence  */
        n_coll = droplet_coalesce(collector, target, P_red, We_coll, B_coll, d_coll);
        break;

    case COLLISION_OUTCOME_BOUNCING:
        /*  bouncing  */
        /*  Message("bouncing\n");  */
        z = -1.0;
        n_coll = droplet_collide(collector, target, P_red, We_coll, B_coll, d_coll, z, collision_outcome);
        break;

    case COLLISION_OUTCOME_STRETCHING:
        /*  stretching  */
        /*  Message("stretching\n");  */
        z = (B_coll - B_stretching)/(1.0 - B_stretching);
        n_coll = droplet_collide(collector, target, P_red, We_coll, B_coll, d_coll, z, collision_outcome);
        break;

    case COLLISION_OUTCOME_REFLEXIVE:
        /*  reflexive  */
        /*  Message("reflexive\n");  */
        z = -sqrt(1.0 - We_reflexive/We_coll);
        n_coll = droplet_collide(collector, target, P_red, We_coll, B_coll, d_coll, z, collision_outcome);
        break;

    case COLLISION_OUTCOME_COUPLED:
        /*  breakup-coupling  */
        /*  Message("breakup-coupling\n");  */
        n_coll = droplet_collide(collector, target, P_red, We_coll, B_coll, d_coll, z, collision_outcome);
        break;

    default:
        /*  no collision  */
        n_coll = 0.0;
        break;
    }
    revert_tp_to_p(collector, collector_p);
    return n_coll;
#endif
}





real droplet_collide(Tracked_Particle* collector, Particle* target, real P_red, const real We_coll, const real B_coll, const real d_coll, const real z, const collision_outcome_t collision_outcome)
{
    /*  parameter z
    //  z =  1  gazing collision        droplets maintain direction
    //  z =  0  inelastic collision     droplets go into the same direction
    //  z = -1  elastic bouncing        droplets maintain relative velocity
    */	
		
	
    /*  breakup-collision coupling  */
    khtab_collision_handler
        (collector, target, NULL, P_red, B_coll, collision_outcome);
    rdiwakar_collision_handler
        (collector, target, NULL, P_red, B_coll, collision_outcome);

        
    real M1 = m1;
    real M2 = m2;    
    /*  multiple collisions  */
    if(n1 > n2)
    {   /*  2 with multiple droplets of 1  */
        M1 *= P_red;                
    }
    else
    {   /*  1 with multiple droplets of 2  */
        M2 *= P_red;
    }
    /*  reduced mass  */
    const real M = fmin(n1,n2) * M1*M2/(M1+M2);

                
    /*  effective collision velocity  */
    real u_rel[ND_3];
    real x_rel[ND_3];    

    v_rel(u2, u1, u_rel);
    v_rel(x2, x1, x_rel);

    if(P_BESSEL(collector) > 0.0)
    {
        /*  velocity fluctuations  */       
        real u_bar[ND_3];
        multiply_matrix_col_vector(J, x_rel, u_bar);
        v_rel(u_rel, u_bar, u_rel);     
    
        /*  unbiased fluctuations  */
        scale(u_rel, 1.0/P_BESSEL(collector), u_rel);       
    }   


    /*  transformed velocities  */    
    real u_aligned[ND_3];    
    
    real E[ND_3][ND_3];
    real cos_psi = 1.0;
    real sin_psi = 0.0;

    
    /*  transform coordinate system  */    
	if(z<0.0)	
    {
        create_orthonormal(u_rel, E);
        /*  rotate into relative velocity direction  */
		multiply_row_vector_matrix(u_rel, E, u_aligned);
		
        /*  rotate coordinate system into collision plane  */
        sin_psi = B_coll;
        cos_psi = sqrt(1.0 - POW2(sin_psi));
                
        rotate_z(u_aligned, cos_psi, sin_psi, u_rel);
    }

	
    const real k = v_sqr(u_rel) * 0.5*M*(1.0-POW2(z));
    /*  momentum transfer  */   
    int i;
    for(i=0; i<ND_3; ++i)
    {
        u_rel[i] *= M*(1.0-(i==0?z:fabs(z)));
    }
    const real p = v_abs(u_rel);


    /*  restore coordinate system  */
    if(z<0.0)
    {
        /*  revert coordinate system from collision plane  */
        rotate_z_inv(u_rel, cos_psi, sin_psi, u_aligned);

        u_aligned[1] /= P_red;
        u_aligned[2] /= P_red;
        
        /*  revert from relative velocity direction  */
        multiply_matrix_col_vector(E, u_aligned, u_rel);
    }


    for(i=0; i<ND_3; ++i)
    {
        u1[i] += u_rel[i]/m1/n1;
        u2[i] -= u_rel[i]/m2/n2;
    }       
	

    /*  report non-permanent collision  */
    statistics_add_event(P_red*fmin(n1,n2),collision_outcome,p,k,We_coll,B_coll,d_coll);
    

#if defined(COLLISION_FRAGMENT_NUMBER) && defined(COLLISION_FRAGMENT_UNIFORMITY)
    if((collision_outcome == COLLISION_OUTCOME_STRETCHING || collision_outcome == COLLISION_OUTCOME_REFLEXIVE) && d_coll > COLLISION_FRAGMENT_UNIFORMITY)    
    {
        const real n_frag = COLLISION_FRAGMENT_NUMBER*P_red*fmin(n1,n2);		
		/*  parcel 1  */
		m1 *= n1/(n1+n_frag/2.0);
		n1 += n_frag/2.0;
		/*  parcel 2  */
		m2 *= n2/(n2+n_frag/2.0);
		n2 += n_frag/2.0;
        /*  adapt diameters   */
        P_MASS_TO_DIAM(p1);
        P_MASS_TO_DIAM(p2);
    }
#endif    
    /*  return number of collisions  */
    return P_red*fmin(n1,n2);
}





real droplet_coalesce(Tracked_Particle* collector, Particle* target, real P_red, const real We_coll, const real B_coll, const real d_coll)
{
	const real P_lim = fmax(n1,n2)/fmin(n1,n2);
	
	if (P_red > P_lim)
		P_red = P_lim;
		
			
    real p;
    real M1 = m1;
    real M2 = m2;   
	/*  multiple collisions    */
    if(n1 >= n2)
    {   /*  coalescing mass    */           
        M1 *= P_red;    
        /*  momentum transfer  */
        p = n2*M1*v_abs(u1);        
    }
    else
    {   /*  coalescing mass    */
        M2 *= P_red;
        /*  momentum transfer  */
        p = n1*M2*v_abs(u2);        
    }
    /*  reduced mass  */
    const real M = fmin(n1,n2) * M1*M2/(M1+M2); 


    /*  effective collision velocity  */
    real u_rel[ND_3];
    real x_rel[ND_3];    
    
    v_rel(u2, u1, u_rel);
    v_rel(x2, x1, x_rel);

    if(P_BESSEL(collector) > 0.0)
    {
        /*  velocity fluctuations  */       
        real u_bar[ND_3];
        multiply_matrix_col_vector(J, x_rel, u_bar);
        v_rel(u_rel, u_bar, u_rel);     
    
        /*  unbiased fluctuations  */
        scale(u_rel, 1.0/P_BESSEL(collector), u_rel);       
    }       
    /*  kinetic energy dissipation   */
    const real k = 0.5*M*v_sqr(u_rel);
    
    
    int i;
    /*  coalesced droplets  */  
    real n_coalesced = fmin(n1,n2);

    /*  coalesced droplet generation */
    real g_coalesced = d1>d2?g1:g2;

    /*  coalesced mass and density   */
    real m_coalesced = M1 + M2;
    real r_coalesced = m_coalesced/(M1/r1 + M2/r2);

    /*  coalesced composition  */
    /*  coalesced temperature  */
    real s_coalesced[VAPOR_MAX_SPECIES];
    real T_coalesced;

    real c1 = 0.0;
    real c2 = 0.0;
    if(P_TYPE(collector) == DPM_TYPE_MULTICOMPONENT)
    {
        /*  coalesced composition  */
        for(i=0; i<TP_N_COMPONENTS(collector); ++i)
            /*  species conservation  */
            s_coalesced[i] = (M1*s1[i] + M2*s2[i]) / (m_coalesced);

        /*  multi-component specific heat  */
        for(i=0; i<TP_N_COMPONENTS(collector); ++i)
        {   /*  mass-weighted mixing law   */
            c1 += s1[i]*MATERIAL_PROP_POLYNOMIAL(collector->injection->material->component[i],PROP_Cp,T1);
            c2 += s2[i]*MATERIAL_PROP_POLYNOMIAL(collector->injection->material->component[i],PROP_Cp,T2);
        }
    }
    else
    {
        /*  single-component specific heat  */
        c1 = MATERIAL_PROP_POLYNOMIAL(collector->injection->material,PROP_Cp,T1);
        c2 = MATERIAL_PROP_POLYNOMIAL(collector->injection->material,PROP_Cp,T2);
    }    
    T_coalesced = (M1*c1*T1 + M2*c2*T2) / (M1*c1+M2*c2);


    /*  coalesced momentum  */
    real u_coalesced[ND_3];
    real x_coalesced[ND_3];     

    for(i=0; i<ND_3; ++i)       
    {   /*  momentum conservation  */
        u_coalesced[i] = (M1*u1[i] + M2*u2[i]) / (m_coalesced);        
        /*  center of mass  */      
        x_coalesced[i] = (M1*x1[i] + M2*x2[i]) / (m_coalesced);     
    }


    /*  update targets  */
    if(n1 > n2)
    {	
		/*  update parcel 1  */
		n1 -= n_coalesced*P_red;
		
        /*  update parcel 2  */
        g2 = g_coalesced;
        m2 = m_coalesced;
        r2 = r_coalesced;
        P_MASS_TO_DIAM(p2);

        if(P_TYPE(collector) == DPM_TYPE_MULTICOMPONENT)
        {
            for(i=0; i<TP_N_COMPONENTS(collector); ++i)
                s2[i] = s_coalesced[i];
        }
        T2 = T_coalesced;
        
        for(i=0; i<ND_3; ++i)
        {
            u2[i] = u_coalesced[i];
#ifdef COLLISION_COALESCE_POSITION
            x2[i] = x_coalesced[i];
#endif
        }

        /*  breakup-collision coupling  */
        khtab_collision_handler
            (collector, target, user2, P_red, 0.0, COLLISION_OUTCOME_COALESCENCE);
        rdiwakar_collision_handler
            (collector, target, user2, P_red, 0.0, COLLISION_OUTCOME_COALESCENCE);        
    }
	else
    { 
		/*  update parcel 2  */
		n2 -= n_coalesced*P_red;

        /*  update parcle 1  */
        g1 = g_coalesced;
        m1 = m_coalesced;
        r1 = r_coalesced;
        P_MASS_TO_DIAM(p1);

        if(P_TYPE(collector) == DPM_TYPE_MULTICOMPONENT)
        {
            for(i=0; i<TP_N_COMPONENTS(collector); ++i)
                s1[i] = s_coalesced[i];
        }
        T1 = T_coalesced;

        for(i=0; i<ND_3; ++i)
        {
            u1[i] = u_coalesced[i];
#ifdef COLLISION_COALESCE_POSITION
            x1[i] = x_coalesced[i];
#endif
        }       
        
        /*  breakup-collision coupling  */
        khtab_collision_handler
            (collector, target, user1, P_red, 0.0, COLLISION_OUTCOME_COALESCENCE);
        rdiwakar_collision_handler
            (collector, target, user1, P_red, 0.0, COLLISION_OUTCOME_COALESCENCE);
    }  

    
    /*  report permanent coalescence  */
    statistics_add_event(P_red*n_coalesced,COLLISION_OUTCOME_COALESCENCE,p,k,We_coll,B_coll,d_coll);
     /*  return number of collisions  */
    return P_red*n_coalesced;
}





bool allow_collision(const Particle* target)
{
    /*  prevent collision  */
    if(is_sleeping_p(target))
        /*  no collision with nozzle flow  */
        return false;
        
    if(n2 > 0.0 && m2 > 0.0 && d2 > 0.0)        
        return true;
    else
        /*  no collision of empty parcels  */
        return false;
}





void parcel_diameter_estimator(Particle** ptrs_minor, Particle** ptrs_major, const int n_minor, const int n_major)
{
#define collector ptrs_minor[h]
#define target ptrs_major[i]   
    int h;
    int i;
    int j;
    int k;    
    
#ifdef __SSE__  
    float_v4sf D_sq;            
    
    float* z = memalign(sizeof(v4sf),sizeof(float)*(n_major+4));
    float* w = memalign(sizeof(v4sf),sizeof(float)*(n_major+4));    
#else
#warning "SSE extensions disabled"
    float* z = malloc(sizeof(float)*(n_major+1));
    float* w = malloc(sizeof(float)*(n_major+1)); 
#endif
        
    
    

	
    for(h=0; h<n_minor; ++h)
    {
        /*  INITIALIZATION  */
        parcel_diameter_initialize_p(collector);        
        if(n_major < 11)
            continue;

        

        
        
        /*  ISOTROPIC PARCEL DIAMETER ESTIMATOR  */        
        
        /*  distances  */
        for(i=0; i<n_major; i++)
        {   /*  distances  */
            z[i] = v_rel_sqr(P_POS(collector),P_POS(ptrs_major[i]));
        }
        
        /*  parcel diameter  */        
        real D = P_PDIAM(collector);
        real F;
             
#ifdef PARCEL_DIAMETER_ITER
        int l;
        for(l=0; l<PARCEL_DIAMETER_ITER; l++)
#endif
        {
            /*  weights  */
#ifdef __SSE__
            D_sq.f[0] = D_sq.f[1] = D_sq.f[2] = D_sq.f[3] = -POW2(D);
            for(i=0; i<n_major; i+=4)
            {
                *(v4sf*)(w+i) = _mm_exp_ps(_mm_div_ps(*(v4sf*)(z+i),D_sq.v));
            }
#else
            for(i=0; i<n_major; ++i)
            {
                w[i] = exp(-z[i]/POW2(D));
            }
#endif
            
            /*  moments  */            
            float w1 = 0.0;
            for(i=0; i<n_major; i++)
            {   
                w1 += w[i];         
            }
            if(w1 <= 0.0)
                break;
                        
            /*  update parcel diameter  */                                            
            D *= F = fminmax(0.5, pow(6.0*M_SQRT_PI/w1,1.0/3.0), 2.0);
#ifdef PARCEL_DIAMETER_TOL        
            if(fabs(F-1.0) < PARCEL_DIAMETER_TOL)
                break;
#endif        
        }
        /*  limit parcel diameter   */      
        P_PDIAM(collector) = D;             

            
            
            

        /*  ANISOTROPIC PARCEL DIAMETER ESTIMATOR  */

        /*  weights  */
#ifdef __SSE__
        D_sq.f[0] = D_sq.f[1] = D_sq.f[2] = D_sq.f[3] = -POW2(D);
        for(i=0; i<n_major; i+=4)
        {
            *(v4sf*)(w+i) = _mm_exp_ps(_mm_div_ps(*(v4sf*)(z+i),D_sq.v));
        }
#else
        for(i=0; i<n_major; ++i)
        {
            w[i] = exp(-z[i]/POW2(D));
        }
#endif
        

        /*  moments  */ 
        float w1 = 0.0;
        float w2 = 0.0;
        for(i=0; i<n_major; i++)
        {   
            w1 += w[i];
            w2 += w[i]*w[i];
        }
        if(w1 <= 0.0)
            continue;

            
        /*  weighted neighbors  */
        const real w_mean
            = P_WMEAN(collector) = w2/w1;                                       
        const real W
            = P_NEIGHBORS(collector) = w1/w_mean;

            
        /*  weighted mean  */   
        real X[ND_3];
        real U[ND_3];
        
        memset(X,0,sizeof(real)*ND_3);    
        memset(U,0,sizeof(real)*ND_3);    

        for(i=0; i<n_major; i++)
        {
            for(j=0; j<ND_3; j++)
            {   
                X[j] += w[i]*P_POS(ptrs_major[i])[j];                  
                U[j] += w[i]*P_VEL(ptrs_major[i])[j];
            }
        }   
        scale(X,1.0/w1,X);          
        scale(U,1.0/w1,U);          

        
        /*  weighted covariance matrices  */
        real covxx[ND_3][ND_3];
        memset(covxx,0,sizeof(real)*ND_3*ND_3);
        real covux[ND_3][ND_3];
        memset(covux,0,sizeof(real)*ND_3*ND_3);
        real covwt[ND_3][ND_3]; 
        memset(covwt,0,sizeof(real)*ND_3*ND_3);

        for(i=0; i<n_major; i++)
        {
            real x[ND_3];
            real u[ND_3];
                            
#ifdef PARCEL_CENTERED_GRADIENT
            v_rel(P_POS(ptrs_major[i]),P_POS(collector),x);
            v_rel(P_VEL(ptrs_major[i]),P_VEL(collector),u);
#else
            v_rel(P_POS(ptrs_major[i]),X,x);
            v_rel(P_VEL(ptrs_major[i]),U,u);
#endif  
            for(j=0; j<ND_3; ++j)
            {
                for(k=0; k<ND_3; ++k)
                {
                    /*  covariance matrix  */
                    covxx[j][k]
                        += w[i] * x[j] * x[k];
                }
                for(k=0; k<ND_3; ++k)
                {
                    /*  correlation matrix */			
                    covux[j][k]
                        += w[i] * u[j] * x[k];
                }
            }

            
#ifdef PARCEL_CENTERED_ELLIPSOID
            v_rel(P_POS(ptrs_major[i]),P_POS(collector),x);
#else
            v_rel(P_POS(ptrs_major[i]),X,x);
#endif  
            for(j=0; j<ND_3; ++j)
            {
                for(k=0; k<ND_3; ++k)
                {
                    /*  covariance matrix  */
                    covwt[j][k]
                        += w[i] * x[j] * x[k];
                }
            }
        }

        for(j=0; j<ND_3; ++j)
        {
            for(k=0; k<ND_3; ++k)
            {
                /*  normalize matrices  */
                covwt[j][k] *= 1.0/(w1-w_mean);                     
            }
        }   


#ifndef COLLISION_DISABLE_GRADIENT  
        /*  jacobian matrix  */                 
        real covinv[ND_3][ND_3];    
            
        /*  regular matrix?  */
        if(det(covxx) != 0.0)
        {   /*  eliminate jacobian  */
            invert_matrix(covxx,covinv);
            multiply_matrices(covux,covinv,J);          
              
            /*  bessel correction   */
            if(W > 4.0)
#ifndef COLLISION_DISABLE_BESSEL    
                P_BESSEL(collector) = sqrt((W-4.0)/(W-1.0));
#else
                P_BESSEL(collector) = 1.0;
#endif      
            else    
                P_BESSEL(collector) = 0.0;
        }
        else
        {   /*  zero velocity gradient  */
            memset(J,0,sizeof(real)*ND_3*ND_3);
            /*  zero bessel correction  */
            P_BESSEL(collector) = 0.0;
        }
#endif  
            
        
#ifndef COLLISION_DISABLE_VOI
        /*  diagonalize covariance matrix  */
        if(eigen(covwt, EV, ev))
        {
            Message("Eigenproblem not converged. det(EV) = %e\n",det(EV));
        }       
        /*  reduced diameter  */
        real D_red = 1.0;

        for(j=0; j<ND_3; ++j)
        {
            /*  half axis  */
            real H = fmax(POW2(d1/D1),2.0*ev[j]/POW2(D1));  
            /*  reduced diameter  */
            D_red *= H;
        }
        P_PDRED(collector) = P_PDIAM(collector) * powf(D_red, 1.0/6.0);
#else
        P_PDRED(collector) = P_PDIAM(collector);    
        /*  isotropic eigenvectors  */
        create_identity(EV);
        /*  isotropic eigenvalues   */
        identity_vector(ev);
        scale(ev,POW2(P_PDIAM(collector))/2.0,ev);  
#endif    
    }
    
    
    /*  free memory  */
    free(z);
    free(w);    

    
    return;
#undef target
#undef collector
}





real parcel_diameter_initialize(Tracked_Particle* collector)
{		
    if(P_PDIAM(collector) <= 0.0)
    {
        /*  initial diameter  */
        P_PDRED(collector) = P_PDIAM(collector) = P_DIAM(collector) * pow(P_N(collector)/PARCEL_DIAMETER_INIT_VFR,1.0/3.0);
        /*  initial eigenvectors  */
        create_identity(EV);
        /*  initial eigenvalues   */
        identity_vector(ev);
        scale(ev,POW2(P_PDIAM(collector))/2.0,ev);
        /*  initial jacobian      */
        memset(J,0,sizeof(real)*ND_3*ND_3);
    }
    return(P_PDIAM(collector));
}

real parcel_diameter_initialize_p(Particle* collector)
{
    if(P_PDIAM(collector) <= 0.0)
    {
        /*  initial diameter  */
        P_PDRED(collector) = P_PDIAM(collector) = P_DIAM(collector) * pow(P_N(collector)/PARCEL_DIAMETER_INIT_VFR,1.0/3.0);
        /*  initial eigenvectors  */
        create_identity(EV);
        /*  initial eigenvalues   */
        identity_vector(ev);
        scale(ev,POW2(P_PDIAM(collector))/2.0,ev);
        /*  initial jacobian      */
        memset(J,0,sizeof(real)*ND_3*ND_3);
    }
    return(P_DIAM(collector));
}





real collision_efficiency(const Particle* collector, const Particle* target)
{
    /*  displacement and relative velocity vectors  */
    real u_rel[ND_3];
    real x_rel[ND_3];    

    v_rel(u2, u1, u_rel);
    v_rel(x2, x1, x_rel);
    
    if(P_BESSEL(collector) > 0.0)
    {
        /*  velocity fluctuations  */       
        real u_bar[ND_3];
        multiply_matrix_col_vector(J, x_rel, u_bar);
        v_rel(u_rel, u_bar, u_rel);     
    
        /*  unbiased fluctuations  */
        scale(u_rel, 1.0/P_BESSEL(collector), u_rel);       
    }   

	
	/*  stokes number  */
	real St;	
	if(d1 > d2)
	{	/*  small droplet stokes number  */
		St = 1.0/18.0 * v_abs(u_rel) * r2/P_FLOW_VISCOSITY(p1)*POW2(d2)/d1;
	}
	else
	{	/*  small droplet stokes number  */
		St = 1.0/18.0 * v_abs(u_rel) * r1/P_FLOW_VISCOSITY(p2)*POW2(d1)/d2;
	}

    
	/*	collision efficiency  */
    real E = pow2((St)/(St+1.0));
    
	return E>=0.0?E:1.0;
}





real collision_probability(const Particle* collector, const Particle* target, const real P_max)
{
    if(n1 <= 0.0 || n2 <= 0.0 || collector == target)
    {   /*  collision probability set zero  */
        return 0.0;
    }
    
    const real d  = 0.5*(d1+d2);
    const real D  = 0.5*(D1+D2);
    const real dt = solver_par.flow_time_step;
            
            
    /*  principal axes  */
    real MA[ND_3][ND_3];
    real f_vol = 1.0;

#ifndef COLLISION_DISABLE_VOI    
    real diag[ND_3];
    
    int i;
    for(i=0; i<ND_3; ++i)
    {
        /*  half axis  */
        real H = fmax(POW2(d/D),2.0*ev[i]/POW2(D1));
        /*  concentration diagonal  */
        f_vol *= diag[i] = 1.0/H;
    }
    multiply_matrix_diag_matrix_t(EV, diag, MA);
    f_vol = sqrt(f_vol);
#else
    create_identity(MA);
    f_vol = 1.0;
#endif
    
    /*  displacement and relative velocity vectors  */
    real u_rel[ND_3];
    real x_rel[ND_3];    

    v_rel(u2, u1, u_rel);
    v_rel(x2, x1, x_rel);
    

    if(P_BESSEL(collector) > 0.0)
    {
        /*  velocity fluctuations  */       
        real u_bar[ND_3];
        multiply_matrix_col_vector(J, x_rel, u_bar);
        v_rel(u_rel, u_bar, u_rel);     
    
        /*  unbiased fluctuations  */
        scale(u_rel, 1.0/P_BESSEL(collector), u_rel);       
    }   
	
    
    /*  displacement and relative velocity scalars  */  
    const real x_sq = distance_norm(x_rel, MA, x_rel);
    const real u_sq = distance_norm(u_rel, MA, u_rel);
    

    /*  axial displacement   */
    const real z0 = -distance_norm(u_rel, MA, x_rel)/sqrt(u_sq);    
    /*  radial displacement  */
    const real r_sq = x_sq - POW2(z0);
    
    
#define COLLISION_DISCRETIZATION 1.0
    /*  axial displacement z1 > z0  */
    real z1 = z0 - sqrt(u_sq)*dt * (COLLISION_DISCRETIZATION-1.0);
#ifdef COLLISION_DETERMINISTIC
    if(z1 < 0.0)
        return 0.0;
#endif
    
    /*  axial displacement z2 < z0  */
    real z2 = z0 - sqrt(u_sq)*dt * (COLLISION_DISCRETIZATION);
#ifdef COLLISION_DETERMINISTIC
    if(z2 < 0.0)
        z2 = 0.0;
#endif


    /*  anisotropy correction  */
    const real fz_sq = u_sq/v_sqr(u_rel);
    const real fr_sq = f_vol/sqrt(fz_sq);

    
    /*  collision probability CSR term  */
    const real CSR = fr_sq * fmax(n1,n2) * POW2(d/D);
    /*  collision probability DMD term  */
    const real DMD = exp(-r_sq/POW2(D));
    
    
#ifdef COLLISION_IGNORE_PROBABILITY
    if(CSR*DMD/P_max < COLLISION_IGNORE_PROBABILITY)
    {   /*  collision probability set zero  */
        return 0.0;
    }
#endif
    

    const real DRD = 0.5*erf(z1/D) - 0.5*erf(z2/D);
    {   /*  collision probability DRD term  */
#ifdef COLLISION_DETERMINISTIC
        return CSR*DMD*DRD * 2.0;
#else
        return CSR*DMD*DRD;
#endif
    }
}





void collision_algorithm(Particle** ptrs_minor, Particle** ptrs_major, const int n_minor, const int n_major)
{
    Particle* collector;
    Particle* target;

    if((n_minor <= 0) || (n_major < n_minor))
        return;
    
    /*  estimate maximum collision probability  */    
    real P_max = 0.0;
#ifdef COLLISION_NTC_ESTIMATOR
    real P_est;
    for(j=0; j<n_major; j++)
    {
        if(P_PDIAM(ptrs_major[j]) > 0.0)
        {            
            /*  estimate collision probability  */
            P_est = P_N(ptrs_major[j])*P_WMEAN(ptrs_major[j]) * POW2(P_DIAM(ptrs_major[j])/P_PDRED(ptrs_major[j]));
            /*  maximize collision probability  */
            P_max = fmax(P_max,P_est);
        }
    }   
    /*  adjust maximum collision probability  */
    P_max = P_max * COLLISION_NTC_ESTIMATOR;    
#else
    P_max = 1.0;
#endif  

    /*  number of tries  */
	int n_tries = floor(0.5 * P_max * (real)n_minor * (real)n_major);
    int i;    
    /*  ntc algorithm    */
    for(i=n_tries; i>0; --i)
    {
        int i_minor = (int)floor(uniform_rand()*(real)n_minor);
        int i_major = (int)floor(uniform_rand()*(real)n_major);

        /*  pick parcel from major cell  */
        collector = ptrs_minor[i_minor];
        target    = ptrs_major[i_major];


        /*  valid pair?  */
        if(!allow_collision(ptrs_minor[i_minor]) || !allow_collision(ptrs_major[i_major]))
            continue;

        
        /*  collision probability  */
#ifdef COLLISION_EFFICIENCY        
        real P = collision_probability(ptrs_minor[i_minor], ptrs_major[i_major], P_max) * collision_efficiency(ptrs_minor[i_minor], ptrs_major[i_major]);
#else
        real P = collision_probability(ptrs_minor[i_minor], ptrs_major[i_major], P_max);
#endif    
        /*  droplets to collide?   */
        if(P <= 0.0)
            continue;
        real P_lim = fmax(n1,n2)/fmin(n1,n2);               
        real P_red = P/P_max;
        
#ifdef COLLISION_NTC_READJUST
        /*  readjust collision frequency  */
        while(P_red > P_lim)
        {
            i *= 2;             /*  increase number of tries        */
            P_red *= 0.5;       /*  decrease collision probability  */
            P_max *= 2.0;       /*  increase maximum frequency      */
            /*  report readjustment */
            statistics_add_readjust();          
        }
#endif

#ifdef COLLISION_NTC_TRUNCATE
    	/*  truncate collision frequency  */
        if (P_red > P_lim)
        {       
            P_red = P_lim;            
            /*  report readjustment */
            statistics_add_readjust();                      
        }
#endif

        
        /*  collision?  */          
        if(P_red > 1.0)
        {   /*  call collision model: mult. collisions  */
            collision_model(ptrs_minor[i_minor], ptrs_major[i_major], P_red);
        }
        else
        if(P_red > uniform_rand())
        {   /*  call collision model: single collision  */
            collision_model(ptrs_minor[i_minor], ptrs_major[i_major], 1.0);
        }           
    }
    /*  (end) loop over candidates  */


    /*  update statistics   */
    statistics_add_tries(n_minor,n_major,n_tries);
    return;
}





void create_tp_from_p(Tracked_Particle* tp, Particle* p)
{   /*  copy state  */
    P_STATE_TO_STATE(p, tp);
	/*  copy part type  */
    P_TYPE(tp)
        = P_TYPE(p);
    /*  copy part ID  */
    P_ID(tp)
        = P_ID(p);
    /*  copy cell ID  */
    P_COPY_CELL_ID(p, tp);
    /*  copy parcel size  */
    P_N(tp)
        = P_N(p);
    /*  link user vars  */
    P_USER(tp)
        = P_USER(p);
    /*  link injection  */
    P_INJECTION(tp)
        = Get_dpm_injections();
    /*  assign time step  */
    P_DT(tp)
        = solver_par.flow_time_step;
        
    /*  link multicomponent state  */
    if(P_TYPE(p) == DPM_TYPE_MULTICOMPONENT)
    {    
        TP_COMPONENT(tp)
            = P_COMPONENT(p);
        TP_N_COMPONENTS(tp)
            = P_N_COMPONENTS(p);        
    }
    return;
}



void revert_tp_to_p(Tracked_Particle* tp, Particle* p)
{   /*  copy state  */
    P_STATE_TO_STATE(tp, p);
    /*  copy parcel  */
    P_N(p)
        = P_N(tp);
}



int linked_list_to_ptrs(Particle* linked_list, Particle** ptrs, const int n_max)
{   /*  item counter  */
    int n=0;

    /*  loop linked list  */
    Particle* item;
    loop(item, linked_list)
    {			
        if(n >= n_max)
        {   /*  prevent list overflow  */
            Message("Warning (linked_list_to_ptrs): particle list overflow!");
            break;
        }
        ptrs[n++] = item;   /*  add item to ptrs  */
    }
    ptrs[n] = NULL;         /*  terminate list  */
    return n;
}





DEFINE_DPM_SPRAY_COLLIDE(collision, collector, target_list)
{
    return;
}





DEFINE_ON_DEMAND(collision_nonadapt_now)
{	
    srand(time(NULL));
    
#if !RP_HOST
    Injection* injection;
    Particle** collector_ptrs;
    Particle** target_ptrs_minor;
    Particle** target_ptrs_major;

    /*  get dpm injection  */
    injection = Get_dpm_injections();

#if PARALLEL
    /*  migrate particles  */
    nodes_to_nodezero_particle_list(&injection->p, injection);
#endif

    /*  start profiling  */
    statistics_reset();

#if PARALLEL
    if(myid == 0)
#endif
    {   /*  collector list  */
        int n_particles = count_p_list(injection->p);
        int n_collectors;
        int n_targets_minor;
        int n_targets_major;

        /*  alloc particle lists  */
        collector_ptrs
            = malloc((n_particles+1)*sizeof(Particle*));
        target_ptrs_minor
            = malloc((n_particles+1)*sizeof(Particle*));
        target_ptrs_major
            = malloc((n_particles+1)*sizeof(Particle*));
        /*  convert lists  */
        n_collectors = n_particles = linked_list_to_ptrs(injection->p, collector_ptrs, n_particles);


        dualcell_build(collector_ptrs, n_particles);
        

        while(n_collectors > 0)
        {   /*  collector index   */
            int i = (int)floor(uniform_rand()*(real)n_collectors--);

            /*  dual-cell search  */
            dualcell_search(collector_ptrs[i], target_ptrs_minor, target_ptrs_major, &n_targets_minor, &n_targets_major, n_particles);

            /*  parcel diameter estimator  */
            parcel_diameter_estimator
                (target_ptrs_minor, target_ptrs_major, n_targets_minor, n_targets_major); 		
            
            /*  collision algorithm  */
            collision_algorithm
                (target_ptrs_minor, target_ptrs_major, n_targets_minor, n_targets_major);

            /*  dual-cell revert  */
            dualcell_revert_periodicity(target_ptrs_major);

            /*  reorder collector list  */
            collector_ptrs[i] = collector_ptrs[n_collectors];
        }


        dualcell_destroy();

        /*  free particle lists  */
        free(target_ptrs_minor);
        free(target_ptrs_major);
        free(collector_ptrs);
    }
#if PARALLEL
    else
    {   /*  free collectors on other nodes  */
        Free_Particles(injection);
    }
#endif

    /*  print statistics  */
    statistics_print();

    /*  relocate particles  */
	Relocate_Unsteady_Particles(true);
#endif
}


DEFINE_EXECUTE_AT_END(collision_nonadapt)
{
#if !RP_HOST
        collision_nonadapt_now();
#endif
}






DEFINE_ON_DEMAND(collision_adaptive_now)
{
    srand(time(NULL));

#if !RP_HOST
    Injection* injection;
    Particle** target_ptrs_minor;
    Particle** target_ptrs_major;
    Particle** ptrs;

    /*  get dpm injection  */
    injection = Get_dpm_injections();

#if PARALLEL
    /*  migrate particles  */
    nodes_to_nodezero_particle_list(&injection->p, injection);
#endif

    /*  start profiling  */
    statistics_reset();

#if PARALLEL
    if(myid == 0)
#endif
    {   /*  collector list  */      
        int n_ptrs = count_p_list(injection->p);        
        int n_targets_minor;
        int n_targets_major;

        /*  alloc particle lists  */
        ptrs = malloc((n_ptrs+1)*sizeof(Particle*));
        /*  convert lists   */
        n_ptrs = linked_list_to_ptrs(injection->p, ptrs, n_ptrs);
                                
        /*  sort particles  */
        int n_cells = adaptive_build(ptrs, n_ptrs);

        int i;
        for(i=0; i<n_cells; ++i)        
        {   /*  search adapt. mesh   */
            adaptive_search(i, &target_ptrs_minor, &target_ptrs_major, &n_targets_minor, &n_targets_major);
            
            /*  parcel diameter estimator  */
            parcel_diameter_estimator
                (target_ptrs_minor, target_ptrs_major, n_targets_minor, n_targets_major); 		
        }
                        
        for(i=0; i<n_cells; ++i)        
        {   /*  search adapt. mesh   */
            adaptive_search(i, &target_ptrs_minor, &target_ptrs_major, &n_targets_minor, &n_targets_major);
            
            /*  collision algorithm  */
            collision_algorithm
                (target_ptrs_minor, target_ptrs_major, n_targets_minor, n_targets_major); 		
        }
        adaptive_destroy();
        
        /*  free particle lists  */
        free(ptrs);
    }
#if PARALLEL
    else
    {   /*  free collectors on other nodes  */
        Free_Particles(injection);
    }
#endif

    /*  print statistics  */
    statistics_print();

    /*  relocate particles  */
	Relocate_Unsteady_Particles(true);
#endif
}


DEFINE_EXECUTE_AT_END(collision_adaptive)
{
#if !RP_HOST
        collision_adaptive_now();
#endif
}
