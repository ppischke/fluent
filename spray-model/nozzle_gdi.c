#include "udf.h"
#include "dpm.h"
#include "math.h"
#include "math_ext.h"



#define DOMAIN_ANGLE    M_PI/4.0

#define TIME_LIFT_OPEN1     0e-6
#define TIME_LIFT_OPEN2    80e-6
#define TIME_LIFT_CLOSE1 -130e-6
#define TIME_LIFT_CLOSE2   -0e-6

#define TIME_VEL_OPEN1      0e-6
#define TIME_VEL_OPEN2    170e-6
#define TIME_VEL_OPEN3    300e-6
#define TIME_VEL_CLOSE1   -30e-6
#define TIME_VEL_CLOSE2    -0e-6





real nozzle_gdi_scale_lift
    (const Particle *particle, const Injection *injection, const real s);
real nozzle_gdi_scale_velocity
    (const Particle *particle, const Injection *injection, const real s);



real nozzle_gdi_scale_lift(const Particle *particle, const Injection *injection, const real s)
{
    /*  particle time  */
    const real time
        = particle->state.time - injection->unsteady_start + s;

    /*  injection duration  */
    const real duration
        = injection->unsteady_stop - injection->unsteady_start;

    /*  injection timing  */
    const real time_open1  = TIME_LIFT_OPEN1;
    const real time_open2  = TIME_LIFT_OPEN2;
    const real time_close1 = TIME_LIFT_CLOSE1 + duration;
    const real time_close2 = TIME_LIFT_CLOSE2 + duration;

    /*  SCALE NEEDLE LIFT  */
    real f_lift;

    if(time < time_open1)
    {   /*  delay  */
        f_lift = 0.40;
    }
    else
    if(time < time_open2)
    {   /*  lift needle  */
        f_lift = 0.40 + 0.60*(time - time_open1)/(time_open2 - time_open1);
    }
    else
    if(time < time_close1)
    {   /*  open  */
        f_lift = 1.00;
    }
    else
    if(time < time_close2)
    {   /*  close needle  */
        f_lift = 0.01 + 0.99*(time - time_close2)/(time_close1 - time_close2);
    }
    else
    {   /*  squeeze out  */
        f_lift = 0.01;
    }
    return f_lift;
}



real nozzle_gdi_scale_velocity(const Particle *particle, const Injection *injection, const real s)
{
    /*  particle time  */
    const real time
        = particle->state.time - injection->unsteady_start + s;

    /*  injection duration  */
    const real duration
        = injection->unsteady_stop - injection->unsteady_start;

    /*  injection timing  */
    const real time_open1  = TIME_VEL_OPEN1;
    const real time_open2  = TIME_VEL_OPEN2;
    const real time_open3  = TIME_VEL_OPEN3;
    const real time_close1 = TIME_VEL_CLOSE1 + duration;
    const real time_close2 = TIME_VEL_CLOSE2 + duration;


    /*  SCALE NEEDLE VELOCITY  */
    real f_vel;

    if(time < time_open1)
    {   /*  delay  */
        f_vel = 0.50;
    }
    else
    if(time < time_open2)
    {   /*  lift needle  */
        f_vel = 0.50 + 0.20*(time - time_open1)/(time_open2 - time_open1);
    }
    else
    if(time < time_open3)
    {   /*  lift needle  */
        f_vel = 0.70 + 0.30*(time - time_open2)/(time_open3 - time_open2);
    }
    else
    if(time < time_close1)
    {   /*  open  */
        f_vel = 1.00;
    }
    else
    if(time < time_close2)
    {   /*  close needle  */
        f_vel = 0.80 + 0.20*(time - time_close2)/(time_close1 - time_close2);
    }
    else
    {   /*  squeeze out  */
        f_vel = 0.80;
    }
    return f_vel;
}



DEFINE_DPM_INJECTION_INIT(nozzle_gdi, injection)
{
#if !RPHOST
    if(!injection->p_init)
        return;
    
    const real div = (real)count_p_list(injection->p_init);
    
    Particle* particle; 
    loop(particle, injection->p_init)
    {
        /*  stagger temporally  */
        const real s = uniform_rand()*dpm_par.time_step;
        
        /*  scale lift  */
        real f_lift = nozzle_gdi_scale_lift(particle, injection, s);
        /*  scale velocity  */
        real f_vel  = nozzle_gdi_scale_velocity(particle, injection, s);
        /*  scale flowrate  */
        real f_rate = f_vel*f_lift;
        /*  scale mass  */
        real f_mass = pow3(f_lift);
        /*  scale number in parcel  */
        real f_num  = f_rate/f_mass;

        /*  scale parcel  */        
        particle->flow_rate *= f_rate;        
        particle->number_in_parcel *= f_num;

        /*  scale particle  */
        particle->state.diam = particle->init_state.diam * f_lift;
        particle->state.mass = particle->init_state.mass * f_mass;

        /*  stagger tangentially  */
        real phi = (0.5-uniform_rand()) * DOMAIN_ANGLE/div;             

        rotate_z(particle->init_state.V,
                    cos(phi), sin(phi),
                        particle->state.V);                     

		rotate_z(particle->init_state.pos,
                    cos(phi), sin(phi),
                        particle->state.pos);

        /*  stagger radially  */						
        int i;
        for(i=0; i<ND_3; ++i)
        {   
            particle->state.V[i] *= f_vel;
            particle->state.pos[i] += particle->state.V[i]*s;              
        }           

		
        /*  status message  */
        Message("Injection scaling:   needle lift/velocity/flow rate/number in parcel:   %e / %e / %e / %e\n", f_lift, f_vel, f_rate, f_num);	
    }
#endif
    return;
}
