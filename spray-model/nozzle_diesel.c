#include "udf.h"
#include "dpm.h"
#include "math.h"
#include "math_ext.h"



#define TIME_VEL_OPEN1      0e-6
#define TIME_VEL_OPEN2    250e-6
#define TIME_VEL_OPEN3    350e-6
#define TIME_VEL_CLOSE1  -200e-6
#define TIME_VEL_CLOSE2    -0e-6





real nozzle_diesel_scale_velocity
    (const Particle *particle, const Injection *injection, const real s);





real nozzle_diesel_scale_velocity(const Particle *particle, const Injection *injection, const real s)
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
        f_vel = 0.40;
    }
    else
    if(time < time_open2)
    {   /*  lift needle  */
        f_vel = 0.40 + 0.50*(time - time_open1)/(time_open2 - time_open1);
    }
    else
    if(time < time_open3)
    {   /*  lift needle  */
        f_vel = 0.90 + 0.10*(time - time_open2)/(time_open3 - time_open2);
    }
    else
    if(time < time_close1)
    {   /*  open  */
        f_vel = 1.00;
     }
     else
     if(time < time_close2)
     {   /*  close needle  */
          f_vel = 0.40 + 0.60*(time - time_close2)/(time_close1 - time_close2);
     }
     else
     {    /*    squeeze out  */
          f_vel = 0.40;
     }
    return f_vel;
}



DEFINE_DPM_INJECTION_INIT(nozzle_diesel, injection)
{
#if !RP_HOST
    if(!injection->p_init)
        return;

    Particle* particle; 
    loop(particle, injection->p_init)
    {   
        /*  stagger temporally  */
        const real s = uniform_rand()*dpm_par.time_step;

        /*  scale lift  */
        real f_lift = 1.0;      
        /*  scale velocity  */
        real f_vel  = nozzle_diesel_scale_velocity(particle, injection, s);
        /*  scale flowrate  */
        real f_rate = f_vel;
        /*  scale number in parcel  */
        real f_num  = f_rate;

        /*  scale parcel  */        
        particle->flow_rate *= f_rate;        
        particle->number_in_parcel *= f_num;

		/*  stagger  */
		memcpy(particle->state.V, particle->init_state.V, sizeof(real)*ND_3);
		memcpy(particle->state.pos, particle->init_state.pos, sizeof(real)*ND_3);

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
