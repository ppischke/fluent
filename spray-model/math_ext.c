#include "udf.h"
#include "dpm.h"
#include "math.h"
#include "math_ext.h"
#include "math_sse.h"





/*  truncated uniform random  */
real uniform_rand()
{
    return ((real)rand()) / ((real)(RAND_MAX) + 1.0);
}



real uniform_rand_lim(const real min, const real max)
{
    return min + (max-min) * uniform_rand();
}



/*  truncated gauss random  */
real gauss_rand_lim(const real sigma)
{
    real r;    
    do
    {   /*  sample gaussian  */
        r = gauss_rand();
    }   /*  reject gaussian  */
    while(fabs(r) > sigma);
    /*  return gauss random  */
    return r;
}



/*  truncated gauss random vector  */
real gauss_rand_vector(real* v, const real sigma)
{
    real r;
    do
    {   /*  sample vector  */
        int i;
        for(i=0; i<ND_3; ++i)
        {
            v[i] = gauss_rand();
        }
        r = v_abs(v);       
    }   /*  reject vector  */
    while(r > sigma);
    /*  return vector absolute  */
    return r;
}



/*	uniform random vector  */
real uniform_rand_vector(real* v)
{
	real r = pow(uniform_rand(), 1.0/3.0);
	
    int i;
    for(i=0; i<ND_3; ++i)
    {
        v[i] = gauss_rand();
    }
	scale_to(v, r, v);
    /*  return vector absolute  */
    return r;
}




/*  rosin rammler distribution  */
real rosin_r_rand(const real d32, const real q, const real uniform)
{    
    /*  invert cumulative probability function  */
    return pow(-log(1.0-uniform), 1.0/q) * d32*tgamma(1.0+2.0/q)/tgamma(1.0+3.0/q);
}

real rosin_r_diam(const real d32, const real q)
{
    real d_max = rosin_r_rand(d32, q, DISTRIB_MAX);    
    real d;
    do
    {   /*  sample diameter  */
        d = rosin_r_rand(d32, q, uniform_rand_lim(DISTRIB_MIN,DISTRIB_MAX));
    }   /*  reject diameter  */
    while(uniform_rand() > pow3(d/d_max));
    
    return d;
}



/*  chi-square distribution  */
real chi_sqr_rand(const real d32, const real uniform)
{            
    /*  dimensionless diam   */
    real x_diam = 6.0 * sqrt(uniform);
    real x_prev;
        
    /*  iterative solution   */
    do
    {   x_prev = x_diam;
        /*  invert cumulative probability function  */    
        x_diam = log((1.0 + x_diam + pow2(x_diam)/2.0 + pow3(x_diam)/6.0) / (1.0-uniform));
    }
    while(fabs(x_diam-x_prev)/x_diam > M_ITER_ACCURACY);

    /*  return random diam  */
    return x_diam*d32/6.0;
}

real chi_sqr_diam(const real d32)
{
    real d_max = chi_sqr_rand(d32, DISTRIB_MAX);    
    real d;
    do
    {   /*  sample diameter  */
        d = chi_sqr_rand(d32, uniform_rand_lim(DISTRIB_MIN,DISTRIB_MAX));
    }   /*  reject diameter  */
    while(uniform_rand() > pow3(d/d_max));
    
    return d;
}





/*  FAST POWER FUNCTIONS  */
real pow2(const real a)
{
    return a*a;
}

real pow3(const real a)
{
    return a*a*a;
}

real pow4(const real a)
{
    return pow2(pow2(a));
}

real pow5(const real a)
{
    return pow2(pow2(a)) * a;
}

real pow6(const real a)
{
    return pow2(pow3(a));
}

real sgn(const real a)
{
    if(a > 0.0)
        return  1.0;
    if(a < 0.0)
        return -1.0;
    /*  pass through nan/inf  */
    return a;
}

real fminmax(const real xmin, const real x, const real xmax)
{
    return fmax(fmin(x,xmax),xmin);
}





/*  VECTOR OPERATIONS  */
real v_sqr(const real v[ND_3])
{
    real v_sq = 0.0;

    int i;
    for(i=0; i<ND_3; ++i)
        v_sq += v[i]*v[i];

    return v_sq;
}


real v_abs(const real v[ND_3])
{
    return sqrt(v_sqr(v));
}


real v_rel_sqr(const real v[ND_3], const real w[ND_3])
{
    real v_rel_sqr = 0.0;

    int i;
    for(i=0; i<ND_3; ++i)
        v_rel_sqr += pow2(v[i]-w[i]);

    return v_rel_sqr;
}


real v_rel_abs(const real v[ND_3], const real w[ND_3])
{
    return sqrt(v_rel_sqr(v,w));
}


void v_rel(const real v[ND_3], const real w[ND_3], real x[ND_3])
{
    int i;
    for(i=0; i<ND_3; ++i)
        x[i] = v[i]-w[i];

    return;
}


void v_sum(const real v[ND_3], const real w[ND_3], real x[ND_3])
{
    int i;
    for(i=0; i<ND_3; ++i)
        x[i] = v[i]+w[i];

    return;
}


real distance_norm(const real x[ND_3], real norm[ND_3][ND_3], const real y[ND_3])
{
    real v = 0.0;
    int i;
    for(i=0; i<ND_3; ++i)
    {
        real u = 0.0;
        int k;
        for(k=0; k<ND_3; ++k)
        {
            u += x[k]*norm[k][i];
        }
        v += u*y[i];
    }
    return v;
}


void scale(const real v[ND_3], const real s, real w[ND_3])
{
    int i;
    for(i=0; i<ND_3; ++i)
        w[i] = v[i]*s;

    return;
}

void scale_to(const real v[ND_3], const real l, real x[ND_3])
{
    real s = v_abs(v);
    scale(v, l/s, x);

    return;
}

real normalize(const real v[ND_3], real x[ND_3])
{
    real s = v_abs(v);
    scale(v, 1.0/s, x);

    return s;
}


/*  3D VECTOR OPERATIONS  */
real dot_product(const real v[ND_3], const real w[ND_3])
{
    real dp = 0.0;

    int i;
    for(i=0; i<ND_3; ++i)
        dp += v[i]*w[i];

    return dp;
}


void cross_product(const real v[ND_3], const real w[ND_3], real cp[ND_3])
{
    int i;
    int j;
    int k;

    for(i=0; i<ND_3; ++i)
    {
        j = (i+1)%ND_3;
        k = (i+2)%ND_3;

        cp[i] = v[j]*w[k] - w[j]*v[k];
    }
    return;
}


void rotate_2d(const real v[ND_3], const real c, const real s, real x[ND_3])
{
    x[0] = v[0]*c - v[1]*s;
    x[1] = v[0]*s + v[1]*c;
    x[2] = v[2];
    return;
}


void rotate_3d(const real v[ND_3], const real c, const real s, real x[ND_3], const int i)
{
    int j = (i+1)%ND_3;
    int k = (i+2)%ND_3;

    x[i] = v[i];
    x[j] = v[j]*c - v[k]*s;
    x[k] = v[j]*s + v[k]*c;
    return;
}


void diagonalize(real A[ND_3][ND_3], real U[ND_3][ND_3], real AD[ND_3][ND_3])
{
    real AU[ND_3][ND_3];

    multiply_matrices(A,U,AU);

    int i;
    for(i=0; i<ND_3; ++i)
    {
        int j;
        for(j=i; j<ND_3; ++j)
        {
            AD[i][j] = 0.0;

            int k;
            for(k=0; k<ND_3; ++k)
            {
                AD[i][j] += U[k][i]*AU[k][j];
            }
            /*  make use of symmetry  */
            AD[j][i] = AD[i][j];
        }
    }
    return;
}


/*  WAS: v_matrix  */
void multiply_matrix_col_vector(real M[ND_3][ND_3], const real cv[ND_3], real a[ND_3])
{   /*
    if(v == w)
        return;
    */
    int i;
    for(i=0; i<ND_3; ++i)
    {
        a[i] = 0.0;

        int k;
        for(k=0; k<ND_3; ++k)
        {
            a[i] += M[i][k]*cv[k];
        }
    }
    return;
}


/*  WAS: v_matrix_inv  */
void multiply_row_vector_matrix(const real rv[ND_3], real M[ND_3][ND_3], real a[ND_3])
{   /*
    if(v == w)
        return;
    */
    int i;
    for(i=0; i<ND_3; ++i)
    {
        a[i] = 0.0;

        int k;
        for(k=0; k<ND_3; ++k)
        {
            a[i] += rv[k]*M[k][i];
        }
    }
    return;
}


void multiply_matrix_diag_matrix_t(real M[ND_3][ND_3], const real diag[ND_3], real A[ND_3][ND_3])
{
    int i;
    int j;
    int k;

    for(i=0; i<ND_3; ++i)
    {
        for(j=i; j<ND_3; ++j)
        {
            A[i][j] = 0.0;

            for(k=0; k<ND_3; ++k)
            {
                A[i][j] += M[i][k] * diag[k] * M[j][k];
            }
            /*  make use of symmetry  */
            A[j][i] = A[i][j];
        }
    }
}


void multiply_matrices(real M[ND_3][ND_3], real N[ND_3][ND_3], real A[ND_3][ND_3])
{
    int i;
    int j;
    int k;

    for(i=0; i<ND_3; ++i)
    {
        for(j=0; j<ND_3; ++j)
        {
            A[i][j] = 0.0;

            for(k=0; k<ND_3; ++k)
            {
                A[i][j] += M[i][k]*N[k][j];
            }
        }
    }
}


void transpose_matrix(real M[ND_3][ND_3], real MT[ND_3][ND_3])
{
    int i;
    for(i=0; i<ND_3; ++i)
    {
        int j;
        for(j=0; j<ND_3; ++j)
        {
            MT[i][j] = M[j][i];
        }
    }
}


void invert_matrix(real M[ND_3][ND_3], real MI[ND_3][ND_3])
{
    real d = det(M);

    if(d == 0.0)
    {
        /*  matrix not invertible  */
        memset(MI,0,sizeof(real)*ND_3*ND_3);
        return;
    }

    int i;
    int j;

    for(i=0; i<ND_3; ++i)
    {
        for(j=0; j<ND_3; ++j)
        {
            /*  sub-determinant  */
            MI[j][i] = ( M[(i+1)%3][(j+1)%3] * M[(i+2)%3][(j+2)%3]
                       - M[(i+2)%3][(j+1)%3] * M[(i+1)%3][(j+2)%3] ) / d;
        }
    }
}


void matrix_mean(real A[ND_3][ND_3], real B[ND_3][ND_3], real AB[ND_3][ND_3])
{
    int i;
    int j;
    
    for(i=0; i<ND_3; ++i)
    {
        for(j=0; j<ND_3; ++j)
        {
            AB[i][j] = 0.5*(A[i][j] + B[i][j]);
        }
    }
    return;
}


void identity_vector(real v[ND_3])
{
    int i;
    for(i=0; i<ND_3; ++i)
    {
        v[i] = 1.0;
    }
    return;
}


void create_identity(real A[ND_3][ND_3])
{
    memset(A,0,sizeof(real)*ND_3*ND_3);

    int i;
    for(i=0; i<ND_3; ++i)
    {        
        A[i][i] = 1.0;
    }
    return;
}


void create_orthonormal(const real x[ND_3], real A[ND_3][ND_3])
{
    real xy[ND_3];
    /*  vector 2 from randomization  */
    int i;
    for(i=0; i<ND_3; ++i)
    {
        xy[i] = gauss_rand();
    }   
    /*  normal base  */
    create_orthonormal_xy(x, xy, A);
    return;
}


void create_orthonormal_xy(const real x[ND_3], const real xy[ND_3], real A[ND_3][ND_3])
{
    real E[ND_3][ND_3];

    /*  vector 1 from normalization  */
    normalize(x, E[0]);
    /*  vector 3 from cross_product  */
    cross_product(E[0],  xy,  E[2]);
    /*  vector 2 from cross product  */
    cross_product(E[2], E[0], E[1]);

    /*  normal base  */
    int i;
    for(i=1; i<ND_3; ++i)
    {
        normalize(E[i],E[i]);
    }
    transpose_matrix(E, A); 
    return;
}


void create_random_base(real A[ND_3][ND_3])
{       
    int i;
    for(i=0; i<ND_3; ++i)
    {
        A[0][i] = gauss_rand();
        A[1][i] = gauss_rand();
    }
    /*  random base  */
    cross_product(A[0],A[1],A[2]);
    cross_product(A[1],A[2],A[0]);
    
    /*  normal base  */
    for(i=0; i<ND_3; ++i)
    {
        normalize(A[i],A[i]);
    }
}


real det(real A[ND_3][ND_3])
{
    real p;     /*  positive diagonal  */
    real n;     /*  negative diagonal  */
    real d = 0.0;

    int i;
    for(i=0; i<ND_3; ++i)
    {
        p = n = 1.0;

        /*  multiply diagonals  */
        int k;
        for(k=0; k<ND_3; ++k)
        {
            p *= A[k][(i+k)%ND_3];
            n *= A[k][(i-k+ND_3)%ND_3];
        }
        /*  determinant  */
        d += (p-n);
    }
    return d;
}
