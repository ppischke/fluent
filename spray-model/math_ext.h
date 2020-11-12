#ifndef __MATH_EXT_H
#define __MATH_EXT_H

#include "udf.h"
#include "random.h"
#include "math.h"
#include "math_sse.h"
#include "dsyevh3.h"
#include "dsyevj3.h"
#include "dsyevq3.h"



#ifndef __cplusplus
/*  c-language bool  */
typedef enum bool_enum
{
    false = 0,
    true = -1
}
bool;
#endif



#define M_ITER_ACCURACY  1.0e-4



#define M_SQRT_PI   1.772453850905516027298167483341

#define M_EXP_M1    0.367879441171442321595523770161    /*  exp(-1.0)  */
#define M_EXP_M2    0.135335283236612691893999494972    /*  exp(-2.0)  */



/*  glibc function declarations  */
double erf(double);
double tgamma(double);
double fmax(double, double);
double fmin(double, double);
double round(double);

float  erff(float);
float  tgammaf(float);
float  fmaxf(float, float);
float  fminf(float, float);
float  roundf(float);

#ifndef isnan
bool   isnan(double);
bool   isinf(double);
bool   isfinite(double);
bool   isnanf(float);
bool   isinff(float);
bool   isfinitef(double);
#endif



#if !RP_DOUBLE
#define erf(x)    erff(x)
#define tgamma(x) tgammaf(x)
#define fmax(x,y) fmaxf(x,y)
#define fmin(x,y) fminf(x,y)
#define round(x)  roundf(x)
#define isnan(x)  isnanf(x)
#define isinf(x)  isinff(x)
#endif



/*  distributions  */
#define DISTRIB_SIGMA                   3.0
#define DISTRIB_MIN                     0.0002
#define DISTRIB_MAX                     0.9998

real uniform_rand();
real uniform_rand_lim(const real, const real);
real uniform_rand_vector(real*);

#define gauss_rand gauss_random
real gauss_rand_lim(const real);
real gauss_rand_vector(real*, const real);

/*  distribution types  */
typedef enum distribution_type_enum
{    
    DIST_CHI_SQR,    
    DIST_ROSIN_R
}
distribution_type_t;

real chi_sqr_rand
    (const real, const real uniform);
real chi_sqr_diam
    (const real);
real rosin_r_rand
    (const real, const real, const real uniform);
real rosin_r_diam
    (const real, const real);
    


/*  fast power functions  */
real pow2(const real);
real pow3(const real);
real pow4(const real);
real pow5(const real);
real pow6(const real);
real sgn (const real);
real fminmax(const real, const real, const real);



/*  fast power macros  */
#define POW1(x) ((x))
#define POW2(x) ((x)*(x))
#define POW3(x) ((x)*(x)*(x))
#define POW4(x) ((x)*(x)*(x)*(x))
#define POW5(x) ((x)*(x)*(x)*(x)*(x))



/*  vector operations  */
real v_sqr
    (const real v[ND_3]);
real v_abs
    (const real v[ND_3]);
real v_rel_sqr
    (const real v[ND_3], const real w[ND_3]);
real v_rel_abs
    (const real v[ND_3], const real w[ND_3]);
void v_rel
    (const real v[ND_3], const real w[ND_3], real x[ND_3]);
void v_sum
    (const real v[ND_3], const real w[ND_3], real x[ND_3]);


real dot_product        /*  dot product    */
    (const real v[ND_3], const real w[ND_3]);
void cross_product      /*  cross product  */
    (const real v[ND_3], const real w[ND_3], real x[ND_3]);
void scale              /*  scale vector to length  */
    (const real v[ND_3], const real s, real x[ND_3]);
void scale_to           /*  scale vector to length  */
    (const real v[ND_3], const real l, real x[ND_3]);
real normalize          /*  scale vector to unity   */
    (const real v[ND_3], real x[ND_3]);


void create_random_base
    (real M[ND_3][ND_3]);
void create_orthonormal
    (const real x[ND_3], real M[ND_3][ND_3]);
void create_orthonormal_xy
    (const real x[ND_3], const real xy[ND_3], real M[ND_3][ND_3]);
void create_identity
    (real M[ND_3][ND_3]);
void identity_vector
    (real v[ND_3]);

void diagonalize
    (real A[ND_3][ND_3], real U[ND_3][ND_3], real AD[ND_3][ND_3]);
void invert_matrix
    (real M[ND_3][ND_3], real MI[ND_3][ND_3]);
void transpose_matrix
    (real M[ND_3][ND_3], real MT[ND_3][ND_3]);
void multiply_matrices
    (real M[ND_3][ND_3], real N[ND_3][ND_3], real A[ND_3][ND_3]);
void multiply_matrix_diag_matrix_t
    (real M[ND_3][ND_3], const real diag[ND_3], real A[ND_3][ND_3]);
void multiply_matrix_col_vector
    (real M[ND_3][ND_3], const real cv[ND_3], real a[ND_3]);
void multiply_row_vector_matrix
    (const real rv[ND_3], real M[ND_3][ND_3], real a[ND_3]);
void matrix_mean
    (real A[ND_3][ND_3], real B[ND_3][ND_3], real AB[ND_3][ND_3]);

#define gauss(A,X,B) Gauss_Elimination(ND_3,ND_3,(real*)A,(real*)X,(real*)B,0)
    
real det
    (real M[ND_3][ND_3]);
real distance_norm
    (const real x[ND_3], real norm[ND_3][ND_3], const real y[ND_3]);

void rotate_3d
    (const real*, const real, const real, real*, const int);
void rotate_2d
    (const real*, const real, const real, real*);

#define rotate_x(v,c,s,w) rotate_3d(v,c,s,w,0)
#define rotate_y(v,c,s,w) rotate_3d(v,c,s,w,1)
#define rotate_z(v,c,s,w) rotate_3d(v,c,s,w,2)

#define rotate_x_inv(v,c,s,w) rotate_3d(v,c,-s,w,0)
#define rotate_y_inv(v,c,s,w) rotate_3d(v,c,-s,w,1)
#define rotate_z_inv(v,c,s,w) rotate_3d(v,c,-s,w,2)

#define eigen dsyevh3
#endif
