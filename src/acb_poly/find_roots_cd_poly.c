/*
    Copyright (C) 2025 Guillaume Moroz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include <stdlib.h>
#include <math.h>
#include <string.h> /* for memcpy */

#if _MSC_VER
#pragma float_control(precise, off)
#else
#pragma GCC optimize("Ofast,unroll-loops")
#endif

/*****************************************
 * Intermediate variables data structure *
 * ***************************************/

/* This field contains 13 double and 1 slong array pointers */
typedef struct {
    /* arrays storing the input with real and imaginary parts separated */
    double *z_r;
    double *z_i;
    double *p_r;
    double *p_i;
    /* arrays for the intermediate numerator (vp) and denominator (dk) used
     * in the Weierstrass-Durand-Kerner iterator */
    double *vp_r;
    double *vp_i;
    double *wdk_r;
    double *wdk_i;
    /* arrays for the permuted points such that intermediate products do not overflow */
    double *pz_r;
    double *pz_i;
    /* arrays for the reciprocal polynomial to evaluate on the inverse points of magnitude >= 1 */
    double *rp_r;
    double *rp_i;
} intermediate_variables;

static void _intermediate_variables_init(intermediate_variables * v, slong n)
{
    double * mem_v = flint_malloc(12*n*sizeof(double));

    /* Set the memory addresses */
    v->z_r   = mem_v+0*n;
    v->z_i   = mem_v+1*n;
    v->p_r   = mem_v+2*n;
    v->p_i   = mem_v+3*n;
    v->vp_r  = mem_v+4*n;
    v->vp_i  = mem_v+5*n;
    v->wdk_r = mem_v+6*n;
    v->wdk_i = mem_v+7*n;
    v->pz_r  = mem_v+8*n;
    v->pz_i  = mem_v+9*n;
    v->rp_r  = mem_v+10*n;
    v->rp_i  = mem_v+11*n;
}

static void _intermediate_variables_clear(intermediate_variables * v)
{
    flint_free(v->z_r);
}

/*********************
 * Utility functions *
 *********************/

static void _d_swap(double* a, double* b)
{
    double t = *a;
    *a = *b;
    *b = t;
}

/* Reorder the z : the values at index < n-c satisfy |z| <= 1, and at index >= n-c they satisfy |z|>1 */
static slong cd_poly_partition_pivot(double* z_r, double* z_i, slong n)
{
    slong i, c = 0;
    for(i=0; i<n; i++) {
        if(z_r[i]*z_r[i] + z_i[i]*z_i[i] > 1) {
            c++;
        } else if(c>0) {
            _d_swap(z_r + i-c, z_r + i);
            _d_swap(z_i + i-c, z_i + i);
        }
    }
    return n-c;
}

/* Order that maximizes the chance of concluding with the first comparison */
static int _compare(const void * first, const void * second)
{
    double f, s, irr = sqrt(2), epsilon=ldexp(1,-25);
    int res;
    f = (*(double*)first)  * irr + (*((double*)first+1));
    s = (*(double*)second) * irr + (*((double*)second+1));
    if(f < s - epsilon) {
        res = -1;
    } else if (f > s + epsilon) {
        res = 1;
    } else {
        f = (*(double*)first)  - irr * (*((double*)first+1));
        s = (*(double*)second) - irr * (*((double*)second+1));
        res = (f > s + epsilon) - (f < s - epsilon);
    }
    return res;
}

/* Remove the leading and trailing zero coeffs of the input polynomial */
static slong _trim_zeros(double * z, slong * fz, const double * p, slong n)
{
    slong i, j;

    /* trailing terms with zero coefficients add roots at infinity*/
    for(i=n-1; i>=1 && p[2*i] == 0.0 && p[2*i+1] == 0.0; i--) {
        z[2*(i-1)] = ldexp(1,1023);
        z[2*(i-1)+1] = ldexp(1,1023);
    }

    /* leading terms with zero coefficients add roots at zero*/
    for(j=0; j<i && p[2*j] == 0.0 && p[2*j+1] == 0.0; j++) {
        z[2*j] = 0.0;
        z[2*j+1] = 0.0;
    }
    *fz = j;

    return i-j+1;
}

static void _vector_inverse(double* iz_r, double* iz_i, const double* z_r, const double* z_i,
                           slong n_start, slong n_end)
{
    for(slong i=n_start; i<n_end; i++) {
        double s, t;
        /* use Smith's algorithm to reduce the risk of overflow */
        if(fabs(z_r[i]) >= fabs(z_i[i]) ) {
            s = z_i[i]/z_r[i];
            t = 1/(z_r[i] + z_i[i] * s);
            iz_r[i] = t;
            iz_i[i] = -s*t;
        } else {
            s = z_r[i]/z_i[i];
            t = 1/(z_r[i] * s + z_i[i]);
            iz_r[i] = s*t;
            iz_i[i] = -t;
        }
    }
}

static inline void divide(double* res_e, double* res_f, double a, double b, double c, double d)
{
    double r, den;
    /* use improved Smith's algorithm to reduce the risk of overflow */
    if(fabs(c) >= fabs(d) ) {
        r = d/c;
        den = c + r*d;
        //if(r != 0) {
            *res_e = (a+b*r)/den;
            *res_f = (b-a*r)/den;
        //} else {
        //    *res_e = (a+d*(b/c))/den;
        //    *res_f = (b-d*(a/c))/den;
        //}
    } else {
        r = c/d;
        den = d + r*c;
        //if(r != 0) {
            *res_e = (a*r+b)/den;
            *res_f = (b*r-a)/den;
        //} else {
        //    *res_e = (c*(a/d)+b)/den;
        //    *res_f = (c*(b/d)-a)/den;
        //}
    }
}


/* Initialize the members of the double_field */
static void _set_poly(intermediate_variables * v, const double * p, slong n)
{
    slong i;

    /* Unzip the real, imaginary parts and numerical valuation of the input polynomial coefficients */
    for(i=0; i<n; i++) {
        v->p_r[i] = p[2*i];
        v->p_i[i] = p[2*i+1];
        v->rp_r[i] = p[2*(n-1-i)];
        v->rp_i[i] = p[2*(n-1-i)+1];
    }
}


/***************************************************
 * Initial values using Newton Polygon computation *
 ***************************************************/

static int _is_clockwise(double a, double b, double x, double y, double u, double v)
{
    return v*(x-a) + y*(a-u) + b*(u-x) <= 0;
}

/* Computes the lower convex hull of the points (i, y[i]) */
static slong lower_convex_hull(slong* result, const double* y, slong n)
{
    slong i, k=0;
    for (i=0; i<n; i++) {
        while(k>1 && _is_clockwise(result[k-2], y[result[k-2]],
                                              result[k-1], y[result[k-1]],
                                              i, y[i])) {
            k--;
        }
        result[k] = i;
        k++;
    }
    return k;
}

static void _initial_values_from_newton_polygon(double* z_r, double* z_i, const double* p_r, const double* p_i,
                                                const slong* np, slong m)
{
    slong i, j, q, p = 0;
    double d, invd, magp, argp, magq, argq, magz, theta,
           pi = acos(-1),
           epsilon = ldexp(1,-50);
    for(i=1; i<m; i++) {
        q = np[i];
        d = q-p;
        invd = 1/d;
        magp = hypot(p_r[p], p_i[p]);
        argp = atan2(p_i[p], p_r[p]);
        magq = hypot(p_r[q], p_i[q]);
        argq = atan2(p_i[q], p_r[q]);
        magz = pow(magp,invd)/pow(magq, invd); //two pow to avoid magp/magq overflow
        if(magz == 0.0) {
            magz = ldexp(1, -500/(q-p));
        }
        for(j=0; j<q-p; j++) {
            /* theta = arg of solutions of x^d = -p[p]/p[q]
             *       = arg(-p[p]/p[q])/d + 2 pi j/d */
            theta = (2*pi*(j+0.5)+argp-argq)/d;
            //theta = 2*pi*j/d;
            /* epsilon is used to avoid initial values being symmetrical to the
             * real axis (see Section 4 of Aberth's paper) */
            z_r[p+j] = magz * (cos(theta)+(1+j)*epsilon) ;
            z_i[p+j] = magz * (sin(theta)+(1+j)*epsilon);
        }
        p = q;
    }
}

#define LOG_ROUNDING 0

void cd_poly_roots_initial_values(double * z_r, double * z_i, const double * p_r, const double * p_i, slong n, const double * z0, slong d)
{
    double epsilon=ldexp(1,-1000);
    double * malp;
    slong * np;
    slong i=0, j, n_np;
    int k=LOG_ROUNDING;

    if(z0 != NULL) {
        for(i=0, j=0; i<n-1 && j<d; j++) {
            if (z0[2*j] != 0 || z0[2*j+1] != 0) {
                z_r[i] = z0[2*j];
                z_i[i] = z0[2*j+1];
                i++;
            }
        }
    }

    if(z0 == NULL || i < n-1) {
        malp = flint_malloc(n*sizeof(double));
        np = flint_malloc(n*sizeof(slong));
        /* Unzip the real, imaginary parts and numerical valuation of the input polynomial coefficients */
        for(i=0; i<n; i++) {
            /* We avoid infinite values so that the code is compatible with
             * unsafe math optimizations */
            malp[i] = -log(hypot(p_r[i],p_i[i])+epsilon);
            /* Rounding the valuation helps avoiding nearby initial radii */
            malp[i] = ldexp(round(ldexp(malp[i],k)),-k);
        }

        /* Computes the Newton polygon as indices in np of length n_np */
        n_np = lower_convex_hull(np, malp, n);

        /* Computes a first estimation of the roots */
        _initial_values_from_newton_polygon(z_r, z_i, p_r, p_i, np, n_np);

        flint_free(malp);
        flint_free(np);
    }
}

#undef LOG_ROUNDING

/**********************************************
 * Weierstrass-Durand-Kerner update functions *
 **********************************************/


/* simd parameters for speeding up Horner and Weirstrass weight computation
 *   Rbits size of registers
 *   Nr number of registers */
#ifdef __AVX512CD__
#define Rbits 512
#define Nr 32
#elif __AVX__
#define Rbits 256
#define Nr 16
#elif __SSE__
#define Rbits 128
#define Nr 16
#else
#define Rbits 64
#define Nr 16
#endif
#define Rbytes (Rbits/8)
#define Nd (Rbytes/8)

#define CDHBlock (Nd*Nr/4)
/* Horner evaluation */
void cd_poly_horner(double* results_r, double* results_i,
                    const double* values_r, const double* values_i, slong n_start, slong n_end,
                    const double* coefficients_r, const double* coefficients_i, slong n)
{
    int p, k;
    slong i, j;
    for(i=n_start; i < n_end; i+=CDHBlock) {
        double x[CDHBlock] = {0};
        double y[CDHBlock] = {0};
        double u[CDHBlock] = {0};
        double v[CDHBlock] = {0};
        p = (i+CDHBlock<=n_end) ? CDHBlock : ((n_end-n_start)%CDHBlock);
        memcpy(u, values_r+i, p*sizeof(double));
        memcpy(v, values_i+i, p*sizeof(double));
        for(j=0; j<n; j++) {
            double a = coefficients_r[n-1-j];
            double b = coefficients_i[n-1-j];
            for(k=0; k<CDHBlock; k++) {
                double s,t;
                s = a + u[k]*x[k] - v[k]*y[k];
                t = b + u[k]*y[k] + v[k]*x[k];
                x[k] = s;
                y[k] = t;
            }
        }
        memcpy(results_r+i, x, p*sizeof(double));
        memcpy(results_i+i, y, p*sizeof(double));
    }
}
#undef CDHBlock

#define CDWBlock (Nd*Nr/4)
/* Weights for the Durand-Kerner or Weierstrass iteration
 * The multiplications are done according to the input order */
void cd_poly_weierstrass_distinct_orders(double* restrict results_r, double* restrict results_i,
                                         double lc_r, double lc_i,
                                         const double* restrict col_values_r, const double* restrict col_values_i, slong d,
                                         const double* restrict row_values_r, const double* restrict row_values_i, slong n_start, slong n_end)
{
    int p, k;
    slong i, j;
    for(i=n_start; i<n_end; i++){
        results_r[i] = 1;
        results_i[i] = 0;
    }
    for(i=n_start; i < n_end; i+=CDWBlock) {
        double xr[CDWBlock]  = {0};
        double yr[CDWBlock]  = {0};
        double u[CDWBlock]   = {0};
        double v[CDWBlock]   = {0};

        p = (i+CDWBlock<=n_end) ? CDWBlock : ((n_end-n_start) % CDWBlock);
        memcpy(xr, results_r+i, p*sizeof(double));
        memcpy(yr, results_i+i, p*sizeof(double));
        memcpy(u, row_values_r+i, p*sizeof(double));
        memcpy(v, row_values_i+i, p*sizeof(double));
        for(j = 0; j < d; j++){
            double a = col_values_r[j];
            double b = col_values_i[j];
            for(k=0; k<CDWBlock; k++) {
                double e, f, s, t;
                // compute the difference e + i f = (u + i v) - (a + i b)
                e = u[k] - a;
                f = v[k] - b;
                // when u + i v == a + i b, let e + i f = 1 instead
                e = ((e==0) & (f==0))? 1 : e;
                // accumulate the multiplications by e + i f = (u + i v) - (a + i b)
                s = e*xr[k] - f*yr[k];
                t = e*yr[k] + f*xr[k];
                xr[k] = s;
                yr[k] = t;
            }
        }
        memcpy(results_r+i, xr, p*sizeof(double));
        memcpy(results_i+i, yr, p*sizeof(double));
    }
    for(i=n_start; i<n_end; i++) {
        double s, t;
        s = lc_r*results_r[i] - lc_i*results_i[i];
        t = lc_r*results_i[i] + lc_i*results_r[i];
        results_r[i] = s;
        results_i[i] = t;
    }
}
#undef CDWBlock

#define CDWBlock (Nd*Nr*4)
/* Weights for the Durand-Kerner or Weierstrass iteration
 * Faster by a small margin (less than 10%), but cannot be adapted to a permutation */
void cd_poly_weierstrass(double* restrict results_r, double* restrict results_i,
                         double lc_r, double lc_i,
                         const double* values_r, const double* values_i,
                         slong n_start, slong n_end, slong d)
{
    int p, q, k, m;
    slong i, j;
    for(i=n_start; i<n_end; i++) {
        results_r[i] = lc_r;
        results_i[i] = lc_i;
    }
    for(i=n_start; i < n_end; i+=CDWBlock) {
        double xr[CDWBlock]  = {0};
        double yr[CDWBlock]  = {0};
        double xc[CDWBlock]  = {0};
        double yc[CDWBlock]  = {0};
        double u[CDWBlock]   = {0};
        double v[CDWBlock]   = {0};
        double a[2*CDWBlock] = {0};
        double b[2*CDWBlock] = {0};

        p = (i+CDWBlock<=n_end) ? CDWBlock : ((n_end-n_start) % CDWBlock);
        memcpy(xr, results_r+i, p*sizeof(double));
        memcpy(yr, results_i+i, p*sizeof(double));
        memcpy(u, values_r+i, p*sizeof(double));
        memcpy(v, values_i+i, p*sizeof(double));
        /* Symmetric part of the matrix with the diagonal */
        memcpy(a, u, p*sizeof(double));
        memcpy(b, v, p*sizeof(double));
        memcpy(a+p, u, p*sizeof(double));
        memcpy(b+p, v, p*sizeof(double));
        for(m=1; m<p; m++) {
            for(k=0; k<p; k++) {
                double e, f, s,t;
                e = u[k] - a[k+m];
                f = v[k] - b[k+m];
                s = e*xr[k] - f*yr[k];
                t = e*yr[k] + f*xr[k];
                xr[k] = s;
                yr[k] = t;
            }
        }
        /* Symmetric part of the matrix without the diagonal */
        for(j=i+p; j<n_end; j+=CDWBlock){
            /* In this loop, p is always CDWBlock */
            q = (j+CDWBlock<=n_end) ? CDWBlock : ((n_end-i-p) % CDWBlock);
            memcpy(a, values_r+j, q*sizeof(double));
            memcpy(b, values_i+j, q*sizeof(double));
            for(m=0; m<q; m++) {
                for(k=0; k<CDWBlock; k++) {
                    double e, f, s,t;
                    e = u[k]-a[m];
                    f = v[k]-b[m];
                    s = e*xr[k] - f*yr[k];
                    t = e*yr[k] + f*xr[k];
                    xr[k] = s;
                    yr[k] = t;
                }
            }
            memcpy(xc, results_r+j, q*sizeof(double));
            memcpy(yc, results_i+j, q*sizeof(double));
            for(k=0; k<CDWBlock; k++) {
                for(m=0; m<q; m++) {
                    double e, f, s,t;
                    e = a[m]-u[k];
                    f = b[m]-v[k];
                    s = e*xc[m] - f*yc[m];
                    t = e*yc[m] + f*xc[m];
                    xc[m] = s;
                    yc[m] = t;
                }
            }
            memcpy(results_r+j, xc, q*sizeof(double));
            memcpy(results_i+j, yc, q*sizeof(double));
        }
        /* Non symmetric part of the matrix */
        slong start = (n_start != 0) ? 0 : n_end;
        slong end = (n_end == d) ? n_start : d ;
        for(j = start; j < end; j += CDWBlock){
            q = (j+CDWBlock<=end) ? CDWBlock : ((end-start) % CDWBlock);
            memcpy(a, values_r+j, q*sizeof(double));
            memcpy(b, values_i+j, q*sizeof(double));
            for(m=0; m<q; m++) {
                for(k=0; k<p; k++) {
                    double e, f, s,t;
                    e = u[k]-a[m];
                    f = v[k]-b[m];
                    s = e*xr[k] - f*yr[k];
                    t = e*yr[k] + f*xr[k];
                    xr[k] = s;
                    yr[k] = t;
                }
            }
        }
        memcpy(results_r+i, xr, p*sizeof(double));
        memcpy(results_i+i, yr, p*sizeof(double));
    }
}
#undef CDWBlock

#undef Nd
#undef Rbytes
#undef Nr
#undef Rbits

double cd_poly_wdk_update(double* z_r, double* z_i,
                             const double* vp_r, const double* vp_i,
                             const double* wdk_r, const double* wdk_i,
                             slong n_start, slong n_end, double stepsize_bound)
{
    slong i;
    double f_r, f_i, mag_ratio, damping, maxstep = 0;
    for(i=n_start; i<n_end; i++) {
        if((wdk_r[i] != 0) || (wdk_i[i] != 0)) {
            divide(&f_r, &f_i, vp_r[i], vp_i[i], wdk_r[i], wdk_i[i]);
            mag_ratio =  hypot(f_r, f_i) / hypot(z_r[i], z_i[i]);
            damping =  fmin(1, stepsize_bound/hypot(f_r, f_i));
            z_r[i] = z_r[i] - damping * f_r;
            z_i[i] = z_i[i] - damping * f_i;
            maxstep = fmax(maxstep, mag_ratio);
            #if 0
            flint_printf("ratio=%.3e, |f|=%.3e, |w|=%.3e, |f/w|=%.3e, |z|=%.3e\n",
                         mag_ratio, hypot(vp_r[i], vp_i[i]), hypot(wdk_r[i], wdk_i[i]),
                         hypot(f_r, f_i), hypot(z_r[i], z_i[i]));
            #endif
        }
    }
    return maxstep;
}

static void cd_poly_no_overflow_permute(double * pz_r, double * pz_i, const double * z_r, const double * z_i, double lc_r, double lc_i, slong n_piv, slong d)
{
    slong i, j=0, k=n_piv;
    double prod_r=lc_r,
           prod_i=lc_i;
    for(i=0; i<d; i++) {
        double s, t;
        if((j >= n_piv) || (k<d && hypot(prod_r, prod_i) < 1)) {
            pz_r[i] = z_r[k];
            pz_i[i] = z_i[k];
            k++;
        } else {
            pz_r[i] = z_r[j];
            pz_i[i] = z_i[j];
            j++;
        }
        s = prod_r*pz_r[i] - prod_i*pz_i[i];
        t = prod_r*pz_i[i] + prod_i*pz_r[i];
        prod_r = s;
        prod_i = t;
    }
    //flint_printf("%.3e\n", hypot(prod_r, prod_i));
}

static double _refine_roots(intermediate_variables * v, slong n, double stepsize_bound)
{
    slong d, n_piv;
    double lc_r, lc_i,
           maxstep, imaxstep;
    d = n-1;
    n_piv = cd_poly_partition_pivot(v->z_r, v->z_i, d);
    /* Dominant computation time */
    /* Direct case */
    lc_r = v->p_r[n-1];
    lc_i = v->p_i[n-1];
    cd_poly_no_overflow_permute(v->pz_r, v->pz_i, v->z_r, v->z_i, 1, 0, n_piv, d);
    //cd_poly_no_overflow_permute(v->pz_r, v->pz_i, v->z_r, v->z_i, lc_r, lc_i, n_piv, d);
    cd_poly_horner(v->vp_r, v->vp_i, v->z_r, v->z_i, 0, n_piv, v->p_r, v->p_i, n);
    cd_poly_weierstrass_distinct_orders(v->wdk_r, v->wdk_i, lc_r, lc_i, v->pz_r, v->pz_i, d, v->z_r, v->z_i, 0, n_piv);
    maxstep = cd_poly_wdk_update(v->z_r, v->z_i, v->vp_r, v->vp_i, v->wdk_r, v->wdk_i, 0, n_piv, stepsize_bound);
    /* Inverse case */
    //divide(&lc_r, &lc_i, 1, 0, v->p_r[0], v->p_i[0]);
    //cd_poly_no_overflow_permute(v->pz_r, v->pz_i, v->z_r, v->z_i, lc_r, lc_i, n_piv, d);
    _vector_inverse(v->z_r, v->z_i, v->z_r, v->z_i, n_piv, d);
    _vector_inverse(v->pz_r, v->pz_i, v->pz_r, v->pz_i, 0, d);
    lc_r = v->p_r[0];
    lc_i = v->p_i[0];
    cd_poly_horner(v->vp_r, v->vp_i, v->z_r, v->z_i, n_piv, d, v->rp_r, v->rp_i, n);
    cd_poly_weierstrass_distinct_orders(v->wdk_r, v->wdk_i, lc_r, lc_i, v->pz_r, v->pz_i, d, v->z_r, v->z_i, n_piv, d);
    imaxstep = cd_poly_wdk_update(v->z_r, v->z_i, v->vp_r, v->vp_i, v->wdk_r, v->wdk_i, n_piv, d, stepsize_bound);
    _vector_inverse(v->z_r, v->z_i, v->z_r, v->z_i, n_piv, d);
    #if 0
    flint_printf("direct %ld, correction %.3e, inverse %ld, correction %.3e\n", n_piv, maxstep, d-n_piv, imaxstep);
    #endif
    return fmax(maxstep, imaxstep);
}

/* One step refine function */
double cd_poly_refine_roots(double * z, const double * p, slong n, double stepsize_bound)
{
    slong i;
    double lc_r, lc_i, maxstep;
    intermediate_variables v[1];
    _intermediate_variables_init(v, n);
    _set_poly(v, p, n);
    for(i=0; i<n-1; i++) {
        v->z_r[i] = z[2*i];
        v->z_i[i] = z[2*i+1];
    }
    lc_r = v->p_r[n-1];
    lc_i = v->p_i[n-1];
    cd_poly_horner(v->vp_r, v->vp_i, v->z_r, v->z_i, 0, n-1, v->p_r, v->p_i, n);
    cd_poly_weierstrass(v->wdk_r, v->wdk_i, lc_r, lc_i, v->z_r, v->z_i, 0, n-1, n-1);
    maxstep = cd_poly_wdk_update(v->z_r, v->z_i, v->vp_r, v->vp_i, v->wdk_r, v->wdk_i, 0, n-1, stepsize_bound);
    for(i=0; i<n-1; i++) {
        z[2*i]   = v->z_r[i];
        z[2*i+1] = v->z_i[i];
    }
    _intermediate_variables_clear(v);
    return maxstep;
}

/* One step refine function using pivot to avoid overflow.
 * The order of the refined roots is not necessarily preserved */
double cd_poly_refine_roots_with_pivot(double * z, const double * p, slong n, double stepsize_bound)
{
    slong i;
    double maxstep;
    intermediate_variables v[1];
    _intermediate_variables_init(v, n);
    _set_poly(v, p, n);
    for(i=0; i<n-1; i++) {
        v->z_r[i] = z[2*i];
        v->z_i[i] = z[2*i+1];
    }
    maxstep = _refine_roots(v, n, stepsize_bound);
    for(i=0; i<n-1; i++) {
        z[2*i]   = v->z_r[i];
        z[2*i+1] = v->z_i[i];
    }
    _intermediate_variables_clear(v);
    return maxstep;
}

#define MAX_STEPSIZE_BOUND 0x1p-3
#define MIN_STEPSIZE_BOUND 0x1p-10

/* Full resolution function */
double cd_poly_find_roots(double * z, const double * p, const double * z0, slong n, slong num_iter, double reltol)
{
    slong nt, fz, i;
    intermediate_variables v[1];
    double stepsize_bound,
           correction=0,
           old_correction=MAX_STEPSIZE_BOUND;

    nt = _trim_zeros(z, &fz, p, n);
    if(nt > 1) {
        _intermediate_variables_init(v, nt);
        _set_poly(v, p + 2*fz, nt);
        cd_poly_roots_initial_values(v->z_r, v->z_i, v->p_r, v->p_i, nt, z0, n-1);
        /* Main solve function */
        stepsize_bound = MIN_STEPSIZE_BOUND;
        for(i=0; i<num_iter; i++) {
            #if 0
            flint_printf("step %ld: ", i);
            #endif
            correction = _refine_roots(v, nt, stepsize_bound);
            stepsize_bound *= (correction > old_correction)? 0.5 : 2;
            stepsize_bound = fmax(fmin(stepsize_bound, MAX_STEPSIZE_BOUND), MIN_STEPSIZE_BOUND);
            if((nt-2)*correction < reltol || correction < ldexp(1,-50)) {
                break;
            }
            old_correction = correction;
        }
        for(i=0; i<nt-1; i++) {
            z[2*(i+fz)] = v->z_r[i];
            z[2*(i+fz)+1] = v->z_i[i];
        }
        qsort(z, n-1, 2*sizeof(double), _compare);

        _intermediate_variables_clear(v);
    }
    return fmax(correction, 0x1p-50);
}

#undef MIN_STEPSIZE_BOUND
#undef MAX_STEPSIZE_BOUND


