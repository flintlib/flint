#include "find_roots_double.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h> /* for memcpy */

/* Remove the leading and trailing zero coeffs of the input polynomial */
static slong trim_zeros(double * z, slong * fz, const double * p, slong n)
{
    slong i, j;

    /* trailing terms with zero coefficients add roots at infinity*/
    for(i=n-1; i>=1 && p[2*i] == 0.0 && p[2*i+1] == 0.0; i--) {
        z[2*(i-1)] = INFINITY;
        z[2*(i-1)+1] = 0.0;
    }

    /* leading terms with zero coefficients add roots at zero*/
    for(j=0; j<i && p[2*j] == 0.0 && p[2*j+1] == 0.0; j++) {
        z[2*j] = 0.0;
        z[2*j+1] = 0.0;
    }
    *fz = j;

    return i-j+1;
}
                  
static int is_not_counter_clockwise(double a, double b, double x, double y, double u, double v)
{
    return v*(x-a) + y*(a-u) + b*(u-x) <= 0;
}

/* Computes the lower convex hull of the points (i, y[i]) */
static slong lower_convex_hull(slong* result, const double* y, slong n)
{
    slong i, k=0;
    for (i=0; i<n; i++) {
        while(k>1 && is_not_counter_clockwise(result[k-2], y[result[k-2]],
                                              result[k-1], y[result[k-1]],
                                              i, y[i])) {
            k--;
        }
        result[k] = i;
        k++;
    }
    return k;
}

static void initial_values(double* z_r, double* z_i, const double* p_r, const double* p_i,
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
        magz = pow(magp/magq, invd);
        if(magz == 0.0) {
            magz = ldexp(1, -500/(q-p));
        }
        for(j=0; j<q-p; j++) {
            /* theta = arg of solutions of x^d = -p[p]/p[q] 
             *       = arg(-p[p]/p[q])/d + 2 pi j/d */
            theta = (2*pi*(j+0.5)+argp-argq)/d;
            /* epsilon is used to avoid initial values being symmetrical to the
             * real axis (see Section 4 of Aberth's paper) */
            z_r[p+j] = magz * (cos(theta)+(1+j)*epsilon) ;
            z_i[p+j] = magz * (sin(theta)+(1+j)*epsilon);
        }
        p = q;
    }

}

static void d_swap(double* a, double* b)
{
    double t = *a;
    *a = *b;
    *b = t;
}

static void vector_inverse(double* iz_r, double* iz_i, const double* z_r, const double* z_i,
                           slong n_start, slong n_end)
{
    for(slong i=n_start; i<n_end; i++) {
        double s, t;
        /* these formula reduces the risk of overflow */
        s = (z_r[i] == 0) ? 0 : 1/(z_r[i] + z_i[i]*(z_i[i]/z_r[i]));
        t = (z_i[i] == 0) ? 0 : -1/(z_r[i]*(z_r[i]/z_i[i]) + z_i[i]); 
        s = (s==0 && t==0) ? ldexp(1,500) : s;
        iz_r[i] = s;
        iz_i[i] = t;
    }
}

/* This field contains 13 array pointers, and 4 of them are twice larger
 * than the others */
typedef struct {
    /* arrays storing the input with real and imaginary parts separated */
    double *z_r;
    double *z_i;
    double *p_r;
    double *p_i;
    /* for the newton polygon, -log(|coeff|) of p */
    double *malp;
    /* arrays for the intermediate numerator (vp) and denominator (dk) used
     * in the Weierstrass-Durand-Kerner iterator */
    double *vp_r;
    double *vp_i;
    double *wdk_r;
    double *wdk_i;
    /* arrays for the inverses of points with magnitude >=1 */
    double *iz_r;
    double *iz_i;
    double *rp_r;
    double *rp_i;
} double_field;

#define LOG_ROUNDING 3

/* Initialize the members of the double_field */
static void _double_cpoly_init(double_field * v, slong* np, double* mem_v, const double * p, slong fz, slong n)
{
    double pivot, epsilon=ldexp(1,-1000);
    slong i, n_np;
    int k=LOG_ROUNDING;

    /* Set the memory addresses */
    /* z_r and z_i will hold 2(n-1) double each */
    v->z_r   = mem_v+0*n; 
    v->z_i   = mem_v+2*n;
    v->p_r   = mem_v+4*n;
    v->p_i   = mem_v+5*n;
    v->malp  = mem_v+6*n;
    v->vp_r  = mem_v+7*n;
    v->vp_i  = mem_v+8*n;
    v->wdk_r  = mem_v+9*n;
    v->wdk_i  = mem_v+10*n;
    /* iz_r and iz_i will hold 2(n-1) double each */
    v->iz_r  = mem_v+11*n;
    v->iz_i  = mem_v+13*n;
    v->rp_r  = mem_v+15*n;
    v->rp_i  = mem_v+16*n;

    /* Unzip the real, imaginary parts and numerical valuation of the input polynomial coefficients */
    for(i=0; i<n; i++) {
        v->p_r[i] = p[2*(fz+i)];
        v->p_i[i] = p[2*(fz+i)+1];
        v->rp_r[i] = p[2*(fz+n-1-i)];   
        v->rp_i[i] = p[2*(fz+n-1-i)+1]; 
        /* We avoid infinite values so that the code is compatible with
         * unsafe math optimizations */
        v->malp[i] = -log(hypot(v->p_r[i],v->p_i[i])+epsilon);
        /* rounding the valuation helps avoiding nearby initial radii */
        v->malp[i] = ldexp(round(ldexp(v->malp[i],k)),-k);
    }

    /* Computes the Newton polygon as indices in np of length n_np */
    n_np = lower_convex_hull(np, v->malp, n);
    
    /* Computes a first estimation of the roots */
    initial_values(v->z_r, v->z_i, v->p_r, v->p_i, np, n_np);
}                         

/* Reorder the z : the values at index < n-c satisfy |z| <= 1, and at index >= n-c they satisfy |z|>1 */
slong double_cpoly_partition_pivot(double* z_r, double* z_i, slong n)
{
    slong i, c = 0;
    for(i=0; i<n; i++) {
        if(z_r[i]*z_r[i] + z_i[i]*z_i[i] > 1) {
            c++;
        } else if(c>0) {
            d_swap(z_r + i-c, z_r + i);
            d_swap(z_i + i-c, z_i + i);
        }
    }             
    return n-c;
}

static int compare(const void * first, const void * second)
{
    double magl, argl, magr, argr, compmag, comparg,epsilon=ldexp(1,-25);
    double left[2];
    double right[2];
    int res;
    memcpy(left, first, 2*sizeof(double));
    memcpy(right, second, 2*sizeof(double));
    magl = hypot(left[0], left[1]);
    magr = hypot(right[0], right[1]);
    compmag = magl/magr - 1;
    if(compmag < -epsilon){
        res = -1;
    } else if(compmag > epsilon) {
        res = 1;
    } else {
        argl = atan2(left[1], left[0]);
        argr = atan2(right[1], right[0]);
        /* returns the -1, 0 or 1 whether argl is lower, equal or greater
         * than argr */
        res = (argr < argl) - (argl < argr);
    }
    return res;
}


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
#define DBlock (Nd*Nr/2)
#define CSBlock (Nd*Nr/2)
#define CDHBlock (Nd*Nr/4)
#define CDWBlock (Nd*Nr)

/* Horner evaluation */
void double_cpoly_horner(double* results_r, double* results_i,
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
                #pragma STDC FP_CONTRACT ON
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

/* Weights for the Durand-Kerner or Weierstrass iteration */
/* twice_values_r and twice_values_i have size 2n,
 * the second half of each array is a copy of the first half */
void double_cpoly_weierstrass_a(double* restrict results_r, double* restrict results_i,
                              double lc_r, double lc_i,
                              const double* twice_values_r, const double* twice_values_i,
                              slong n_start, slong n_end, slong d)
{
    int p, k;
    slong i, j;
    for(i=n_start; i<n_end; i++){
        results_r[i] = lc_r;
        results_i[i] = lc_i;
    }
    for(j=1; j<d; j++) {
        for(i=n_start; i < n_end; i++) {
            double q, r, s, t;
            q = twice_values_r[i] - twice_values_r[i+j];
            r = twice_values_i[i] - twice_values_i[i+j];
            #pragma STDC FP_CONTRACT ON
            s = q*results_r[i] - r*results_i[i];
            t = q*results_i[i] + r*results_r[i];
            results_r[i] = s;
            results_i[i] = t;
        }
    }
}

/* Weights for the Durand-Kerner or Weierstrass iteration
 * twice_values_r and twice_values_i have size 2n,
 * the second half of each array is a copy of the first half
 * A variant of this version could be usefull for the Borsh Supan weights*/                                                         
void double_cpoly_weierstrass_b(double* restrict results_r, double* restrict results_i,
                              double lc_r, double lc_i,
                              const double* twice_values_r, const double* twice_values_i,
                              slong n_start, slong n_end, slong d)
{
    int p, k;
    slong i, j;
    for(i=n_start; i<n_end; i++){
        results_r[i] = lc_r;
        results_i[i] = lc_i;
    }
    for(j=1; j<n_end-n_start; j++) {
        for(i=n_start; i < n_end-j; i++) {
            double q, r, s, t;
            q = twice_values_r[i] - twice_values_r[i+j];
            r = twice_values_i[i] - twice_values_i[i+j];
            #pragma STDC FP_CONTRACT ON
            s = q*results_r[i] - r*results_i[i];
            t = q*results_i[i] + r*results_r[i];
            results_r[i] = s;
            results_i[i] = t;
            s = -q*results_r[i+j] + r*results_i[i+j];
            t = -q*results_i[i+j] - r*results_r[i+j];
            results_r[i+j] = s;
            results_i[i+j] = t;
        }
    }
    for(j=n_end % d; ((j<n_start) || (j>= n_end)) && (n_start < n_end); j = (j + 1) % d) {
        for(i=n_start; i < n_end; i++) {
            double q, r, s, t;
            q = twice_values_r[i] - twice_values_r[j];
            r = twice_values_i[i] - twice_values_i[j];
            #pragma STDC FP_CONTRACT ON
            s = q*results_r[i] - r*results_i[i];
            t = q*results_i[i] + r*results_r[i];
            results_r[i] = s;
            results_i[i] = t;
        }
    }
}

/* Weights for the Durand-Kerner or Weierstrass iteration */
/* twice_values_r and twice_values_i have size 2n,
 * the second half of each array is a copy of the first half */
void double_cpoly_weierstrass(double* restrict results_r, double* restrict results_i,
                              double lc_r, double lc_i,
                              const double* twice_values_r, const double* twice_values_i,
                              slong n_start, slong n_end, slong d)
{
    int p, q, k, m;
    slong i, j;
    for(i=n_start; i<n_end; i++){
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
        memcpy(u, twice_values_r+i, p*sizeof(double));
        memcpy(v, twice_values_i+i, p*sizeof(double));
        /* Symmetric part of the matrix with the diagonal */
        memcpy(a, u, p*sizeof(double));
        memcpy(b, v, p*sizeof(double));
        memcpy(a+p, u, p*sizeof(double));
        memcpy(b+p, v, p*sizeof(double));
        for(m=1; m<p; m++) {
            for(k=0; k<p; k++) {
                #pragma STDC FP_CONTRACT ON
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
            memcpy(a, twice_values_r+j, q*sizeof(double));
            memcpy(b, twice_values_i+j, q*sizeof(double));
            for(m=0; m<q; m++) {
                for(k=0; k<CDWBlock; k++) {
                    #pragma STDC FP_CONTRACT ON
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
                    #pragma STDC FP_CONTRACT ON
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
            memcpy(a, twice_values_r+j, q*sizeof(double));
            memcpy(b, twice_values_i+j, q*sizeof(double));
            for(m=0; m<q; m++) {
                for(k=0; k<p; k++) {
                    #pragma STDC FP_CONTRACT ON
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

static void double_cpoly_wdk_update(double* z_r, double* z_i,
                             const double* vp_r, const double* vp_i,
                             const double* wdk_r, const double* wdk_i,
                             slong n_start, slong n_end)
{
    slong i;
    double mag_vp, arg_vp, mag_wdk, arg_wdk, mag_z, ratio,
           pi = acos(-1);
    for(i=n_start; i<n_end; i++) {
        if((wdk_r[i] != 0) || (wdk_i[i] != 0)) {
            mag_vp = hypot(vp_r[i], vp_i[i]);
            arg_vp = atan2(vp_i[i], vp_r[i]);
            mag_wdk = hypot(wdk_r[i], wdk_i[i]);
            arg_wdk = atan2(wdk_i[i], wdk_r[i]);
            mag_z = hypot(z_r[i], z_i[i]);
            /* The min expression avoids spurious large step */
            ratio = fmin(mag_vp / mag_wdk, 0.5*mag_z);
            z_r[i] = z_r[i] - ( ratio * cos(arg_vp - arg_wdk) );
            z_i[i] = z_i[i] - ( ratio * sin(arg_vp - arg_wdk) );
        }
    }
}

static void _double_cpoly_refine_roots(double_field * v, slong n)
{
    slong d, n_piv;
    double lc_r, lc_i;
    d = n-1;
    n_piv = double_cpoly_partition_pivot(v->z_r, v->z_i, d);
    vector_inverse(v->iz_r, v->iz_i, v->z_r, v->z_i, 0, d);
    /* the vectors are repeated for the function `double_cpoly_weierstrass` */
    memcpy(v->z_r+d, v->z_r, d*sizeof(double));
    memcpy(v->z_i+d, v->z_i, d*sizeof(double));
    memcpy(v->iz_r+d, v->iz_r, d*sizeof(double));
    memcpy(v->iz_i+d, v->iz_i, d*sizeof(double));
    /* Dominant computation time */
    /* Direct case */
    lc_r = v->p_r[n-1];
    lc_i = v->p_i[n-1];
    double_cpoly_horner(v->vp_r, v->vp_i, v->z_r, v->z_i, 0, n_piv, v->p_r, v->p_i, n);
    double_cpoly_weierstrass(v->wdk_r, v->wdk_i, lc_r, lc_i, v->z_r, v->z_i, 0, n_piv, d);
    double_cpoly_wdk_update(v->z_r, v->z_i, v->vp_r, v->vp_i, v->wdk_r, v->wdk_i, 0, n_piv);
    /* Inverse case */
    lc_r = v->p_r[0];
    lc_i = v->p_i[0];
    double_cpoly_horner(v->vp_r, v->vp_i, v->iz_r, v->iz_i, n_piv, d, v->rp_r, v->rp_i, n);
    double_cpoly_weierstrass(v->wdk_r, v->wdk_i, lc_r, lc_i, v->iz_r, v->iz_i, n_piv, d, d);
    double_cpoly_wdk_update(v->iz_r, v->iz_i, v->vp_r, v->vp_i, v->wdk_r, v->wdk_i, n_piv, d);
    vector_inverse(v->z_r, v->z_i, v->iz_r, v->iz_i, n_piv, d);
}

void double_cpoly_find_roots(double * z, const double * p, slong n, slong max_iter, int verbose)
{
    slong nt, fz, i, status;
    double_field v[1];
    double * mem_v;
    slong * np;

    nt = trim_zeros(z, &fz, p, n); 
    if(nt > 1) {
        mem_v = flint_malloc(17*nt*sizeof(double));
        np = flint_malloc(nt*sizeof(slong));
        _double_cpoly_init(v, np, mem_v, p, fz, nt);
        /* Main solve function */
        for(i=0; i<max_iter; i++) {
            _double_cpoly_refine_roots(v, nt);
        }
        for(i=0; i<nt-1; i++) {
            z[2*(i+fz)] = v->z_r[i];
            z[2*(i+fz)+1] = v->z_i[i];
        }
        qsort(z, n-1, 2*sizeof(double), compare);

        flint_free(np);
        flint_free(mem_v);
    }
}



