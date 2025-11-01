#include "find_roots_double.h"
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
        if(isfinite(y[i])) {
            while(k>1 && is_not_counter_clockwise(result[k-2], y[result[k-2]],
                                                  result[k-1], y[result[k-1]],
                                                  i, y[i])) {
                k--;
            }
            result[k] = i;
            k++;
        }
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
            z_r[p+j] = magz * cos(theta);
            /* epsilon is used to avoid initial values being symmetrical to the
             * real axis (see Section 4 of Aberth's paper) */
            z_i[p+j] = magz * (sin(theta) + epsilon);
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

/* Initialize the members of the double_field */
static void _double_cpoly_init(double_field * v, slong* np, double* mem_v, const double * p, slong fz, slong n)
{
    double pivot;
    slong i, n_np;

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
        v->malp[i] = -log(hypot(v->p_r[i],v->p_i[i]));
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
#define CDWBlock (Nd*Nr/6)

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
void double_cpoly_weierstrass(double* results_r, double* results_i,
                              double lc_r, double lc_i,
                              const double* twice_values_r, const double* twice_values_i,
                              slong n_start, slong n_end, slong d)
{
    int p, k;
    slong i, j;
    for(i=n_start; i < n_end; i+=CDWBlock) {
        double x[CDWBlock] = {0};
        double y[CDWBlock] = {0};
        double u[CDWBlock] = {0};
        double v[CDWBlock] = {0};
        double a[CDWBlock] = {0};
        double b[CDWBlock] = {0};

        slong p = (i+CDWBlock<=n_end) ? CDWBlock : ((n_end-n_start) % CDWBlock);
        for(int j=0; j<CDWBlock; j++){
            x[j] = lc_r;
            y[j] = lc_i;
        }
        memcpy(u, twice_values_r+i, p*sizeof(double));
        memcpy(v, twice_values_i+i, p*sizeof(double));
        /* starts product after diagonal and continues periodically horizontally
         * until the next diagonal \:::\:
         *                         :\:::\ */
        for(j=1; j<d; j++) {
            memcpy(a, twice_values_r+i+j, p*sizeof(double));
            memcpy(b, twice_values_i+i+j, p*sizeof(double));
            for(k=0; k<CDWBlock; k++) {
                #pragma STDC FP_CONTRACT ON
                double q, r, s,t;
                q = u[k]-a[k];
                r = v[k]-b[k];
                s = q*x[k] - r*y[k];
                t = q*y[k] + r*x[k];
                x[k] = s;
                y[k] = t;
            }
        }
        memcpy(results_r+i, x, p*sizeof(double));
        memcpy(results_i+i, y, p*sizeof(double));
    }
}

static void double_cpoly_wdk_update(double* z_r, double* z_i,
                             const double* vp_r, const double* vp_i,
                             const double* wdk_r, const double* wdk_i,
                             slong n_start, slong n_end)
{
    slong i;
    double mag_vp, arg_vp, mag_wdk, arg_wdk, pi = acos(-1);
    for(i=n_start; i<n_end; i++) {
        if(isfinite(vp_r[i]) && isfinite(vp_i[i]) && ((wdk_r[i]) != 0 || (wdk_i[i] != 0))) {
            mag_vp = hypot(vp_r[i], vp_i[i]);
            arg_vp = atan2(vp_i[i], vp_r[i]);
            mag_wdk = hypot(wdk_r[i], wdk_i[i]);
            arg_wdk = atan2(wdk_i[i], wdk_r[i]);
            z_r[i] = z_r[i] - ( mag_vp / mag_wdk * cos(arg_vp - arg_wdk) );
            z_i[i] = z_i[i] - ( mag_vp / mag_wdk * sin(arg_vp - arg_wdk) );
        }
    }
}

static void double_cpoly_refine_roots(double_field * v, slong n)
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
            double_cpoly_refine_roots(v, nt);
        }
        for(i=0; i<nt-1; i++) {
            z[2*(i+fz)] = v->z_r[i];
            z[2*(i+fz)+1] = v->z_i[i];
        }

        flint_free(np);
        flint_free(mem_v);
    }
}



