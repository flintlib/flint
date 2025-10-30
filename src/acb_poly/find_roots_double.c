#include "ulong_extras.h"

typedef struct {
    /* arrays storing the input with real and imaginary parts separated */
    double *z_r;
    double *z_i;
    double *p_r;
    double *p_i;
    /* arrays for the intermediate numerator (vp) and denominator (dk) used
     * in the Durand-Kerner iterator */
    double *vp_r;
    double *vp_i;
    double *dk_r;
    double *dk_i;
    /* arrays for the inverses of points with magnitude >=1 */
    double *iz_r;
    double *iz_i;
    double *Rp_r;
    double *Rp_i;
} double_field;


static void reciprocal(double* Rp_r, double* Rp_i, const double* p_r, const double* p_i, slong n)
{
    for(int i=0; i<n; i++) {
        Rp_r[i] = p_r[n-1-i];
        Rp_i[i] = p_i[n-1-i];
    }
}

/* Remove the leading and trailing zero coeffs of the input polynomial */
static slong trim_zeros(double* z, slong* fz, slong* bz, const double* p, slong n)
{
    int i;
    double bc = INFINITY,
           br = 0.0;

    for(i=0; i<n-1 && p[2*i] == 0.0 && p[2*i+1] == 0.0 && in_rads[i] == 0.0; i++) {
        z[2*i] = 0.0;
        z[2*i+1] = 0.0;
        rad[i] = 0.0;
    }
    *fz = i;

    if(*fz == n-1 && p[2*(n-1)] == 0.0 && p[2*(n-1)+1] == 0.0 && in_rads[n-1] == 0.0) {
        bc = 0;
        br = INFINITY;
        *fz = 1;
    }

    for(i=0; i<n-1 && p[2*(n-1-i)] == 0.0 && p[2*(n-1-i)+1] == 0 && in_rads[n-1-i] == 0.0; i++) {
        z[2*(n-2-i)] = bc;
        z[2*(n-2-i)+1] = 0.0;
        rad[n-2-i] = br;
    }
    *bz = i;
    return n-*fz-*bz;
}
                  

static void _double_poly_init(double_field* v, slong n)
{
    reciprocal(v->Rp_r, v->Rp_i, n, v->p_r, v->p_i, slong n);
}


int double_poly_find_roots_iter()
{
}



