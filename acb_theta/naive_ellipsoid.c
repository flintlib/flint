
#include "acb_theta.h"

void
acb_theta_naive_ellipsoid(acb_theta_eld_t E, arf_t eps, ulong ab, int all,
        int unif, slong ord, acb_srcptr z, const acb_mat_t tau, slong prec)
{  
    arb_t pi, temp;
    arf_t R2, bound;
    slong scl = -1;
    arb_mat_t im;
    arb_mat_t cho;
    arb_mat_t imz;
    arb_t normz;
    arb_ptr offset;
    slong g = acb_mat_nrows(tau);
    slong eld_prec = ACB_THETA_ELD_DEFAULT_PREC;
    int res;
    slong k;

    arb_init(pi);
    arb_init(temp);
    arf_init(R2);
    arf_init(bound);
    arb_mat_init(im, g, g);
    arb_mat_init(cho, g, g);
    arb_mat_init(imz, g, 1);
    arb_init(normz);
    offset = _arb_vec_init(g);
  
    arf_one(eps);
    arf_mul_2exp_si(eps, eps, -prec + ACB_THETA_NAIVE_EPS_2EXP);

    if (all)
    {
        ab = 0;
        scl = -2;
    }

    acb_mat_get_imag(im, tau);
    arb_const_pi(pi, prec);
    arb_mat_scalar_mul_arb(cho, im, pi, prec);
  
    res = arb_mat_cho(cho, cho, eld_prec);
    if (!res)
    {
        eld_prec = prec;
        arb_mat_cho(cho, cho, eld_prec);
    }
    if (!res)
    {
        flint_printf("acb_theta_naive_ellipsoid: Error (imaginary part is not positive definite)\n");
        fflush(stdout);
        flint_abort();
    }
    
    arb_mat_transpose(cho, cho);
    acb_theta_naive_radius(R2, cho, ord, eps, eld_prec);
  
    if (unif) /* any offset less than 1/2 */
    {
        flint_printf("(acb_theta_naive_ellipsoid) Not implemented\n");
        flint_abort();
    }
    
    else /* set offset in terms of z */
    {
        for (k = 0; k < g; k++)
	{
            arb_set(arb_mat_entry(imz, k, 0), acb_imagref(&z[k]));	  
	}
        arb_mat_inv(im, im, eld_prec);
        arb_mat_mul(imz, im, imz, prec);
        arb_mat_mul(imz, cho, imz, prec);
        for (k = 0; k < g; k++)
	{
            arb_set(&offset[k], arb_mat_entry(imz, k, 0));
	}
    }
    
    /* exponential error factor in terms of z */
    for (k = 0; k < g; k++)
    {
        arb_set(arb_mat_entry(imz, k, 0), acb_imagref(&z[k]));	  
    }
    arb_mat_mul(imz, im, imz, prec);
    arb_zero(normz);
    for (k = 0; k < g; k++)
    {
        arb_mul(temp, arb_mat_entry(imz, k, 0), acb_imagref(&z[k]), prec);
        arb_add(normz, normz, temp, prec);
    }
    arb_mul(normz, normz, pi, prec);
    arb_exp(normz, normz, prec);
    arb_get_ubound_arf(bound, normz, prec);
    arf_mul(eps, eps, bound, prec, ARF_RND_CEIL);
    
    /* Fill ellipsoid */
    arb_mat_scalar_mul_2exp_si(cho, cho, scl);
    acb_theta_eld_fill(E, cho, R2, offset, NULL, ab >> g, eld_prec);
    
    arb_clear(pi);
    arb_clear(temp);
    arf_clear(R2);
    arf_clear(bound);
    arb_mat_clear(im);
    arb_mat_clear(cho);  
    arb_clear(normz);
    arb_mat_clear(imz);
    _arb_vec_clear(offset, g);
}
