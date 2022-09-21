
#include "acb_theta.h"

void
acb_theta_agm_conv_rate(arf_t r, arf_t e, acb_srcptr a, slong g, slong prec)
{
    acb_t diff;
    arb_t temp;
    arb_t max;
    arb_t eps;
    arb_t res;
    slong n = 1<<g;
    slong k;

    acb_init(diff);
    arb_init(temp);
    arb_init(max);
    arb_init(eps);
    arb_init(res);

    arb_pos_inf(eps);
    arb_zero(max);
    for (k = 0; k < n; k++)
    {
        acb_abs(temp, &a[k], prec);
        arb_max(max, max, temp, prec);
        acb_sub(diff, &a[k], &a[0], prec);
        acb_abs(temp, diff, prec);
        arb_min(eps, eps, temp, prec);
    }
    acb_abs(temp, &a[0], prec);
    arb_div(eps, eps, temp, prec);

    /* Abort if not eps < 1/8, to be safe */
    arb_one(temp);
    arb_div_si(temp, temp, 8, prec);
    if (!arb_lt(eps, temp))
    {
        flint_printf("agm_conv_rate: Error (quadratic convergence not reached)");
        arb_printd(eps, 10); flint_printf("\n");
        fflush(stdout);
        flint_abort();
    }
    
    /* Get eta = 1/8 + 1/12*eps^3/(1-eps) */
    arb_sub_si(res, eps, 1, prec);
    arb_pow_si(temp, eps, 3, prec);
    arb_mul(res, res, temp, prec);
    arb_div_si(res, res, -12, prec);

    arb_one(temp);
    arb_mul_2exp_si(temp, temp, -3);
    arb_add(res, res, temp, prec);

    /* Get res = (1 + eta*(2+eps)^2)/(2(1-eps)) */
    arb_add_si(temp, eps, 2, prec);
    arb_sqr(temp, temp, prec);
    arb_mul(res, res, temp, prec);
    arb_add_si(res, res, 1, prec);

    arb_sub_si(temp, eps, 1, prec);
    arb_mul_si(temp, temp, -2, prec);
    arb_div(res, res, temp, prec);

    /* Replace eps by res*eps, res by M/res */
    arb_mul(eps, eps, res, prec);
    arb_div(res, max, res, prec);

    arb_get_ubound_arf(r, res, prec);
    arb_get_ubound_arf(e, eps, prec);
    
    acb_clear(diff);
    arb_clear(temp);
    arb_clear(max);
    arb_clear(eps);
    arb_clear(res);    
}
