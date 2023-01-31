
#include "acb_theta.h"

void
acb_theta_agm_conv_rate(arf_t c, arf_t r, const arf_t eps, slong prec)
{
    arb_t eps_arb;
    arb_t temp;
    arb_t res;

    arb_init(temp);
    arb_init(eps_arb);
    arb_init(res);

    arb_set_arf(eps_arb, eps);

    /* Get eta = 1/8 + 1/12*eps^3/(1-eps) */
    arb_sub_si(res, eps_arb, 1, prec);
    arb_pow_ui(temp, eps_arb, 3, prec);
    arb_mul(res, res, temp, prec);
    arb_div_si(res, res, -12, prec);

    arb_one(temp);
    arb_mul_2exp_si(temp, temp, -3);
    arb_add(res, res, temp, prec);

    /* Get res = (1 + eta*(2+eps)^2)/(2(1-eps)) */
    arb_add_si(temp, eps_arb, 2, prec);
    arb_sqr(temp, temp, prec);
    arb_mul(res, res, temp, prec);
    arb_add_si(res, res, 1, prec);

    arb_sub_si(temp, eps_arb, 1, prec);
    arb_mul_si(temp, temp, -2, prec);
    arb_div(res, res, temp, prec);

    /* Replace eps by res*eps, res by 1/res */
    arb_mul(eps_arb, eps_arb, res, prec);
    arb_inv(res, res, prec);

    /* Abort if not eps < 1 */
    arb_set_si(temp, 1);
    if (!arb_lt(eps_arb, temp))
    {
        flint_printf
            ("agm_conv_rate: Error (quadratic convergence not reached)");
        arb_printd(eps_arb, 10);
        flint_printf("\n");
        fflush(stdout);
        flint_abort();
    }

    arb_get_ubound_arf(c, res, prec);
    arb_get_ubound_arf(r, eps_arb, prec);

    arb_clear(temp);
    arb_clear(eps_arb);
    arb_clear(res);
}
