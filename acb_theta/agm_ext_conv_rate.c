
#include "acb_theta.h"

void
acb_theta_agm_ext_conv_rate(arf_t c, arf_t r, arf_t e, acb_srcptr a, slong g,
        slong prec)
{
    arb_t Mu, mu, Ms, ms;
    arb_t temp;
    arb_t res;
    slong n = 1<<g;
    slong k;

    arb_init(Mu);
    arb_init(mu);
    arb_init(Ms);
    arb_init(ms);
    arb_init(temp);
    arb_init(res);
    
    /* Get convergence rate of regular Borchardt */
    acb_theta_agm_conv_rate(r, e, &a[n], g, prec);

    /* Get maxima, minima */
    arb_zero(Mu);
    arb_zero(Ms);
    arb_pos_inf(mu);
    arb_pos_inf(ms);
    for (k = 0; k < n; k++)
    {
        acb_abs(temp, &a[k], prec);
        arb_max(Mu, Mu, temp, prec);
        arb_min(mu, mu, temp, prec);
        acb_abs(temp, &a[k+n], prec);
        arb_max(Ms, Ms, temp, prec);
        arb_min(ms, ms, temp, prec);
    }

    /* Get lambda s.t. |u_0^(n+1) - v_0^n t_0^n| <= x_n:= lambda e^(2^(n-1)) */
    arb_div(res, Mu, ms, prec);
    arb_sqrt(res, res, prec);
    arb_mul_arf(res, res, e, prec);

    arb_div(temp, Ms, mu, prec);
    arb_sqrt(temp, temp, prec);
    arb_add(res, res, temp, prec);
    
    arb_mul_arf(res, res, r, prec);
    arb_mul_2exp_si(res, res, -1);
    arb_mul(res, res, Ms, prec);

    /* Get lambda' s.t. x_n^2 + 2 x_n v_0^(n) t_0^(n) <= lambda' e^(2^(n-1)) */
    arb_mul(temp, Mu, Ms, prec);
    arb_sqrt(temp, temp, prec);
    arb_mul_2exp_si(temp, temp, 1);
    arb_addmul_arf(temp, res, r, prec);
    arb_mul(res, res, temp, prec);

    /* Get lambda'' s.t. |q_{n+1} - 1| <= lambda'' e^(2^(n-1)) */
    arb_set_arf(temp, e);
    arb_sub_si(temp, temp, 1, prec);
    arb_neg(temp, temp);
    arb_inv(temp, temp, prec);
    arb_mul_arf(temp, temp, e, prec);
    arb_mul_arf(temp, temp, r, prec);
    arb_mul(temp, temp, Ms, prec);
    arb_mul(temp, temp, Mu, prec);
    arb_add(res, res, temp, prec);

    arb_mul(temp, mu, ms, prec);
    arb_div(res, res, temp, prec);

    arb_get_ubound_arf(c, res, prec);
    
    arb_clear(Mu);
    arb_clear(mu);
    arb_clear(Ms);
    arb_clear(ms);
    arb_clear(temp);
    arb_clear(res);    
}
