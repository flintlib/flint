
#include "acb_theta.h"

void
acb_theta_agm_ext_conv_rate(arf_t c1, arf_t c2, arf_t r, const arf_t eps,
                            const arf_t m, const arf_t M, slong prec)
{
    arb_t M_arb, m_arb;
    arb_t temp;
    arb_t res;

    arb_init(M_arb);
    arb_init(m_arb);
    arb_init(temp);
    arb_init(res);

    arb_set_arf(M_arb, M);
    arb_set_arf(m_arb, m);

    /* Get convergence rate of regular Borchardt */
    acb_theta_agm_conv_rate(c1, r, eps, prec);

    /* Get lambda s.t. |u_0^(n+1) - v_0^n t_0^n| <= x_n:= lambda e^(2^(n-1)) */
    arb_div(res, M_arb, m_arb, prec);
    arb_sqrt(res, res, prec);
    arb_one(temp);
    arb_add_arf(temp, temp, r, prec);
    arb_mul(res, res, temp, prec);

    arb_mul_arf(res, res, c1, prec);
    arb_mul_2exp_si(res, res, -1);
    arb_mul(res, res, M_arb, prec);

    /* Get lambda' s.t. x_n^2 + 2 x_n v_0^(n) t_0^(n) <= lambda' e^(2^(n-1)) */
    arb_mul_2exp_si(temp, M_arb, 1);
    arb_addmul_arf(temp, res, c1, prec);
    arb_mul(res, res, temp, prec);

    /* Get lambda'' s.t. |q_{n+1} - 1| <= lambda'' e^(2^(n-1)) */
    arb_set_arf(temp, r);
    arb_sub_si(temp, temp, 1, prec);
    arb_neg(temp, temp);
    arb_inv(temp, temp, prec);
    arb_mul_arf(temp, temp, r, prec);
    arb_mul_arf(temp, temp, c1, prec);
    arb_mul(temp, temp, M_arb, prec);
    arb_mul(temp, temp, M_arb, prec);
    arb_add(res, res, temp, prec);

    arb_sqr(temp, m_arb, prec);
    arb_div(res, res, temp, prec);

    arb_get_ubound_arf(c2, res, prec);

    arb_clear(M_arb);
    arb_clear(m_arb);
    arb_clear(temp);
    arb_clear(res);
}
