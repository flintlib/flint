
#include "acb_theta.h"

/* Number of good steps and final relative error for input which is at
   relative distance at most 1/20th */
/* Therefore relative error after k steps is 10/7 * (7eps/2)^(2^k) * (1+eps)
   for eps=1/20 */

slong acb_theta_agm_nb_good_steps(arf_t rel_err, slong g, slong prec)
{
    arb_t eps, target, t;
    arf_t u;
    fmpz_t exp;
    slong nb;
  
    slong lowprec = ACB_THETA_AGM_LOWPREC;

    arb_init(eps);
    arb_init(target);
    arb_init(t);
    arf_init(u);
    fmpz_init(exp);

    arb_one(eps);
    arb_div_si(eps, eps, 20, lowprec);
    arb_one(target);
    arb_mul_2exp_si(target, target, -prec);

    arb_add_si(t, eps, 1, lowprec);
    arb_div(target, target, t, lowprec);
    arb_set_si(t, 10);
    arb_div_si(t, t, 7, lowprec);
    arb_div(target, target, t, lowprec);

    arb_mul_si(eps, eps, 7, lowprec);
    arb_mul_2exp_si(eps, eps, -1);

    /* Now solve for eps^(2^k) <= target */
    arb_log(target, target, lowprec);
    arb_log(t, eps, lowprec);
    arb_div(target, target, t, lowprec);
    arb_get_ubound_arf(u, target, lowprec);
    
    if (!arf_is_finite(u))
    {
        flint_printf("agm_nb_good_steps: Error (infinite value)\n");
        fflush(stdout);
        flint_abort();
    }

    arf_frexp(u, exp, u);
    arf_one(rel_err);
    arf_mul_2exp_si(rel_err, rel_err, -prec);
    nb = fmpz_get_si(exp);

    arb_clear(eps);
    arb_clear(target);
    arb_clear(t);
    arf_clear(u);
    fmpz_clear(exp);
    return nb;
}
