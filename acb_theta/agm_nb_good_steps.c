
#include "acb_theta.h"

slong
acb_theta_agm_nb_good_steps(const arf_t r, const arf_t e, slong prec)
{
    arb_t x;
    arb_t temp;
    arf_t b;
    fmpz_t exp;
    slong lowprec = ACB_THETA_AGM_LOWPREC;
    slong nb;

    arb_init(x);
    arb_init(temp);
    arf_init(b);
    fmpz_init(exp);
    
    /* Solve for r * e^(2^k) * (1+re)/(1-re)^2 <= target */
    arb_one(x);
    arb_mul_2exp_si(x, x, -prec);
    arb_div_arf(x, x, r, lowprec);

    arb_set_arf(temp, e);
    arb_mul_arf(temp, temp, r, lowprec);
    arb_sub_si(temp, temp, 1, lowprec);
    arb_sqr(temp, temp, lowprec);
    arb_mul(x, x, temp, lowprec);
    
    arb_set_arf(temp, e);
    arb_mul_arf(temp, temp, r, lowprec);
    arb_add_si(temp, temp, 1, lowprec);
    arb_div(x, x, temp, lowprec);

    arb_log(x, x, lowprec);
    arb_set_arf(temp, e);
    arb_log(temp, temp, lowprec);
    arb_div(x, x, temp, lowprec);
    arb_get_ubound_arf(b, x, lowprec);
    
    if (!arf_is_finite(b))
    {
        flint_printf("agm_nb_good_steps: Error (infinite value)\n");
        fflush(stdout);
        flint_abort();
    }

    arf_frexp(b, exp, b);
    nb = fmpz_get_si(exp);

    flint_printf("(agm_nb_good_steps) Make %wd good steps\n", nb);
    
    arb_clear(x);
    arb_clear(temp);
    arf_clear(b);
    fmpz_clear(exp);
    return nb;
}
