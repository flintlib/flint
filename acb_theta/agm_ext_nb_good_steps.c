
#include "acb_theta.h"

slong
acb_theta_agm_ext_nb_good_steps(const arf_t c, const arf_t r, const arf_t e,
        slong g, slong prec)
{
    arb_t x;
    arb_t temp;
    arf_t b;
    fmpz_t exp;
    slong nb1;
    slong lowprec = ACB_THETA_AGM_LOWPREC;
    
    arb_init(x);
    arb_init(temp);
    arf_init(b);
    fmpz_init(exp);
    
    nb1 = acb_theta_agm_good_steps(r, e, lowprec);
    nb1 = FLINT_MAX(1, nb1);

    /* e must be at most 1/2 */
    arb_set_arf(x, e);
    arb_mul_si(x, x, 2, lowprec);
    arb_sub_si(x, x, 1, lowprec);
    if (!arb_is_negative(x))
    {
        flint_printf("acb_theta_agm_ext_nb_good_steps: Error");
        flint_printf(" (quadratic convergence not reached)\n");
        arf_printd(e, 10); flint_printf("\n");
        fflush(stdout);
        flint_abort();
    }

    /* Solve 4*c*e^(2^(k-1)) <= target */
    arb_one(x);
    arb_mul_2exp_si(x, x, -prec - 2);
    arb_div_arf(x, x, c, lowprec);
    
    arb_log(x, x, lowprec);
    arb_set_arf(temp, e);    
    arb_log(temp, temp, lowprec);
    arb_div(x, x, temp, lowprec);
    arb_get_ubound_arf(b, x, lowprec);
    
    if (!arf_is_finite(b))
    {
        flint_printf("agm_ext_nb_good_steps: Error (infinite value)\n");
        fflush(stdout);
        flint_abort();
    }

    arf_frexp(b, exp, b);
    nb = fmpz_get_si(exp);

    flint_printf("agm_ext_nb_good_steps: Make %wd good steps\n", nb);
    
    arb_clear(x);
    arb_clear(temp);
    arf_clear(b);
    fmpz_clear(exp);
    return nb;    
}
