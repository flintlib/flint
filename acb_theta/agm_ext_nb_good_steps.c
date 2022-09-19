
#include "acb_theta.h"

slong acb_theta_agm_ext_nb_good_steps(arf_t rel_err, slong g, slong prec)
{
    arb_t B;
    arf_t bound;
    slong n;
    slong lowprec = ACB_THETA_AGM_LOWPREC;

    arb_init(B);
    arf_init(bound);
    
    arb_set_si(B, 21);
    arb_div_si(B, B, 19, lowprec); /* cf nb_bad_steps */
    arb_pow_ui(B, B, 3, lowprec);
    arb_mul_si(B, B, 5, lowprec);
    arb_log_base_ui(B, B, 2, lowprec);

    arb_add_si(B, B, prec + 1, lowprec);
    arb_log_base_ui(B, B, 2, lowprec);
    arb_add_si(B, B, 2, lowprec);

    arb_get_ubound_arf(bound, B, lowprec);

    if (!arf_is_finite(bound) || arf_cmp_si(bound, WORD_MAX) > 0)
    {
        flint_printf("agm_ext_nb_good_steps: Error (cannot convert to integer)\n");
        arf_printd(bound, 30); flint_printf("\n");
    }

    n = arf_get_si(bound, ARF_RND_CEIL);
    arf_one(rel_err);
    arf_mul_2exp_si(rel_err, rel_err, -prec);

    arb_clear(B);
    arf_clear(bound);
    return n;
}
