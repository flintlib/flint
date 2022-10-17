
#include "acb_theta.h"

slong acb_theta_balance_lowprec(slong g, slong prec)
{
    arb_t x;
    slong lowprec;
    
    arb_init(x);
    
    arb_set_si(x, prec);
    arb_root_ui(x, x, g+1, prec);
    arb_mul_si(x, x, ACB_THETA_BALANCE_LOWPREC_MUL, prec);
    lowprec = arf_get_si(arb_midref(x), ARF_RND_CEIL);

    arb_clear(x);
    return lowprec;
}
