
#include "acb_theta.h"

void acb_theta_agm_rel_dist(arb_t eps, acb_srcptr a, slong nb, slong lowprec,
        slong prec)
{
    arb_t abs;
    slong k;
    
    arb_init(abs);

    acb_theta_agm_abs_dist(eps, a, nb, lowprec, prec);
    acb_abs(abs, &a[0], lowprec);
    arb_div(eps, eps, abs, lowprec);
    
    arb_clear(abs);
}
