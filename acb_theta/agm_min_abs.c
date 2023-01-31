
#include "acb_theta.h"

void
acb_theta_agm_min_abs(arb_t min, acb_srcptr a, slong nb, slong prec)
{
    arb_t abs;
    slong k;

    arb_init(abs);

    arb_pos_inf(min);
    for (k = 0; k < nb; k++)
    {
        acb_abs(abs, &a[k], prec);
        arb_min(min, min, abs, prec);
    }

    arb_clear(abs);
}
