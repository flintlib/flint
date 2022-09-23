
#include "acb_theta.h"

void acb_theta_agm_max_abs(arb_t max, acb_srcptr a, slong nb, slong prec)
{
    arb_t abs;
    slong k;

    arb_init(abs);

    arb_zero(max);
    for (k = 0; k < nb; k++)
    {
        acb_abs(abs, &a[k], prec);
        arb_max(max, max, abs, prec);
    }

    arb_clear(abs);
}
