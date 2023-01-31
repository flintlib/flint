
#include "acb_theta.h"

void
acb_theta_agm_ext_step_good(acb_ptr r, acb_srcptr a, slong g, slong prec)
{
    slong k;

    for (k = 0; k < (1 << (g + 1)); k++)
    {
        acb_sqrt(&r[k], &a[k], prec);
    }
    acb_theta_agm_ext_step_sqrt(r, r, g, prec);
}
