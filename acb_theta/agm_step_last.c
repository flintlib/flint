
#include "acb_theta.h"

void acb_theta_agm_step_last(acb_t r, acb_srcptr a, slong g, slong prec)
{
    slong k;
    slong n = 1<<g;

    acb_zero(r);
    for (k = 0; k < n; k++)
    {
        acb_add(r, r, &a[k], prec);
    }
    acb_mul_2exp_si(r, r, -g);
}
