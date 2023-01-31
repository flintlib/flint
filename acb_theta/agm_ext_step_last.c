
#include "acb_theta.h"

void
acb_theta_agm_ext_step_last(acb_t r, const acb_t s, acb_srcptr a, slong g,
                            slong prec)
{
    slong n = 1 << g;
    slong k;
    acb_t res;
    acb_t temp;

    acb_init(res);
    acb_init(temp);

    for (k = 0; k < n; k++)
    {
        acb_sqrt(temp, &a[k], prec);
        acb_add(res, res, temp, prec);
    }
    acb_mul_2exp_si(res, res, -g);
    acb_sqrt(temp, s, prec);
    acb_div(res, res, temp, prec);

    acb_set(r, res);
    acb_clear(res);
    acb_clear(temp);
}
