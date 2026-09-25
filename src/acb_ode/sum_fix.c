#include "acb_types.h"
#include "acb.h"
#include "acb_ode.h"


void
_acb_ode_sum_fix(acb_ode_sum_struct * sum)
{
    acb_t inv, invpow;
    acb_init(inv);
    acb_init(invpow);

    for (slong j = 0; j < sum->npts; j++)
    {
        if (!sum->shifted_sums[j])
            continue;

        acb_inv(inv, sum->pts + j, sum->sum_wp);

        acb_one(invpow);
        for (slong i = 1; i < sum->nder; i++)
        {
            acb_mul(invpow, invpow, inv, sum->sum_wp);
            for (slong m = 0; m < sum->nsols; m++)
            {
                for (slong k = 0; k < sum->sol[m].nlogs; k++)
                {
                    acb_ptr c = acb_ode_sol_sum_ptr(sum->sol + m, j, k, i);
                    acb_mul(c, c, invpow, sum->sum_wp);
                }
            }
        }

        sum->shifted_sums[j] = 0;
    }

    /* XXX in WANT_SERIES mode, set series lengths as required */

    acb_clear(invpow);
    acb_clear(inv);
}
