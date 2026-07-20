#include "acb_types.h"
#include "acb.h"
#include "acb_ode.h"

/* reads series[low-n0:mid-n0], writes series[mid-n0:high-n0] */
void
_acb_ode_sum_update_residuals(acb_ode_sum_struct * sum,
                             slong low, slong mid, slong high)
{
    /* flint_printf("apply_diffop: n0=%wd low=%wd mid=%wd high=%wd\n"); */

    if (high - mid >= 7 * sum->dop_len)  /* threshold from ore_algebra */
    {
        for (slong m = 0; m < sum->nsols; m++)
            _acb_ode_apply_diffop_polmul(
                    sum->sol[m].series, mid - sum->n0,
                    sum->dop, sum->dop_len,
                    sum->group->leader, low,
                    sum->sol[m].series, low - sum->n0, mid - low,
                    sum->sol[m].nlogs,
                    mid - low, high - mid,
                    sum->wp);
    }
    else
    {
        slong nlogs = acb_ode_sum_max_nlogs(sum);

        slong weights_len = (high - mid) * nlogs * (mid - low);
        acb_ptr weights = _acb_vec_init(weights_len);

        _acb_ode_apply_diffop_basecase_weights(
                weights, sum->dop, sum->dop_len, sum->group->leader,
                low, mid - low, nlogs, mid - low, high - mid,
                sum->wp);

        for (slong m = 0; m < sum->nsols; m++)
            _acb_ode_apply_diffop_basecase_precomp(
                    sum->sol[m].series, mid - sum->n0,
                    weights, nlogs, sum->sol[m].series, low - sum->n0, mid - low,
                    sum->sol[m].nlogs, mid - low, high - mid, sum->wp);

        _acb_vec_clear(weights, weights_len);
    }

    /* todo: also do this in rigorous mode but keep the radii for linear error
       analysis */
    if (sum->flags & ACB_ODE_APPROX)
    {
        for (slong m = 0; m < sum->nsols; m++)
            for (slong k = 0; k < sum->sol[m].nlogs; k++)
                for (slong n = mid; n < high; n++)
                {
                    acb_ptr c = sum->sol[m].series[k].coeffs + (n - sum->n0);
                    mag_zero(arb_radref(acb_realref(c)));
                    mag_zero(arb_radref(acb_imagref(c)));
                }
    }
}
