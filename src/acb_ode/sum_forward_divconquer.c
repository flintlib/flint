#include "acb_ode.h"

/* consumes residual coefficients n..high-1, writes the corresponding solution
   coefficients, does not update the residuals */
static void
divconquer_block(acb_ode_sum_struct * sum, slong high)
{
    slong low = sum->n;
    FLINT_ASSERT(low <= high);

    if (high == low)
        return;
    else if (high == low + 1)
    {
        _acb_ode_sum_forward_1(sum);
        return;
    }

    slong mid = (high + low)/2;
    divconquer_block(sum, mid);
    _acb_ode_sum_update_residuals(sum, low, mid, high);
    divconquer_block(sum, high);
}

int
_acb_ode_sum_forward_divconquer(acb_ode_sum_struct * sum, slong high,
                                slong block_len, slong stride,
                                acb_ode_bound_t bound,
                                acb_ode_group_bound_t gbound, slong prec)
{
    FLINT_ASSERT((bound == NULL) == (gbound == NULL));

    while (sum->n < high)
    {
        slong low = sum->n;
        slong mid = FLINT_MIN(low + block_len, high);

        /* at the moment we only test convergence when a full residual is
           naturally available, but in principle it would make sense to do so
           between the two recursive calls in sum_divconquer, or even in/just
           after sum_forward_1 (and that would be useful for operators of high degree) */

        if (high >= low + block_len &&
            /* enough previous terms for the heuristic convergence check */
            low >= sum->n0 + sum->dop_clen - 1)
        {
            if (bound != NULL)
            {
                FLINT_ASSERT(block_len >= sum->dop_clen - 1);
                int done = _acb_ode_sum_done(sum, stride, bound, gbound, prec);
                if (done)
                    return done;
            }

            /* gradually decrease the working precision when it appears wasteful */
            slong prec_wanted = 0;
            for (slong m = 0; m < sum->nsols; m++)
                prec_wanted = FLINT_MAX(prec_wanted, sum->sol[m].cvest->prec_wanted);
            if (prec_wanted < sum->wp)
                sum->wp -= (sum->wp - prec_wanted)/2;
        }

        divconquer_block(sum, mid);

        if (high > mid)
        {
            slong next = FLINT_MIN(mid + block_len, high);
            if(sum->flags & ACB_ODE_WANT_SERIES)
                acb_ode_sum_fit_length(sum, next);
            else
            {
                FLINT_ASSERT(mid - sum->n0 >= block_len);
                acb_ode_sum_discard_head(sum, mid - block_len);
            }

            _acb_ode_sum_update_residuals(sum, low, mid, next);
        }
    }
    return 0;
}

