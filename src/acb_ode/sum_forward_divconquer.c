#include "acb_types.h"
#include "acb.h"
#include "acb_poly.h"
#include "acb_ode.h"


/* todo: document how solutions are represented */

/* XXX better decouple rhs from sol */


/* reads series[low-n0:mid-n0], writes series[mid-n0:high-n0] */
static void
apply_diffop(acb_ode_sum_struct * sum, slong low, slong mid, slong high)
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
        {
            for (slong n = mid; n < high; n++)
            {
                acb_ptr c = sum->sol[m].series->coeffs + (n - sum->n0);
                mag_zero(arb_radref(acb_realref(c)));
                mag_zero(arb_radref(acb_imagref(c)));
            }
        }
    }
}


static void
discard_block(acb_ode_sum_struct * sum, slong block_len)
{
    for (slong m = 0; m < sum->nsols; m++)
    {
        for (slong k = 0; k < sum->sol[m].nlogs; k++)
        {
            acb_poly_struct * f = sum->sol[m].series + k;
            /* the sage version uses slightly different sizes near max_terms */
            _acb_poly_shift_right(f->coeffs, f->coeffs, 2 * block_len, block_len);
            _acb_vec_zero(f->coeffs + block_len, block_len);  /* XXX useful? */
        }
    }
    sum->n0 += block_len;
}


// assumes all sol have enough space for 2*block_len coefficients
int
_acb_ode_sum_forward_divconquer(acb_ode_sum_struct * sum, slong high,
                                slong block_len, slong stride, slong prec)
{
    slong low = sum->n;
    FLINT_ASSERT(low <= high);

    if (high == low)
        return 0;
    else if (high == low + 1) {
        _acb_ode_sum_forward_1(sum);
        return 0;
    }

    slong mid = high <= low + block_len ? (high + low)/2 : low + block_len;

    int done = _acb_ode_sum_forward_divconquer(sum, mid, block_len, stride, prec);

    /*
       acb_poly_t tmp;
       acb_poly_init(tmp);
       flint_printf("divconquer: block_len=%wd n0=%wd low=%wd mid=%wd high=%wd\n", block_len, n0, low, mid, high);
       for (slong m = 0; m < sum->nsols; m++)
       {
       acb_poly_set(tmp, sum->sol[m].series);
       acb_poly_truncate(tmp, mid - n0);
       flint_printf("series[%wd]=%{acb_poly}\n", m, tmp);
       }
       acb_poly_clear(tmp);
       */

    if (done)
        return 1;

    slong resid_len = FLINT_MIN(high - mid, block_len);
    apply_diffop(sum, low, mid, mid + resid_len);

    /*
       for (slong m = 0; m < sum->nsols; m++)
       flint_printf("series+res[%wd]=%{acb_poly}\n", m, sum->sol[m].series);
       for (slong j = 0; j < sum->npts; j++)
       for (slong m = 0; m < sum->nsols; m++)
       flint_printf("sums[%wd][%wd]=%{acb_poly}\n", m, j, sum->sol[m].sums + j);
       */

    /* Doing the convergence test here allows the call to apply_diffop to be
       shared, but otherwise it would make sense as well to do it in or just
       after sum_forward_1. */

    if (mid % stride == 0)
    {
        /* TODO we currently rely on the choice of block_len and stride to
           ensure that sum_done has access to (1) enough terms for a robust
           quick convergence check and (2) the full residual */
        FLINT_ASSERT(sum->n >= sum->n0 + sum->dop_clen - 1);  /* (1) */
        FLINT_ASSERT(resid_len >= sum->dop_clen - 1);         /* (2) */
        if (acb_ode_sum_done(sum, stride, prec))
            return 1;

        /* gradually decrease the working precision when it appears wasteful */
        slong prec_wanted = 0;
        for (slong m = 0; m < sum->nsols; m++)
            prec_wanted = FLINT_MAX(prec_wanted, sum->sol[m].cvest->prec_wanted);
        if (prec_wanted < sum->wp)
            sum->wp -= (sum->wp - prec_wanted)/2;
    }

    if (mid - low >= block_len)
    {
        /* extend or shift buffer
           XXX TBI once we have a better representation of vectors of
               polynomials */
        if(sum->flags & ACB_ODE_WANT_SERIES)
            /* support returning series coefficients even when the number of
               terms is not known in advance */
            for (slong m = 0; m < sum->nsols; m++)
                acb_ode_sol_fit_length(sum->sol + m, high + block_len);
        else
            discard_block(sum, block_len);
    }

    /* XXX terminal recursion --> loop for -O0? */
    return _acb_ode_sum_forward_divconquer(sum, high, block_len, stride, prec);
}
