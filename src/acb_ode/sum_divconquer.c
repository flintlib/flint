#include "acb.h"
#include "acb_ode.h"

static void
acb_ode_sum_reset(acb_ode_sum_t sum)
{
    sum->n = 0;
    sum->n0 = 0;

    for (slong m = 0; m < sum->nsols; m++)
    {
        acb_ode_sol_zero(sum->sol + m);

        sum->sol[m].cvest->accuracy = sum->wp;
        sum->sol[m].cvest->loss_rate = sum->dop_len * sum->dop_clen;
    }

    for (slong i = 0; i < sum->npts; i++)
        acb_one(sum->pows + i);

    mag_one(sum->magpow);
}

/* XXX make it possible to customize the working precision--how exactly? */
void
acb_ode_sum_divconquer(acb_ode_sum_t sum, slong nterms, slong prec)
{
    if (nterms < 0)
        nterms = WORD_MAX;

    acb_ode_sum_precompute(sum);
    acb_ode_sum_reset(sum);

    /* todo: take into account nterms and the APPROX flag */
    sum->sum_wp = acb_ode_choose_prec(&sum->wp, sum->dop, sum->dop_len,
                                         sum->mag, sum->cvrad, prec);

    slong block_len = FLINT_MAX(1, sum->dop_clen - 1);
    slong stride = block_len;  /* todo: less when block_len is large */

    for (slong m = 0; m < sum->nsols; m++)
        acb_ode_sol_fit_length(sum->sol + m, 2 * block_len);

    _acb_ode_sum_forward_divconquer(sum, nterms, block_len, stride, prec);
    _acb_ode_sum_fix(sum);
}
