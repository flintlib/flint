#include "acb_types.h"
#include "acb.h"
#include "acb_poly.h"
#include "acb_holonomic.h"
#include "fmpz_vec.h"


typedef acb_holonomic_sum_context_struct * ctx_ptr;


void
acb_holonomic_sum_context_init(ctx_ptr ctx, slong dop_len, slong npts,
                               slong nsols, slong nder)
{
    slong dop_order = dop_len - 1;

    ctx->dop = _acb_poly_vec_init(dop_len);
    ctx->dop_len = dop_len;

    acb_init(ctx->expo);

    acb_poly_init(ctx->ind);

    ctx->sing_shifts = flint_calloc(2 * dop_order, sizeof(slong));
    for (slong s = 0; s < dop_order; s++)
        ctx->sing_shifts[s].n = -1;

    ctx->sol = flint_malloc(nsols * sizeof(acb_holonomic_sol_struct));
    /* using dop_order as a bound
     * - for max possible log prec (will need updating to support inhomogeneous
     *   equations),
     * - for number of initial value positions */
    for (slong m = 0; m < nsols; m++)
        acb_holonomic_sol_init(ctx->sol + m, dop_order, dop_order, npts);
    ctx->nsols = nsols;

    ctx->pts = _acb_vec_init(npts);
    ctx->pows = _acb_vec_init(npts * nder);
    ctx->shifted_sums = flint_malloc(npts);
    ctx->npts = npts;
    ctx->nder = nder;

    ctx->bounds_prec = MAG_BITS;  /* default value */

    ctx->binom_n = _fmpz_vec_init(nder);

    ctx->flags = 0;
    ctx->have_precomputed = 0;
}


void
acb_holonomic_sum_context_clear(ctx_ptr ctx)
{
    _acb_poly_vec_clear(ctx->dop, ctx->dop_len);

    acb_clear(ctx->expo);
    acb_poly_clear(ctx->ind);

    flint_free(ctx->sing_shifts);

    for (slong m = 0; m < ctx->nsols; m++)
        acb_holonomic_sol_clear(ctx->sol + m);
    flint_free(ctx->sol);

    _acb_vec_clear(ctx->pts, ctx->npts);
    _acb_vec_clear(ctx->pows, ctx->npts * ctx->nder);
    flint_free(ctx->shifted_sums);

    _fmpz_vec_clear(ctx->binom_n, ctx->nder);
}


void
acb_holonomic_sum_ordinary(ctx_ptr ctx)
{
    for (slong n = 0; n < ctx->dop_len - 1; n++)
    {
        ctx->sing_shifts[n].n = n;
        ctx->sing_shifts[n].mult = 1;
    }
}


void
acb_holonomic_sum_mum(ctx_ptr ctx)
{
    if (ctx->dop_len <= 1)
        return;
    ctx->sing_shifts[0].n = 0;
    ctx->sing_shifts[0].n = ctx->dop_len - 1;
}


/* XXX in the singular case, it is probably more useful to take the solutions of
 * highest level only... */
void
acb_holonomic_sum_canonical_basis(ctx_ptr ctx)
{
    for (int m = 0; m < ctx->nsols; m++)
        acb_holonomic_sol_unit_ini(ctx->sol + m, m, ctx->sing_shifts);
}


void
_acb_holonomic_sum_precompute(ctx_ptr ctx)
{
    if (ctx->have_precomputed)
        return;

    slong dop_len = ctx->dop_len;

    ctx->dop_clen = 0;
    for (slong i = 0; i < dop_len; i++)
        ctx->dop_clen = FLINT_MAX(ctx->dop_clen, ctx->dop[i].length);

    ctx->block_size = FLINT_MAX(1, ctx->dop_clen - 1);

    acb_poly_fit_length(ctx->ind, dop_len);
    _acb_poly_set_length(ctx->ind, dop_len);
    for (slong i = 0; i < dop_len; i++)
        acb_poly_get_coeff_acb((ctx->ind)->coeffs + i, ctx->dop + i, 0);
    _acb_poly_normalise(ctx->ind);

    /* We can trade some full-width muls for muls by small integers by
     * multiplying everyone by ξ^n), but this optimization is only
     * valid when ξ does not contain 0. */
    for (slong j = 0; j < ctx->npts; j++)  /* XXX move? */
        ctx->shifted_sums[j] = !acb_contains_zero(ctx->pts + j);

    /* TODO error bounds etc. */

    ctx->have_precomputed = 1;
}


void
_acb_holonomic_sum_reset(ctx_ptr ctx)
{
    for (slong m = 0; m < ctx->nsols; m++)
    {
        acb_holonomic_sol_reset(ctx->sol + m);
        acb_holonomic_sol_fit_length(ctx->sol + m, 2 * ctx->block_size,
                                     ctx->nder);
    }

    for (slong i = 0; i < ctx->npts; i++)
        acb_one(ctx->pows + i);
}
