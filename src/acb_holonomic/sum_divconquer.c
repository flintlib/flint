#include "acb_types.h"
#include "acb.h"
#include "acb_poly.h"
#include "acb_holonomic.h"
#include "fmpz.h"
#include "fmpz_vec.h"


static void
acb_poly_shift_series(acb_poly_struct * g, const acb_poly_struct * f,
                      acb_srcptr a, slong n, slong len, slong prec)
{
    acb_poly_fit_length(g, len);

    if (len == 1 && acb_is_zero(a))  /* covers ordinary points */
    {
        acb_ptr c = g->coeffs;
        acb_zero(c);

        for (slong i = acb_poly_degree(f); i >= 0; i--)
        {
            acb_mul_si(c, c, n, prec);
            acb_add(c, c, f->coeffs + i, prec);
        }

        _acb_poly_set_length(g, 1);
        _acb_poly_normalise(g);
    }
    else
    {
        acb_t b;
        acb_init(b);

        acb_add_si(b, a, n, prec);
        acb_poly_taylor_shift(g, f, b, prec);  /* wasteful */
        acb_poly_truncate(g, len);

        acb_clear(b);
    }
}


static slong
max_nlogs(acb_holonomic_sum_context_struct * ctx)
{
    slong nlogs = 0;
    for (slong m = 0; m < ctx->nsols; m++)
        nlogs = FLINT_MAX(nlogs, ctx->sol[m].nlogs);
    return nlogs;
}


static slong
max_nlogs_xn(acb_holonomic_sum_context_struct * ctx, slong off)
{
    for (slong nlogs = max_nlogs(ctx); nlogs > 0; nlogs--)
    {
        for (slong m = 0; m < ctx->nsols; m++)
        {
            if (nlogs <= ctx->sol[m].nlogs)
            {
                acb_poly_struct * p = ctx->sol[m].series + nlogs - 1;
                if (!acb_is_zero(p->coeffs + off))
                    return nlogs;
            }
        }
    }
    return 0;
}


/* XXX better decouple rhs from sol */


static void
apply_diffop(acb_holonomic_sum_context_struct * ctx,
             slong base, slong low, slong mid, slong high)
{
    /* flint_printf("apply_diffop: base=%wd low=%wd mid=%wd high=%wd\n"); */

    if (high - mid >= 7 * ctx->dop_len)  /* threshold from ore_algebra */
    {
        for (slong m = 0; m < ctx->nsols; m++)
            _acb_holonomic_apply_diffop_polmul(
                    ctx->sol[m].series, mid - base,
                    ctx->dop, ctx->dop_len,
                    ctx->expo, low,
                    ctx->sol[m].series, low - base, mid - low,
                    ctx->sol[m].nlogs,
                    mid - low, high - mid,
                    ctx->prec);
    }
    else
    {
        slong nlogs = max_nlogs(ctx);

        slong weights_len = (high - mid) * nlogs * (mid - low);
        acb_ptr weights = _acb_vec_init(weights_len);

        _acb_holonomic_apply_diffop_basecase_weights(
                weights, ctx->dop, ctx->dop_len, ctx->expo,
                low, mid - low, nlogs, mid - low, high - mid,
                ctx->prec);

        for (slong m = 0; m < ctx->nsols; m++)
            _acb_holonomic_apply_diffop_basecase_precomp(
                    ctx->sol[m].series, mid - base,
                    weights, nlogs, ctx->sol[m].series, low - base, mid - low,
                    ctx->sol[m].nlogs, mid - low, high - mid, ctx->prec);

        _acb_vec_clear(weights, weights_len);
    }
}


static void
next_binom_n(fmpz * binom_n, slong n, slong len)
{
    if (n == 0)
    {
        _fmpz_vec_zero(binom_n, len);
        fmpz_one(binom_n);
    }
    else
    {
        for (slong i = len - 1; i > 0; i--)
            fmpz_add(binom_n + i, binom_n + i, binom_n + i - 1);
    }
}


void
_acb_holonomic_sum_forward_1(acb_holonomic_sum_context_struct * ctx,
                             slong base, slong n)
{
    slong ini_i, mult, len;
    acb_poly_t ind_n;

    acb_poly_init(ind_n);

    FLINT_ASSERT(n >= base);

    /* find the position and multiplicity of the current shift as an indicial
     * root */

    mult = 0;
    for (slong i = 0; i < ctx->dop_len - 1; i++)
    {
        /* flint_printf("i=%wd shift=%wd mult=%wd\n", i, ctx->sing_shifts[i].n, ctx->sing_shifts[i].mult); */
        if (ctx->sing_shifts[i].n == n)
        {
            ini_i = i;
            mult = ctx->sing_shifts[i].mult;
            break;
        }
    }

    /* find the maximum log(x)-length of the corresponding coefficient of x on
     * the rhs */

    len = max_nlogs_xn(ctx, n - base);

    /* evaluate the indicial polynomial */

    acb_poly_shift_series(ind_n, ctx->ind, ctx->expo, n, len + mult, ctx->prec);

    next_binom_n(ctx->binom_n, n, ctx->nder);

    /* flint_printf("fwd1: n=%wd mult=%wd len=%wd ind_n=%{acb_poly}\n", n, mult, len, ind_n); */

    /* advance solutions */

    for (slong m = 0; m < ctx->nsols; m++)
    {
        acb_holonomic_sol_struct * sol = ctx->sol + m;
        _acb_holonomic_sol_add_term(
                sol, base, n,
                sol->series, sol->nlogs, mult, ini_i, ind_n,
                ctx->pows, ctx->shifted_sums, ctx->npts, ctx->binom_n, ctx->nder,
                ctx->prec, ctx->sums_prec);
    }

    /* compute the next power of each evaluation point */
    /* XXX wasteful near the end */

    for (slong j = 0; j < ctx->npts; j++)
    {
        acb_mul(ctx->pows + ((n + 1) % ctx->nder) * ctx->npts + j,
                ctx->pows + (n       % ctx->nder) * ctx->npts + j,
                ctx->pts + j,
                ctx->sums_prec);
    }

    acb_poly_clear(ind_n);
}


static void
discard_block(acb_holonomic_sum_context_struct * ctx, slong size)
{
    for (slong m = 0; m < ctx->nsols; m++)
    {
        for (slong k = 0; k < ctx->sol[m].nlogs; k++)
        {
            acb_poly_struct * f = ctx->sol[m].series + k;
            /* NOTE the sage version uses slightly different sizes near
             * max_terms */
            _acb_poly_shift_right(f->coeffs, f->coeffs, 2 * size, size);
            _acb_vec_zero(f->coeffs + size, size);  /* XXX useful? */
        }
    }
}


void
_acb_holonomic_sum_divconquer(acb_holonomic_sum_context_struct * ctx,
                              slong base, slong low, slong high)
{
    slong mid;

    slong bs = ctx->block_size;

    FLINT_ASSERT(base <= low);
    FLINT_ASSERT(low <= high);

    if (high == low)
        return;
    else if (high == low + 1) {
        _acb_holonomic_sum_forward_1(ctx, base, low);
        return;
    }

    if (high <= low + bs)
        mid = (high + low)/2;
    else
        mid = low + bs;

    _acb_holonomic_sum_divconquer(ctx, base, low, mid);

    /*
    acb_poly_t tmp;
    acb_poly_init(tmp);
    flint_printf("divconquer: bs=%wd base=%wd low=%wd mid=%wd high=%wd\n", bs, base, low, mid, high);
    for (slong m = 0; m < ctx->nsols; m++)
    {
        acb_poly_set(tmp, ctx->sol[m].series);
        acb_poly_truncate(tmp, mid - base);
        flint_printf("series[%wd]=%{acb_poly}\n", m, tmp);
    }
    acb_poly_clear(tmp);
    */

    slong resid_len = FLINT_MIN(high - mid, bs);

    apply_diffop(ctx, base, low, mid, mid + resid_len);

    /*
    for (slong m = 0; m < ctx->nsols; m++)
        flint_printf("series+res[%wd]=%{acb_poly}\n", m, ctx->sol[m].series);
    for (slong j = 0; j < ctx->npts; j++)
        for (slong m = 0; m < ctx->nsols; m++)
            flint_printf("sums[%wd][%wd]=%{acb_poly}\n", m, j, ctx->sol[m].sums + j);
    */

    /* TODO convergence check (requires the residual!) */

    /* TODO adjust working precisions */

    if (mid - low >= bs)
    {
        if(ctx->flags & ACB_HOLONOMIC_WANT_SERIES)
        {
            for (slong m = 0; m < ctx->nsols; m++)
                acb_holonomic_sol_fit_length(ctx->sol + m, high + bs, 0);
        }
        else
        {
            discard_block(ctx, bs);
            base += bs;
        }
    }

    /* XXX in WANT_SERIES mode, set the correct lengths before terminating? */

    /* XXX We sometimes do not compute the residual corresponding to the last
     * block. This may need to change when adding support for error bounds. */
    _acb_holonomic_sum_divconquer(ctx, base, mid, high);
}


void
_acb_holonomic_sum_fix(acb_holonomic_sum_context_struct * ctx)
{
    acb_t inv, invpow;
    acb_init(inv);
    acb_init(invpow);

    for (slong j = 0; j < ctx->npts; j++)
    {
        if (!ctx->shifted_sums[j])
            continue;

        acb_inv(inv, ctx->pts + j, ctx->sums_prec);

        acb_one(invpow);
        for (slong i = 1; i < ctx->nder; i++)
        {
            acb_mul(invpow, invpow, inv, ctx->sums_prec);
            for (slong m = 0; m < ctx->nsols; m++)
            {
                for (slong k = 0; k < ctx->sol[m].nlogs; k++)
                {
                    acb_poly_struct * f = acb_holonomic_sol_sum_ptr(
                            ctx->sol + m, j, k);
                    acb_ptr c = f->coeffs + i;
                    acb_mul(c, c, invpow, ctx->sums_prec);
                }
            }
        }
    }

    acb_clear(invpow);
    acb_clear(inv);
}


/* XXX temporary, to be generalized */
void
acb_holonomic_sum_divconquer(
        acb_holonomic_sum_context_struct * ctx,
        slong nterms)
{
    _acb_holonomic_sum_precompute(ctx);
    _acb_holonomic_sum_reset(ctx);
    _acb_holonomic_sum_divconquer(ctx, 0, 0, nterms);
    _acb_holonomic_sum_fix(ctx);
}
