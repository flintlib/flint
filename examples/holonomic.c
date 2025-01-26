#include "acb_types.h"
#include "acb_poly.h"
#include "acb_mat.h"
#include "fmpq_types.h"
#include "fmpq_poly.h"
#include "acb_holonomic.h"


void
_acb_holonomic_sum_swap_mat_ordinary(
        acb_mat_t mat,
        const acb_holonomic_sum_context_struct * ctx, slong s)
{
    for (slong j = 0; j < ctx->nsols; j++)
    {
        FLINT_ASSERT(ctx->sol[j].nlogs <= 1);
        if (ctx->sol[j].nlogs < 1)
            continue;
        for (slong i = 0; i < ctx->nder; i++)
            acb_swap(acb_mat_entry(mat, i, j),
                     (ctx->sol[j].sums + s)->coeffs + i);
    }
}

/* TODO
 * - series
 * - singular
 */

void
ordinary(void)
{
    acb_holonomic_sum_context_t ctx;

    /* (x^2 + 1)*Dx^3 + 7*x
     * = (x^2 + 1)*Tx^3 + (-3*x^2 - 3)*Tx^2 + (2*x^2 + 2)*Tx + 7*x^4 */
    acb_holonomic_sum_context_init(ctx, 4, 2, 3, 3);
    acb_poly_set_coeff_si(ctx->dop + 3, 2, 1);
    acb_poly_set_coeff_si(ctx->dop + 3, 0, 1);
    acb_poly_set_coeff_si(ctx->dop + 2, 2, -3);
    acb_poly_set_coeff_si(ctx->dop + 2, 0, -3);
    acb_poly_set_coeff_si(ctx->dop + 1, 2, 2);
    acb_poly_set_coeff_si(ctx->dop + 1, 0, 2);
    acb_poly_set_coeff_si(ctx->dop + 0, 4, 7);

    ctx->flags |= ACB_HOLONOMIC_WANT_SERIES;

    acb_holonomic_sum_ordinary(ctx);
    acb_holonomic_sum_canonical_basis(ctx);

    acb_set_d_d(ctx->pts, 0.25, 0.25);
    acb_set_si(ctx->pts + 1, 0);
    mag_set_d(&((ctx->pts + 1)->real.rad), 1e-8);

    ctx->prec = 64;
    ctx->sums_prec = 64;

    acb_holonomic_sum_divconquer(ctx, 20);

    acb_mat_t mat;
    acb_mat_init(mat, ctx->nder, ctx->nsols);

    for (slong i = 0; i < ctx->npts; i++)
    {
        _acb_holonomic_sum_swap_mat_ordinary(mat, ctx, i);
        acb_mat_printd(mat, 8);
        /* flint_printf("%{acb_mat}\n\n", mat); */
    }

    acb_holonomic_sum_context_clear(ctx);
    acb_mat_clear(mat);
}


void
series(void)
{
    acb_holonomic_sum_context_t ctx;

    acb_holonomic_sum_context_init(ctx, 2, 0, 1, 1);
    acb_poly_set_coeff_si(ctx->dop + 1, 0, 1);
    acb_poly_set_coeff_si(ctx->dop + 0, 1, -1);

    acb_holonomic_sum_ordinary(ctx);
    acb_holonomic_sum_canonical_basis(ctx);
    ctx->prec = 64;
    ctx->flags |= ACB_HOLONOMIC_WANT_SERIES;

    slong len = 5;
    acb_holonomic_sum_divconquer(ctx, len);
    /* Clear the high part reserved for the residual. (In this special case, the
     * high part is zero because block_length = 1.) */
    acb_poly_truncate(ctx->sol[0].series, len);

    flint_printf("%{acb_poly}\n", ctx->sol[0].series);

    acb_holonomic_sum_context_clear(ctx);
}


void
bessel_j0(void)
{
    acb_holonomic_sum_context_t ctx;

    slong dop_order = 2;
    acb_holonomic_sum_context_init(ctx, dop_order + 1, 1, dop_order, dop_order);
    acb_poly_set_coeff_si(ctx->dop + 2, 0, 1);
    acb_poly_set_coeff_si(ctx->dop + 0, 2, 1);

    ctx->sing_shifts[0].n = 0;
    ctx->sing_shifts[0].mult = 2;

    acb_holonomic_sum_canonical_basis(ctx);

    acb_set_d(ctx->pts, 0.25);

    ctx->prec = 64;
    ctx->sums_prec = 64;

    acb_holonomic_sum_divconquer(ctx, 10);

    for (slong m = 0; m < ctx->nsols; m++)
    {
        acb_holonomic_sol_struct * sol = ctx->sol + m;
        /* series in x */
        flint_printf("f%wd =", m);
        for (slong k = 0; k < sol->nlogs; k++)
        {
            flint_printf(" + (%{acb_poly})*log(%{acb} + x)^%wd/%wd",
                         acb_holonomic_sol_sum_ptr(sol, 0, k),
                         ctx->pts, k, k);
        }
        flint_printf("\n");
    }

    acb_holonomic_sum_context_clear(ctx);
}


void
whittaker_m(void)  /* non-integer exponent, no logs */
{
    acb_t kappa, mu2, half;
    acb_init(kappa);
    acb_init(mu2);
    acb_init(half);

    acb_set_si(kappa, 2);
    acb_set_si(mu2, 3);
    acb_set_d(half, .5);

    slong prec = 64;
    slong dop_order = 2;
    slong len = dop_order;

    acb_holonomic_sum_context_t ctx;

    acb_holonomic_sum_context_init(ctx, dop_order + 1, 1, 1, len);

    acb_poly_set_coeff_si(ctx->dop + 2, 0, 4);
    acb_poly_set_coeff_si(ctx->dop + 1, 0, -4);
    acb_poly_set_coeff_si(ctx->dop + 0, 2, -1);
    acb_mul_si((ctx->dop + 0)->coeffs + 1, kappa, 4, prec);
    acb_poly_set_coeff_si(ctx->dop + 0, 0, 1);
    acb_addmul_si((ctx->dop + 0)->coeffs + 0, mu2, -4, prec);

    acb_sqrt(ctx->expo, mu2, prec);
    acb_add(ctx->expo, ctx->expo, half, prec);  /* sub for other expo */

    ctx->sing_shifts[0].n = 0;
    ctx->sing_shifts[0].mult = 1;

    acb_holonomic_sum_canonical_basis(ctx);

    acb_set_d(ctx->pts, 1.4242);

    ctx->prec = prec;
    ctx->sums_prec = prec;

    acb_holonomic_sum_divconquer(ctx, prec);

    flint_printf("(%{acb} + x)^(%{acb}) * (%{acb_poly})\n",
                 ctx->pts, ctx->expo,
                 acb_holonomic_sol_sum_ptr(ctx->sol, 0, 0));

    acb_poly_t tmp;
    acb_poly_init(tmp);
    acb_poly_set_coeff_si(tmp, 1, 1);
    acb_poly_set_coeff_acb(tmp, 0, ctx->pts);
    acb_poly_pow_acb_series(tmp, tmp, ctx->expo, len, prec);
    acb_poly_mullow(tmp, tmp, acb_holonomic_sol_sum_ptr(ctx->sol, 0, 0), len, prec);

    flint_printf("M(%{acb} + x) = %{acb_poly} + O(x^%wd)\n", ctx->pts, tmp, len);

    acb_poly_clear(tmp);

    acb_holonomic_sum_context_clear(ctx);
    acb_clear(kappa);
    acb_clear(mu2);
    acb_clear(half);
}


/* TODO once we have more support for extracting the results: Apéry */


int
main(void)
{
    ordinary();
    series();
    bessel_j0();
    whittaker_m();

    flint_cleanup_master();
}
