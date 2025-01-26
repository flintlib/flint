#include "test_helpers.h"
#include "acb_types.h"
#include "acb_poly.h"
#include "fmpq_types.h"
#include "fmpq_poly.h"
#include "acb_holonomic.h"

TEST_FUNCTION_START(acb_holonomic_sum_divconquer, state)
{

    for (slong iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        /* classical functions -- well for now just atan */

        acb_holonomic_sum_context_t ctx;

        slong dop_order = 2;
        acb_holonomic_sum_context_init(ctx, dop_order + 1,
                                       n_randint(state, 5),
                                       dop_order,
                                       1 + n_randint(state, 3));

        acb_poly_set_coeff_si(ctx->dop + 2, 2, 1);
        acb_poly_set_coeff_si(ctx->dop + 2, 0, 1);
        acb_poly_set_coeff_si(ctx->dop + 1, 2, 1);
        acb_poly_set_coeff_si(ctx->dop + 1, 0, -1);

        acb_holonomic_sum_ordinary(ctx);
        acb_holonomic_sum_canonical_basis(ctx);

        ctx->prec = n_randint(state, 1000);
        ctx->sums_prec = ctx->prec;

        for (slong i = 0; i < ctx->npts; i++)
        {
            acb_ptr a = ctx->pts + i;
            acb_randtest_precise(a, state, ctx->prec, 0);
            acb_div_si(a, a, 4, ctx->prec);
        }

        ctx->flags = n_randint(state, 2);

        slong nterms = ctx->prec;
        acb_holonomic_sum_divconquer(ctx, nterms);

        acb_poly_struct * ref = _acb_poly_vec_init(dop_order);

        acb_poly_one(ref);

        for (slong i = 0; i < ctx->npts; i++)
        {
            acb_poly_zero(ref + 1);
            acb_poly_set_coeff_acb(ref + 1, 0, ctx->pts + i);
            acb_poly_set_coeff_si(ref + 1, 1, 1);
            acb_poly_atan_series(ref + 1, ref + 1, ctx->nder, ctx->prec);

            for (slong m = 0; m < dop_order; m++)
            {
                if (!acb_poly_overlaps(ref + m, acb_holonomic_sol_sum_ptr(ctx->sol + m, i, 0)))
                {
                    flint_printf("FAIL\n\n");
                    flint_printf("i = %wd, m = %wd, a = %{acb}\n", i, m, ctx->pts + i);
                    for (slong j = 0; j < 2; j++)
                        flint_printf("sum[%wd] = %{acb_poly}\n", j,
                                     acb_holonomic_sol_sum_ptr(ctx->sol + m, i, 0));
                    flint_printf("ref = %{acb_poly}\n", ref + m);
                    flint_abort();
                }
            }
        }

        if (ctx->flags & ACB_HOLONOMIC_WANT_SERIES)
        {
            for (int m = 0; m < dop_order; m++)  /* XXX tmp */
                acb_poly_truncate(ctx->sol[m].series, nterms);

            acb_poly_zero(ref + 1);
            acb_poly_set_coeff_si(ref + 1, 1, 1);
            acb_poly_atan_series(ref + 1, ref + 1, nterms, ctx->prec);

            for (slong m = 0; m < dop_order; m++)
            {
                if (!acb_poly_overlaps(ref + m, ctx->sol[m].series))
                {
                    flint_printf("FAIL\n\n");
                    flint_printf("m = %wd, nterms = %wd\n", m, nterms);
                    flint_printf("series = %{acb_poly}\n", ctx->sol[m].series);
                    flint_printf("ref = %{acb_poly}\n", ref);
                    flint_abort();
                }
            }

        }

        _acb_poly_vec_clear(ref, dop_order);

        acb_holonomic_sum_context_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
