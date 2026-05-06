#include "test_helpers.h"
#include "acb_types.h"
#include "acb_poly.h"
#include "acb_ode.h"

// XXX very very incomplete and too involved

TEST_FUNCTION_START(acb_ode_sum_divconquer, state)
{

    for (slong iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        /* classical functions -- well for now just atan */

        acb_ode_sum_t sum;

        slong dop_order = 2;
        acb_ode_sum_init(sum, dop_order + 1,
                                 n_randint(state, 5),
                                 dop_order,
                                 1 + n_randint(state, 3));

        /* sum_set_diffop */
        acb_poly_set_coeff_si(sum->dop + 2, 2, 1);
        acb_poly_set_coeff_si(sum->dop + 2, 0, 1);
        acb_poly_set_coeff_si(sum->dop + 1, 2, 1);
        acb_poly_set_coeff_si(sum->dop + 1, 0, -1);
        mag_one(sum->cvrad);

        acb_ode_sum_set_ordinary(sum);
        acb_ode_sum_set_ini_echelon(sum);

        slong prec = n_randint(state, 1000);
        /* XXX at the moment nterms must be a multiple of the degrees due to
           incomplete termination logic in sum_divconquer in other cases */
        slong nterms = n_randint(state, 500) * 2;

        for (slong i = 0; i < sum->npts; i++)
        {
            acb_ptr a = sum->pts + i;
            acb_randtest_precise(a, state, prec, 0);
            acb_div_si(a, a, 4, prec);
        }

        /* note that sum->bound is null */
        sum->flags = ACB_ODE_APPROX
                   | (n_randint(state, 2) * ACB_ODE_WANT_SERIES);

        acb_ode_sum_divconquer(sum, nterms, prec);

        acb_poly_struct * ref = _acb_poly_vec_init(dop_order);
        acb_poly_one(ref);  /* atan: first sol=1 */

        if (prec <= 12)
            goto cleanup;

        if (nterms >= prec)  /* todo: test rigorous version */
        {
            for (slong i = 0; i < sum->npts; i++)
            {
                acb_poly_zero(ref + 1);  /* 2nd sol = atan */
                acb_poly_set_coeff_acb(ref + 1, 0, sum->pts + i);
                acb_poly_set_coeff_si(ref + 1, 1, 1);
                /* APPROX => compute the ref to lower precision */
                acb_poly_atan_series(ref + 1, ref + 1, sum->nder, prec - 8);

                for (slong m = 0; m < dop_order; m++)
                {
                    if (!_acb_poly_overlaps(
                            acb_ode_sol_sum_ptr(sum->sol + m, i, 0, 0),
                            sum->sol->nder,
                            (ref + m)->coeffs,
                            (ref + m)->length))
                    {
                        flint_printf("FAIL (nonrigorous partial sum)\n\n");
                        flint_printf("prec = %wd, nterms = %wd\n", prec, nterms);
                        flint_printf("i = %wd, m = %wd, a = %{acb}\n", i, m, sum->pts + i);
                        for (slong j = 0; j < 2; j++)
                            flint_printf("sum[%wd] = %{acb*}\n", j, acb_ode_sol_sum_ptr(sum->sol + m, i, 0, 0), sum->sol->nder);

                        flint_printf("ref = %{acb_poly}\n", ref + m);
                        flint_abort();
                    }
                }
            }
        }

        if (sum->flags & ACB_ODE_WANT_SERIES)
        {
            if (sum->npts == 0 && sum->n != nterms)
                flint_printf("FAIL (wrong number of terms)\n");

            acb_poly_zero(ref + 1);
            acb_poly_set_coeff_si(ref + 1, 1, 1);
            /* APPROX => compute the ref to lower precision */
            acb_poly_atan_series(ref + 1, ref + 1, nterms, prec - 8);

            for (slong m = 0; m < dop_order; m++)
            {
                /* in general, only the first sum->n terms are correct */
                for (slong n = 0; n < sum->n; n++)
                {
                    if (n >= (ref + m)->length)
                        break;  // for now
                    if (!acb_overlaps((sum->sol[m].series)->coeffs + n,
                                      (ref + m)->coeffs + n))
                    {
                        flint_printf("FAIL (series)\n\n");
                        flint_printf("m = %wd, nterms = %wd, sum->n = %wd, n=%wd\n", m, nterms, sum->n, n);
                        flint_printf("series = %{acb_poly}\n", sum->sol[m].series);
                        flint_printf("ref = %{acb_poly}\n", ref + m);
                        flint_abort();
                    }
                }
            }
        }

cleanup:
        _acb_poly_vec_clear(ref, dop_order);
        acb_ode_sum_clear(sum);
    }

    TEST_FUNCTION_END(state);
}
