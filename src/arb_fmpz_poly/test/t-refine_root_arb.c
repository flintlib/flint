/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_vec.h"
#include "fmpq_vec.h"
#include "fmpz_poly.h"
#include "arb.h"
#include "arb_fmpz_poly.h"

/* Check that res is accurate to prec bits and that poly changes sign
   between the endpoints of res. */
static void
_check_refined_root(const arb_t res, const arb_t initial, const fmpz_poly_t poly, slong prec)
{
    arb_t a, b, fa, fb;
    slong wp;
    int sa = 0, sb = 0;

    if (arb_rel_accuracy_bits(res) < prec || !arb_contains(initial, res))
    {
        flint_printf("FAIL (accuracy)\n");
        flint_printf("prec = %wd\n", prec);
        flint_printf("poly = "); fmpz_poly_print(poly); flint_printf("\n");
        flint_printf("initial = "); arb_printd(initial, 30); flint_printf("\n");
        flint_printf("res = "); arb_printd(res, 30); flint_printf("\n");
        flint_abort();
    }

    arb_init(a);
    arb_init(b);
    arb_init(fa);
    arb_init(fb);

    arb_get_lbound_arf(arb_midref(a), res, ARF_PREC_EXACT);
    arb_get_ubound_arf(arb_midref(b), res, ARF_PREC_EXACT);

    for (wp = 64; wp < 100000; wp *= 2)
    {
        arb_fmpz_poly_evaluate_arb(fa, poly, a, wp);
        arb_fmpz_poly_evaluate_arb(fb, poly, b, wp);

        if (!arb_contains_zero(fa) && !arb_contains_zero(fb))
        {
            sa = arf_sgn(arb_midref(fa));
            sb = arf_sgn(arb_midref(fb));
            break;
        }

        /* exact zero at an endpoint is fine if res is exact */
        if (arb_is_exact(res) && (arb_is_zero(fa) || arb_is_zero(fb)))
        {
            sa = 1;
            sb = -1;
            break;
        }
    }

    if (sa == sb)
    {
        flint_printf("FAIL (sign change)\n");
        flint_printf("prec = %wd\n", prec);
        flint_printf("poly = "); fmpz_poly_print(poly); flint_printf("\n");
        flint_printf("initial = "); arb_printd(initial, 30); flint_printf("\n");
        flint_printf("res = "); arb_printd(res, 30); flint_printf("\n");
        flint_abort();
    }

    arb_clear(a);
    arb_clear(b);
    arb_clear(fa);
    arb_clear(fb);
}

TEST_FUNCTION_START(arb_fmpz_poly_refine_root_arb, state)
{
    slong iter;

    for (iter = 0; iter < 300 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpz_poly_t f, g;
        fmpq * ex;
        fmpz * iv;
        slong * iv_exp;
        slong i, j, n, n_ex, n_iv, deg, prec;
        arb_t x, y;
        fmpz_t t;

        fmpz_poly_init(f);
        fmpz_poly_init(g);
        arb_init(x);
        arb_init(y);
        fmpz_init(t);

        n = 1 + n_randint(state, 30);
        prec = 2 + n_randint(state, 1000);

        switch (n_randint(state, 4))
        {
            case 0:
                /* c * prod (x - j) + d: roots extremely close to integers */
                fmpz_poly_one(f);
                for (j = 1; j <= n; j++)
                {
                    fmpz_poly_zero(g);
                    fmpz_poly_set_coeff_si(g, 0, -j);
                    fmpz_poly_set_coeff_si(g, 1, 1);
                    fmpz_poly_mul(f, f, g);
                }
                fmpz_poly_scalar_mul_ui(f, f, 1 + n_randint(state, 3));
                fmpz_poly_get_coeff_fmpz(t, f, 0);
                fmpz_add_si(t, t, (slong) n_randint(state, 5) - 2);
                fmpz_poly_set_coeff_fmpz(f, 0, t);
                break;
            case 1:
                /* x^2 - eps: two roots close to 0 */
                fmpz_poly_zero(f);
                fmpz_poly_set_coeff_ui(f, 2, 1);
                fmpz_one(t);
                fmpz_mul_2exp(t, t, n_randint(state, 300));
                fmpz_poly_set_coeff_fmpz(f, 2, t);
                fmpz_poly_set_coeff_si(f, 0, -1 - (slong) n_randint(state, 3));
                fmpz_poly_zero(g);
                fmpz_poly_set_coeff_si(g, 0, -1 - (slong) n_randint(state, 100));
                fmpz_poly_set_coeff_si(g, 1, 1);
                fmpz_poly_pow(g, g, n_randint(state, 4));
                fmpz_poly_mul(f, f, g);
                fmpz_poly_get_coeff_fmpz(t, f, 0);
                fmpz_add_si(t, t, (slong) n_randint(state, 5) - 2);
                fmpz_poly_set_coeff_fmpz(f, 0, t);
                break;
            case 2:
                /* random roots near dyadic numbers */
                fmpz_poly_one(f);
                for (j = 0; j < FLINT_MIN(n, 10); j++)
                {
                    fmpz_poly_zero(g);
                    fmpz_poly_set_coeff_si(g, 0, (slong) n_randint(state, 64) - 32);
                    fmpz_poly_set_coeff_ui(g, 1, UWORD(1) << n_randint(state, 8));
                    fmpz_poly_mul(f, f, g);
                }
                fmpz_poly_get_coeff_fmpz(t, f, 0);
                fmpz_add_si(t, t, (slong) n_randint(state, 3) - 1);
                fmpz_poly_set_coeff_fmpz(f, 0, t);
                break;
            default:
                fmpz_poly_randtest(f, state, 1 + n, 1 + n_randint(state, 200));
                break;
        }

        fmpz_poly_squarefree_part(f, f);
        deg = fmpz_poly_degree(f);

        if (deg >= 1)
        {
            ex = _fmpq_vec_init(deg);
            iv = _fmpz_vec_init(deg);
            iv_exp = flint_malloc(sizeof(slong) * deg);

            fmpz_poly_isolate_real_roots(ex, &n_ex, iv, iv_exp, &n_iv, f);

            /* remove the exact roots so that the closed isolating intervals
               contain exactly one root */
            if (n_ex != 0)
            {
                fmpz_poly_product_roots_fmpq_vec(g, ex, n_ex);
                fmpz_poly_divexact(f, f, g);
            }

            for (i = 0; i < n_iv; i++)
            {
                fmpz_mul_2exp(t, iv + i, 1);
                fmpz_add_ui(t, t, 1);
                arb_set_fmpz(x, t);
                mag_one(arb_radref(x));
                arb_mul_2exp_si(x, x, iv_exp[i] - 1);

                arb_fmpz_poly_refine_root_arb(y, f, x, prec);
                _check_refined_root(y, x, f, prec);

                /* aliasing */
                arb_fmpz_poly_refine_root_arb(x, f, x, prec);
                if (!arb_overlaps(x, y))
                {
                    flint_printf("FAIL (aliasing)\n");
                    flint_abort();
                }
            }

            _fmpq_vec_clear(ex, deg);
            _fmpz_vec_clear(iv, deg);
            flint_free(iv_exp);
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        arb_clear(x);
        arb_clear(y);
        fmpz_clear(t);
    }

    TEST_FUNCTION_END(state);
}
