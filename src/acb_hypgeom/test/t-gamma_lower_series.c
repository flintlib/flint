/*
    Copyright (C) 2016 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_poly.h"
#include "acb_hypgeom.h"

TEST_FUNCTION_START(acb_hypgeom_gamma_lower_series, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n1, n2, bits1, bits2, bits3;
        acb_poly_t S, A, B, C;
        acb_t s, t, c;
        int regularized;

        regularized = n_randint(state, 3);

        bits1 = 2 + n_randint(state, 200);
        bits2 = 2 + n_randint(state, 200);
        bits3 = 2 + n_randint(state, 200);

        m = 1 + n_randint(state, 10);
        n1 = 1 + n_randint(state, 10);
        n2 = 1 + n_randint(state, 10);

        acb_poly_init(S);
        acb_poly_init(A);
        acb_poly_init(B);
        acb_poly_init(C);
        acb_init(s);
        acb_init(t);
        acb_init(c);

        acb_poly_randtest(S, state, m, bits1, 3);
        acb_poly_randtest(A, state, m, bits1, 3);
        acb_poly_randtest(B, state, m, bits1, 3);

        if (n_randint(state, 10))
        {
            acb_randtest(s, state, bits1, 3);
        }
        else
        {
            fmpz_t k;
            fmpz_init(k);
            fmpz_randtest(k, state, 100);
            acb_set_fmpz(s, k);
            fmpz_clear(k);
        }

        acb_hypgeom_gamma_lower_series(A, s, S, regularized, n1, bits2);
        acb_hypgeom_gamma_lower_series(B, s, S, regularized, n2, bits3);

        acb_poly_set(C, A);
        acb_poly_truncate(C, FLINT_MIN(n1, n2));
        acb_poly_truncate(B, FLINT_MIN(n1, n2));

        if (!acb_poly_overlaps(B, C))
        {
            flint_printf("FAIL (consistency)\n\n");
            flint_printf("regularized = %d\n\n", regularized);
            flint_printf("S = "); acb_poly_printd(S, 15); flint_printf("\n\n");
            flint_printf("A = "); acb_poly_printd(A, 15); flint_printf("\n\n");
            flint_printf("B = "); acb_poly_printd(B, 15); flint_printf("\n\n");
            flint_printf("C = "); acb_poly_printd(C, 15); flint_printf("\n\n");
            flint_abort();
        }

        /* f(h(x)) = exp(-h(x)) h(x)^(s-1) */
        acb_poly_neg(C, S);
        acb_poly_exp_series(C, C, n1, bits2);
        acb_sub_ui(t, s, 1, bits2);
        acb_poly_pow_acb_series(B, S, t, n1, bits2);
        acb_poly_mullow(C, C, B, n1, bits2);

        if (regularized == 0)
        {
            /* integral(f(h(x)) h'(x))' = f(h(x)) h'(x) */
            acb_poly_derivative(B, S, bits2);
            acb_poly_mullow(C, C, B, n1, bits2);

            acb_poly_truncate(C, n1 - 1);
        }
        else if (regularized == 1)
        {
            /* (integral(f(h(x)) h'(x)) / c)' = (f(h(x)) h'(x)) / c */
            acb_poly_derivative(B, S, bits2);
            acb_poly_mullow(C, C, B, n1, bits2);

            acb_gamma(c, s, bits2);
            acb_poly_scalar_div(C, C, c, bits2);

            acb_poly_truncate(C, n1 - 1);
        }
        else if (regularized == 2)
        {
            if (acb_is_int(s) &&
                arf_sgn(arb_midref(acb_realref(s))) <= 0 &&
                arf_cmpabs_2exp_si(arb_midref(acb_realref(s)), 30) < 0)
            {
                /* (h(x)^-s)' = -s h(x)^-(s+1) h'(x) */
                acb_poly_derivative(C, S, bits2);

                acb_add_ui(t, s, 1, bits2);
                acb_neg(t, t);
                acb_poly_pow_acb_series(B, S, t, n1, bits2);

                acb_poly_mullow(C, C, B, n1, bits2);

                acb_poly_scalar_mul(C, C, s, bits2);
                acb_poly_neg(C, C);

                acb_poly_truncate(C, n1 - 1);
            }
            else
            {
                /* (h(x)^-s integral(f(h(x)) h'(x)) / c)' =
                 * h(x)^-(s+1) (h(x) f(h(x)) - s integral(f(h(x)) h'(x))) h'(x) / c */
                acb_poly_t D;
                acb_poly_init(D);

                acb_poly_derivative(B, S, bits2);
                acb_poly_mullow(D, C, B, n1, bits2);
                acb_poly_integral(D, D, bits2);
                acb_poly_scalar_mul(D, D, s, bits2);
                acb_poly_mullow(C, C, S, n1, bits2);
                acb_poly_sub(D, C, D, bits2);

                acb_add_ui(t, s, 1, bits2);
                acb_neg(t, t);
                acb_poly_pow_acb_series(B, S, t, n1, bits2);

                acb_poly_mullow(C, D, B, n1, bits2);

                acb_poly_derivative(B, S, bits2);
                acb_poly_mullow(C, C, B, n1, bits2);

                acb_gamma(c, s, bits2);
                acb_poly_scalar_div(C, C, c, bits2);

                acb_poly_truncate(C, n1 - 1);

                acb_poly_clear(D);
            }
        }

        acb_poly_derivative(B, A, bits2);

        if (!acb_poly_overlaps(B, C))
        {
            flint_printf("FAIL (derivative)\n\n");
            flint_printf("regularized = %d\n\n", regularized);
            flint_printf("S = "); acb_poly_printd(S, 15); flint_printf("\n\n");
            flint_printf("A = "); acb_poly_printd(A, 15); flint_printf("\n\n");
            flint_printf("B = "); acb_poly_printd(B, 15); flint_printf("\n\n");
            flint_printf("C = "); acb_poly_printd(C, 15); flint_printf("\n\n");
            flint_abort();
        }

        acb_hypgeom_gamma_lower_series(S, s, S, regularized, n1, bits2);

        if (!acb_poly_overlaps(A, S))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_printf("regularized = %d\n\n", regularized);
            flint_printf("S = "); acb_poly_printd(S, 15); flint_printf("\n\n");
            flint_printf("A = "); acb_poly_printd(A, 15); flint_printf("\n\n");
            flint_abort();
        }

        acb_poly_clear(S);
        acb_poly_clear(A);
        acb_poly_clear(B);
        acb_poly_clear(C);
        acb_clear(s);
        acb_clear(t);
        acb_clear(c);
    }

    TEST_FUNCTION_END(state);
}
