/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "test_helpers.h"
#include "ulong_extras.h"
#include "padic.h"
#include "padic_poly.h"

TEST_FUNCTION_START(padic_poly_get_set_coeff_padic, state)
{
    int i, result;

    /* Random coefficients: get(set(c)) agrees with c at the lower
       precision, and the polynomial stays normalised and canonical */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        padic_ctx_t ctx;
        fmpz_t p;
        slong N, N2, n;
        padic_poly_t f;
        padic_t c, d, e;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;
        N2 = N - 10 + n_randint(state, 21);
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_poly_init2(f, 0, N);
        padic_init2(c, N2);
        padic_init2(d, N2);
        padic_init2(e, FLINT_MIN(N, N2));

        padic_poly_randtest(f, state, n_randint(state, 10), ctx);
        n = n_randint(state, 12);
        if (n_randint(state, 4) == 0)
            padic_zero(c);
        else
            padic_randtest(c, state, ctx);

        padic_poly_set_coeff_padic(f, n, c, ctx);
        padic_poly_get_coeff_padic(d, f, n, ctx);
        padic_set(e, c, ctx);

        result = (padic_prec(d) == N2)
              && padic_poly_is_canonical(f, ctx)
              && padic_poly_is_reduced(f, ctx)
              && (f->length == 0 || !fmpz_is_zero(f->coeffs + f->length - 1));

        /* d = c mod p^min(N, N2) */
        padic_prec(d) = FLINT_MIN(N, N2);
        padic_reduce(d, ctx);
        result = result && padic_equal(d, e);

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("N = %wd, N2 = %wd, n = %wd\n", N, N2, n);
            flint_printf("f = "), padic_poly_debug(f), flint_printf("\n\n");
            flint_printf("c = "), padic_print(c, ctx), flint_printf("\n\n");
            flint_printf("d = "), padic_print(d, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        padic_poly_clear(f);
        padic_clear(c);
        padic_clear(d);
        padic_clear(e);

        padic_ctx_clear(ctx);
        fmpz_clear(p);
    }

    /* Deterministic cases */
    {
        padic_ctx_t ctx;
        fmpz_t p;
        padic_poly_t f;
        padic_t c;

        fmpz_init_set_ui(p, 7);
        padic_ctx_init(ctx, p, 0, 30, PADIC_TERSE);
        padic_poly_init2(f, 0, 20);
        padic_init2(c, 20);

        /* get_coeff keeps the precision of c */
        padic_poly_set_si(f, 1000, ctx);
        padic_prec(c) = 5;
        padic_poly_get_coeff_padic(c, f, 0, ctx);
        if (padic_prec(c) != 5 || fmpz_cmp_ui(padic_unit(c), 1000 % 16807) != 0)
        {
            flint_printf("FAIL (get_coeff precision)\n");
            fflush(stdout);
            flint_abort();
        }

        /* setting the leading coefficient to zero normalises */
        padic_prec(c) = 20;
        padic_poly_zero(f);
        padic_set_si(c, 3, ctx);
        padic_poly_set_coeff_padic(f, 0, c, ctx);
        padic_poly_set_coeff_padic(f, 2, c, ctx);
        padic_zero(c);
        padic_poly_set_coeff_padic(f, 2, c, ctx);
        if (padic_poly_length(f) != 1)
        {
            flint_printf("FAIL (set_coeff normalisation)\n");
            fflush(stdout);
            flint_abort();
        }

        /* printing 3x has no constant term (assume tmpfile() is broken
           on windows) */
#if !defined(_MSC_VER) && !defined(__MINGW32__)
        {
            FILE * tmp = tmpfile();
            char buf[64] = { 0 };

            if (tmp != NULL)
            {
                padic_poly_zero(f);
                padic_set_si(c, 3, ctx);
                padic_poly_set_coeff_padic(f, 1, c, ctx);
                padic_poly_fprint_pretty(tmp, f, "x", ctx);
                rewind(tmp);
                if (fgets(buf, sizeof(buf), tmp) == NULL)
                    buf[0] = '\0';
                fclose(tmp);

                if (strcmp(buf, "(3)*x") != 0)
                {
                    flint_printf("FAIL (print_pretty)\n");
                    flint_printf("%s\n", buf);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }
#endif

        padic_poly_clear(f);
        padic_clear(c);
        padic_ctx_clear(ctx);
        fmpz_clear(p);
    }

    TEST_FUNCTION_END(state);
}
