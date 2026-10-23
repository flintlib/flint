/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include <string.h>
#include "arb.h"
#include "decimal.h"
#include "gr.h"

TEST_FUNCTION_START(arb_get_str, state)
{
    slong iter;

    /* just test no crashing... */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x;
        char * s;
        slong n;

        arb_init(x);

        arb_randtest_special(x, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));

        n = 1 + n_randint(state, 300);

        s = arb_get_str(x, n, (n_randint(state, 2) * ARB_STR_MORE)
                            |  (n_randint(state, 2) * ARB_STR_NO_RADIUS)
                            | (ARB_STR_CONDENSE * n_randint(state, 50)));

        flint_free(s);
        arb_clear(x);
    }

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, y;
        char * s;
        slong n, prec;
        int conversion_error;

        arb_init(x);
        arb_init(y);

        arb_randtest_special(x, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        arb_randtest_special(y, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));

        n = 1 + n_randint(state, 300);
        prec = 2 + n_randint(state, 1000);

        s = arb_get_str(x, n, n_randint(state, 2) * ARB_STR_MORE);
        conversion_error = arb_set_str(y, s, prec);

        if (conversion_error || !arb_contains(y, x))
        {
            flint_printf("FAIL (roundtrip)  iter = %wd\n", iter);
            flint_printf("x = "); arb_printd(x, 50); flint_printf("\n\n");
            flint_printf("s = %s", s); flint_printf("\n\n");
            flint_printf("y = "); arb_printd(y, 50); flint_printf("\n\n");
            flint_abort();
        }

        flint_free(s);
        arb_clear(x);
        arb_clear(y);
    }

    /* tightness: exact values, tiny radii, one-ulp rule */
    {
        arb_t x;
        fmpz_t t;
        char * s;
        slong iter2;

        arb_init(x);
        fmpz_init(t);

        fmpz_ui_pow_ui(t, 10, 100);
        arb_set_fmpz(x, t);
        s = arb_get_str(x, 31, 0);
        if (strcmp(s, "1.000000000000000000000000000000e+100"))
        {
            flint_printf("FAIL (exact 10^100)\n%s\n", s);
            flint_abort();
        }
        flint_free(s);

        arb_one(x);
        arb_add_error_2exp_si(x, -1000);
        s = arb_get_str(x, 30, 0);
        if (strcmp(s, "[1.00000000000000000000000000000 +/- 9.34e-302]"))
        {
            flint_printf("FAIL (tiny radius)\n%s\n", s);
            flint_abort();
        }
        flint_free(s);

        arb_set_d(x, 0.125);
        s = arb_get_str(x, 2, 0);
        if (strcmp(s, "[0.13 +/- 5.00e-3]"))
        {
            flint_printf("FAIL (rounding)\n%s\n", s);
            flint_abort();
        }
        flint_free(s);

        /* d 10^k with huge k: exact, and the printed value is exact */
        for (iter2 = 0; iter2 < 20 * flint_test_multiplier(); iter2++)
        {
            fmpz_t d;
            char * s2;
            fmpz_init(d);
            fmpz_randtest_not_zero(d, state, 40);
            fmpz_ui_pow_ui(t, 10, 20000 + n_randint(state, 20000));
            fmpz_mul(t, t, d);
            arb_set_fmpz(x, t);
            s = arb_get_str(x, 20, 0);
            if (strchr(s, '+') != NULL && strstr(s, "+/-") != NULL)
            {
                flint_printf("FAIL (exact huge)\n%s\n", s);
                flint_abort();
            }
            s2 = arb_get_str(x, 20, ARB_STR_NO_RADIUS);
            if (strcmp(s, s2))
            {
                flint_printf("FAIL (exact huge no radius)\n%s\n%s\n", s, s2);
                flint_abort();
            }
            flint_free(s);
            flint_free(s2);
            fmpz_clear(d);
        }

        /* the printed radius is at most one ulp and covers the true error */
        for (iter2 = 0; iter2 < 1000 * flint_test_multiplier(); iter2++)
        {
            arb_t y;
            slong n;
            char * p;

            arb_init(y);
            arb_randtest(x, state, 1 + n_randint(state, 300), 1 + n_randint(state, 200));
            if (n_randint(state, 2))
                mag_zero(arb_radref(x));
            n = 1 + n_randint(state, 40);
            s = arb_get_str(x, n, 0);
            p = strstr(s, "+/-");

            if (arb_set_str(y, s, 400) || !arb_contains(y, x))
            {
                flint_printf("FAIL (containment)\n%s\n", s);
                flint_abort();
            }

            if (p != NULL && s[1] != '+')
            {
                /* [mid +/- rad]: check rad <= 10^(E - k + 1) where the
                   midpoint has k significant digits and exponent E */
                gr_ctx_t ctx;
                decfloat_t r, ulp;
                char * q;
                char * epos;
                slong k = 0, zeros = 0, intdigits = 0;
                fmpz_t E, one;
                int seen_nonzero = 0, seen_point = 0;

                fmpz_init(E);
                fmpz_init(one);
                fmpz_one(one);

                gr_ctx_init_decfloat(ctx, 10, 0);
                decfloat_init(r, ctx);
                decfloat_init(ulp, ctx);

                epos = NULL;
                for (q = s + 1; q < p; q++)
                {
                    if (*q == 'e') { epos = q; break; }
                    if (*q == '.') { seen_point = 1; continue; }
                    if (*q < '0' || *q > '9') continue;
                    if (!seen_point) intdigits++;
                    if (*q != '0' || seen_nonzero) { seen_nonzero = 1; k++; }
                    else if (seen_point) zeros++;
                }

                if (epos != NULL)
                {
                    q = flint_malloc(p - epos);
                    strncpy(q, epos + 1 + (epos[1] == '+'), p - epos - 2);
                    q[p - epos - 2 - (epos[1] == '+')] = '\0';
                    if (fmpz_set_str(E, q, 10))
                        flint_abort();
                    flint_free(q);
                }
                else if (s[1] == '0' || (s[1] == '-' && s[2] == '0'))
                    fmpz_set_si(E, -1 - zeros);
                else
                    fmpz_set_si(E, intdigits - 1);

                /* parse the radius (up to the bracket) */
                q = flint_malloc(strlen(p));
                strcpy(q, p + 4);
                q[strlen(q) - 1] = '\0';
                if (decfloat_set_str(r, q, ctx) != GR_SUCCESS)
                    flint_abort();
                flint_free(q);
                fmpz_sub_ui(E, E, k - 1);
                GR_MUST_SUCCEED(decfloat_set_fmpz_10exp_fmpz(ulp, one, E, ctx));

                if (_decfloat_cmp(r, ulp, ctx) > 0)
                {
                    flint_printf("FAIL (one ulp)\n%s\nk = %wd E - k + 1 = %{fmpz}\n", s, k, E);
                    flint_abort();
                }

                fmpz_clear(E);
                fmpz_clear(one);
                decfloat_clear(r, ctx);
                decfloat_clear(ulp, ctx);
                gr_ctx_clear(ctx);
            }
            flint_free(s);
            arb_clear(y);
        }

        arb_clear(x);
        fmpz_clear(t);
    }

    /* test ARB_STR_NO_RADIUS */
    {
        arb_t x;
        char * s;

        arb_init(x);

        arb_set_str(x, "3.1415926535897932", 53);
        s = arb_get_str(x, 10, ARB_STR_NO_RADIUS);
        if (strcmp(s, "3.141592654"))
        {
            flint_printf("FAIL (ARB_STR_NO_RADIUS)\n");
            flint_printf("%s\n", s);
            flint_abort();
        }
        flint_free(s);

        arb_set_str(x, "+/- 3.45e-10", 53);
        s = arb_get_str(x, 10, ARB_STR_NO_RADIUS);
        if (strcmp(s, "0e-9"))
        {
            flint_printf("FAIL (ARB_STR_NO_RADIUS)\n");
            flint_printf("%s\n", s);
            flint_abort();
        }
        flint_free(s);

        arb_set_str(x, "+/- 3.45e10", 53);
        s = arb_get_str(x, 10, ARB_STR_NO_RADIUS);
        if (strcmp(s, "0e+11"))
        {
            flint_printf("FAIL (ARB_STR_NO_RADIUS)\n");
            flint_printf("%s\n", s);
            flint_abort();
        }
        flint_free(s);

        arb_set_str(x, "5e10 +/- 6e10", 53);
        s = arb_get_str(x, 10, ARB_STR_NO_RADIUS);
        if (strcmp(s, "0e+12"))
        {
            flint_printf("FAIL (ARB_STR_NO_RADIUS)\n");
            flint_printf("%s\n", s);
            flint_abort();
        }
        flint_free(s);

        arb_set_str(x, "5e-100000000000000000002 +/- 5e-100000000000000000002", 53);
        s = arb_get_str(x, 10, ARB_STR_NO_RADIUS);
        if (strcmp(s, "0e-100000000000000000000"))
        {
            flint_printf("FAIL (ARB_STR_NO_RADIUS)\n");
            flint_printf("%s\n", s);
            flint_abort();
        }
        flint_free(s);

        arb_clear(x);
    }

    TEST_FUNCTION_END(state);
}
