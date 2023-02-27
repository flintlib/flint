/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"


int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("add_naive....");
    fflush(stdout);

    flint_randinit(state);

    /* test exact addition: (x + y) + z == x + (y + z) */
    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        slong bits, res1, res2, res3, res4;
        fmpr_t x, y, z, t, u;

        bits = 2 + n_randint(state, 200);

        fmpr_init(x);
        fmpr_init(y);
        fmpr_init(z);
        fmpr_init(t);
        fmpr_init(u);

        fmpr_randtest_special(x, state, bits, 10);
        fmpr_randtest_special(y, state, bits, 10);
        fmpr_randtest_special(z, state, bits, 10);
        fmpr_randtest_special(t, state, bits, 10);
        fmpr_randtest_special(u, state, bits, 10);

        res1 = fmpr_add_naive(t, x, y, FMPR_PREC_EXACT, FMPR_RND_DOWN);
        res2 = fmpr_add_naive(t, t, z, FMPR_PREC_EXACT, FMPR_RND_DOWN);
        res3 = fmpr_add_naive(u, y, z, FMPR_PREC_EXACT, FMPR_RND_DOWN);
        res4 = fmpr_add_naive(u, x, u, FMPR_PREC_EXACT, FMPR_RND_DOWN);

        if (!fmpr_equal(t, u) ||
            res1 != FMPR_RESULT_EXACT || res2 != FMPR_RESULT_EXACT ||
            res3 != FMPR_RESULT_EXACT || res4 != FMPR_RESULT_EXACT)
        {
            flint_printf("FAIL\n\n");
            flint_printf("bits = %wd\n", bits);
            flint_printf("x = "); fmpr_print(x); flint_printf("\n\n");
            flint_printf("y = "); fmpr_print(y); flint_printf("\n\n");
            flint_printf("z = "); fmpr_print(z); flint_printf("\n\n");
            flint_printf("t = "); fmpr_print(t); flint_printf("\n\n");
            flint_printf("u = "); fmpr_print(u); flint_printf("\n\n");
            flint_abort();
        }

        fmpr_clear(x);
        fmpr_clear(y);
        fmpr_clear(z);
        fmpr_clear(t);
        fmpr_clear(u);
    }

    /* compare rounding with mpfr */
    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        slong bits, res;
        int mpfr_res;
        fmpr_t x, y, z, w;
        mpfr_t X, Y, Z;

        bits = 2 + n_randint(state, 200);

        fmpr_init(x);
        fmpr_init(y);
        fmpr_init(z);
        fmpr_init(w);

        mpfr_init2(X, bits + 100);
        mpfr_init2(Y, bits + 100);
        mpfr_init2(Z, bits);

        fmpr_randtest_special(x, state, bits + n_randint(state, 100), 10);
        fmpr_randtest_special(y, state, bits + n_randint(state, 100), 10);
        fmpr_randtest_special(z, state, bits, 10);

        fmpr_get_mpfr(X, x, MPFR_RNDN);
        fmpr_get_mpfr(Y, y, MPFR_RNDN);

        switch (n_randint(state, 4))
        {
            case 0:
                mpfr_res = mpfr_add(Z, X, Y, MPFR_RNDZ);
                res = fmpr_add_naive(z, x, y, bits, FMPR_RND_DOWN);
                break;
            case 1:
                mpfr_res = mpfr_add(Z, X, Y, MPFR_RNDA);
                res = fmpr_add_naive(z, x, y, bits, FMPR_RND_UP);
                break;
            case 2:
                mpfr_res = mpfr_add(Z, X, Y, MPFR_RNDD);
                res = fmpr_add_naive(z, x, y, bits, FMPR_RND_FLOOR);
                break;
            default:
                mpfr_res = mpfr_add(Z, X, Y, MPFR_RNDU);
                res = fmpr_add_naive(z, x, y, bits, FMPR_RND_CEIL);
                break;
        }

        fmpr_set_mpfr(w, Z);

        if (!fmpr_equal(z, w) || (res == FMPR_RESULT_EXACT) != (mpfr_res == 0)
            || !fmpr_check_ulp(z, res, bits))
        {
            flint_printf("FAIL\n\n");
            flint_printf("bits = %wd\n", bits);
            flint_printf("x = "); fmpr_print(x); flint_printf("\n\n");
            flint_printf("y = "); fmpr_print(y); flint_printf("\n\n");
            flint_printf("z = "); fmpr_print(z); flint_printf("\n\n");
            flint_printf("w = "); fmpr_print(w); flint_printf("\n\n");
            flint_printf("returned: %wd, %d\n", res, mpfr_res);
            flint_abort();
        }

        /* check error bound */
        if (!fmpr_is_nan(z) && !fmpr_is_inf(z))
        {
            fmpr_t z_exact, error_bound, true_error;

            fmpr_init(z_exact);
            fmpr_init(error_bound);
            fmpr_init(true_error);

            fmpr_set_error_result(error_bound, z, res);
            fmpr_add_naive(z_exact, x, y, FMPR_PREC_EXACT, FMPR_RND_DOWN);

            fmpr_sub(true_error, z, z_exact, FMPR_PREC_EXACT, FMPR_RND_DOWN);
            fmpr_abs(true_error, true_error);

            if (fmpr_is_zero(error_bound) != fmpr_is_zero(true_error) ||
                fmpr_cmp(true_error, error_bound) > 0)
            {
                flint_printf("FAIL: error bound\n\n");
                flint_printf("bits = %wd\n", bits);
                flint_printf("x = "); fmpr_print(x); flint_printf("\n\n");
                flint_printf("y = "); fmpr_print(y); flint_printf("\n\n");
                flint_printf("z = "); fmpr_print(z); flint_printf("\n\n");
                flint_printf("z_exact = "); fmpr_print(z_exact); flint_printf("\n\n");
                flint_printf("true_error = "); fmpr_print(true_error); flint_printf("\n\n");
                flint_printf("error_bound = "); fmpr_print(error_bound); flint_printf("\n\n");
                flint_abort();
            }

            fmpr_clear(z_exact);
            fmpr_clear(error_bound);
            fmpr_clear(true_error);
        }

        fmpr_clear(x);
        fmpr_clear(y);
        fmpr_clear(z);
        fmpr_clear(w);

        mpfr_clear(X);
        mpfr_clear(Y);
        mpfr_clear(Z);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
