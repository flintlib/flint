/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Martin Lee

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_powmod_x_ui_preinv, state)
{
    int i, result;

    /* No aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, res1, res2, t, f, finv;
        mp_limb_t n;
        ulong exp;

        n = n_randtest_prime(state, 0);
        exp = n_randlimb(state) % 32;

        nmod_poly_init(f, n);
        nmod_poly_init(finv, n);
        nmod_poly_init(res1, n);
        nmod_poly_init(res2, n);
        nmod_poly_init(t, n);

        do {
            nmod_poly_randtest(f, state, n_randint(state, 50));
        } while (nmod_poly_is_zero(f));

        nmod_poly_reverse(finv, f, f->length);
        nmod_poly_inv_series(finv, finv, f->length);

        nmod_poly_init2 (a, n, f->length-1);
        nmod_poly_set_coeff_ui (a, 1, 1);

        nmod_poly_powmod_x_ui_preinv(res1, exp, f, finv);
        nmod_poly_powmod_ui_binexp_preinv (res2, a, exp, f, finv);

        result = (nmod_poly_equal(res1, res2));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("exp: %wu\n\n", exp);
            flint_printf("a:\n"); nmod_poly_print(a), flint_printf("\n\n");
            flint_printf("f:\n"); nmod_poly_print(f), flint_printf("\n\n");
            flint_printf("res1:\n"); nmod_poly_print(res1), flint_printf("\n\n");
            flint_printf("res2:\n"); nmod_poly_print(res2), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(f);
        nmod_poly_clear(finv);
        nmod_poly_clear(res1);
        nmod_poly_clear(res2);
        nmod_poly_clear(t);
    }

    /* Aliasing of res and f */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        nmod_poly_t res1, t, f, finv;
        mp_limb_t n;
        ulong exp;

        n = n_randtest_prime(state, 0);
        exp = n_randlimb(state) % 32;

        nmod_poly_init(f, n);
        nmod_poly_init(finv, n);
        nmod_poly_init(res1, n);
        nmod_poly_init(t, n);

        do {
            nmod_poly_randtest(f, state, n_randint(state, 50));
        } while (nmod_poly_is_zero(f));

        nmod_poly_reverse(finv, f, f->length);
        nmod_poly_inv_series(finv, finv, f->length);

        nmod_poly_powmod_x_ui_preinv(res1, exp, f, finv);
        nmod_poly_powmod_x_ui_preinv(f, exp, f, finv);

        result = (nmod_poly_equal(res1, f));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("exp: %wu\n\n", exp);
            flint_printf("f:\n"); nmod_poly_print(f), flint_printf("\n\n");
            flint_printf("res1:\n"); nmod_poly_print(res1), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(f);
        nmod_poly_clear(finv);
        nmod_poly_clear(res1);
        nmod_poly_clear(t);
    }

    /* Aliasing of res and finv */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        nmod_poly_t res1, t, f, finv;
        mp_limb_t n;
        ulong exp;

        n = n_randtest_prime(state, 0);
        exp = n_randlimb(state) % 32;

        nmod_poly_init(f, n);
        nmod_poly_init(finv, n);
        nmod_poly_init(res1, n);
        nmod_poly_init(t, n);

        do {
            nmod_poly_randtest(f, state, n_randint(state, 50));
        } while (nmod_poly_is_zero(f));

        nmod_poly_reverse(finv, f, f->length);
        nmod_poly_inv_series(finv, finv, f->length);

        nmod_poly_powmod_x_ui_preinv(res1, exp, f, finv);
        nmod_poly_powmod_x_ui_preinv(finv, exp, f, finv);

        result = (nmod_poly_equal(res1, finv));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("exp: %wu\n\n", exp);
            flint_printf("f:\n"); nmod_poly_print(f), flint_printf("\n\n");
            flint_printf("finv:\n"); nmod_poly_print(finv), flint_printf("\n\n");
            flint_printf("res1:\n"); nmod_poly_print(res1), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(f);
        nmod_poly_clear(finv);
        nmod_poly_clear(res1);
        nmod_poly_clear(t);
    }

    TEST_FUNCTION_END(state);
}
