/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Martin Lee

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"
#include "fmpz.h"

TEST_FUNCTION_START(nmod_poly_powmod_fmpz_binexp_preinv, state)
{
    int i, result;

    /* Aliasing of res and a */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, res1, t, f, finv;
        mp_limb_t n;
        fmpz_t exp;

        fmpz_init(exp);

        n = n_randtest_prime(state, 0);
        fmpz_randtest_unsigned(exp, state, n_randint(state, 100) + 1);

        nmod_poly_init(a, n);
        nmod_poly_init(f, n);
        nmod_poly_init(finv, n);
        nmod_poly_init(res1, n);
        nmod_poly_init(t, n);

        nmod_poly_randtest(a, state, n_randint(state, 50));
        do {
            nmod_poly_randtest(f, state, n_randint(state, 50));
        } while (nmod_poly_is_zero(f));

        nmod_poly_reverse(finv, f, f->length);
        nmod_poly_inv_series(finv, finv, f->length);

        nmod_poly_powmod_fmpz_binexp_preinv(res1, a, exp, f, finv);
        nmod_poly_powmod_fmpz_binexp_preinv(a, a, exp, f, finv);

        result = (nmod_poly_equal(res1, a));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("exp: "); fmpz_print(exp); flint_printf("\n\n");
            flint_printf("a:\n"); nmod_poly_print(a), flint_printf("\n\n");
            flint_printf("f:\n"); nmod_poly_print(f), flint_printf("\n\n");
            flint_printf("res1:\n"); nmod_poly_print(res1), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(exp);
        nmod_poly_clear(a);
        nmod_poly_clear(f);
        nmod_poly_clear(finv);
        nmod_poly_clear(res1);
        nmod_poly_clear(t);
    }

    /* Aliasing of res and f */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, res1, t, f, finv;
        mp_limb_t n;
        fmpz_t exp;

        fmpz_init(exp);

        n = n_randtest_prime(state, 0);
        fmpz_randtest_unsigned(exp, state, n_randint(state, 100) + 1);

        nmod_poly_init(a, n);
        nmod_poly_init(f, n);
        nmod_poly_init(finv, n);
        nmod_poly_init(res1, n);
        nmod_poly_init(t, n);

        nmod_poly_randtest(a, state, n_randint(state, 50));
        do {
            nmod_poly_randtest(f, state, n_randint(state, 50));
        } while (nmod_poly_is_zero(f));

        nmod_poly_reverse(finv, f, f->length);
        nmod_poly_inv_series(finv, finv, f->length);

        nmod_poly_powmod_fmpz_binexp_preinv(res1, a, exp, f, finv);
        nmod_poly_powmod_fmpz_binexp_preinv(f, a, exp, f, finv);

        result = (nmod_poly_equal(res1, f));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("exp: "); fmpz_print(exp); flint_printf("\n\n");
            flint_printf("a:\n"); nmod_poly_print(a), flint_printf("\n\n");
            flint_printf("f:\n"); nmod_poly_print(f), flint_printf("\n\n");
            flint_printf("res1:\n"); nmod_poly_print(res1), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(exp);
        nmod_poly_clear(a);
        nmod_poly_clear(f);
        nmod_poly_clear(finv);
        nmod_poly_clear(res1);
        nmod_poly_clear(t);
    }

    /* Aliasing of res and finv */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, res1, t, f, finv;
        mp_limb_t n;
        fmpz_t exp;

        fmpz_init(exp);

        n = n_randtest_prime(state, 0);
        fmpz_randtest_unsigned(exp, state, n_randint(state, 100) + 1);

        nmod_poly_init(a, n);
        nmod_poly_init(f, n);
        nmod_poly_init(finv, n);
        nmod_poly_init(res1, n);
        nmod_poly_init(t, n);

        nmod_poly_randtest(a, state, n_randint(state, 50));
        do {
            nmod_poly_randtest(f, state, n_randint(state, 50));
        } while (nmod_poly_is_zero(f));

        nmod_poly_reverse(finv, f, f->length);
        nmod_poly_inv_series(finv, finv, f->length);

        nmod_poly_powmod_fmpz_binexp_preinv(res1, a, exp, f, finv);
        nmod_poly_powmod_fmpz_binexp_preinv(finv, a, exp, f, finv);

        result = (nmod_poly_equal(res1, finv));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("exp: "); fmpz_print(exp); flint_printf("\n\n");
            flint_printf("a:\n"); nmod_poly_print(a), flint_printf("\n\n");
            flint_printf("f:\n"); nmod_poly_print(f), flint_printf("\n\n");
            flint_printf("finv:\n"); nmod_poly_print(finv), flint_printf("\n\n");
            flint_printf("res1:\n"); nmod_poly_print(res1), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(exp);
        nmod_poly_clear(a);
        nmod_poly_clear(f);
        nmod_poly_clear(finv);
        nmod_poly_clear(res1);
        nmod_poly_clear(t);
    }

    /* No aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, res1, res2, t, f, finv;
        mp_limb_t n;
        fmpz_t exp;

        fmpz_init(exp);

        n = n_randtest_prime(state, 0);
        fmpz_randtest_unsigned(exp, state, n_randint(state, 100) + 1);;

        nmod_poly_init(a, n);
        nmod_poly_init(f, n);
        nmod_poly_init(finv, n);
        nmod_poly_init(res1, n);
        nmod_poly_init(res2, n);
        nmod_poly_init(t, n);

        nmod_poly_randtest(a, state, n_randint(state, 50));
        do {
            nmod_poly_randtest(f, state, n_randint(state, 50));
        } while (nmod_poly_is_zero(f));

        nmod_poly_reverse(finv, f, f->length);
        nmod_poly_inv_series(finv, finv, f->length);

        nmod_poly_powmod_fmpz_binexp_preinv(res1, a, exp, f, finv);

        nmod_poly_powmod_fmpz_binexp(res2, a, exp, f);

        result = (nmod_poly_equal(res1, res2));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("exp: "); fmpz_print(exp); flint_printf("\n\n");
            flint_printf("a:\n"); nmod_poly_print(a), flint_printf("\n\n");
            flint_printf("f:\n"); nmod_poly_print(f), flint_printf("\n\n");
            flint_printf("res1:\n"); nmod_poly_print(res1), flint_printf("\n\n");
            flint_printf("res2:\n"); nmod_poly_print(res2), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(exp);
        nmod_poly_clear(a);
        nmod_poly_clear(f);
        nmod_poly_clear(finv);
        nmod_poly_clear(res1);
        nmod_poly_clear(res2);
        nmod_poly_clear(t);
    }

    TEST_FUNCTION_END(state);
}
