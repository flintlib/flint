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

TEST_FUNCTION_START(nmod_poly_mulmod_preinv, state)
{
    int i, result;

    /* Aliasing res and a */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, res, t, f, finv;

        mp_limb_t n = n_randtest_prime(state, 0);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(f, n);
        nmod_poly_init(finv, n);
        nmod_poly_init(res, n);
        nmod_poly_init(t, n);

        nmod_poly_randtest(a, state, n_randint(state, 50));
        nmod_poly_randtest(b, state, n_randint(state, 50));
        do {
            nmod_poly_randtest(f, state, n_randint(state, 50));
        } while (nmod_poly_is_zero(f));
        if (a->length >= f->length)
          nmod_poly_rem (a, a, f);
        if (b->length >= f->length)
          nmod_poly_rem (b, b, f);

        nmod_poly_reverse(finv, f, f->length);
        nmod_poly_inv_series(finv, finv, f->length);

        nmod_poly_mulmod_preinv(res, a, b, f, finv);
        nmod_poly_mulmod_preinv(a, a, b, f, finv);

        result = (nmod_poly_equal(res, a));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a:\n"); nmod_poly_print(a), flint_printf("\n\n");
            flint_printf("b:\n"); nmod_poly_print(b), flint_printf("\n\n");
            flint_printf("f:\n"); nmod_poly_print(f), flint_printf("\n\n");
            flint_printf("res1:\n"); nmod_poly_print(res), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(f);
        nmod_poly_clear(finv);
        nmod_poly_clear(res);
        nmod_poly_clear(t);
    }

    /* Aliasing res and b */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, res, t, f, finv;

        mp_limb_t n = n_randtest_prime(state, 0);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(f, n);
        nmod_poly_init(finv, n);
        nmod_poly_init(res, n);
        nmod_poly_init(t, n);

        nmod_poly_randtest(a, state, n_randint(state, 50));
        nmod_poly_randtest(b, state, n_randint(state, 50));
        do {
            nmod_poly_randtest(f, state, n_randint(state, 50));
        } while (nmod_poly_is_zero(f));
        if (a->length >= f->length)
          nmod_poly_rem (a, a, f);
        if (b->length >= f->length)
          nmod_poly_rem (b, b, f);

        nmod_poly_reverse(finv, f, f->length);
        nmod_poly_inv_series(finv, finv, f->length);

        nmod_poly_mulmod_preinv(res, a, b, f, finv);
        nmod_poly_mulmod_preinv(b, a, b, f, finv);

        result = (nmod_poly_equal(res, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a:\n"); nmod_poly_print(a), flint_printf("\n\n");
            flint_printf("b:\n"); nmod_poly_print(b), flint_printf("\n\n");
            flint_printf("f:\n"); nmod_poly_print(f), flint_printf("\n\n");
            flint_printf("res1:\n"); nmod_poly_print(res), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(f);
        nmod_poly_clear(finv);
        nmod_poly_clear(res);
        nmod_poly_clear(t);
    }

    /* Aliasing res and f */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, res, t, f, finv;

        mp_limb_t n = n_randtest_prime(state, 0);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(f, n);
        nmod_poly_init(finv, n);
        nmod_poly_init(res, n);
        nmod_poly_init(t, n);

        nmod_poly_randtest(a, state, n_randint(state, 50));
        nmod_poly_randtest(b, state, n_randint(state, 50));
        do {
            nmod_poly_randtest(f, state, n_randint(state, 50));
        } while (nmod_poly_is_zero(f));
        if (a->length >= f->length)
          nmod_poly_rem (a, a, f);
        if (b->length >= f->length)
          nmod_poly_rem (b, b, f);

        nmod_poly_reverse(finv, f, f->length);
        nmod_poly_inv_series(finv, finv, f->length);

        nmod_poly_mulmod_preinv(res, a, b, f, finv);
        nmod_poly_mulmod_preinv(f, a, b, f, finv);

        result = (nmod_poly_equal(res, f));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a:\n"); nmod_poly_print(a), flint_printf("\n\n");
            flint_printf("b:\n"); nmod_poly_print(b), flint_printf("\n\n");
            flint_printf("f:\n"); nmod_poly_print(f), flint_printf("\n\n");
            flint_printf("res1:\n"); nmod_poly_print(res), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(f);
        nmod_poly_clear(finv);
        nmod_poly_clear(res);
        nmod_poly_clear(t);
    }

    /* Aliasing res and finv */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, res, t, f, finv;

        mp_limb_t n = n_randtest_prime(state, 0);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(f, n);
        nmod_poly_init(finv, n);
        nmod_poly_init(res, n);
        nmod_poly_init(t, n);

        nmod_poly_randtest(a, state, n_randint(state, 50));
        nmod_poly_randtest(b, state, n_randint(state, 50));
        do {
            nmod_poly_randtest(f, state, n_randint(state, 50));
        } while (nmod_poly_is_zero(f));
        if (a->length >= f->length)
          nmod_poly_rem (a, a, f);
        if (b->length >= f->length)
          nmod_poly_rem (b, b, f);

        nmod_poly_reverse(finv, f, f->length);
        nmod_poly_inv_series(finv, finv, f->length);

        nmod_poly_mulmod_preinv(res, a, b, f, finv);
        nmod_poly_mulmod_preinv(finv, a, b, f, finv);

        result = (nmod_poly_equal(res, finv));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a:\n"); nmod_poly_print(a), flint_printf("\n\n");
            flint_printf("b:\n"); nmod_poly_print(b), flint_printf("\n\n");
            flint_printf("f:\n"); nmod_poly_print(f), flint_printf("\n\n");
            flint_printf("finv:\n"); nmod_poly_print(finv), flint_printf("\n\n");
            flint_printf("res1:\n"); nmod_poly_print(res), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(f);
        nmod_poly_clear(finv);
        nmod_poly_clear(res);
        nmod_poly_clear(t);
    }

    /* No aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, res1, res2, t, f, finv;

        mp_limb_t n = n_randtest_prime(state, 0);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(f, n);
        nmod_poly_init(finv, n);
        nmod_poly_init(res1, n);
        nmod_poly_init(res2, n);
        nmod_poly_init(t, n);

        nmod_poly_randtest(a, state, n_randint(state, 50));
        nmod_poly_randtest(b, state, n_randint(state, 50));
        do {
            nmod_poly_randtest(f, state, n_randint(state, 50));
        } while (nmod_poly_is_zero(f));
        if (a->length >= f->length)
          nmod_poly_rem (a, a, f);
        if (b->length >= f->length)
          nmod_poly_rem (b, b, f);

        nmod_poly_reverse(finv, f, f->length);
        nmod_poly_inv_series(finv, finv, f->length);

        nmod_poly_mulmod_preinv(res1, a, b, f, finv);

        nmod_poly_mul(res2, a, b);
        nmod_poly_divrem(t, res2, res2, f);

        result = (nmod_poly_equal(res1, res2));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a:\n"); nmod_poly_print(a), flint_printf("\n\n");
            flint_printf("b:\n"); nmod_poly_print(b), flint_printf("\n\n");
            flint_printf("f:\n"); nmod_poly_print(f), flint_printf("\n\n");
            flint_printf("res1:\n"); nmod_poly_print(res1), flint_printf("\n\n");
            flint_printf("res2:\n"); nmod_poly_print(res2), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(f);
        nmod_poly_clear(finv);
        nmod_poly_clear(res1);
        nmod_poly_clear(res2);
        nmod_poly_clear(t);
    }

    TEST_FUNCTION_END(state);
}
