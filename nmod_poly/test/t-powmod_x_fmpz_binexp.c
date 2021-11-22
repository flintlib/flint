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

#undef ulong
#define ulong ulongxx/* interferes with system includes */

#include <stdlib.h>
#include <stdio.h>

#undef ulong

#include <gmp.h>

#define ulong mp_limb_t

#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);
    
    flint_printf("powmod_x_fmpz_preinv....");
    fflush(stdout);

    /* No aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, res1, res2, t, f, finv;
        mp_limb_t n;
        fmpz_t exp;

        fmpz_init(exp);

        n = n_randtest_prime(state, 0);
        fmpz_randtest_unsigned(exp, state, n_randint(state, 100) + 1);

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

        nmod_poly_init2(a, n, f->length-1);
        nmod_poly_set_coeff_ui(a, 1, 1);

        nmod_poly_powmod_x_fmpz_preinv(res1, exp, f, finv);
        nmod_poly_powmod_fmpz_binexp_preinv(res2, a, exp, f, finv);

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

    /* Aliasing of res and f */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        nmod_poly_t res1, t, f, finv;
        mp_limb_t n;
        fmpz_t exp;

        fmpz_init(exp);

        n = n_randtest_prime(state, 0);
        fmpz_randtest_unsigned(exp, state, n_randint(state, 100) + 1);

        nmod_poly_init(f, n);
        nmod_poly_init(finv, n);
        nmod_poly_init(res1, n);
        nmod_poly_init(t, n);

        do {
            nmod_poly_randtest(f, state, n_randint(state, 50));
        } while (nmod_poly_is_zero(f));

        nmod_poly_reverse(finv, f, f->length);
        nmod_poly_inv_series(finv, finv, f->length);

        nmod_poly_powmod_x_fmpz_preinv(res1, exp, f, finv);
        nmod_poly_powmod_x_fmpz_preinv(f, exp, f, finv);

        result = (nmod_poly_equal(res1, f));
        if (!result)
        {
            flint_printf("FAIL: aliasing1\n");
            flint_printf("exp: "); fmpz_print(exp); flint_printf("\n\n");
            flint_printf("f:\n"); nmod_poly_print(f), flint_printf("\n\n");
            flint_printf("res1:\n"); nmod_poly_print(res1), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(exp);
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
        fmpz_t exp;

        fmpz_init(exp);

        n = n_randtest_prime(state, 0);
        fmpz_randtest_unsigned(exp, state, n_randint(state, 100) + 1);

        nmod_poly_init(f, n);
        nmod_poly_init(finv, n);
        nmod_poly_init(res1, n);
        nmod_poly_init(t, n);

        do {
            nmod_poly_randtest(f, state, n_randint(state, 50));
        } while (nmod_poly_is_zero(f));

        nmod_poly_reverse(finv, f, f->length);
        nmod_poly_inv_series(finv, finv, f->length);

        nmod_poly_powmod_x_fmpz_preinv(res1, exp, f, finv);
        nmod_poly_powmod_x_fmpz_preinv(finv, exp, f, finv);

        result = (nmod_poly_equal(res1, finv));
        if (!result)
        {
            flint_printf("FAIL: aliasing2\n");
            flint_printf("exp: "); fmpz_print(exp); flint_printf("\n\n");
            flint_printf("f:\n"); nmod_poly_print(f), flint_printf("\n\n");
            flint_printf("finv:\n"); nmod_poly_print(finv), flint_printf("\n\n");
            flint_printf("res1:\n"); nmod_poly_print(res1), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(exp);
        nmod_poly_clear(f);
        nmod_poly_clear(finv);
        nmod_poly_clear(res1);
        nmod_poly_clear(t);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
