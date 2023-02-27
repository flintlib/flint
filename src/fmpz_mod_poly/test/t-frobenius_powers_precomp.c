/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
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
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    fmpz_mod_ctx_t ctx;
    FLINT_TEST_INIT(state);

    flint_printf("frobenius_powers_precomp....");
    fflush(stdout);

    fmpz_mod_ctx_init_ui(ctx, 2);

    /* Aliasing of res and f */
    for (i = 0; i < 25 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t res, f, finv;
        fmpz_mod_poly_frobenius_powers_2exp_t pow;
        fmpz_t p;
        ulong exp;
        
        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, p);
        exp = n_randint(state, 50);
        
        fmpz_mod_poly_init(f, ctx);
        fmpz_mod_poly_init(finv, ctx);
        fmpz_mod_poly_init(res, ctx);

        fmpz_mod_poly_randtest_not_zero(f, state, n_randint(state, 50) + 1, ctx);

        fmpz_mod_poly_reverse (finv, f, f->length, ctx);
        fmpz_mod_poly_inv_series_newton (finv, finv, f->length, ctx);

        fmpz_mod_poly_frobenius_powers_2exp_precomp(pow, f, finv, exp, ctx);

        fmpz_mod_poly_frobenius_power(res, pow, f, exp, ctx);
        fmpz_mod_poly_frobenius_power(f, pow, f, exp, ctx);

        result = (fmpz_mod_poly_equal(res, f, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("f:\n"); fmpz_mod_poly_print(f, ctx), flint_printf("\n\n");
            flint_printf("res:\n"); fmpz_mod_poly_print(res, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(p);
        fmpz_mod_poly_clear(f, ctx);
        fmpz_mod_poly_clear(finv, ctx);
        fmpz_mod_poly_clear(res, ctx);
        fmpz_mod_poly_frobenius_powers_2exp_clear(pow, ctx);
    }

    /* Compare powers_2exp and powers */
    for (i = 0; i < 25 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t res, f, finv;
        fmpz_mod_poly_frobenius_powers_t pow;
        fmpz_mod_poly_frobenius_powers_2exp_t pow2;
        fmpz_t p;
        ulong exp, exp2;
        
        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, p);
        exp = n_randint(state, 50) + 1;
        exp2 = n_randint(state, exp);

        fmpz_mod_poly_init(f, ctx);
        fmpz_mod_poly_init(finv, ctx);
        fmpz_mod_poly_init(res, ctx);

        fmpz_mod_poly_randtest_not_zero(f, state, n_randint(state, 50) + 1, ctx);

        fmpz_mod_poly_reverse (finv, f, f->length, ctx);
        fmpz_mod_poly_inv_series_newton (finv, finv, f->length, ctx);

        fmpz_mod_poly_frobenius_powers_precomp(pow, f, finv, exp, ctx);
        fmpz_mod_poly_frobenius_powers_2exp_precomp(pow2, f, finv, exp, ctx);
        
        fmpz_mod_poly_frobenius_power(res, pow2, f, exp2, ctx);
        
        result = (fmpz_mod_poly_equal(res, pow->pow + exp2, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("pow->pow + exp2:\n"); fmpz_mod_poly_print(pow->pow + exp2, ctx), flint_printf("\n\n");
            flint_printf("res:\n"); fmpz_mod_poly_print(res, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(p);
        fmpz_mod_poly_clear(f, ctx);
        fmpz_mod_poly_clear(finv, ctx);
        fmpz_mod_poly_clear(res, ctx);
        fmpz_mod_poly_frobenius_powers_clear(pow, ctx);
        fmpz_mod_poly_frobenius_powers_2exp_clear(pow2, ctx);
    }

    fmpz_mod_ctx_clear(ctx);
    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
