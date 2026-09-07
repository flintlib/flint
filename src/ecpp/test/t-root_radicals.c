/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.1

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "ecpp.h"

TEST_FUNCTION_START(ecpp_root_radicals, state)
{
    slong iter;

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        fmpz_t n, r, x;
        fmpz_mod_ctx_t ctx;
        fmpz_mod_poly_t f, lin;
        slong d = 2 + n_randint(state, 3), i;
        int ok;

        fmpz_init(n); fmpz_init(r); fmpz_init(x);
        do
            fmpz_randprime(n, state, 20 + n_randint(state, 200), 0);
        while (d >= 3 && fmpz_fdiv_ui(n, 9) % 3 != 1 && fmpz_fdiv_ui(n, 9) != 2 && fmpz_fdiv_ui(n, 9) != 5);
        fmpz_mod_ctx_init(ctx, n);
        fmpz_mod_poly_init(f, ctx);
        fmpz_mod_poly_init(lin, ctx);

        fmpz_mod_poly_one(f, ctx);
        for (i = 0; i < d; i++)
        {
            fmpz_randm(r, state, n);
            if (n_randint(state, 4) == 0)
                fmpz_set_ui(r, n_randint(state, 3));    /* repeated / small roots */
            fmpz_mod_poly_set_coeff_ui(lin, 1, 1, ctx);
            fmpz_mod_neg(r, r, ctx);
            fmpz_mod_poly_set_coeff_fmpz(lin, 0, r, ctx);
            fmpz_mod_poly_mul(f, f, lin, ctx);
        }

        ok = ecpp_root_radicals(x, f, state, ctx);
        if (!ok)
        {
            flint_printf("FAIL: no root found, degree %wd\n", d);
            flint_printf("n = "); fmpz_print(n); flint_printf("\n");
            fmpz_mod_poly_print(f, ctx); flint_printf("\n");
            flint_abort();
        }
        fmpz_mod_poly_evaluate_fmpz(r, f, x, ctx);
        if (!fmpz_is_zero(r))
        {
            flint_printf("FAIL: not a root, degree %wd\n", d);
            flint_abort();
        }

        fmpz_mod_poly_clear(f, ctx);
        fmpz_mod_poly_clear(lin, ctx);
        fmpz_mod_ctx_clear(ctx);
        fmpz_clear(n); fmpz_clear(r); fmpz_clear(x);
    }

    TEST_FUNCTION_END(state);
}
