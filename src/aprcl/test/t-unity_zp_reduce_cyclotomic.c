/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "aprcl.h"

TEST_FUNCTION_START(aprcl_unity_zp_reduce_cyclotomic, state)
{
    int i, j;
    fmpz_mod_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctx, 2);

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t cyclo_poly;
        ulong p, exp;
        fmpz_t n;
        unity_zp f, g;

        p = n_randprime(state, 2 + n_randint(state, 2), 0);
        exp = n_randint(state, 4);
        while (exp == 0)
            exp = n_randint(state, 4);

        fmpz_init(n);
        fmpz_randtest_unsigned(n, state, 200);
        while (fmpz_cmp_ui(n, 2) < 0)
            fmpz_randtest_unsigned(n, state, 200);

        fmpz_mod_ctx_set_modulus(ctx, n);

        unity_zp_init(f, p, exp, n);
        unity_zp_init(g, p, exp, n);

        for (j = 0; j < 100; j++)
        {
            ulong ind;
            fmpz_t val;

            fmpz_init(val);

            ind = n_randint(state, n_pow(p, exp));

            fmpz_randtest_unsigned(val, state, 200);

            unity_zp_coeff_set_fmpz(f, ind, val);

            fmpz_clear(val);
        }

        unity_zp_reduce_cyclotomic(g, f);

        fmpz_mod_poly_init(cyclo_poly, ctx);
        for (j = 0; j < p; j++)
            fmpz_mod_poly_set_coeff_ui(cyclo_poly, j * n_pow(p, exp - 1), 1, ctx);
        fmpz_mod_poly_rem(f->poly, f->poly, cyclo_poly, ctx);

        if (unity_zp_equal(f, g) == 0)
        {
            flint_printf("FAIL\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(n);
        fmpz_mod_poly_clear(cyclo_poly, ctx);
        unity_zp_clear(f);
        unity_zp_clear(g);
    }

    fmpz_mod_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
}
