/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_mod.h"

TEST_FUNCTION_START(fmpz_mod_pow_ui, state)
{
    flint_bitcnt_t max_modulus_bits = 200;
    slong i, j;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p; /* p not nec. prime */
        fmpz_t a, b, u, v, w, x;
        ulong e, f, t;
        fmpz_mod_ctx_t fpctx;

        fmpz_init_set_ui(p, 2);
        fmpz_mod_ctx_init(fpctx, p);

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(u);
        fmpz_init(v);
        fmpz_init(w);
        fmpz_init(x);

        for (j = 0; j < 100; j++)
        {
            if (n_randint(state, 1000) == 0)
                fmpz_randtest_unsigned(p, state, 100000);
            else
                fmpz_randtest_unsigned(p, state, max_modulus_bits);

            fmpz_add_ui(p, p, 1);
            fmpz_mod_ctx_set_modulus(fpctx, p);

            fmpz_randtest_mod(a, state, p);
            fmpz_randtest_mod(b, state, p);
            e = n_randlimb(state);
            f = n_randlimb(state);
            t = e + f;
            if (t < e || t < f)
            {
                continue;
            }

            /* check a^(e + f) =  a^e * a^f */
            fmpz_mod_pow_ui(u, a, e, fpctx);
            fmpz_mod_pow_ui(v, a, f, fpctx);
            fmpz_mod_pow_ui(w, a, t, fpctx);
            fmpz_mod_mul(x, u, v, fpctx);
            if (!fmpz_equal(w, x))
            {
                printf("FAIL\n");
                flint_printf("check a^(e + f) =  a^e * a^f\n"
                                            "i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            /* check a^(e + f) =  a^e * a^f with aliasing */
            fmpz_mod_pow_ui(u, a, e, fpctx);
            fmpz_mod_pow_ui(v, a, f, fpctx);
            fmpz_mod_pow_ui(a, a, t, fpctx);
            fmpz_mod_mul(x, u, v, fpctx);
            if (!fmpz_equal(a, x))
            {
                printf("FAIL\n");
                flint_printf("check a^(e + f) =  a^e * a^f with aliasing\n"
                                            "i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            /* check 0^0 = 1 */
            fmpz_zero(a);
            fmpz_mod_pow_ui(u, a, 0, fpctx);
            if (!fmpz_mod_is_one(u, fpctx))
            {
                printf("FAIL\n");
                flint_printf("check 0^0 = 1\n"
                                            "i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            /* check (a*b)^e =  a^e * b^e with aliasing */
            fmpz_mod_mul(u, a, b, fpctx);
            fmpz_mod_pow_ui(a, a, e, fpctx);
            fmpz_mod_pow_ui(b, b, e, fpctx);
            fmpz_mod_pow_ui(u, u, e, fpctx);
            fmpz_mod_mul(a, a, b, fpctx);
            if (!fmpz_equal(a, u))
            {
                printf("FAIL\n");
                flint_printf("check (a*b)^e =  a^e * b^e with aliasing\n"
                                            "i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_clear(p);
        fmpz_mod_ctx_clear(fpctx);

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(u);
        fmpz_clear(v);
        fmpz_clear(w);
        fmpz_clear(x);
    }

    TEST_FUNCTION_END(state);
}
