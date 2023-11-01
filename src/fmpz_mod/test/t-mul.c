/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mod.h"

TEST_FUNCTION_START(fmpz_mod_mul, state)
{
    flint_bitcnt_t max_bits = 200, max_bits_huge = 500 * FLINT_BITS;
    slong i, j, k;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t p;   /* p not nec prime */
        fmpz_t a, b, c, d, ab, bc;
        fmpz_mod_ctx_t fpctx;

        fmpz_init_set_ui(p, 2);
        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(d);
        fmpz_init(ab);
        fmpz_init(bc);
        fmpz_mod_ctx_init(fpctx, p);

        for (j = 0; j < 10; j++)
        {
            if (j == 0)
            {
                /* exact powers of 2 */
                fmpz_one(p);
                fmpz_mul_2exp(p, p, n_randint(state, max_bits));
            }
            else if (j == 1)
            {
                /* sum of two powers of 2 */
                fmpz_one(p);
                fmpz_mul_2exp(p, p, n_randint(state, max_bits/2));
                fmpz_add_ui(p, p, 1);
                fmpz_mul_2exp(p, p, n_randint(state, max_bits/2));
            }
            else if (j == 2)
            {
                /* one less than a power of 2 */
                fmpz_one(p);
                fmpz_mul_2exp(p, p, 1 + n_randint(state, max_bits));
                fmpz_sub_ui(p, p, 1);
            }
            else if (j == 3 && n_randint(state, 30) == 0)
            {
                fmpz_randtest_unsigned(p, state, max_bits_huge);
                fmpz_add_ui(p, p, 1);
            }
            else
            {
                fmpz_randtest_unsigned(p, state, max_bits);
                fmpz_add_ui(p, p, 1);
            }

            fmpz_mod_ctx_set_modulus(fpctx, p);

            for (k = 0; k < 10; k++)
            {
                if (k == 0)
                {
                    fmpz_sub_ui(a, p, 1);
                    fmpz_sub_ui(b, p, 1);
                    fmpz_sub_ui(c, p, 1);
                    fmpz_sub_ui(d, p, 1);
                    fmpz_sub_ui(ab, p, 1);
                    fmpz_sub_ui(bc, p, 1);
                }
                else
                {
                    fmpz_randtest_mod(a, state, p);
                    fmpz_randtest_mod(b, state, p);
                    fmpz_randtest_mod(c, state, p);
                    fmpz_randtest_mod(d, state, p);
                    fmpz_randtest_mod(ab, state, p);
                    fmpz_randtest_mod(bc, state, p);
                }

                fmpz_mul(ab, a, b);
                fmpz_mod(ab, ab, p);
                fmpz_mul(bc, b, c);
                fmpz_mod(bc, bc, p);

                fmpz_mod_mul(d, a, b, fpctx);
                if (!fmpz_equal(d, ab))
                {
                    printf("FAIL\n");
                    flint_printf("i = %wd, j = %wd, k = %wd\n", i, j, k);
                    fflush(stdout);
                    flint_abort();
                }

                fmpz_mod_mul(a, a, b, fpctx);
                if (!fmpz_equal(a, ab))
                {
                    printf("FAIL\ncheck aliasing first");
                    flint_printf("i = %wd, j = %wd, k = %wd\n", i, j, k);
                    fflush(stdout);
                    flint_abort();
                }

                fmpz_mod_mul(c, b, c, fpctx);
                if (!fmpz_equal(c, bc))
                {
                    printf("FAIL\ncheck aliasing second");
                    flint_printf("i = %wd, j = %wd, k = %wd\n", i, j, k);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        fmpz_mod_ctx_clear(fpctx);
        fmpz_clear(p);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(d);
        fmpz_clear(ab);
        fmpz_clear(bc);
    }

    TEST_FUNCTION_END(state);
}
