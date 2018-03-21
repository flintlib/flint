/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, j, success;
    FLINT_TEST_INIT(state);

    flint_printf("gcd_prs....");
    fflush(stdout);


    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        int ok;
        fmpz_t ac, bc, d;
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, c, g, t;
        ordering_t ord;
        slong nvars, len1, len2, len3, exp_bound1, exp_bound2, exp_bound3;
        slong coeff_bits;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 3) + 1;

        fmpz_mpoly_ctx_init(ctx, nvars, ord);

        fmpz_init(ac);
        fmpz_init(bc);
        fmpz_init(d);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(c, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(t, ctx);

        for (j = 0; j < 4; j++)
        {
            len1 = n_randint(state, 18);
            len2 = n_randint(state, 18);
            len3 = n_randint(state, 14);
            exp_bound1 = n_randint(state, 4) + 2;
            exp_bound2 = n_randint(state, 4) + 2;
            exp_bound3 = n_randint(state, 3) + 2;
            coeff_bits = n_randint(state, 3) + 1;
            fmpz_mpoly_randtest_bound(a, state, len1, coeff_bits, exp_bound1, ctx);
            fmpz_mpoly_randtest_bound(b, state, len2, coeff_bits, exp_bound2, ctx);
            fmpz_mpoly_randtest_bound(c, state, len3, coeff_bits, exp_bound3, ctx);

            fmpz_mpoly_mul_johnson(a, a, c, ctx);
            fmpz_mpoly_mul_johnson(b, b, c, ctx);

            success = fmpz_mpoly_gcd_prs(g, a, b, ctx);
            if (!success)
                continue;

            fmpz_mpoly_assert_canonical(g, ctx);

            if (fmpz_mpoly_is_zero(g, ctx))
            {
                if (!fmpz_mpoly_is_zero(a, ctx) || !fmpz_mpoly_is_zero(b, ctx))
                {
                    flint_printf("FAIL\ngcd is zero but both inputs are not\n"
                                                     "i: %wd  j: %wd\n", i, j);
                    flint_abort();
                }
                continue;
            }

            ok = fmpz_mpoly_divides_monagan_pearce(a, a, g, ctx);
            ok = ok && fmpz_mpoly_divides_monagan_pearce(b, b, g, ctx);
            if (!ok)
            {
                flint_printf("FAIL\ngcd doesn't divide both inputs\n"
                                                     "i: %wd  j: %wd\n", i, j);
                flint_abort();
            }

            if (!fmpz_mpoly_gcd_is_unit(a, b, ctx))
                flint_printf("FAIL\ncofactors are not relatively prime\n"
                                                     "i: %wd  j: %wd\n", i, j);

        }

        fmpz_clear(ac);
        fmpz_clear(bc);
        fmpz_clear(d);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(c, ctx);
        fmpz_mpoly_clear(g, ctx);       
        fmpz_mpoly_clear(t, ctx);
    }

    FLINT_TEST_CLEANUP(state);

    printf("PASS\n");
    return 0;
}

