/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"
#include "radix.h"

TEST_FUNCTION_START(radix_divmod_bn, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        radix_t radix;
        nn_ptr a, b, qd, qc, remd, remc;
        slong an, bn, n, i;
        int rd, rc;

        radix_init_randtest(radix, state);

        an = 1 + n_randint(state, 60);
        bn = 1 + n_randint(state, 60);
        n = 1 + n_randint(state, 80);

        a = flint_malloc(an * sizeof(ulong));
        b = flint_malloc(bn * sizeof(ulong));
        qd = flint_malloc(n * sizeof(ulong));
        qc = flint_malloc(n * sizeof(ulong));
        remd = flint_malloc(bn * sizeof(ulong));
        remc = flint_malloc(bn * sizeof(ulong));

        radix_randtest_limbs(a, state, an, radix);
        radix_randtest_limbs(b, state, bn, radix);

        rd = radix_divmod_bn(qd, remd, a, an, b, bn, n, radix);
        rc = radix_divmod_bn_classical(qc, remc, a, an, b, bn, n, radix);

        /* The dispatcher forwards to classical or Karp-Markstein, which are
           proven to agree; either way the result must equal classical. */
        if (rd != rc)
        {
            flint_printf("FAIL: dispatcher return %d != classical %d\n", rd, rc);
            flint_printf("bn=%wd n=%wd\n", bn, n);
            flint_abort();
        }

        if (rd)
        {
            for (i = 0; i < n; i++)
            {
                if (qd[i] != qc[i])
                {
                    flint_printf("FAIL: dispatcher quotient differs at limb %wd\n", i);
                    flint_printf("bn=%wd n=%wd\n", bn, n);
                    flint_abort();
                }
            }
            for (i = 0; i < bn; i++)
            {
                if (remd[i] != remc[i])
                {
                    flint_printf("FAIL: dispatcher remainder differs at limb %wd\n", i);
                    flint_abort();
                }
            }
        }

        radix_clear(radix);
        flint_free(a);
        flint_free(b);
        flint_free(qd);
        flint_free(qc);
        flint_free(remd);
        flint_free(remc);
    }

    TEST_FUNCTION_END(state);
}
