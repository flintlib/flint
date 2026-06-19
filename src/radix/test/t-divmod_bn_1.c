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

TEST_FUNCTION_START(radix_divmod_bn_1, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        radix_t radix;
        nn_ptr a, q, qk, prod, chk;
        ulong b, rem, remk;
        slong an, n, i;
        int ret, retk;

        radix_init_randtest(radix, state);

        an = 1 + n_randint(state, 50);
        n = 1 + n_randint(state, 90);

        a = flint_malloc(an * sizeof(ulong));
        q = flint_malloc(n * sizeof(ulong));
        qk = flint_malloc(n * sizeof(ulong));
        prod = flint_malloc((n + 1) * sizeof(ulong));
        chk = flint_malloc((n + 1) * sizeof(ulong));

        radix_randtest_limbs(a, state, an, radix);
        {
            nn_ptr bb = flint_malloc(sizeof(ulong));
            radix_randtest_limbs(bb, state, 1, radix);
            b = bb[0];
            flint_free(bb);
        }

        ret = radix_divmod_bn_1(q, &rem, a, an, b, n, radix);

        /* Cross-check against the independent Karp-Markstein routine with
           bn = 1 (which does not delegate to radix_divmod_bn_1). */
        {
            ulong bb = b;
            retk = radix_divmod_bn_karp_markstein(qk, &remk, a, an, &bb, 1, n, radix);
        }

        if (ret != retk)
        {
            flint_printf("FAIL: return %d != Karp-Markstein %d\n", ret, retk);
            flint_abort();
        }

        if (ret == 0)
        {
            if (n_gcd(b, LIMB_RADIX(radix)) == 1)
            {
                flint_printf("FAIL: reported non-invertible but b is a unit\n");
                flint_abort();
            }
            goto cleanup;
        }

        for (i = 0; i < n; i++)
        {
            if (q[i] != qk[i])
            {
                flint_printf("FAIL: quotient differs from Karp-Markstein at limb %wd\n", i);
                flint_printf("radix %wu ^ %u, an=%wd n=%wd\n", DIGIT_RADIX(radix), radix->exp, an, n);
                flint_abort();
            }
        }

        if (rem != remk)
        {
            flint_printf("FAIL: remainder differs from Karp-Markstein\n");
            flint_abort();
        }

        /* Direct Hensel property: q*b == a (mod B^n). */
        {
            ulong bb = b;
            radix_mulmid(prod, q, n, &bb, 1, 0, n + 1, radix);
        }
        for (i = 0; i < n; i++)
        {
            ulong ai = (i < an) ? a[i] : 0;
            if (prod[i] != ai)
            {
                flint_printf("FAIL: q*b != a (mod B^n) at limb %wd\n", i);
                flint_abort();
            }
        }

        /* Remainder relation: a == q*b + B^n*rem (mod B^{n+1}). */
        flint_mpn_copyi(chk, prod, n + 1);
        {
            ulong rr = rem;
            radix_add(chk + n, chk + n, 1, &rr, 1, radix);
        }
        for (i = 0; i < n + 1; i++)
        {
            ulong ai = (i < an) ? a[i] : 0;
            if (chk[i] != ai)
            {
                flint_printf("FAIL: remainder relation fails at limb %wd\n", i);
                flint_abort();
            }
        }

        /* rem == NULL accepted; quotient unchanged. */
        {
            nn_ptr q2 = flint_malloc(n * sizeof(ulong));
            int ret2 = radix_divmod_bn_1(q2, NULL, a, an, b, n, radix);
            if (ret2 != ret)
            {
                flint_printf("FAIL: rem==NULL changed the return value\n");
                flint_abort();
            }
            for (i = 0; i < n; i++)
            {
                if (q2[i] != q[i])
                {
                    flint_printf("FAIL: rem==NULL changed the quotient at limb %wd\n", i);
                    flint_abort();
                }
            }
            flint_free(q2);
        }

        /* Aliasing: q may alias a. */
        {
            slong copyn = FLINT_MAX(an, n);
            nn_ptr qa = flint_malloc(copyn * sizeof(ulong));
            int reta;

            flint_mpn_zero(qa, copyn);
            flint_mpn_copyi(qa, a, an);
            reta = radix_divmod_bn_1(qa, NULL, qa, an, b, n, radix);

            if (reta != ret)
            {
                flint_printf("FAIL: aliasing q==a changed the return value\n");
                flint_abort();
            }
            for (i = 0; i < n; i++)
            {
                if (qa[i] != q[i])
                {
                    flint_printf("FAIL: aliased q differs at limb %wd\n", i);
                    flint_abort();
                }
            }
            flint_free(qa);
        }

cleanup:
        radix_clear(radix);
        flint_free(a);
        flint_free(q);
        flint_free(qk);
        flint_free(prod);
        flint_free(chk);
    }

    TEST_FUNCTION_END(state);
}
