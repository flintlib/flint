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

TEST_FUNCTION_START(radix_divmod_bn_karp_markstein, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        radix_t radix;
        nn_ptr a, b, q, qc, rem, remc;
        slong an, bn, n, i;
        int ret, retc;

        radix_init_randtest(radix, state);

        an = 1 + n_randint(state, 60);
        bn = 1 + n_randint(state, 60);   /* large bn is the Karp-Markstein regime */
        n = 1 + n_randint(state, 80);

        a = flint_malloc(an * sizeof(ulong));
        b = flint_malloc(bn * sizeof(ulong));
        q = flint_malloc(n * sizeof(ulong));
        qc = flint_malloc(n * sizeof(ulong));
        rem = flint_malloc(bn * sizeof(ulong));
        remc = flint_malloc(bn * sizeof(ulong));

        radix_randtest_limbs(a, state, an, radix);
        radix_randtest_limbs(b, state, bn, radix);

        ret = radix_divmod_bn_karp_markstein(q, rem, a, an, b, bn, n, radix);
        retc = radix_divmod_bn_classical(qc, remc, a, an, b, bn, n, radix);

        /* Must agree with the independent classical implementation on the
           invertibility verdict, the quotient, and the remainder. */
        if (ret != retc)
        {
            flint_printf("FAIL: return %d != classical %d\n", ret, retc);
            flint_printf("radix %wu ^ %u\n", DIGIT_RADIX(radix), radix->exp);
            flint_printf("an=%wd bn=%wd n=%wd\n", an, bn, n);
            flint_abort();
        }

        if (ret == 0)
        {
            if (n_gcd(b[0], LIMB_RADIX(radix)) == 1)
            {
                flint_printf("FAIL: reported non-invertible but b[0] is a unit\n");
                flint_abort();
            }
            goto cleanup;
        }

        for (i = 0; i < n; i++)
        {
            if (q[i] != qc[i])
            {
                flint_printf("FAIL: quotient differs from classical at limb %wd\n", i);
                flint_printf("radix %wu ^ %u\n", DIGIT_RADIX(radix), radix->exp);
                flint_printf("an=%wd bn=%wd n=%wd\n", an, bn, n);
                flint_printf("a = %{ulong*}\n", a, an);
                flint_printf("b = %{ulong*}\n", b, bn);
                flint_abort();
            }
        }

        for (i = 0; i < bn; i++)
        {
            if (rem[i] != remc[i])
            {
                flint_printf("FAIL: remainder differs from classical at limb %wd\n", i);
                flint_printf("radix %wu ^ %u\n", DIGIT_RADIX(radix), radix->exp);
                flint_printf("an=%wd bn=%wd n=%wd\n", an, bn, n);
                flint_abort();
            }
        }

        /* rem == NULL accepted; quotient unchanged. */
        {
            nn_ptr q2 = flint_malloc(n * sizeof(ulong));
            int ret2 = radix_divmod_bn_karp_markstein(q2, NULL, a, an, b, bn, n, radix);
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
            reta = radix_divmod_bn_karp_markstein(qa, NULL, qa, an, b, bn, n, radix);

            if (reta != ret)
            {
                flint_printf("FAIL: aliasing q==a changed the return value\n");
                flint_printf("an=%wd bn=%wd n=%wd\n", an, bn, n);
                flint_abort();
            }
            for (i = 0; i < n; i++)
            {
                if (qa[i] != q[i])
                {
                    flint_printf("FAIL: aliased q differs at limb %wd\n", i);
                    flint_printf("an=%wd bn=%wd n=%wd\n", an, bn, n);
                    flint_abort();
                }
            }
            flint_free(qa);
        }

cleanup:
        radix_clear(radix);
        flint_free(a);
        flint_free(b);
        flint_free(q);
        flint_free(qc);
        flint_free(rem);
        flint_free(remc);
    }

    TEST_FUNCTION_END(state);
}
