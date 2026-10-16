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

TEST_FUNCTION_START(radix_divmod_bn_classical, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        radix_t radix;
        nn_ptr a, b, q, rem, prod, chk;
        slong an, bn, n, remn, i;
        int ret, exact_ref, exact_caller;

        radix_init_randtest(radix, state);

        an = 1 + n_randint(state, 50);
        bn = 1 + n_randint(state, 30);
        n = 1 + n_randint(state, 60);

        a = flint_malloc(an * sizeof(ulong));
        b = flint_malloc(bn * sizeof(ulong));
        q = flint_malloc(n * sizeof(ulong));
        rem = flint_malloc(bn * sizeof(ulong));
        prod = flint_malloc((n + bn) * sizeof(ulong));
        chk = flint_malloc((n + bn) * sizeof(ulong));

        radix_randtest_limbs(a, state, an, radix);
        radix_randtest_limbs(b, state, bn, radix);

        ret = radix_divmod_bn_classical(q, rem, a, an, b, bn, n, radix);

        if (ret == 0)
        {
            /* must be a genuine non-unit leading limb */
            if (n_gcd(b[0], LIMB_RADIX(radix)) == 1)
            {
                flint_printf("FAIL: reported non-invertible but b[0] is a unit\n");
                flint_printf("radix %wu ^ %u = %wu\n", DIGIT_RADIX(radix), radix->exp, LIMB_RADIX(radix));
                flint_printf("b = %{ulong*}\n", b, bn);
                flint_abort();
            }
            goto cleanup;
        }

        if (n_gcd(b[0], LIMB_RADIX(radix)) != 1)
        {
            flint_printf("FAIL: should have reported non-invertible\n");
            flint_printf("radix %wu ^ %u = %wu\n", DIGIT_RADIX(radix), radix->exp, LIMB_RADIX(radix));
            flint_printf("b = %{ulong*}\n", b, bn);
            flint_abort();
        }

        /* prod = q * b (length n + bn) */
        radix_mulmid(prod, q, n, b, bn, 0, n + bn, radix);

        /* (1) Hensel property: q * b == a (mod B^n) */
        for (i = 0; i < n; i++)
        {
            ulong av = (i < an) ? a[i] : 0;
            if (prod[i] != av)
            {
                flint_printf("FAIL: q*b != a (mod B^n) at limb %wd\n", i);
                flint_printf("radix %wu ^ %u = %wu\n", DIGIT_RADIX(radix), radix->exp, LIMB_RADIX(radix));
                flint_printf("an=%wd bn=%wd n=%wd\n", an, bn, n);
                flint_printf("a = %{ulong*}\n", a, an);
                flint_printf("b = %{ulong*}\n", b, bn);
                flint_printf("q = %{ulong*}\n", q, n);
                flint_abort();
            }
        }

        /* (2) the documented caller-side exactness recipe (two zero tests)
               must agree with the reference test prod == a. */
        exact_ref = 1;
        {
            slong cmpn = FLINT_MAX(n + bn, an);
            for (i = 0; i < cmpn; i++)
            {
                ulong av = (i < an) ? a[i] : 0;
                ulong pv = (i < n + bn) ? prod[i] : 0;
                if (av != pv)
                {
                    exact_ref = 0;
                    break;
                }
            }
        }

        {
            slong hi = n + bn;
            int rz = flint_mpn_zero_p(rem, bn);
            int az = (an <= hi) ? 1 : flint_mpn_zero_p(a + hi, an - hi);
            exact_caller = rz && az;
        }

        if (exact_caller != exact_ref)
        {
            flint_printf("FAIL: caller exactness %d != reference %d\n", exact_caller, exact_ref);
            flint_printf("radix %wu ^ %u = %wu\n", DIGIT_RADIX(radix), radix->exp, LIMB_RADIX(radix));
            flint_printf("an=%wd bn=%wd n=%wd\n", an, bn, n);
            flint_printf("a = %{ulong*}\n", a, an);
            flint_printf("b = %{ulong*}\n", b, bn);
            flint_printf("q = %{ulong*}\n", q, n);
            flint_printf("rem = %{ulong*}\n", rem, bn);
            flint_abort();
        }

        /* (3) resume remainder: a == q*b + B^n * rem  (mod B^{n+bn}).
               Normalise rem here (the routine returns it unnormalised). */
        remn = bn;
        while (remn > 0 && rem[remn - 1] == 0)
            remn--;

        flint_mpn_copyi(chk, prod, n + bn);
        if (remn > 0)
            radix_add(chk + n, chk + n, bn, rem, remn, radix);

        for (i = 0; i < n + bn; i++)
        {
            ulong av = (i < an) ? a[i] : 0;
            if (chk[i] != av)
            {
                flint_printf("FAIL: a != q*b + B^n*rem (mod B^{n+bn}) at limb %wd\n", i);
                flint_printf("radix %wu ^ %u = %wu\n", DIGIT_RADIX(radix), radix->exp, LIMB_RADIX(radix));
                flint_printf("an=%wd bn=%wd n=%wd\n", an, bn, n);
                flint_printf("a = %{ulong*}\n", a, an);
                flint_printf("b = %{ulong*}\n", b, bn);
                flint_printf("q = %{ulong*}\n", q, n);
                flint_printf("rem = %{ulong*}\n", rem, bn);
                flint_abort();
            }
        }

        /* (4) rem == NULL accepted; same return value and same quotient. */
        {
            nn_ptr q2 = flint_malloc(n * sizeof(ulong));
            int ret2 = radix_divmod_bn_classical(q2, NULL, a, an, b, bn, n, radix);
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

        /* (5) aliasing: q may alias a (in-place). The result must match. */
        {
            slong copyn = FLINT_MAX(an, n);
            nn_ptr qa = flint_malloc(copyn * sizeof(ulong));
            int reta;

            flint_mpn_zero(qa, copyn);
            flint_mpn_copyi(qa, a, an);
            reta = radix_divmod_bn_classical(qa, NULL, qa, an, b, bn, n, radix);

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
        flint_free(rem);
        flint_free(prod);
        flint_free(chk);
    }

    TEST_FUNCTION_END(state);
}
