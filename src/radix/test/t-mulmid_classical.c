/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "radix.h"

TEST_FUNCTION_START(radix_mulmid_classical, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        radix_t radix;
        nn_ptr a, b, c, d;
        slong an, bn, lo, hi;
        int squaring;

        radix_init_randtest(radix, state);

        an = 1 + n_randint(state, 30);
        bn = 1 + n_randint(state, 30);

        squaring = n_randint(state, 2);
        if (squaring)
            an = bn;

        lo = n_randint(state, an + bn);
        hi = lo + 1 + n_randint(state, an + bn - lo);

        a = flint_malloc(an * sizeof(ulong));
        b = flint_malloc(bn * sizeof(ulong));
        c = flint_malloc((hi - lo) * sizeof(ulong));
        d = flint_malloc((hi - lo) * sizeof(ulong));

        radix_randtest_limbs(a, state, an, radix);
        radix_randtest_limbs(b, state, bn, radix);
        radix_randtest_limbs(c, state, hi - lo, radix);
        radix_randtest_limbs(d, state, hi - lo, radix);

        if (squaring)
        {
            radix_mulmid_classical(c, a, an, a, an, lo, hi, radix);
            radix_mulmid_naive(d, a, an, a, an, lo, hi, radix);
        }
        else
        {
            radix_mulmid_classical(c, a, an, b, bn, lo, hi, radix);
            radix_mulmid_naive(d, a, an, b, bn, lo, hi, radix);
        }

        if (mpn_cmp(c, d, hi - lo) != 0)
        {
            flint_printf("FAIL: mul\n");
            flint_printf("an = %wd, bn = %wd, lo = %wd, hi = %wd, squaring = %d\n", an, bn, lo, hi, squaring);
            flint_printf("%{ulong*}\n", a, an);
            flint_printf("%{ulong*}\n", b, bn);
            flint_printf("%{ulong*}\n", c, hi - lo);
            flint_printf("%{ulong*}\n", d, hi - lo);
            flint_abort();
        }

        radix_clear(radix);

        flint_free(a);
        flint_free(b);
        flint_free(c);
        flint_free(d);
    }

    /*
        Two-limb column sums at their largest values (limbs B - 1, sometimes
        B - 2), for B of nu > FLINT_BITS / 2 bits, around the length of the
        shorter operand where the Barrett step stops applying
        (2^(FLINT_BITS - 2 - nu) limbs); products and squarings, full
        products and windows
    */
    {
        flint_bitcnt_t nu;

        for (nu = FLINT_BITS / 2 + 1; nu <= FLINT_BITS; nu++)
        {
            int kind;

            for (kind = 0; kind < 3; kind++)
            {
                radix_t radix;
                ulong B;
                slong lens[6];
                int nlens = 0, t;

                if (kind == 0)
                    B = (nu == FLINT_BITS) ? UWORD_MAX : (UWORD(1) << nu) - 1;
                else if (kind == 1)
                    B = (UWORD(1) << (nu - 1)) + 1;
                else
                    B = n_randbits(state, nu);

                radix_init(radix, B, 1);

                lens[nlens++] = 2;
                lens[nlens++] = 3;
                if (nu + 2 <= FLINT_BITS && FLINT_BITS - 2 - nu <= 5)
                {
                    lens[nlens++] = WORD(1) << (FLINT_BITS - 2 - nu);
                    lens[nlens++] = (WORD(1) << (FLINT_BITS - 2 - nu)) + 1;
                }
                lens[nlens++] = 2 + n_randint(state, 40);

                for (t = 0; t < nlens; t++)
                {
                    nn_ptr a, b, c, d;
                    slong an, bn, lo, hi, i;
                    int squaring;

                    for (squaring = 0; squaring < 2; squaring++)
                    {
                        bn = lens[t];
                        an = squaring ? bn : bn + n_randint(state, 8);

                        if (n_randint(state, 2))
                        {
                            lo = 0;
                            hi = an + bn;
                        }
                        else
                        {
                            lo = n_randint(state, an + bn);
                            hi = lo + 1 + n_randint(state, an + bn - lo);
                        }

                        a = flint_malloc(an * sizeof(ulong));
                        b = flint_malloc(bn * sizeof(ulong));
                        c = flint_malloc((hi - lo) * sizeof(ulong));
                        d = flint_malloc((hi - lo) * sizeof(ulong));

                        for (i = 0; i < an; i++)
                            a[i] = n_randint(state, 8) ? B - 1 : B - 2;
                        for (i = 0; i < bn; i++)
                            b[i] = n_randint(state, 8) ? B - 1 : B - 2;

                        if (squaring)
                        {
                            radix_mulmid_classical(c, a, an, a, an, lo, hi, radix);
                            radix_mulmid_naive(d, a, an, a, an, lo, hi, radix);
                        }
                        else
                        {
                            radix_mulmid_classical(c, a, an, b, bn, lo, hi, radix);
                            radix_mulmid_naive(d, a, an, b, bn, lo, hi, radix);
                        }

                        if (mpn_cmp(c, d, hi - lo) != 0)
                        {
                            flint_printf("FAIL: largest column sums\n");
                            flint_printf("B = %wu, an = %wd, bn = %wd, lo = %wd, hi = %wd, squaring = %d\n",
                                         B, an, bn, lo, hi, squaring);
                            flint_printf("%{ulong*}\n", c, hi - lo);
                            flint_printf("%{ulong*}\n", d, hi - lo);
                            flint_abort();
                        }

                        flint_free(a);
                        flint_free(b);
                        flint_free(c);
                        flint_free(d);
                    }
                }

                radix_clear(radix);
            }
        }
    }

    TEST_FUNCTION_END(state);
}
