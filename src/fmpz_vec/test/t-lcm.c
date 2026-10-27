/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"

TEST_FUNCTION_START(fmpz_vec_lcm, state)
{
    int i, result;

    /* Check that lcm(a f) = abs(a) lcm(f) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, c, d;
        fmpz *f;
        slong len = n_randint(state, 100);

        fmpz_init(a);
        fmpz_init(c);
        fmpz_init(d);
        f = _fmpz_vec_init(len);
        _fmpz_vec_randtest(f, state, len, 200);
        fmpz_randtest(a, state, 100);

        _fmpz_vec_lcm(c, f, len);

        if (len == 0)
        {
            result = fmpz_is_one(c);
        }
        else
        {
            _fmpz_vec_scalar_mul_fmpz(f, f, len, a);
            fmpz_abs(a, a);
            fmpz_mul(c, a, c);
            _fmpz_vec_lcm(d, f, len);
            result = (fmpz_equal(c, d));
        }

        if (!result)
        {
            fmpz_print(c), flint_printf("\n\n");
            fmpz_print(d), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(c);
        fmpz_clear(d);
        _fmpz_vec_clear(f, len);
    }

    /* Compare with chained fmpz_lcm on denominator-like shapes, also
       strided and with an initial value, possibly aliased */
    for (i = 0; i < 3000 * flint_test_multiplier(); i++)
    {
        fmpz_t c, d, D, inp;
        fmpz * f;
        slong len, j, stride, shape;
        flint_bitcnt_t bits;

        len = n_randint(state, 3) == 0 ? n_randint(state, 300) : n_randint(state, 30);
        stride = 1 + n_randint(state, 3);
        bits = n_randint(state, 2) ? 1 + n_randint(state, 70) : 1 + n_randint(state, 400);
        shape = n_randint(state, 7);

        fmpz_init(c);
        fmpz_init(d);
        fmpz_init(D);
        fmpz_init(inp);
        f = _fmpz_vec_init(len * stride + 1);

        fmpz_randtest_not_zero(D, state, bits);
        for (j = 0; j < len; j++)
        {
            fmpz * e = f + j * stride;

            switch (shape)
            {
                case 0: /* small */
                    fmpz_set_si(e, n_randint(state, 300) - 150);
                    break;
                case 1: /* equal up to sign, or 1 */
                    if (n_randint(state, 3))
                        fmpz_set(e, D);
                    else
                        fmpz_one(e);
                    if (n_randint(state, 2))
                        fmpz_neg(e, e);
                    break;
                case 2: /* divisors of D */
                    fmpz_gcd(e, D, f + (j / 2) * stride);
                    if (fmpz_is_zero(e))
                        fmpz_set_ui(e, n_randint(state, 10) + 1);
                    break;
                case 3: /* growing, like factorials */
                    if (j == 0)
                        fmpz_set_ui(e, 1 + n_randint(state, 5));
                    else
                        fmpz_mul_ui(e, f + (j - 1) * stride, 1 + n_randint(state, 20));
                    break;
                case 4: /* random, mixed sizes */
                    fmpz_randtest(e, state, n_randint(state, 2) ? 30 : bits);
                    if (n_randint(state, 4) != 0 && fmpz_is_zero(e))
                        fmpz_one(e);
                    break;
                case 5: /* products of small primes, forcing a multi-limb lcm */
                    fmpz_set_ui(e, n_nth_prime(1 + n_randint(state, 200)));
                    fmpz_mul_ui(e, e, n_randint(state, 4) + 1);
                    break;
                default: /* large equal entries and small divisors */
                    if (n_randint(state, 2))
                        fmpz_set(e, D);
                    else
                        fmpz_set_ui(e, 1 + n_randint(state, 50));
                    break;
            }

            /* rarely a zero */
            if (n_randint(state, 200) == 0)
                fmpz_zero(e);

            /* the gaps of the stride hold garbage */
            if (stride > 1)
                fmpz_randtest(e + 1, state, 100);
        }

        fmpz_one(c);
        for (j = 0; j < len; j++)
            fmpz_lcm(c, c, f + j * stride);

        if (stride == 1 && n_randint(state, 2))
            _fmpz_vec_lcm(d, f, len);
        else
            _fmpz_vec_lcm_strided(d, f, stride, len, NULL);

        if (!fmpz_equal(c, d))
        {
            flint_printf("FAIL (shapes): len = %wd, stride = %wd, shape = %wd\n", len, stride, shape);
            fmpz_print(c), flint_printf("\n\n");
            fmpz_print(d), flint_printf("\n\n");
            flint_abort();
        }

        /* with an initial value, possibly aliased */
        if (n_randint(state, 2))
            fmpz_randtest(inp, state, n_randint(state, 2) ? 20 : 200);
        else
            fmpz_set(inp, D);
        fmpz_lcm(c, c, inp);
        if (n_randint(state, 2))
        {
            _fmpz_vec_lcm_strided(d, f, stride, len, inp);
        }
        else
        {
            fmpz_set(d, inp);
            _fmpz_vec_lcm_strided(d, f, stride, len, d);
        }

        if (!fmpz_equal(c, d))
        {
            flint_printf("FAIL (inp): len = %wd, stride = %wd, shape = %wd\n", len, stride, shape);
            fmpz_print(inp), flint_printf("\n\n");
            fmpz_print(c), flint_printf("\n\n");
            fmpz_print(d), flint_printf("\n\n");
            flint_abort();
        }

        fmpz_clear(c);
        fmpz_clear(d);
        fmpz_clear(D);
        fmpz_clear(inp);
        _fmpz_vec_clear(f, len * stride + 1);
    }

    TEST_FUNCTION_END(state);
}
