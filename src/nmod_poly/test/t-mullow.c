/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_mullow, state)
{
    int i, result;

    /* Compare with truncated product of a and b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c;
        slong trunc;
        ulong n = n_randtest_not_zero(state);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        trunc = n_randint(state, 50);
        nmod_poly_randtest(b, state, trunc);
        nmod_poly_randtest(c, state, trunc);

        nmod_poly_mullow(a, b, c, trunc);
        nmod_poly_mul(b, b, c);
        nmod_poly_truncate(b, trunc);

        result = (nmod_poly_equal(a, b));
        if (!result)
        {
            flint_printf(":\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
    }

    /* Compare with the classical algorithm with a short operand (rows of
       scalar products in _nmod_poly_mullow), moduli of 30 bits or more and
       powers of 2 */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c, d;
        slong trunc;
        ulong n = n_randint(state, 4) == 0
                ? UWORD(1) << (1 + n_randint(state, FLINT_BITS - 1))
                : n_randtest_bits(state, 30 + n_randint(state, FLINT_BITS - 29));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_init(d, n);
        nmod_poly_randtest(b, state, 1 + n_randint(state, 300));
        nmod_poly_randtest(c, state, 1 + n_randint(state, 18));

        if (n_randint(state, 2))
            nmod_poly_swap(b, c);

        trunc = n_randint(state, b->length + c->length + 1);

        nmod_poly_mullow(a, b, c, trunc);
        nmod_poly_mullow_classical(d, b, c, trunc);

        result = (nmod_poly_equal(a, d));
        if (!result)
        {
            flint_printf("FAIL (short operand):\n");
            flint_printf("n = %wu, lengths %wd, %wd, trunc = %wd\n",
                         n, b->length, c->length, trunc);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(d);
    }

    /* Compare with the classical algorithm for moduli of any size and a
       shorter operand of up to 600 coefficients (IFMA kernel, then
       fft_small) */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c, d;
        slong trunc;
        ulong n = FLINT_MAX(n_randtest_bits(state, 1 + n_randint(state, FLINT_BITS)), 1);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_init(d, n);
        nmod_poly_randtest(b, state, 1 + n_randint(state, n_randint(state, 4) ? 40 : 600));
        nmod_poly_randtest(c, state, b->length + n_randint(state, n_randint(state, 8) ? 300 : 3000));

        if (n_randint(state, 2))
            nmod_poly_swap(b, c);
        if (n_randint(state, 5) == 0)
            nmod_poly_set(c, b);

        trunc = n_randint(state, b->length + c->length + 1);

        nmod_poly_mullow(a, b, c, trunc);
        nmod_poly_mullow_classical(d, b, c, trunc);

        result = (nmod_poly_equal(a, d));
        if (!result)
        {
            flint_printf("FAIL (moduli of any size):\n");
            flint_printf("n = %wu, lengths %wd, %wd, trunc = %wd\n",
                         n, b->length, c->length, trunc);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(d);
    }

    TEST_FUNCTION_END(state);
}
