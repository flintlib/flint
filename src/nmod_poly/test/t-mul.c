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
#if FLINT_HAVE_FFT_SMALL
# include "fft_small.h"
#endif

TEST_FUNCTION_START(nmod_poly_mul, state)
{
    int i, result;

    /* Check aliasing of a and b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c;

        ulong n = n_randtest_not_zero(state);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_randtest(b, state, n_randint(state, 50));
        nmod_poly_randtest(c, state, n_randint(state, 50));

        nmod_poly_mul(a, b, c);
        nmod_poly_mul(b, b, c);

        result = (nmod_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
    }

    /* Check aliasing of a and c */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c;

        ulong n = n_randtest_not_zero(state);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_randtest(b, state, n_randint(state, 50));
        nmod_poly_randtest(c, state, n_randint(state, 50));

        nmod_poly_mul(a, b, c);
        nmod_poly_mul(c, b, c);

        result = (nmod_poly_equal(a, c));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(c), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
    }

    /* Check (b*c)+(b*d) = b*(c+d) */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a1, a2, b, c, d;

        ulong n = n_randtest_not_zero(state);

        nmod_poly_init(a1, n);
        nmod_poly_init(a2, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_init(d, n);
        nmod_poly_randtest(b, state, n_randint(state, 100));
        nmod_poly_randtest(c, state, n_randint(state, 100));
        nmod_poly_randtest(d, state, n_randint(state, 100));

        nmod_poly_mul(a1, b, c);
        nmod_poly_mul(a2, b, d);
        nmod_poly_add(a1, a1, a2);

        nmod_poly_add(c, c, d);
        nmod_poly_mul(a2, b, c);

        result = (nmod_poly_equal(a1, a2));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a1), flint_printf("\n\n");
            nmod_poly_print(a2), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a1);
        nmod_poly_clear(a2);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(d);
    }

#if FLINT_HAVE_FFT_SMALL
    /* Check against the classical algorithm modulo the primes of the
       fft_small context, for which the dispatcher has its own cutoff,
       with lengths on both sides of it, squarings included */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c, d;

        ulong n = get_default_mpn_ctx()->ffts[n_randint(state, MPN_CTX_NCRTS)].mod.n;

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_init(d, n);
        nmod_poly_randtest(b, state, 1 + n_randint(state, n_randint(state, 2) ? 60 : 1000));
        nmod_poly_randtest(c, state, 1 + n_randint(state, n_randint(state, 2) ? 60 : 1000));

        if (n_randint(state, 5) == 0)
            nmod_poly_set(c, b);

        nmod_poly_mul(a, b, c);
        nmod_poly_mul_classical(d, b, c);

        result = (nmod_poly_equal(a, d));
        if (!result)
        {
            flint_printf("FAIL (fft_small prime):\n");
            flint_printf("n = %wu, lengths %wd, %wd\n", n, b->length, c->length);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(d);
    }
#endif

    TEST_FUNCTION_END(state);
}
