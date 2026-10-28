/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "thread_support.h"
#include "fmpz.h"
#include "fmpz_vec.h"

TEST_FUNCTION_START(fmpz_vec_content, state)
{
    int i, result;

    /* Check that content(a f) = abs(a) content(f) */
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

        _fmpz_vec_content(c, f, len);
        _fmpz_vec_scalar_mul_fmpz(f, f, len, a);
        fmpz_abs(a, a);
        fmpz_mul(c, a, c);
        _fmpz_vec_content(d, f, len);

        result = (fmpz_equal(c, d));
        if (!result)
        {
            flint_printf("FAIL:\n");
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

    /* Compare with chained gcds on vectors of various shapes, also for
       _fmpz_vec_content_chained and with aliasing */
    for (i = 0; i < 3000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, c, d, inp;
        fmpz *f;
        slong len, j;
        flint_bitcnt_t bits, abits;
        int shape;

        len = n_randint(state, 3) == 0 ? n_randint(state, 1000) : n_randint(state, 40);
        bits = n_randint(state, 3) == 0 ? 2 + n_randint(state, 600) : 2 + n_randint(state, 130);
        abits = n_randint(state, 200);
        shape = n_randint(state, 6);

        fmpz_init(a);
        fmpz_init(c);
        fmpz_init(d);
        fmpz_init(inp);
        f = _fmpz_vec_init(len);

        fmpz_randtest_not_zero(a, state, abits + 1);
        if (n_randint(state, 4) == 0)
            fmpz_mul_2exp(a, a, n_randint(state, 150));

        for (j = 0; j < len; j++)
        {
            flint_bitcnt_t b = bits;

            if (shape == 1)         /* growing sizes */
                b = 1 + bits * (j + 1) / len;
            else if (shape == 2)    /* shrinking sizes */
                b = 1 + bits * (len - j) / len;
            else if (shape == 3)    /* mixed */
                b = n_randint(state, 2) ? 1 + n_randint(state, 40) : bits;

            fmpz_randtest(f + j, state, b);
            if (shape == 4 && n_randint(state, 4) != 0)
                fmpz_zero(f + j);
            if (shape == 5)
                fmpz_set_si(f + j, n_randint(state, 7) - 3);
            fmpz_mul(f + j, f + j, a);
        }

        /* break the content at a random place */
        if (len != 0 && n_randint(state, 4) == 0)
            fmpz_randtest(f + n_randint(state, len), state, bits);

        fmpz_zero(c);
        for (j = 0; j < len; j++)
            fmpz_gcd(c, c, f + j);

        if (n_randint(state, 2))
            _fmpz_vec_content(d, f, len);
        else
        {
            /* aliased with an entry */
            fmpz * g = _fmpz_vec_init(len + 1);
            _fmpz_vec_set(g, f, len);
            _fmpz_vec_content(g + len, g, len);
            fmpz_set(d, g + len);
            if (len != 0)
            {
                slong k = n_randint(state, len);
                fmpz_t e;
                fmpz_init(e);
                _fmpz_vec_content(e, f, len);
                fmpz_set(g + k, f + k);
                _fmpz_vec_content(g + k, g, len);
                if (!fmpz_equal(g + k, e))
                    TEST_FUNCTION_FAIL("aliased entry: len = %wd\n", len);
                fmpz_clear(e);
            }
            _fmpz_vec_clear(g, len + 1);
        }

        if (!fmpz_equal(c, d))
        {
            flint_printf("FAIL (shapes): len = %wd, bits = %wu, shape = %d\n", len, bits, shape);
            fmpz_print(c), flint_printf("\n\n");
            fmpz_print(d), flint_printf("\n\n");
            flint_abort();
        }

        /* chained, with inp possibly aliased with the result */
        fmpz_randtest(inp, state, n_randint(state, 2) ? 30 : 300);
        if (n_randint(state, 2))
            fmpz_mul(inp, inp, a);
        fmpz_gcd(c, c, inp);
        if (n_randint(state, 2))
        {
            _fmpz_vec_content_chained(d, f, len, inp);
        }
        else
        {
            fmpz_set(d, inp);
            _fmpz_vec_content_chained(d, f, len, d);
        }

        if (!fmpz_equal(c, d))
        {
            flint_printf("FAIL (chained): len = %wd, bits = %wu, shape = %d\n", len, bits, shape);
            fmpz_print(inp), flint_printf("\n\n");
            fmpz_print(c), flint_printf("\n\n");
            fmpz_print(d), flint_printf("\n\n");
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(c);
        fmpz_clear(d);
        fmpz_clear(inp);
        _fmpz_vec_clear(f, len);
    }

    /* Multithreaded: long vectors with large entries */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpz_t a, c, d;
        fmpz *f;
        slong len, j;
        flint_bitcnt_t bits;

        flint_set_num_threads(1 + n_randint(state, 4));

        len = 64 + n_randint(state, 400);
        bits = n_randint(state, 2) ? 64 * (50 + n_randint(state, 400)) : 30 + n_randint(state, 30);
        if (bits < 64)
            len = 60000 + n_randint(state, 10000);

        fmpz_init(a);
        fmpz_init(c);
        fmpz_init(d);
        f = _fmpz_vec_init(len);

        fmpz_randtest_not_zero(a, state, n_randint(state, 2) ? 20 : 2000);
        for (j = 0; j < len; j++)
        {
            fmpz_randtest(f + j, state, bits);
            fmpz_mul(f + j, f + j, a);
        }
        if (n_randint(state, 3) == 0)
            fmpz_add_ui(f + n_randint(state, len), f + n_randint(state, len), 1);

        fmpz_zero(c);
        for (j = 0; j < len; j++)
            fmpz_gcd(c, c, f + j);

        _fmpz_vec_content(d, f, len);

        if (!fmpz_equal(c, d))
        {
            flint_printf("FAIL (threads): len = %wd, bits = %wu\n", len, bits);
            fmpz_print(c), flint_printf("\n\n");
            fmpz_print(d), flint_printf("\n\n");
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(c);
        fmpz_clear(d);
        _fmpz_vec_clear(f, len);
    }

    flint_set_num_threads(1);

    TEST_FUNCTION_END(state);
}
