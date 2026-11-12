/*
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_mullow_classical, state)
{
    int i, result;

    /* Check aliasing of a and b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c;
        slong len, trunc;
        ulong n = n_randtest_not_zero(state);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_randtest(b, state, n_randint(state, 50));
        nmod_poly_randtest(c, state, n_randint(state, 50));

        len = b->length + c->length - 1;
        if (len <= 0)
            trunc = 0;
        else
            trunc = n_randint(state, b->length + c->length);

        nmod_poly_mullow_classical(a, b, c, trunc);
        nmod_poly_mullow_classical(b, b, c, trunc);

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
        slong len, trunc;
        ulong n = n_randtest_not_zero(state);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_randtest(b, state, n_randint(state, 50));
        nmod_poly_randtest(c, state, n_randint(state, 50));

        len = b->length + c->length - 1;
        if (len <= 0)
            trunc = 0;
        else
            trunc = n_randint(state, b->length + c->length - 1);

        nmod_poly_mullow_classical(a, b, c, trunc);
        nmod_poly_mullow_classical(c, b, c, trunc);

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

    /* Compare with mul_basecase */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c, d;
        slong len, trunc;
        ulong n = n_randtest_not_zero(state);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_init(d, n);
        nmod_poly_randtest(b, state, n_randint(state, 50));
        nmod_poly_randtest(c, state, n_randint(state, 50));

        len = b->length + c->length - 1;
        if (len <= 0)
            trunc = 0;
        else
            trunc = n_randint(state, b->length + c->length - 1);

        if (n_randint(state, 2))  /* check squaring */
        {
            nmod_poly_mul_classical(a, b, b);
            nmod_poly_truncate(a, trunc);
            nmod_poly_mullow_classical(d, b, b, trunc);
        }
        else
        {
            nmod_poly_mul_classical(a, b, c);
            nmod_poly_truncate(a, trunc);
            nmod_poly_mullow_classical(d, b, c, trunc);
        }

        result = (nmod_poly_equal(a, d));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(d), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(d);
    }

    /* Check the underscore function with the operands in either order */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        nmod_t mod;
        nn_ptr a, b, p, r;
        slong la, lb, trunc;

        nmod_init(&mod, n_randtest_not_zero(state));
        la = 1 + n_randint(state, 50);
        lb = 1 + n_randint(state, 50);
        trunc = 1 + n_randint(state, la + lb - 1);

        a = _nmod_vec_init(la);
        b = _nmod_vec_init(lb);
        p = _nmod_vec_init(la + lb - 1);
        r = _nmod_vec_init(trunc);
        _nmod_vec_randtest(a, state, la, mod);
        _nmod_vec_randtest(b, state, lb, mod);

        if (la >= lb)
            _nmod_poly_mul(p, a, la, b, lb, mod);
        else
            _nmod_poly_mul(p, b, lb, a, la, mod);

        _nmod_poly_mullow_classical(r, a, la, b, lb, trunc, mod);

        result = _nmod_vec_equal(r, p, trunc);
        if (!result)
        {
            flint_printf("FAIL (operands in either order):\n");
            flint_printf("n = %wu, la = %wd, lb = %wd, trunc = %wd\n", mod.n, la, lb, trunc);
            fflush(stdout);
            flint_abort();
        }

        _nmod_vec_clear(a);
        _nmod_vec_clear(b);
        _nmod_vec_clear(p);
        _nmod_vec_clear(r);
    }

    TEST_FUNCTION_END(state);
}
