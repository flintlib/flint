/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ca_vec.h"
#include "ca_poly.h"

static int
_ca_poly_test_equal_repr(const ca_poly_t A, const ca_poly_t B, ca_ctx_t ctx)
{
    slong i;

    if (A->length != B->length)
        return 0;

    for (i = 0; i < A->length; i++)
        if (!ca_equal_repr(A->coeffs + i, B->coeffs + i, ctx))
            return 0;

    return 1;
}

static int
_ca_vec_test_equal_repr(ca_srcptr A, ca_srcptr B, slong len, ca_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        if (!ca_equal_repr(A + i, B + i, ctx))
            return 0;

    return 1;
}

TEST_FUNCTION_START(ca_poly_shift_reverse, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_poly_t A, B, E;
        ca_ptr u, v, w;
        slong i, n, len, alloc;
        int op;

        ca_ctx_init(ctx);
        ca_poly_init(A, ctx);
        ca_poly_init(B, ctx);
        ca_poly_init(E, ctx);

        if (n_randint(state, 4) == 0)
            ca_poly_randtest(A, state, 1 + n_randint(state, 6), 1, 10, ctx);
        else
            ca_poly_randtest_rational(A, state, n_randint(state, 10), 10, ctx);

        ca_poly_randtest_rational(B, state, n_randint(state, 10), 10, ctx);

        n = n_randint(state, 12);
        op = n_randint(state, 3);

        /* Polynomial-level functions, compared with an expected polynomial
           built coefficient by coefficient. */
        ca_poly_zero(E, ctx);

        if (op == 0)
        {
            for (i = 0; i < A->length; i++)
                ca_poly_set_coeff_ca(E, i + n, A->coeffs + i, ctx);
        }
        else if (op == 1)
        {
            for (i = n; i < A->length; i++)
                ca_poly_set_coeff_ca(E, i - n, A->coeffs + i, ctx);
        }
        else
        {
            for (i = 0; i < FLINT_MIN(n, A->length); i++)
                ca_poly_set_coeff_ca(E, n - 1 - i, A->coeffs + i, ctx);
        }

        for (i = 0; i < 2; i++)
        {
            /* i == 0: no aliasing; i == 1: aliasing */
            if (i == 0)
            {
                if (op == 0)
                    ca_poly_shift_left(B, A, n, ctx);
                else if (op == 1)
                    ca_poly_shift_right(B, A, n, ctx);
                else
                    ca_poly_reverse(B, A, n, ctx);
            }
            else
            {
                ca_poly_set(B, A, ctx);

                if (op == 0)
                    ca_poly_shift_left(B, B, n, ctx);
                else if (op == 1)
                    ca_poly_shift_right(B, B, n, ctx);
                else
                    ca_poly_reverse(B, B, n, ctx);
            }

            if (!_ca_poly_test_equal_repr(B, E, ctx))
            {
                flint_printf("FAIL (poly, op = %d, aliasing = %wd, n = %wd)\n", op, i, n);
                flint_printf("A = "); ca_poly_print(A, ctx); flint_printf("\n");
                flint_printf("B = "); ca_poly_print(B, ctx); flint_printf("\n");
                flint_printf("E = "); ca_poly_print(E, ctx); flint_printf("\n");
                flint_abort();
            }
        }

        /* Underscore functions, including aliasing. */
        len = A->length;

        if (len > 0 && (op != 2 || n >= len))
        {
            if (op == 1)
                n = n_randint(state, len);

            alloc = len + n;
            u = _ca_vec_init(alloc, ctx);
            v = _ca_vec_init(alloc, ctx);
            w = _ca_vec_init(alloc, ctx);

            /* expected result in w */
            if (op == 0)
            {
                for (i = 0; i < len; i++)
                    ca_set(w + i + n, A->coeffs + i, ctx);
                for (i = 0; i < n; i++)
                    ca_zero(w + i, ctx);
            }
            else if (op == 1)
            {
                for (i = n; i < len; i++)
                    ca_set(w + i - n, A->coeffs + i, ctx);
            }
            else
            {
                for (i = 0; i < n - len; i++)
                    ca_zero(w + i, ctx);
                for (i = 0; i < len; i++)
                    ca_set(w + n - 1 - i, A->coeffs + i, ctx);
            }

            /* fill output with garbage */
            for (i = 0; i < alloc; i++)
                ca_randtest_rational(v + i, state, 10, ctx);

            if (op == 0)
                _ca_poly_shift_left(v, A->coeffs, len, n, ctx);
            else if (op == 1)
                _ca_poly_shift_right(v, A->coeffs, len, n, ctx);
            else
                _ca_poly_reverse(v, A->coeffs, len, n, ctx);

            _ca_vec_set(u, A->coeffs, len, ctx);
            for (i = len; i < alloc; i++)
                ca_randtest_rational(u + i, state, 10, ctx);

            if (op == 0)
                _ca_poly_shift_left(u, u, len, n, ctx);
            else if (op == 1)
                _ca_poly_shift_right(u, u, len, n, ctx);
            else
                _ca_poly_reverse(u, u, len, n, ctx);

            {
                slong rlen = (op == 0) ? len + n : (op == 1) ? len - n : n;

                if (!_ca_vec_test_equal_repr(v, w, rlen, ctx) ||
                    !_ca_vec_test_equal_repr(u, w, rlen, ctx))
                {
                    flint_printf("FAIL (underscore, op = %d, len = %wd, n = %wd)\n", op, len, n);
                    flint_printf("A = "); ca_poly_print(A, ctx); flint_printf("\n");
                    flint_abort();
                }
            }

            _ca_vec_clear(u, alloc, ctx);
            _ca_vec_clear(v, alloc, ctx);
            _ca_vec_clear(w, alloc, ctx);
        }

        ca_poly_clear(A, ctx);
        ca_poly_clear(B, ctx);
        ca_poly_clear(E, ctx);
        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
