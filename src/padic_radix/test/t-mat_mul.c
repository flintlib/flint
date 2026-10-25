/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "gr.h"
#include "gr_mat.h"
#include "padic_radix.h"

/* random element with valuations and precisions spread over a few limbs,
   zeros (exact and inexact) and signed exact units */
static void
_randelem_mat(padic_radix_t x, flint_rand_t state, gr_ctx_t ctx)
{
    slong e = radix_limb_exponent(PADIC_RADIX_CTX_RADIX(ctx));

    padic_radix_randtest(x, state, ctx);

    switch (n_randint(state, 8))
    {
        case 0:
            GR_IGNORE(padic_radix_zero(x, ctx));
            break;
        case 1:
            radix_integer_zero(&x->u, PADIC_RADIX_CTX_RADIX(ctx));
            x->v = 0;
            x->N = (slong) n_randint(state, 8 * e + 1) - 4 * e;
            GR_IGNORE(_padic_radix_finalize(x, ctx));
            break;
        case 2:
        case 3:
            if (x->u.size != 0)
            {
                slong rel = (x->N == PADIC_RADIX_EXACT) ? 0 : x->N - x->v;
                x->v = (slong) n_randint(state, 6 * e + 1) - 3 * e;
                if (x->N != PADIC_RADIX_EXACT)
                    x->N = x->v + rel;
                if (_padic_radix_finalize(x, ctx) != GR_SUCCESS)
                    GR_IGNORE(padic_radix_zero(x, ctx));
            }
            break;
        case 4:
            GR_IGNORE(padic_radix_set_si(x, (slong) n_randint(state, 41) - 20, ctx));
            break;
        default:
            break;
    }
}

TEST_FUNCTION_START(padic_radix_mat_mul, state)
{
    slong iter;

    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        gr_mat_t A, B, C, D;
        slong m, k, n, i, j;
        int s1, s2;

        gr_ctx_init_padic_radix_randtest(ctx, state, 1 + n_randint(state, 200));

        if (n_randint(state, 8) == 0)
        {
            /* larger shapes: more than one block in the fmpz route */
            m = n_randint(state, 12);
            k = n_randint(state, 60);
            n = n_randint(state, 12);
            if (n_randint(state, 2))
                m = n = k;
        }
        else
        {
            m = n_randint(state, 7);
            k = n_randint(state, 12);
            n = n_randint(state, 7);
        }

        gr_mat_init(A, m, k, ctx);
        gr_mat_init(B, k, n, ctx);
        gr_mat_init(C, m, n, ctx);
        gr_mat_init(D, m, n, ctx);

        for (i = 0; i < m; i++)
            for (j = 0; j < k; j++)
                _randelem_mat(gr_mat_entry_ptr(A, i, j, ctx), state, ctx);
        for (i = 0; i < k; i++)
            for (j = 0; j < n; j++)
                _randelem_mat(gr_mat_entry_ptr(B, i, j, ctx), state, ctx);
        GR_IGNORE(gr_mat_randtest(C, state, ctx));

        /* the packed kernel gives exactly the delayed dot products */
        s1 = padic_radix_mat_mul_packed(C, A, B, ctx);

        if (s1 != GR_UNABLE)
        {
            s2 = GR_SUCCESS;
            for (i = 0; i < m; i++)
                for (j = 0; j < n; j++)
                    s2 |= padic_radix_dot_strided_delayed(gr_mat_entry_ptr(D, i, j, ctx), NULL, 0,
                        gr_mat_entry_srcptr(A, i, 0, ctx), 1,
                        (k == 0) ? NULL : gr_mat_entry_srcptr(B, 0, j, ctx), B->stride,
                        k, ctx);

            if (s1 == GR_SUCCESS && s2 == GR_SUCCESS)
            {
                for (i = 0; i < m; i++)
                    for (j = 0; j < n; j++)
                    {
                        padic_radix_struct * c = gr_mat_entry_ptr(C, i, j, ctx);
                        padic_radix_struct * d = gr_mat_entry_ptr(D, i, j, ctx);

                        if (c->v != d->v || c->N != d->N ||
                            !radix_integer_equal(&c->u, &d->u, PADIC_RADIX_CTX_RADIX(ctx)))
                        {
                            flint_printf("FAIL: packed mat_mul vs delayed dot, entry (%wd, %wd)\n", i, j);
                            gr_ctx_println(ctx);
                            flint_printf("A = "); gr_mat_print(A, ctx); flint_printf("\n");
                            flint_printf("B = "); gr_mat_print(B, ctx); flint_printf("\n");
                            flint_printf("c = "); gr_println(c, ctx);
                            flint_printf("d = "); gr_println(d, ctx);
                            flint_abort();
                        }
                    }
            }
        }

        /* the fmpz_mat_mul route agrees exactly with the packed kernel */
        s1 = padic_radix_mat_mul_packed(C, A, B, ctx);
        s2 = padic_radix_mat_mul_fmpz(D, A, B, ctx);

        if (s1 == GR_SUCCESS && s2 == GR_SUCCESS)
        {
            for (i = 0; i < m; i++)
                for (j = 0; j < n; j++)
                {
                    padic_radix_struct * c = gr_mat_entry_ptr(C, i, j, ctx);
                    padic_radix_struct * d = gr_mat_entry_ptr(D, i, j, ctx);

                    if (c->v != d->v || c->N != d->N ||
                        !radix_integer_equal(&c->u, &d->u, PADIC_RADIX_CTX_RADIX(ctx)))
                    {
                        flint_printf("FAIL: fmpz mat_mul vs packed, entry (%wd, %wd)\n", i, j);
                        gr_ctx_println(ctx);
                        flint_printf("A = "); gr_mat_print(A, ctx); flint_printf("\n");
                        flint_printf("B = "); gr_mat_print(B, ctx); flint_printf("\n");
                        flint_printf("c = "); gr_println(c, ctx);
                        flint_printf("d = "); gr_println(d, ctx);
                        flint_abort();
                    }
                }
        }
        else if (s1 != GR_UNABLE && s2 != GR_UNABLE && s1 != s2)
        {
            flint_printf("FAIL: fmpz mat_mul status %d vs %d\n", s2, s1);
            flint_abort();
        }

        /* aliasing */
        if (m == k && k == n && s2 == GR_SUCCESS)
        {
            gr_mat_t E;
            gr_mat_init(E, m, m, ctx);
            GR_MUST_SUCCEED(gr_mat_set(E, A, ctx));
            s1 = padic_radix_mat_mul_fmpz(E, E, B, ctx);
            if (s1 != GR_SUCCESS || gr_mat_equal(E, D, ctx) == T_FALSE)
            {
                flint_printf("FAIL: fmpz mat_mul aliasing\n");
                flint_abort();
            }
            gr_mat_clear(E, ctx);
        }

        /* aliasing through the method */
        if (m == k && k == n)
        {
            s1 = gr_mat_mul(C, A, B, ctx);
            s2 = gr_mat_mul(A, A, B, ctx);
            if (s1 == GR_SUCCESS && s2 == GR_SUCCESS && gr_mat_equal(A, C, ctx) == T_FALSE)
            {
                flint_printf("FAIL: aliasing\n");
                flint_abort();
            }
        }

        gr_mat_clear(A, ctx);
        gr_mat_clear(B, ctx);
        gr_mat_clear(C, ctx);
        gr_mat_clear(D, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
