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
#include "fmpz_mat.h"
#include "fmpq.h"
#include "fmpq_mat.h"
#include "gr.h"
#include "gr_mat.h"
#include "padic_radix.h"

static slong
_fmpq_val(const fmpq_t x, ulong p)
{
    fmpz_t t, q;
    slong v;
    if (fmpq_is_zero(x))
        return WORD_MAX;
    fmpz_init(t);
    fmpz_init_set_ui(q, p);
    v = fmpz_remove(t, fmpq_numref(x), q) - fmpz_remove(t, fmpq_denref(x), q);
    fmpz_clear(t);
    fmpz_clear(q);
    return v;
}

/* whether the p-adic x agrees with the rational q to its precision */
static int
_agrees(const padic_radix_struct * x, const fmpq_t q, ulong p, gr_ctx_t ctx)
{
    fmpq_t t, y;
    fmpz_t f, pp;
    int ok;

    if (x->N == PADIC_RADIX_EXACT)
        return 1;   /* (not produced here for inexact inputs) */

    fmpq_init(t); fmpq_init(y); fmpz_init(f); fmpz_init_set_ui(pp, p);
    if (x->u.size != 0)
    {
        radix_integer_get_fmpz(f, &x->u, PADIC_RADIX_CTX_RADIX(ctx));
        fmpq_set_fmpz(y, f);
        fmpz_pow_ui(f, pp, FLINT_ABS(x->v));
        if (x->v >= 0)
            fmpq_mul_fmpz(y, y, f);
        else
            fmpq_div_fmpz(y, y, f);
    }
    fmpq_sub(t, y, q);
    ok = fmpq_is_zero(t) || _fmpq_val(t, p) >= x->N;
    fmpq_clear(t); fmpq_clear(y); fmpz_clear(f); fmpz_clear(pp);
    return ok;
}

TEST_FUNCTION_START(padic_radix_solve, state)
{
    slong iter;

    /* A = U diag(p^a) V, B known to O(p^N): correct values and the optimal
       precision N + min_j (min(v(Ainv_ij) + min(0, min_k v(x_k)))) */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        fmpz_mat_t U, V, D, M, Bz;
        fmpq_mat_t Mi, Xq, Bq;
        gr_mat_t A, B, X;
        ulong p;
        slong n, c, N, i, j, k, vx;
        int status, alias, fast = (n_randint(state, 3) == 0);

        p = n_nth_prime(1 + n_randint(state, 6));
        n = 1 + n_randint(state, n_randint(state, 4) ? 8 : 30);
        c = 1 + n_randint(state, 3);
        N = 5 + n_randint(state, 30);
        alias = n_randint(state, 2);

        gr_ctx_init_padic_radix(ctx, p, PADIC_RADIX_PREC_INF, N + 300, fast ? PADIC_RADIX_FAST_LINALG : 0);
        fmpz_mat_init(U, n, n); fmpz_mat_init(V, n, n); fmpz_mat_init(D, n, n);
        fmpz_mat_init(M, n, n); fmpz_mat_init(Bz, n, c);
        fmpq_mat_init(Mi, n, n); fmpq_mat_init(Xq, n, c); fmpq_mat_init(Bq, n, c);
        gr_mat_init(A, n, n, ctx); gr_mat_init(B, n, c, ctx); gr_mat_init(X, n, c, ctx);

        fmpz_mat_randbits(U, state, 6);
        fmpz_mat_randbits(V, state, 6);
        fmpz_mat_randbits(Bz, state, 6);
        for (i = 0; i < n; i++)
            fmpz_ui_pow_ui(fmpz_mat_entry(D, i, i), p, n_randint(state, 4));
        fmpz_mat_mul(M, U, D);
        fmpz_mat_mul(M, M, V);

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                padic_radix_struct * x = gr_mat_entry_ptr(A, i, j, ctx);
                GR_MUST_SUCCEED(padic_radix_exact_set_fmpz(x, fmpz_mat_entry(M, i, j), ctx));
                x->N = N;
                GR_MUST_SUCCEED(_padic_radix_finalize(x, ctx));
            }
            for (j = 0; j < c; j++)
            {
                padic_radix_struct * x = gr_mat_entry_ptr(B, i, j, ctx);
                GR_MUST_SUCCEED(padic_radix_exact_set_fmpz(x, fmpz_mat_entry(Bz, i, j), ctx));
                x->N = N;
                GR_MUST_SUCCEED(_padic_radix_finalize(x, ctx));
            }
        }

        if (alias)
        {
            GR_MUST_SUCCEED(gr_mat_set(X, B, ctx));
            status = padic_radix_mat_nonsingular_solve(X, A, X, ctx);
        }
        else
            status = padic_radix_mat_nonsingular_solve(X, A, B, ctx);

        fmpq_mat_set_fmpz_mat(Mi, M);
        if (fmpq_mat_inv(Mi, Mi) && status == GR_SUCCESS)
        {
            fmpq_mat_set_fmpz_mat(Bq, Bz);
            fmpq_mat_mul(Xq, Mi, Bq);

            for (k = 0; k < c; k++)
            {
                vx = 0;
                for (i = 0; i < n; i++)
                    vx = FLINT_MIN(vx, _fmpq_val(fmpq_mat_entry(Xq, i, k), p));

                for (i = 0; i < n; i++)
                {
                    const padic_radix_struct * x = gr_mat_entry_srcptr(X, i, k, ctx);
                    slong opt = WORD_MAX;

                    for (j = 0; j < n; j++)
                    {
                        slong v = _fmpq_val(fmpq_mat_entry(Mi, i, j), p);
                        if (v != WORD_MAX)
                            opt = FLINT_MIN(opt, N + v + vx);
                    }

                    if (!_agrees(x, fmpq_mat_entry(Xq, i, k), p, ctx))
                    {
                        flint_printf("FAIL: value (%wd, %wd)\n", i, k);
                        fmpz_mat_print_pretty(M); flint_printf("\n");
                        gr_println(x, ctx);
                        flint_abort();
                    }

                    /* (the differential bound is the true precision when
                       the solution is determined to positive precision) */
                    if (!fast && opt != WORD_MAX && opt > 0 && x->N < opt)
                    {
                        flint_printf("FAIL: precision (%wd, %wd): %wd, optimal %wd\n", i, k, x->N, opt);
                        fmpz_mat_print_pretty(M); flint_printf("\n");
                        flint_abort();
                    }
                }
            }
        }
        else if (status == GR_SUCCESS)
        {
            flint_printf("FAIL: singular matrix solved\n");
            flint_abort();
        }

        fmpz_mat_clear(U); fmpz_mat_clear(V); fmpz_mat_clear(D);
        fmpz_mat_clear(M); fmpz_mat_clear(Bz);
        fmpq_mat_clear(Mi); fmpq_mat_clear(Xq); fmpq_mat_clear(Bq);
        gr_mat_clear(A, ctx); gr_mat_clear(B, ctx); gr_mat_clear(X, ctx);
        gr_ctx_clear(ctx);
    }

    /* random contexts and matrices: A X = B holds to the precision */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        gr_mat_t A, B, X, AX;
        slong n, i, j;
        int status;

        gr_ctx_init_padic_radix_randtest(ctx, state, 1 + n_randint(state, 60));
        n = n_randint(state, 8);
        gr_mat_init(A, n, n, ctx); gr_mat_init(B, n, 2, ctx);
        gr_mat_init(X, n, 2, ctx); gr_mat_init(AX, n, 2, ctx);
        GR_IGNORE(gr_mat_randtest(A, state, ctx));
        GR_IGNORE(gr_mat_randtest(B, state, ctx));

        status = padic_radix_mat_nonsingular_solve(X, A, B, ctx);
        if (status == GR_SUCCESS)
        {
            status = gr_mat_mul(AX, A, X, ctx);
            if (status == GR_SUCCESS && gr_mat_equal(AX, B, ctx) == T_FALSE)
            {
                flint_printf("FAIL: A X != B\n");
                gr_ctx_println(ctx);
                gr_mat_print(A, ctx); flint_printf("\n");
                gr_mat_print(B, ctx); flint_printf("\n");
                gr_mat_print(X, ctx); flint_printf("\n");
                flint_abort();
            }
        }
        (void) i; (void) j;

        gr_mat_clear(A, ctx); gr_mat_clear(B, ctx); gr_mat_clear(X, ctx); gr_mat_clear(AX, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
