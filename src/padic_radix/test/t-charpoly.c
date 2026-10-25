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
#include "fmpz_poly.h"
#include "gr.h"
#include "gr_mat.h"
#include "gr_poly.h"
#include "padic_radix.h"

TEST_FUNCTION_START(padic_radix_charpoly, state)
{
    slong iter;

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        gr_mat_t A;
        slong n, i, j, e;
        int s1, s2;
        padic_radix_struct * c1, * c2;
        padic_radix_t d;

        gr_ctx_init_padic_radix_randtest(ctx, state, 1 + n_randint(state, 100));
        e = radix_limb_exponent(PADIC_RADIX_CTX_RADIX(ctx));
        n = n_randint(state, 12);

        gr_mat_init(A, n, n, ctx);
        c1 = gr_heap_init_vec(n + 1, ctx);
        c2 = gr_heap_init_vec(n + 1, ctx);
        padic_radix_init(d, ctx);

        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
            {
                padic_radix_struct * x = gr_mat_entry_ptr(A, i, j, ctx);
                padic_radix_randtest(x, state, ctx);
                if (n_randint(state, 3) == 0)
                    GR_IGNORE(padic_radix_set_si(x, (slong) n_randint(state, 21) - 10, ctx));
                if (x->u.size != 0 && n_randint(state, 3) == 0)
                {
                    slong rel = (x->N == PADIC_RADIX_EXACT) ? 0 : x->N - x->v;
                    x->v = (slong) n_randint(state, 4 * e + 1) - 2 * e;
                    if (x->N != PADIC_RADIX_EXACT)
                        x->N = x->v + rel;
                    if (_padic_radix_finalize(x, ctx) != GR_SUCCESS)
                        GR_IGNORE(padic_radix_zero(x, ctx));
                }
            }

        s1 = _padic_radix_mat_charpoly_berkowitz(c1, A, ctx);
        s2 = _gr_mat_charpoly_berkowitz(c2, A, ctx);

        if (s1 == GR_SUCCESS && s2 == GR_SUCCESS)
            for (i = 0; i <= n; i++)
            {
                slong N = FLINT_MIN(c1[i].N, c2[i].N);
                GR_IGNORE(padic_radix_sub(d, c1 + i, c2 + i, ctx));
                if (d->u.size != 0 && d->v < N)
                {
                    flint_printf("FAIL: charpoly coefficient %wd\n", i);
                    gr_ctx_println(ctx);
                    flint_printf("c1 = "); gr_println(c1 + i, ctx);
                    flint_printf("c2 = "); gr_println(c2 + i, ctx);
                    flint_abort();
                }
            }

        gr_heap_clear_vec(c1, n + 1, ctx);
        gr_heap_clear_vec(c2, n + 1, ctx);
        padic_radix_clear(d, ctx);
        gr_mat_clear(A, ctx);
        gr_ctx_clear(ctx);
    }

    /* consistency of all the algorithms on random (jagged, signed, exact
       or not) matrices: values agree to the common precision */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        gr_mat_t A;
        slong n, i, j, k, m;
        padic_radix_struct * c[4];
        padic_radix_t d;
        int s[4];

        gr_ctx_init_padic_radix_randtest(ctx, state, 1 + n_randint(state, 60));
        n = n_randint(state, 10);
        gr_mat_init(A, n, n, ctx);
        padic_radix_init(d, ctx);
        for (m = 0; m < 4; m++)
            c[m] = gr_heap_init_vec(n + 1, ctx);

        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
            {
                padic_radix_struct * x = gr_mat_entry_ptr(A, i, j, ctx);
                padic_radix_randtest(x, state, ctx);
                if (n_randint(state, 4) == 0)
                    GR_IGNORE(padic_radix_set_si(x, (slong) n_randint(state, 21) - 10, ctx));
            }

        s[0] = _gr_mat_charpoly_berkowitz(c[0], A, ctx);
        s[1] = padic_radix_mat_charpoly(c[1], A, ctx);
        s[2] = padic_radix_mat_charpoly_lift(c[2], A, ctx);
        s[3] = padic_radix_mat_charpoly_hessenberg(c[3], A, ctx);

        for (m = 1; m < 4; m++)
            if (s[0] == GR_SUCCESS && s[m] == GR_SUCCESS)
                for (k = 0; k <= n; k++)
                {
                    GR_IGNORE(padic_radix_sub(d, c[0] + k, c[m] + k, ctx));
                    if (d->u.size != 0 && d->v < FLINT_MIN(c[0][k].N, c[m][k].N))
                    {
                        flint_printf("FAIL: algorithm %wd vs Berkowitz, coefficient %wd\n", m, k);
                        gr_ctx_println(ctx);
                        gr_mat_print(A, ctx); flint_printf("\n");
                        gr_println(c[0] + k, ctx);
                        gr_println(c[m] + k, ctx);
                        flint_abort();
                    }
                }

        for (m = 0; m < 4; m++)
            gr_heap_clear_vec(c[m], n + 1, ctx);
        padic_radix_clear(d, ctx);
        gr_mat_clear(A, ctx);
        gr_ctx_clear(ctx);
    }

    /* matrices U diag(p^a) V known to O(p^N): correct values and the
       optimal precision N + min v(B_k) (when below 2N), where
       adj(xI - A) = sum B_k x^k (with PADIC_RADIX_FAST_LINALG, only
       correct values) */
    for (iter = 0; iter < 30 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        fmpz_mat_t U, V, M, Bk, T;
        fmpz_poly_t cz;
        fmpz_t t, q;
        gr_mat_t A;
        padic_radix_struct * cp;
        padic_radix_t d, z;
        ulong p;
        slong n, N, i, j, k, * opt;
        int fast = (n_randint(state, 3) == 0);

        p = n_nth_prime(1 + n_randint(state, 6));
        n = 1 + n_randint(state, 9);
        N = 2 + n_randint(state, 30);

        gr_ctx_init_padic_radix(ctx, p, PADIC_RADIX_PREC_INF, N + 200, fast ? PADIC_RADIX_FAST_LINALG : 0);
        fmpz_mat_init(U, n, n); fmpz_mat_init(V, n, n); fmpz_mat_init(M, n, n);
        fmpz_mat_init(Bk, n, n); fmpz_mat_init(T, n, n);
        fmpz_poly_init(cz); fmpz_init(t); fmpz_init_set_ui(q, p);
        gr_mat_init(A, n, n, ctx);
        cp = gr_heap_init_vec(n + 1, ctx);
        padic_radix_init(d, ctx);
        padic_radix_init(z, ctx);
        opt = flint_malloc(sizeof(slong) * n);

        fmpz_mat_randbits(U, state, 6);
        fmpz_mat_randbits(V, state, 6);
        fmpz_mat_zero(T);
        for (i = 0; i < n; i++)
            fmpz_ui_pow_ui(fmpz_mat_entry(T, i, i), p, n_randint(state, 4));
        fmpz_mat_mul(M, U, T);
        fmpz_mat_mul(M, M, V);
        fmpz_mat_charpoly(cz, M);

        fmpz_mat_one(Bk);
        for (k = n - 1; k >= 0; k--)
        {
            slong mv = WORD_MAX;
            for (i = 0; i < n; i++)
                for (j = 0; j < n; j++)
                    if (!fmpz_is_zero(fmpz_mat_entry(Bk, i, j)))
                        mv = FLINT_MIN(mv, (slong) fmpz_remove(t, fmpz_mat_entry(Bk, i, j), q));
            opt[k] = (mv == WORD_MAX) ? WORD_MAX : N + mv;
            if (k > 0)
            {
                fmpz_mat_mul(T, M, Bk);
                fmpz_mat_swap(T, Bk);
                for (i = 0; i < n; i++)
                    fmpz_add(fmpz_mat_entry(Bk, i, i), fmpz_mat_entry(Bk, i, i), cz->coeffs + k);
            }
        }

        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
            {
                padic_radix_struct * x = gr_mat_entry_ptr(A, i, j, ctx);
                GR_MUST_SUCCEED(padic_radix_exact_set_fmpz(x, fmpz_mat_entry(M, i, j), ctx));
                x->N = N;
                GR_MUST_SUCCEED(_padic_radix_finalize(x, ctx));
            }

        GR_MUST_SUCCEED(padic_radix_mat_charpoly(cp, A, ctx));

        for (k = 0; k < n; k++)
        {
            GR_MUST_SUCCEED(padic_radix_exact_set_fmpz(d, cz->coeffs + k, ctx));
            GR_MUST_SUCCEED(padic_radix_sub(z, cp + k, d, ctx));
            if (z->u.size != 0)
            {
                flint_printf("FAIL: value of coefficient %wd\n", k);
                fmpz_mat_print_pretty(M); flint_printf("\n");
                gr_println(cp + k, ctx);
                flint_abort();
            }

            if (!fast && opt[k] != WORD_MAX && opt[k] < 2 * N && cp[k].N < opt[k])
            {
                flint_printf("FAIL: precision of coefficient %wd: %wd, optimal %wd\n", k, cp[k].N, opt[k]);
                fmpz_mat_print_pretty(M); flint_printf("\n");
                flint_abort();
            }
        }

        flint_free(opt);
        padic_radix_clear(d, ctx);
        padic_radix_clear(z, ctx);
        gr_heap_clear_vec(cp, n + 1, ctx);
        gr_mat_clear(A, ctx);
        fmpz_mat_clear(U); fmpz_mat_clear(V); fmpz_mat_clear(M);
        fmpz_mat_clear(Bk); fmpz_mat_clear(T);
        fmpz_poly_clear(cz); fmpz_clear(t); fmpz_clear(q);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
