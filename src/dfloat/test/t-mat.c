/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_mat.h"
#include "dfloat.h"

/* The matrix methods of the plain real and complex rings (and their
   strong variants):

   - gr_mat_mul against the exact product (through acb), entrywise:
     |C - AB|_ij <= 2^(-53n + 6) (k + 2) sum_l |A_il| |B_lj|, for
     random shapes (on both sides of the cutoffs and of the SIMD block
     widths), entries with zeros, real values, spread magnitudes, and
     with C aliasing A or B; the generic max-norm test against the
     classical product; the products through arb/acb_mat_approx_mul (used
     for large matrices) normwise against the classical product;
   - gr_mat_lu: P A = L U to within a normwise tolerance, and the
     rank of rank-deficient matrices (zero rows and columns), with and
     without rank_check, as gr_mat_lu_classical finds it;
   - gr_mat_nonsingular_solve_tril / triu (unit or not, with few or
     many right-hand sides, X aliasing B) and gr_mat_nonsingular_solve:
     residuals within a normwise tolerance. */

void _d2__mat_mul_arb(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, int cplx);
void _d3__mat_mul_arb(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, int cplx);
void _d4__mat_mul_arb(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, int cplx);

static void
_t_m_rand_entry(gr_ptr x, int n, int cplx, int kind, flint_rand_t state, gr_ctx_t ctx)
{
    double * d = x;
    int k, part, parts = cplx ? 2 : 1;

    for (part = 0; part < parts; part++)
    {
        double * p = d + part * n;
        int r = n_randint(state, 8);
        for (k = 0; k < n; k++)
            p[k] = 0.0;
        if (kind == 1 || r == 0)
            continue;                       /* zero (part) */
        if (r == 1)
        {
            p[0] = (double) ((slong) n_randint(state, 7) - 3);
            continue;
        }
        p[0] = (1.0 + (double) n_randint(state, 1 << 30) * 0x1p-30) * (n_randint(state, 2) ? 1 : -1);
        if (kind == 2)
            p[0] = ldexp(p[0], (int) n_randint(state, 40) - 20);
        for (k = 1; k < n; k++)
            p[k] = ldexp(p[k - 1], -53 - (int) n_randint(state, 3)) * (n_randint(state, 2) ? 0.75 : -0.625);
    }
}

static void
_t_m_rand(gr_mat_t A, int n, int cplx, flint_rand_t state, gr_ctx_t ctx)
{
    slong i, j;
    int kind = n_randint(state, 3);     /* generic, sparse, spread */
    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            _t_m_rand_entry(GR_MAT_ENTRY(A, i, j, ctx->sizeof_elem), n, cplx,
                (kind == 1 && n_randint(state, 2)) ? 1 : (kind == 1 ? 0 : kind), state, ctx);
}

/* a random well-conditioned matrix: diagonally dominant up to a
   permutation of the rows */
static void
_t_m_rand_wellcond(gr_mat_t A, int n, int cplx, flint_rand_t state, gr_ctx_t ctx)
{
    slong i, sz = ctx->sizeof_elem, len = A->r;
    gr_ptr t;
    _t_m_rand(A, n, cplx, state, ctx);
    GR_TMP_INIT(t, ctx);
    for (i = 0; i < len; i++)
    {
        GR_MUST_SUCCEED(gr_set_si(t, 2 * len + 1, ctx));
        GR_MUST_SUCCEED(gr_add(GR_MAT_ENTRY(A, i, i, sz), GR_MAT_ENTRY(A, i, i, sz), t, ctx));
    }
    for (i = 0; i < len; i++)
        _gr_vec_swap(GR_MAT_ENTRY(A, i, 0, sz), GR_MAT_ENTRY(A, n_randint(state, len), 0, sz), len, ctx);
    GR_TMP_CLEAR(t, ctx);
}

/* max |entries| as a double */
static double
_t_m_norm(const gr_mat_t A, gr_ctx_t ctx, gr_ctx_t CC)
{
    slong i, j;
    double m = 0.0;
    acb_t z;
    arb_t t;
    acb_init(z);
    arb_init(t);
    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
        {
            GR_MUST_SUCCEED(gr_set_other(z, GR_MAT_ENTRY(A, i, j, ctx->sizeof_elem), ctx, CC));
            acb_abs(t, z, 100);
            m = FLINT_MAX(m, arf_get_d(arb_midref(t), ARF_RND_UP));
        }
    acb_clear(z);
    arb_clear(t);
    return m;
}

/* max |A B - C| (exactly, through acb) */
static double
_t_m_residual(const gr_mat_t A, const gr_mat_t B, const gr_mat_t C, gr_ctx_t ctx, gr_ctx_t CC)
{
    slong i, j, l, sz = ctx->sizeof_elem;
    double m = 0.0;
    acb_t s, x, y;
    arb_t t;
    acb_init(s); acb_init(x); acb_init(y); arb_init(t);
    for (i = 0; i < C->r; i++)
        for (j = 0; j < C->c; j++)
        {
            GR_MUST_SUCCEED(gr_set_other(s, GR_MAT_ENTRY(C, i, j, sz), ctx, CC));
            acb_neg(s, s);
            for (l = 0; l < A->c; l++)
            {
                GR_MUST_SUCCEED(gr_set_other(x, GR_MAT_ENTRY(A, i, l, sz), ctx, CC));
                GR_MUST_SUCCEED(gr_set_other(y, GR_MAT_ENTRY(B, l, j, sz), ctx, CC));
                acb_addmul(s, x, y, 2000);
            }
            acb_abs(t, s, 100);
            m = FLINT_MAX(m, arf_get_d(arb_midref(t), ARF_RND_UP));
        }
    acb_clear(s); acb_clear(x); acb_clear(y); arb_clear(t);
    return m;
}

TEST_FUNCTION_START(mat, state)
{
    slong iter;

    DFLOAT_TEST_SKIP_IF_UNSUPPORTED(state);

    /* the products through arb_mat_approx_mul / acb_mat_approx_mul
       (used above large cutoffs), normwise against the classical one */
    for (iter = 0; iter < 8; iter++)
    {
        gr_ctx_t ctx, CC;
        gr_mat_t A, B, C, D;
        slong len = 1 + n_randint(state, 40);
        int n = 2 + iter % 3, cplx = (iter >= 4);
        double err;

        GR_MUST_SUCCEED(gr_ctx_init_dfloat(ctx, n, cplx ? DFLOAT_COMPLEX : 0));
        gr_ctx_init_complex_acb(CC, 400);
        gr_mat_init(A, len, len, ctx);
        gr_mat_init(B, len, len, ctx);
        gr_mat_init(C, len, len, ctx);
        gr_mat_init(D, len, len, ctx);
        _t_m_rand(A, n, cplx, state, ctx);
        _t_m_rand(B, n, cplx, state, ctx);
        switch (n)
        {
            case 2: _d2__mat_mul_arb(C, A, B, cplx); break;
            case 3: _d3__mat_mul_arb(C, A, B, cplx); break;
            default: _d4__mat_mul_arb(C, A, B, cplx); break;
        }
        GR_MUST_SUCCEED(gr_mat_mul_classical(D, A, B, ctx));
        GR_MUST_SUCCEED(gr_mat_sub(C, C, D, ctx));
        err = _t_m_norm(C, ctx, CC);
        if (!(err <= ldexp(len * _t_m_norm(A, ctx, CC) * _t_m_norm(B, ctx, CC), -53 * n + 8)))
        {
            flint_printf("FAIL: mat_mul (arb route) %g\n", err);
            flint_abort();
        }
        gr_mat_clear(A, ctx);
        gr_mat_clear(B, ctx);
        gr_mat_clear(C, ctx);
        gr_mat_clear(D, ctx);
        gr_ctx_clear(ctx);
        gr_ctx_clear(CC);
    }

    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx, CC;
        int n, cplx, strong, which;
        slong sz;

        n = 1 + n_randint(state, DFLOAT_MAX_N);
        cplx = n_randint(state, 2);
        strong = (n_randint(state, 4) == 0) ? DFLOAT_STRONG : 0;
        GR_MUST_SUCCEED(gr_ctx_init_dfloat(ctx, n, (cplx ? DFLOAT_COMPLEX : 0) | strong));
        gr_ctx_init_complex_acb(CC, 53 * n + 100);
        sz = ctx->sizeof_elem;
        which = n_randint(state, 5);

        if (which == 0)
        {
            /* products, entrywise against the exact product */
            slong a, b, c, i, j, l;
            int alias = n_randint(state, 3);
            gr_mat_t A, B, C;
            acb_t s, x, y;
            arb_t e, m, t;

            a = n_randint(state, 4) ? n_randint(state, 20) : n_randint(state, 70);
            b = n_randint(state, 4) ? n_randint(state, 20) : n_randint(state, 70);
            c = n_randint(state, 4) ? n_randint(state, 20) : n_randint(state, 70);
            if (alias == 1) c = b;
            if (alias == 2) a = b;

            gr_mat_init(A, a, b, ctx);
            gr_mat_init(B, b, c, ctx);
            gr_mat_init(C, a, c, ctx);
            _t_m_rand(A, n, cplx, state, ctx);
            _t_m_rand(B, n, cplx, state, ctx);

            if (alias == 1)
            {
                GR_MUST_SUCCEED(gr_mat_set(C, A, ctx));
                GR_MUST_SUCCEED(gr_mat_mul(C, C, B, ctx));
            }
            else if (alias == 2)
            {
                GR_MUST_SUCCEED(gr_mat_set(C, B, ctx));
                GR_MUST_SUCCEED(gr_mat_mul(C, A, C, ctx));
            }
            else
                GR_MUST_SUCCEED(gr_mat_mul(C, A, B, ctx));

            acb_init(s); acb_init(x); acb_init(y);
            arb_init(e); arb_init(m); arb_init(t);
            for (i = 0; i < a; i++)
                for (j = 0; j < c; j++)
                {
                    acb_zero(s);
                    arb_zero(m);
                    for (l = 0; l < b; l++)
                    {
                        GR_MUST_SUCCEED(gr_set_other(x, GR_MAT_ENTRY(A, i, l, sz), ctx, CC));
                        GR_MUST_SUCCEED(gr_set_other(y, GR_MAT_ENTRY(B, l, j, sz), ctx, CC));
                        acb_mul(x, x, y, 2000);
                        acb_add(s, s, x, 2000);
                        acb_abs(t, x, 2000);
                        arb_add(m, m, t, 2000);
                    }
                    GR_MUST_SUCCEED(gr_set_other(x, GR_MAT_ENTRY(C, i, j, sz), ctx, CC));
                    acb_sub(x, x, s, 2000);
                    acb_abs(e, x, 2000);
                    arb_mul_si(m, m, b + 2, 2000);
                    arb_mul_2exp_si(m, m, -53 * n + 6);
                    /* (an absolute allowance for underflow) */
                    arb_set_d(t, 1e-290);
                    arb_add(m, m, t, 2000);
                    if (!arb_le(e, m))
                    {
                        flint_printf("FAIL: mat_mul (n = %d, cplx = %d, %wd x %wd x %wd, alias = %d, entry %wd %wd)\n", n, cplx, a, b, c, alias, i, j);
                        arb_printd(e, 10); flint_printf("\n");
                        arb_printd(m, 10); flint_printf("\n");
                        flint_abort();
                    }
                }
            acb_clear(s); acb_clear(x); acb_clear(y);
            arb_clear(e); arb_clear(m); arb_clear(t);
            gr_mat_clear(A, ctx);
            gr_mat_clear(B, ctx);
            gr_mat_clear(C, ctx);
        }
        else if (which == 1)
        {
            gr_ptr tol = gr_heap_init(ctx);
            GR_MUST_SUCCEED(gr_one(tol, ctx));
            GR_MUST_SUCCEED(gr_mul_2exp_si(tol, tol, -53 * n + 14, ctx));
            gr_mat_test_approx_mul_max_norm((gr_method_mat_binary_op) gr_mat_mul, tol, state, 2, 40, ctx);
            gr_heap_clear(tol, ctx);
        }
        else if (which == 2)
        {
            /* LU: P A = L U */
            slong len = n_randint(state, 4) ? n_randint(state, 24) : n_randint(state, 60), rank, i, j;
            slong * P = flint_malloc(sizeof(slong) * FLINT_MAX(len, 1));
            gr_mat_t A, LU, L, U, PA;
            double nA, res;

            gr_mat_init(A, len, len, ctx);
            gr_mat_init(LU, len, len, ctx);
            gr_mat_init(L, len, len, ctx);
            gr_mat_init(U, len, len, ctx);
            gr_mat_init(PA, len, len, ctx);
            _t_m_rand_wellcond(A, n, cplx, state, ctx);

            if (n_randint(state, 2))
            {
                GR_MUST_SUCCEED(gr_mat_set(LU, A, ctx));
                GR_MUST_SUCCEED(gr_mat_lu(&rank, P, LU, LU, 0, ctx));
            }
            else
                GR_MUST_SUCCEED(gr_mat_lu(&rank, P, LU, A, 0, ctx));

            if (rank != len)
            {
                flint_printf("FAIL: lu rank (n = %d, cplx = %d, len = %wd, rank = %wd)\n", n, cplx, len, rank);
                flint_abort();
            }

            for (i = 0; i < len; i++)
            {
                for (j = 0; j < len; j++)
                {
                    if (j < i)
                        GR_MUST_SUCCEED(gr_set(GR_MAT_ENTRY(L, i, j, sz), GR_MAT_ENTRY(LU, i, j, sz), ctx));
                    else
                        GR_MUST_SUCCEED(gr_set(GR_MAT_ENTRY(U, i, j, sz), GR_MAT_ENTRY(LU, i, j, sz), ctx));
                }
                GR_MUST_SUCCEED(gr_one(GR_MAT_ENTRY(L, i, i, sz), ctx));
                GR_MUST_SUCCEED(_gr_vec_set(GR_MAT_ENTRY(PA, i, 0, sz), GR_MAT_ENTRY(A, P[i], 0, sz), len, ctx));
            }

            nA = _t_m_norm(A, ctx, CC);
            res = _t_m_residual(L, U, PA, ctx, CC);
            if (!(res <= ldexp(nA * (len + 1) * (len + 1), -53 * n + 10)))
            {
                flint_printf("FAIL: lu residual (n = %d, cplx = %d, len = %wd): %g vs %g\n", n, cplx, len, res, nA);
                flint_abort();
            }

            gr_mat_clear(A, ctx);
            gr_mat_clear(LU, ctx);
            gr_mat_clear(L, ctx);
            gr_mat_clear(U, ctx);
            gr_mat_clear(PA, ctx);
            flint_free(P);
        }
        else if (which == 3)
        {
            /* the rank of rank-deficient matrices: r random rows and
               zero rows, some columns zero, shuffled (the elimination
               keeps the zero rows exactly zero); with rank_check, the
               rank is 0 unless the matrix has full rank */
            slong rows = n_randint(state, 40), cols = n_randint(state, 40), r, i, j, rank1, rank2;
            slong * P = flint_malloc(sizeof(slong) * FLINT_MAX(rows, 1));
            int rank_check = n_randint(state, 2);
            gr_mat_t A, LU;

            r = n_randint(state, rows + 1);
            gr_mat_init(A, rows, cols, ctx);
            gr_mat_init(LU, rows, cols, ctx);
            for (i = 0; i < r; i++)
                for (j = 0; j < cols; j++)
                {
                    _t_m_rand_entry(GR_MAT_ENTRY(A, i, j, sz), n, cplx, 0, state, ctx);
                    GR_MUST_SUCCEED(gr_add_si(GR_MAT_ENTRY(A, i, j, sz), GR_MAT_ENTRY(A, i, j, sz), (i == j) ? 5 : 0, ctx));
                }
            for (j = 0; j < cols; j++)
                if (n_randint(state, 8) == 0)
                    for (i = 0; i < rows; i++)
                        GR_MUST_SUCCEED(gr_zero(GR_MAT_ENTRY(A, i, j, sz), ctx));
            for (i = 0; i < rows; i++)
                _gr_vec_swap(GR_MAT_ENTRY(A, i, 0, sz), GR_MAT_ENTRY(A, n_randint(state, rows), 0, sz), cols, ctx);

            GR_MUST_SUCCEED(gr_mat_lu(&rank1, P, LU, A, rank_check, ctx));
            GR_MUST_SUCCEED(gr_mat_lu_classical(&rank2, P, LU, A, rank_check, ctx));
            if (rank1 != rank2)
            {
                flint_printf("FAIL: lu rank (n = %d, cplx = %d, %wd x %wd, r = %wd, rank_check = %d: %wd vs %wd)\n", n, cplx, rows, cols, r, rank_check, rank1, rank2);
                flint_abort();
            }
            gr_mat_clear(A, ctx);
            gr_mat_clear(LU, ctx);
            flint_free(P);
        }
        else
        {
            /* triangular solving and solving */
            slong len = n_randint(state, 4) ? n_randint(state, 24) : n_randint(state, 60);
            slong m = n_randint(state, 3) ? n_randint(state, 4) + 1 : n_randint(state, 30);
            int kind = n_randint(state, 3), unit = n_randint(state, 2), alias = n_randint(state, 2);
            slong i, j;
            gr_mat_t A, B, X;
            double nA, nX, res;

            gr_mat_init(A, len, len, ctx);
            gr_mat_init(B, len, m, ctx);
            gr_mat_init(X, len, m, ctx);
            _t_m_rand_wellcond(A, n, cplx, state, ctx);
            _t_m_rand(B, n, cplx, state, ctx);

            if (kind < 2)
            {
                /* make A triangular, with a dominant diagonal */
                gr_ptr t;
                GR_TMP_INIT(t, ctx);
                for (i = 0; i < len; i++)
                {
                    for (j = 0; j < len; j++)
                        if ((kind == 0 && j > i) || (kind == 1 && j < i))
                            GR_MUST_SUCCEED(gr_zero(GR_MAT_ENTRY(A, i, j, sz), ctx));
                    GR_MUST_SUCCEED(gr_set_si(t, unit ? 1 : 2 * len + 1, ctx));
                    GR_MUST_SUCCEED(gr_set(GR_MAT_ENTRY(A, i, i, sz), t, ctx));
                }
                GR_TMP_CLEAR(t, ctx);
                if (alias)
                    GR_MUST_SUCCEED(gr_mat_set(X, B, ctx));
                if (kind == 0)
                    GR_MUST_SUCCEED(gr_mat_nonsingular_solve_tril(X, A, alias ? X : B, unit, ctx));
                else
                    GR_MUST_SUCCEED(gr_mat_nonsingular_solve_triu(X, A, alias ? X : B, unit, ctx));
            }
            else
            {
                GR_MUST_SUCCEED(gr_mat_nonsingular_solve(X, A, B, ctx));
            }

            nA = _t_m_norm(A, ctx, CC);
            nX = _t_m_norm(X, ctx, CC);
            res = _t_m_residual(A, X, B, ctx, CC);
            if (!(res <= ldexp((nA * nX + 1e-300) * (len + 1) * (len + 1), -53 * n + 10)))
            {
                flint_printf("FAIL: solve residual (n = %d, cplx = %d, kind = %d, unit = %d, len = %wd, m = %wd): %g vs %g %g\n", n, cplx, kind, unit, len, m, res, nA, nX);
                flint_abort();
            }

            gr_mat_clear(A, ctx);
            gr_mat_clear(B, ctx);
            gr_mat_clear(X, ctx);
        }

        gr_ctx_clear(ctx);
        gr_ctx_clear(CC);
    }

    TEST_FUNCTION_END(state);
}
