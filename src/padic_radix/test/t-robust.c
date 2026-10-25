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
#include "gr_vec.h"
#include "gr_mat.h"
#include "padic_radix.h"

/*
    Inputs with huge valuation, precision or exponent differences (up to
    WORD_MAX / 8, or far beyond the precision of the context). In contexts
    with a finite precision, intermediate results must be bounded in terms
    of the working precision (so the operations succeed quickly), and the
    block kernels must give up on (or not be chosen for) frames spread far
    beyond the precision; in contexts without finite precision, results
    beyond the size_limit of the context give GR_UNABLE. A regression shows
    up as an allocation failure (abort) or a timeout.

    (In a context with only an absolute precision a, an exact element of
    valuation v far below a has a working precision of a - v digits, which
    is the user's choice: those inputs are avoided here.)
*/

#define ROBUST_SIZE_LIMIT (WORD(1) << 20)

/* u p^v + O(p^N), u a random nonzero integer (N = EXACT for exact) */
static void
_robust_elem(padic_radix_struct * x, slong v, slong N, int neg, flint_rand_t state, gr_ctx_t ctx)
{
    fmpz_t f;
    fmpz_init(f);
    fmpz_randtest_not_zero(f, state, 30);
    fmpz_abs(f, f);
    if (neg)
        fmpz_neg(f, f);
    GR_MUST_SUCCEED(padic_radix_exact_set_fmpz(x, f, ctx));
    x->v += v;
    x->N = (N == PADIC_RADIX_EXACT) ? N : x->v + N;
    if (_padic_radix_finalize(x, ctx) != GR_SUCCESS)
        GR_MUST_SUCCEED(padic_radix_one(x, ctx));
    fmpz_clear(f);
}

/* the balls of x and y intersect */
static int
_robust_agree(const padic_radix_struct * x, const padic_radix_struct * y, gr_ctx_t ctx)
{
    padic_radix_t t;
    int ok;
    padic_radix_init(t, ctx);
    ok = (padic_radix_sub(t, x, y, ctx) != GR_SUCCESS) || (padic_radix_is_zero(t, ctx) != T_FALSE);
    padic_radix_clear(t, ctx);
    return ok;
}

static int
_robust_identical(const padic_radix_struct * x, const padic_radix_struct * y, gr_ctx_t ctx)
{
    return radix_integer_equal(&x->u, &y->u, PADIC_RADIX_CTX_RADIX(ctx))
        && (x->u.size == 0 || x->v == y->v) && x->N == y->N;
}

TEST_FUNCTION_START(padic_radix_robust, state)
{
    slong iter, big = WORD_MAX / 8;

    /* scalar operations with exponents +-WORD_MAX/8 in all kinds of
       contexts: success or GR_UNABLE */
    for (iter = 0; iter < 30 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        padic_radix_t x, y, z;
        slong rel, pabs, vx, vy, Nx, Ny;
        int k;

        switch (n_randint(state, 4))
        {
            case 0: rel = 20; pabs = PADIC_RADIX_PREC_INF; break;
            case 1: rel = PADIC_RADIX_PREC_INF; pabs = 20; break;
            case 2: rel = PADIC_RADIX_PREC_INF; pabs = PADIC_RADIX_PREC_INF; break;
            default: rel = 30; pabs = 50; break;
        }

        gr_ctx_init_padic_radix(ctx, n_nth_prime(1 + n_randint(state, 5)), rel, pabs, PADIC_RADIX_SIGNED);
        ctx->size_limit = ROBUST_SIZE_LIMIT;
        padic_radix_init(x, ctx);
        padic_radix_init(y, ctx);
        padic_radix_init(z, ctx);

        vx = n_randint(state, 2) ? big : -big;
        vy = n_randint(state, 3) ? 0 : (n_randint(state, 2) ? big / 2 : -big / 2);
        Nx = n_randint(state, 2) ? PADIC_RADIX_EXACT : 20 + n_randint(state, 1000);
        Ny = n_randint(state, 2) ? PADIC_RADIX_EXACT : 20 + n_randint(state, 1000);
        if (rel == PADIC_RADIX_PREC_INF && pabs != PADIC_RADIX_PREC_INF)
        {
            if (Nx == PADIC_RADIX_EXACT)
                vx = big;
            if (Ny == PADIC_RADIX_EXACT)
                vy = FLINT_ABS(vy);
        }
        _robust_elem(x, vx, Nx, n_randint(state, 2), state, ctx);
        _robust_elem(y, vy, Ny, n_randint(state, 2), state, ctx);

        for (k = 0; k < 9; k++)
        {
            int s;
            switch (k)
            {
                case 0: s = padic_radix_add(z, x, y, ctx); break;
                case 1: s = padic_radix_sub(z, y, x, ctx); break;
                case 2: s = padic_radix_mul(z, x, y, ctx); break;
                case 3: s = padic_radix_div(z, x, y, ctx); break;
                case 4: s = padic_radix_inv(z, x, ctx); break;
                case 5: s = padic_radix_sqrt(z, x, ctx); break;
                case 6: s = padic_radix_pow_ui(z, x, 1000003, ctx); break;
                case 7: s = padic_radix_pow_ui(z, y, UWORD(1) << (FLINT_BITS - 4), ctx); break;
                default: s = padic_radix_addmul(z, x, y, ctx); break;
            }
            if (s & GR_TEST_FAIL)
            {
                flint_printf("FAIL (scalar op %d)\n", k);
                flint_abort();
            }
        }

        padic_radix_clear(x, ctx);
        padic_radix_clear(y, ctx);
        padic_radix_clear(z, ctx);
        gr_ctx_clear(ctx);
    }

    /* huge powers in a context with finite precision saturate (they are
       computed with truncated products, whatever the exponent) */
    {
        gr_ctx_t ctx;
        padic_radix_t x, y;
        fmpz_t n;

        gr_ctx_init_padic_radix(ctx, 7, 3, PADIC_RADIX_PREC_INF, 0);
        padic_radix_init(x, ctx);
        padic_radix_init(y, ctx);
        fmpz_init(n);
        fmpz_ui_pow_ui(n, 10, 100);

        GR_MUST_SUCCEED(padic_radix_exact_set_ui(x, 7, ctx));
        GR_MUST_SUCCEED(padic_radix_pow_fmpz(y, x, n, ctx));
        if (y->u.size != 0 || y->N != PADIC_RADIX_ERR_MAX)
        {
            flint_printf("FAIL (7^(10^100) with relative precision 3)\n");
            flint_abort();
        }

        GR_MUST_SUCCEED(padic_radix_exact_set_ui(x, 14, ctx));
        GR_MUST_SUCCEED(padic_radix_pow_fmpz(y, x, n, ctx));
        if (y->u.size != 0 || y->N != PADIC_RADIX_ERR_MAX)
        {
            flint_printf("FAIL (14^(10^100) with relative precision 3)\n");
            flint_abort();
        }

        padic_radix_clear(x, ctx);
        padic_radix_clear(y, ctx);
        fmpz_clear(n);
        gr_ctx_clear(ctx);
    }

    /* exact sums with a huge valuation gap: beyond the size limit of an
       exact context (GR_UNABLE), but the delayed dot product falls back to
       the addmul chain, which gets exact cancellations right */
    {
        gr_ctx_t ctx;
        padic_radix_struct * a, * b;
        padic_radix_t r;

        gr_ctx_init_padic_radix(ctx, n_nth_prime(1 + n_randint(state, 5)), PADIC_RADIX_PREC_INF, PADIC_RADIX_PREC_INF, PADIC_RADIX_SIGNED);
        ctx->size_limit = ROBUST_SIZE_LIMIT;
        a = gr_heap_init_vec(3, ctx);
        b = gr_heap_init_vec(3, ctx);
        padic_radix_init(r, ctx);

        GR_MUST_SUCCEED(padic_radix_exact_set_si(a + 0, 5, ctx));
        a[0].v = big;
        GR_MUST_SUCCEED(padic_radix_exact_set_si(a + 1, -5, ctx));
        a[1].v = big;
        GR_MUST_SUCCEED(padic_radix_exact_set_si(a + 2, 3, ctx));
        for (iter = 0; iter < 3; iter++)
            GR_MUST_SUCCEED(padic_radix_exact_set_si(b + iter, 2, ctx));

        if (padic_radix_add(r, a + 0, a + 2, ctx) != GR_UNABLE)
        {
            flint_printf("FAIL (exact add with a huge gap)\n");
            flint_abort();
        }

        if (padic_radix_dot_strided_delayed(r, NULL, 0, a, 1, b, 1, 3, ctx) != GR_SUCCESS)
        {
            flint_printf("FAIL (exact dot with a huge gap)\n");
            flint_abort();
        }

        GR_MUST_SUCCEED(padic_radix_exact_set_si(a + 0, 6, ctx));
        if (r->N != PADIC_RADIX_EXACT || !_robust_identical(r, a + 0, ctx))
        {
            flint_printf("FAIL (exact dot with a huge gap: value)\n");
            flint_abort();
        }

        gr_heap_clear_vec(a, 3, ctx);
        gr_heap_clear_vec(b, 3, ctx);
        padic_radix_clear(r, ctx);
        gr_ctx_clear(ctx);
    }

    /* polynomial products with valuations spread far beyond the precision:
       the Kronecker substitution / modular routes give up, and the
       dispatcher gives the result of the packed classical kernel */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        padic_radix_struct * a, * b, * c1, * c2;
        slong len1, len2, i, shape, rel = 10 + n_randint(state, 40);
        int s;

        gr_ctx_init_padic_radix(ctx, n_nth_prime(1 + n_randint(state, 5)), rel, PADIC_RADIX_PREC_INF, PADIC_RADIX_SIGNED);
        len1 = 64 + n_randint(state, 200);
        len2 = 64 + n_randint(state, 200);
        shape = n_randint(state, 3);
        a = gr_heap_init_vec(len1, ctx);
        b = gr_heap_init_vec(len2, ctx);
        c1 = gr_heap_init_vec(len1 + len2 - 1, ctx);
        c2 = gr_heap_init_vec(len1 + len2 - 1, ctx);

        for (i = 0; i < len1; i++)
            _robust_elem(a + i, (shape == 0) ? (i % 2) * 1000000 : (shape == 1) ? 1000 * i
                : (slong) n_randint(state, 2000001) - 1000000, 20, 0, state, ctx);
        for (i = 0; i < len2; i++)
            _robust_elem(b + i, (shape == 0) ? ((i + 1) % 2) * 1000000 : (shape == 1) ? 1000 * i
                : (slong) n_randint(state, 2000001) - 1000000, 20, 0, state, ctx);

        s = _padic_radix_poly_mulmid_KS(c1, a, len1, b, len2, 0, len1 + len2 - 1, ctx);
        s |= _padic_radix_poly_mulmid_mod(c1, a, len1, b, len2, 0, len1 + len2 - 1, ctx);
        if (s != GR_UNABLE)
        {
            flint_printf("FAIL (poly frames routes with spread valuations, shape %wd)\n", shape);
            flint_abort();
        }

        GR_MUST_SUCCEED(_padic_radix_poly_mullow(c1, a, len1, b, len2, len1 + len2 - 1, ctx));
        GR_MUST_SUCCEED(_padic_radix_poly_mulmid_packed(c2, a, len1, b, len2, 0, len1 + len2 - 1, ctx));
        for (i = 0; i < len1 + len2 - 1; i++)
            if (!_robust_identical(c1 + i, c2 + i, ctx))
            {
                flint_printf("FAIL (poly dispatcher with spread valuations)\n");
                flint_abort();
            }

        gr_heap_clear_vec(a, len1, ctx);
        gr_heap_clear_vec(b, len2, ctx);
        gr_heap_clear_vec(c1, len1 + len2 - 1, ctx);
        gr_heap_clear_vec(c2, len1 + len2 - 1, ctx);
        gr_ctx_clear(ctx);
    }

    /* a few coefficients of huge valuation (at the ends, or sparse): the
       few outputs whose terms all have huge valuations are computed by dot
       products, the others by the transforms; identical to the packed
       kernel */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        padic_radix_struct * a, * b, * c1, * c2;
        slong len1, len2, i, alg, rel = 10 + n_randint(state, 40), per = 20 + n_randint(state, 40);
        int s;

        gr_ctx_init_padic_radix(ctx, n_nth_prime(1 + n_randint(state, 5)), rel, PADIC_RADIX_PREC_INF, PADIC_RADIX_SIGNED);
        len1 = 64 + n_randint(state, 200);
        len2 = 64 + n_randint(state, 200);
        a = gr_heap_init_vec(len1, ctx);
        b = gr_heap_init_vec(len2, ctx);
        c1 = gr_heap_init_vec(len1 + len2 - 1, ctx);
        c2 = gr_heap_init_vec(len1 + len2 - 1, ctx);

        for (i = 0; i < len1; i++)
            _robust_elem(a + i, (i < 2 || i % per == 0) ? 1000000 : 0, 20,
                n_randint(state, 2), state, ctx);
        for (i = 0; i < len2; i++)
            _robust_elem(b + i, (i < 2 || i > len2 - 3) ? 1000000 : 0,
                n_randint(state, 4) ? 20 : PADIC_RADIX_EXACT, n_randint(state, 2), state, ctx);

        GR_MUST_SUCCEED(_padic_radix_poly_mulmid_packed(c2, a, len1, b, len2, 0, len1 + len2 - 1, ctx));

        for (alg = 0; alg < 3; alg++)
        {
            if (alg == 0)
                s = _padic_radix_poly_mulmid_KS(c1, a, len1, b, len2, 0, len1 + len2 - 1, ctx);
            else if (alg == 1)
                s = _padic_radix_poly_mulmid_mod(c1, a, len1, b, len2, 0, len1 + len2 - 1, ctx);
            else
                s = _padic_radix_poly_mullow(c1, a, len1, b, len2, len1 + len2 - 1, ctx);

            if (s == GR_UNABLE && alg == 1)
                continue;

            if (s != GR_SUCCESS)
            {
                flint_printf("FAIL (poly with a few huge valuations: alg %wd gives up)\n", alg);
                flint_abort();
            }

            for (i = 0; i < len1 + len2 - 1; i++)
                if (!_robust_identical(c1 + i, c2 + i, ctx))
                {
                    flint_printf("FAIL (poly with a few huge valuations: alg %wd, coefficient %wd)\n", alg, i);
                    flint_abort();
                }
        }

        gr_heap_clear_vec(a, len1, ctx);
        gr_heap_clear_vec(b, len2, ctx);
        gr_heap_clear_vec(c1, len1 + len2 - 1, ctx);
        gr_heap_clear_vec(c2, len1 + len2 - 1, ctx);
        gr_ctx_clear(ctx);
    }

    /* matrix products with spread frames (rows of A and columns of B
       with valuations 0 and 10^6 in alternating positions) */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        gr_mat_t A, B, C1, C2;
        slong n = 16 + n_randint(state, 32), i, j;

        gr_ctx_init_padic_radix(ctx, n_nth_prime(1 + n_randint(state, 5)), 20, PADIC_RADIX_PREC_INF, PADIC_RADIX_SIGNED);
        gr_mat_init(A, n, n, ctx);
        gr_mat_init(B, n, n, ctx);
        gr_mat_init(C1, n, n, ctx);
        gr_mat_init(C2, n, n, ctx);

        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
            {
                _robust_elem(gr_mat_entry_ptr(A, i, j, ctx), ((i + j) % 2) * 1000000, 20, 0, state, ctx);
                _robust_elem(gr_mat_entry_ptr(B, i, j, ctx), ((i + j + 1) % 2) * 1000000, 20, 0, state, ctx);
            }

        if (padic_radix_mat_mul_fmpz(C1, A, B, ctx) != GR_UNABLE
            || padic_radix_mat_mul_mpn_mod(C1, A, B, ctx) != GR_UNABLE)
        {
            flint_printf("FAIL (matrix fmpz/mpn_mod routes with spread frames)\n");
            flint_abort();
        }

        GR_MUST_SUCCEED(gr_mat_mul(C1, A, B, ctx));
        GR_MUST_SUCCEED(padic_radix_mat_mul_packed(C2, A, B, ctx));
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                if (!_robust_identical(gr_mat_entry_ptr(C1, i, j, ctx), gr_mat_entry_ptr(C2, i, j, ctx), ctx))
                {
                    flint_printf("FAIL (matrix dispatcher with spread frames)\n");
                    flint_abort();
                }

        gr_mat_clear(A, ctx);
        gr_mat_clear(B, ctx);
        gr_mat_clear(C1, ctx);
        gr_mat_clear(C2, ctx);
        gr_ctx_clear(ctx);
    }

    /* characteristic polynomials: valuations spread beyond the precision
       (the lift reduces the entries), and a huge scaling (the lift gives
       up); the result agrees with Berkowitz's algorithm */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        gr_mat_t A;
        padic_radix_struct * c1, * c2;
        slong n = 2 + n_randint(state, 8), i, j;
        int scaled = n_randint(state, 2);

        gr_ctx_init_padic_radix(ctx, n_nth_prime(1 + n_randint(state, 5)), 30, PADIC_RADIX_PREC_INF, PADIC_RADIX_SIGNED);
        gr_mat_init(A, n, n, ctx);
        c1 = gr_heap_init_vec(n + 1, ctx);
        c2 = gr_heap_init_vec(n + 1, ctx);

        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                _robust_elem(gr_mat_entry_ptr(A, i, j, ctx), n_randint(state, 2) ? 1000000 : 0,
                    20, n_randint(state, 2), state, ctx);
        if (scaled)
            _robust_elem(gr_mat_entry_ptr(A, 0, n - 1, ctx), -1000000, PADIC_RADIX_EXACT, 0, state, ctx);

        if (scaled && padic_radix_mat_charpoly_lift(c1, A, ctx) != GR_UNABLE)
        {
            flint_printf("FAIL (charpoly lift with a huge scaling)\n");
            flint_abort();
        }

        GR_MUST_SUCCEED(padic_radix_mat_charpoly(c1, A, ctx));
        GR_MUST_SUCCEED(_padic_radix_mat_charpoly_berkowitz(c2, A, ctx));
        for (i = 0; i <= n; i++)
            if (!_robust_agree(c1 + i, c2 + i, ctx))
            {
                flint_printf("FAIL (charpoly with spread valuations)\n");
                flint_abort();
            }

        gr_mat_clear(A, ctx);
        gr_heap_clear_vec(c1, n + 1, ctx);
        gr_heap_clear_vec(c2, n + 1, ctx);
        gr_ctx_clear(ctx);
    }

    /* linear algebra and block kernels with exponents +-WORD_MAX/8 and
       +-WORD_MAX/32 (beyond / within the range of the packed kernels) and
       huge precisions: success or GR_UNABLE */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        gr_mat_t A, B, C;
        padic_radix_struct * a, * b, * c;
        gr_ptr d;
        slong n = 2 + n_randint(state, 10), i, j, rel, pabs, vb, s, sb;
        int exact = n_randint(state, 2);

        switch (n_randint(state, 3))
        {
            case 0: rel = 20; pabs = PADIC_RADIX_PREC_INF; break;
            case 1: rel = PADIC_RADIX_PREC_INF; pabs = 20; break;
            default: rel = PADIC_RADIX_PREC_INF; pabs = PADIC_RADIX_PREC_INF; break;
        }
        vb = n_randint(state, 2) ? big : WORD_MAX / 32;
        /* (no exact elements of huge negative valuation with only an
           absolute precision: their working precision is huge) */
        sb = (exact && rel == PADIC_RADIX_PREC_INF && pabs != PADIC_RADIX_PREC_INF) ? 1 : -1;

        gr_ctx_init_padic_radix(ctx, n_nth_prime(1 + n_randint(state, 5)), rel, pabs, PADIC_RADIX_SIGNED);
        ctx->size_limit = ROBUST_SIZE_LIMIT;
        gr_mat_init(A, n, n, ctx);
        gr_mat_init(B, n, n, ctx);
        gr_mat_init(C, n, n, ctx);
        a = gr_heap_init_vec(4 * n, ctx);
        b = gr_heap_init_vec(4 * n, ctx);
        c = gr_heap_init_vec(8 * n, ctx);
        d = gr_heap_init(ctx);

        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
            {
                slong w = n_randint(state, 3);
                w = (w == 0) ? vb : (w == 1) ? sb * vb : 0;
                _robust_elem(gr_mat_entry_ptr(A, i, j, ctx), w, exact ? PADIC_RADIX_EXACT : 20, n_randint(state, 2), state, ctx);
                _robust_elem(gr_mat_entry_ptr(B, i, j, ctx), (sb == 1) ? w : -w, exact ? PADIC_RADIX_EXACT : 1000, n_randint(state, 2), state, ctx);
            }
        for (i = 0; i < 4 * n; i++)
        {
            slong w = n_randint(state, 3);
            w = (w == 0) ? vb : (w == 1) ? sb * vb : 0;
            _robust_elem(a + i, w, exact ? PADIC_RADIX_EXACT : 20, n_randint(state, 2), state, ctx);
            _robust_elem(b + i, (sb == 1) ? w : -w, exact ? PADIC_RADIX_EXACT : 20, n_randint(state, 2), state, ctx);
        }

        s = gr_mat_mul(C, A, B, ctx);
        s |= padic_radix_mat_mul_packed(C, A, B, ctx);
        s |= padic_radix_mat_mul_fmpz(C, A, B, ctx);
        s |= padic_radix_mat_mul_mpn_mod(C, A, B, ctx);
        s |= gr_mat_det(d, A, ctx);
        s |= gr_mat_nonsingular_solve(C, A, B, ctx);
        s |= _gr_mat_charpoly(c, A, ctx);
        s |= padic_radix_mat_charpoly_lift(c, A, ctx);
        s |= _padic_radix_poly_mullow(c, a, 4 * n, b, 4 * n, 8 * n - 1, ctx);
        s |= _padic_radix_poly_mulmid_KS(c, a, 4 * n, b, 4 * n, 0, 8 * n - 1, ctx);
        s |= _padic_radix_poly_mulmid_mod(c, a, 4 * n, b, 4 * n, 0, 8 * n - 1, ctx);
        s |= _padic_radix_poly_mulmid_packed(c, a, 4 * n, b, 4 * n, 0, 8 * n - 1, ctx);
        s |= padic_radix_dot(d, NULL, 0, a, b, 4 * n, ctx);

        if (s & GR_TEST_FAIL)
        {
            flint_printf("FAIL (block operations with huge exponents)\n");
            flint_abort();
        }

        gr_mat_clear(A, ctx);
        gr_mat_clear(B, ctx);
        gr_mat_clear(C, ctx);
        gr_heap_clear_vec(a, 4 * n, ctx);
        gr_heap_clear_vec(b, 4 * n, ctx);
        gr_heap_clear_vec(c, 8 * n, ctx);
        gr_heap_clear(d, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
