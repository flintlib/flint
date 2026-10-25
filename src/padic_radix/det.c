/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mat.h"
#include "perm.h"
#include "gr.h"
#include "gr_mat.h"
#include "padic_radix.h"

#define ENTRY(A, i, j) ((padic_radix_struct *) gr_mat_entry_ptr(A, i, j, ctx))
#define CENTRY(A, i, j) ((const padic_radix_struct *) gr_mat_entry_srcptr(A, i, j, ctx))

/*
    Pivot of minimal valuation in a column. Elimination with such pivots
    involves only multipliers of nonnegative valuation, so that it does not
    amplify the absolute errors of the entries below the pivot row.
*/
int
padic_radix_mat_find_nonzero_pivot(slong * pivot_row, gr_mat_t mat,
    slong start_row, slong end_row, slong column, gr_ctx_t ctx)
{
    slong i, best = -1, bv = WORD_MAX;
    int unknown = 0;

    if (end_row <= start_row)
        return GR_DOMAIN;

    for (i = start_row; i < end_row; i++)
    {
        const padic_radix_struct * x = ENTRY(mat, i, column);

        if (x->u.size != 0)
        {
            if (x->v < bv)
            {
                bv = x->v;
                best = i;
            }
        }
        else if (x->N != PADIC_RADIX_EXACT)
            unknown = 1;
    }

    *pivot_row = best;

    if (best >= 0)
        return GR_SUCCESS;

    return unknown ? GR_UNABLE : GR_DOMAIN;
}

/* a + b, or WORD_MIN when the sum is below -PADIC_RADIX_ERR_MAX (nothing
   known); sums above PADIC_RADIX_ERR_MAX are clamped (still lower bounds) */
static slong
_sat_add(slong a, slong b)
{
    if (a == WORD_MIN || b == WORD_MIN)
        return WORD_MIN;
    a = FLINT_MIN(a, PADIC_RADIX_ERR_MAX);
    b = FLINT_MIN(b, PADIC_RADIX_ERR_MAX);
    a += b;     /* no overflow: both are in [-ERR_MAX, ERR_MAX] */
    if (a < -PADIC_RADIX_ERR_MAX)
        return WORD_MIN;
    return FLINT_MIN(a, PADIC_RADIX_ERR_MAX);
}

/*
    Lower bound for the precision of the determinant of a matrix of zeros
    known to O(p^N_ij): the minimum over permutations of sum N_(i,s(i)) is
    at least both the sum of the row minima and the sum of the column
    minima. Exact zeros count as infinite precision.
*/
static slong
_zero_block_prec(const gr_mat_t A, slong k, slong n, gr_ctx_t ctx)
{
    slong i, j, rs = 0, cs = 0;

    for (i = k; i < n; i++)
    {
        slong m = PADIC_RADIX_EXACT;
        for (j = k; j < n; j++)
            m = FLINT_MIN(m, CENTRY(A, i, j)->N);
        if (m == PADIC_RADIX_EXACT)
            return PADIC_RADIX_EXACT;
        rs = _sat_add(rs, m);
    }

    for (j = k; j < n; j++)
    {
        slong m = PADIC_RADIX_EXACT;
        for (i = k; i < n; i++)
            m = FLINT_MIN(m, CENTRY(A, i, j)->N);
        if (m == PADIC_RADIX_EXACT)
            return PADIC_RADIX_EXACT;
        cs = _sat_add(cs, m);
    }

    return FLINT_MAX(rs, cs);
}

/*
    Gaussian elimination with full pivoting on valuations: the pivot is an
    entry of minimal valuation in the whole remaining submatrix. Then all
    multipliers and all entries of the pivot row have valuation at least
    that of the pivot, the Schur complements keep the absolute precision
    of the input, the pivots have the valuations of the elementary divisors
    a_1 <= ... <= a_n, and the determinant is obtained to precision
    N + a_1 + ... + a_(n-1) for entries known to O(p^N), which is optimal
    (it is that of the differential d det = tr(adj(A) dA)).
*/
int
padic_radix_mat_det_lu_full(padic_radix_t res, const gr_mat_t A0, gr_ctx_t ctx)
{
    slong n = A0->r, i, j, k, r, c;
    gr_mat_t A;
    padic_radix_t inv, f;
    int sign = 0, status = GR_SUCCESS;

    if (n != A0->c)
        return GR_DOMAIN;

    if (n == 0)
        return padic_radix_one(res, ctx);

    gr_mat_init(A, n, n, ctx);
    status |= gr_mat_set(A, A0, ctx);
    padic_radix_init(inv, ctx);
    padic_radix_init(f, ctx);
    status |= padic_radix_one(res, ctx);

    for (k = 0; k < n; k++)
    {
        slong bv = WORD_MAX;

        r = c = -1;
        for (i = k; i < n; i++)
            for (j = k; j < n; j++)
            {
                const padic_radix_struct * x = ENTRY(A, i, j);
                if (x->u.size != 0 && x->v < bv)
                {
                    bv = x->v;
                    r = i;
                    c = j;
                }
            }

        if (r < 0)
        {
            /* the remaining block consists of zeros */
            padic_radix_t z;
            slong zN = _zero_block_prec(A, k, n, ctx);

            if (zN == WORD_MIN)
            {
                status = GR_UNABLE;
                break;
            }

            padic_radix_init(z, ctx);
            z->N = zN;
            status |= _padic_radix_finalize(z, ctx);
            status |= padic_radix_mul(res, res, z, ctx);
            padic_radix_clear(z, ctx);
            break;
        }

        if (r != k)
        {
            _gr_vec_swap(ENTRY(A, r, 0), ENTRY(A, k, 0), n, ctx);
            sign ^= 1;
        }

        if (c != k)
        {
            for (i = 0; i < n; i++)
                padic_radix_swap(ENTRY(A, i, c), ENTRY(A, i, k), ctx);
            sign ^= 1;
        }

        status |= padic_radix_mul(res, res, ENTRY(A, k, k), ctx);

        if (k == n - 1)
            break;

        status |= padic_radix_inv(inv, ENTRY(A, k, k), ctx);

        for (i = k + 1; i < n; i++)
        {
            padic_radix_struct * x = ENTRY(A, i, k);

            if (x->u.size == 0 && x->N == PADIC_RADIX_EXACT)
                continue;

            status |= padic_radix_mul(f, x, inv, ctx);
            status |= _padic_radix_vec_submul_scalar(ENTRY(A, i, k + 1),
                ENTRY(A, k, k + 1), n - k - 1, f, ctx);
        }
    }

    if (status & GR_UNABLE)
    {
        radix_integer_zero(&res->u, PADIC_RADIX_CTX_RADIX(ctx));
        res->v = 0;
        res->N = -PADIC_RADIX_ERR_MAX;
    }
    else if (sign)
        status |= padic_radix_neg(res, res, ctx);

    padic_radix_clear(inv, ctx);
    padic_radix_clear(f, ctx);
    gr_mat_clear(A, ctx);

    return status;
}

/*
    Blocked elimination with full pivoting by valuation levels. At level a
    (the minimal valuation of the current Schur complement S), the matrix
    T = (S / p^a) mod p has some rank r; a rank-revealing LU decomposition of
    T over F_p selects r rows and columns such that S11 = p^a U with U
    invertible modulo p. These r pivots are all of valuation a, as in full
    pivoting, and they are eliminated at once: det(S) = +- det(S11)
    det(S22 - S21 S11^(-1) S12), where det(S11) comes from an LU
    decomposition of S11 (whose pivots are of valuation a) and the Schur
    complement, whose entries have valuation at least a + 1, from a
    triangular solve and a matrix product. As with full pivoting, the
    Schur complements keep the absolute precision of the input, so the
    determinant has the optimal precision; the work is that of a blocked LU
    decomposition (a single level when A is invertible modulo p).
*/
int
padic_radix_mat_det_levels(padic_radix_t res, const gr_mat_t A, gr_ctx_t ctx)
{
    const radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    ulong p = GR_PADIC_RADIX_CTX(ctx)->p;
    slong n = A->r, m, i, j, k, r, a;
    gr_mat_t S;
    int sign = 0, status = GR_SUCCESS;
    nmod_t mod;

    if (n != A->c)
        return GR_DOMAIN;

    if (n == 0)
        return padic_radix_one(res, ctx);

    nmod_init(&mod, p);
    gr_mat_init(S, n, n, ctx);
    status |= gr_mat_set(S, A, ctx);
    status |= padic_radix_one(res, ctx);

    while ((m = S->r) > 0 && status == GR_SUCCESS)
    {
        nmod_mat_t T;
        slong * P, * rows, * cols, * pc, * P2;
        char * used;
        gr_mat_t S11, S12, S21, S22, LU;
        slong rk;

        a = WORD_MAX;
        for (i = 0; i < m; i++)
            for (j = 0; j < m; j++)
            {
                const padic_radix_struct * x = CENTRY(S, i, j);
                if (x->u.size != 0)
                    a = FLINT_MIN(a, x->v);
            }

        if (a == WORD_MAX)
        {
            /* a block of zeros */
            padic_radix_t z;
            slong zN = _zero_block_prec(S, 0, m, ctx);

            if (zN == WORD_MIN)
            {
                status = GR_UNABLE;
                break;
            }

            padic_radix_init(z, ctx);
            z->N = zN;
            status |= _padic_radix_finalize(z, ctx);
            status |= padic_radix_mul(res, res, z, ctx);
            padic_radix_clear(z, ctx);
            break;
        }

        /* rank profile of S / p^a modulo p */
        nmod_mat_init(T, m, m, p);
        for (i = 0; i < m; i++)
            for (j = 0; j < m; j++)
            {
                const padic_radix_struct * x = CENTRY(S, i, j);
                ulong t = 0;
                if (x->u.size != 0 && x->v == a)
                {
                    t = nmod_set_ui(x->u.d[0], mod);
                    if (x->u.size < 0)
                        t = nmod_neg(t, mod);
                }
                nmod_mat_entry(T, i, j) = t;
            }

        P = flint_malloc(sizeof(slong) * (4 * m));
        rows = P + m;
        cols = rows + m;
        pc = cols + m;
        used = flint_calloc(m, 1);

        r = nmod_mat_lu(P, T, 0);

        /* pivot columns: the first nonzero entry of each row of U */
        for (k = 0, j = 0; k < r; k++)
        {
            while (j < m && nmod_mat_entry(T, k, j) == 0)
                j++;
            pc[k] = j++;
        }
        nmod_mat_clear(T);

        for (k = 0; k < m; k++)
            rows[k] = P[k];
        for (k = 0; k < r; k++)
        {
            cols[k] = pc[k];
            used[pc[k]] = 1;
        }
        for (j = 0, k = r; j < m; j++)
            if (!used[j])
                cols[k++] = j;

        sign ^= _perm_parity(rows, m) ^ _perm_parity(cols, m);

        gr_mat_init(S11, r, r, ctx);
        gr_mat_init(S12, r, m - r, ctx);
        gr_mat_init(S21, m - r, r, ctx);
        gr_mat_init(S22, m - r, m - r, ctx);

        for (i = 0; i < m; i++)
            for (j = 0; j < m; j++)
            {
                padic_radix_struct * x = ENTRY(S, rows[i], cols[j]);
                padic_radix_struct * y;
                if (i < r)
                    y = (j < r) ? ENTRY(S11, i, j) : ENTRY(S12, i, j - r);
                else
                    y = (j < r) ? ENTRY(S21, i - r, j) : ENTRY(S22, i - r, j - r);
                padic_radix_swap(x, y, ctx);
            }

        /* det(S11) by LU (the pivots have valuation a) */
        P2 = flint_malloc(sizeof(slong) * FLINT_MAX(r, 1));
        gr_mat_init(LU, r, r, ctx);
        status |= gr_mat_lu(&rk, P2, LU, S11, 1, ctx);
        if (status == GR_SUCCESS && rk == r)
        {
            for (k = 0; k < r; k++)
                status |= padic_radix_mul(res, res, ENTRY(LU, k, k), ctx);
            sign ^= _perm_parity(P2, r);

            /* Schur complement S22 - S21 S11^(-1) S12 */
            if (r < m)
            {
                gr_mat_t W, Q;
                gr_mat_init(W, r, m - r, ctx);
                gr_mat_init(Q, m - r, m - r, ctx);
                status |= gr_mat_nonsingular_solve_lu_precomp(W, P2, LU, S12, ctx);
                status |= gr_mat_mul(Q, S21, W, ctx);
                status |= gr_mat_sub(S22, S22, Q, ctx);
                gr_mat_clear(W, ctx);
                gr_mat_clear(Q, ctx);
            }
        }
        else
            status |= GR_UNABLE;

        gr_mat_clear(LU, ctx);
        flint_free(P2);
        flint_free(P);
        flint_free(used);
        gr_mat_clear(S11, ctx);
        gr_mat_clear(S12, ctx);
        gr_mat_clear(S21, ctx);
        gr_mat_swap(S, S22, ctx);
        gr_mat_clear(S22, ctx);
    }

    gr_mat_clear(S, ctx);

    if (status & GR_UNABLE)
    {
        radix_integer_zero(&res->u, radix);
        res->v = 0;
        res->N = -PADIC_RADIX_ERR_MAX;
        return status;
    }

    if (sign)
        status |= padic_radix_neg(res, res, ctx);

    return status;
}

/*
    Both algorithms give the optimal precision; full pivoting with rank-one
    updates has less overhead for small matrices, the blocked version by
    valuation levels is faster above a size depending on the precision.
*/
int
padic_radix_mat_det(padic_radix_t res, const gr_mat_t A, gr_ctx_t ctx)
{
    slong n = A->r, cut, L;
    int status;

    if (n != A->c)
        return GR_DOMAIN;

    /* ad - bc has optimal precision */
    if (n <= 2)
        return gr_mat_det_cofactor(res, A, ctx);

    if (PADIC_RADIX_CTX_PREC_REL(ctx) != PADIC_RADIX_PREC_INF)
    {
        slong e = PADIC_RADIX_CTX_RADIX(ctx)->exp;
        L = (FLINT_MIN(PADIC_RADIX_CTX_PREC_REL(ctx), WORD_MAX / 4) + e - 1) / e;
    }
    else
    {
        slong i;
        L = 1;
        for (i = 0; i < n; i++)
            L = FLINT_MAX(L, FLINT_ABS(CENTRY(A, i, i)->u.size));
    }

    cut = (L <= 1) ? 40 : (L <= 2) ? 24 : 16;

    if (n > cut)
    {
        status = padic_radix_mat_det_levels(res, A, ctx);
        if (status == GR_SUCCESS)
            return status;
    }

    return padic_radix_mat_det_lu_full(res, A, ctx);
}
