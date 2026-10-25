/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mat.h"
#include "gr.h"
#include "gr_mat.h"
#include "padic_radix.h"

#define ENTRY(A, i, j) ((padic_radix_struct *) gr_mat_entry_ptr(A, i, j, ctx))
#define CENTRY(A, i, j) ((const padic_radix_struct *) gr_mat_entry_srcptr(A, i, j, ctx))

/*
    Solve A X = B by block Gauss-Jordan elimination with full pivoting by
    valuation levels on G = [A | B]. At each level, with a the minimal
    valuation over the remaining rows and columns of A, a rank-revealing LU
    decomposition of (G_rem / p^a) mod p selects r pivot rows and columns
    with S11 = p^a U, U invertible modulo p (all pivots of valuation a, as in
    full pivoting). The pivot rows are replaced by S11^(-1) times themselves
    (LU decomposition of S11 and triangular solves), and the pivot columns
    are eliminated from all the other rows, including the pivot rows of the
    previous levels, by one matrix product. This is the order of operations
    of Gauss-Jordan elimination with full pivoting, which gives the solution
    with the optimal (differential) precision; the work is that of a blocked
    elimination (a single level, i.e. one LU solve, when A is invertible
    modulo p).
*/
static int
_solve_gauss_jordan_levels(gr_mat_t X, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    ulong p = GR_PADIC_RADIX_CTX(ctx)->p;
    slong n = A->r, c = B->c, nc = n + c, i, j, k, r, a, rk, nrem;
    gr_mat_t G;
    slong * rrem, * crem, * prow_of_col;
    nmod_t mod;
    int status = GR_SUCCESS;

    nmod_init(&mod, p);

    /* when A / p^a is invertible modulo p (a the minimal valuation), all
       pivots of the LU decomposition with minimal-valuation column pivots
       have valuation a, which is as good as full pivoting */
    {
        nmod_mat_t T;
        slong aa = WORD_MAX;

        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
            {
                const padic_radix_struct * x = CENTRY(A, i, j);
                if (x->u.size != 0)
                    aa = FLINT_MIN(aa, x->v);
            }

        if (aa != WORD_MAX)
        {
            nmod_mat_init(T, n, n, p);
            for (i = 0; i < n; i++)
                for (j = 0; j < n; j++)
                {
                    const padic_radix_struct * x = CENTRY(A, i, j);
                    ulong t = 0;
                    if (x->u.size != 0 && x->v == aa)
                    {
                        t = nmod_set_ui(x->u.d[0], mod);
                        if (x->u.size < 0)
                            t = nmod_neg(t, mod);
                    }
                    nmod_mat_entry(T, i, j) = t;
                }
            r = nmod_mat_rank(T);
            nmod_mat_clear(T);

            if (r == n)
                return gr_mat_nonsingular_solve_lu(X, A, B, ctx);
        }
    }

    gr_mat_init(G, n, nc, ctx);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
            status |= padic_radix_set(ENTRY(G, i, j), CENTRY(A, i, j), ctx);
        for (j = 0; j < c; j++)
            status |= padic_radix_set(ENTRY(G, i, n + j), CENTRY(B, i, j), ctx);
    }

    rrem = flint_malloc(sizeof(slong) * 3 * n);
    crem = rrem + n;
    prow_of_col = crem + n;
    for (i = 0; i < n; i++)
    {
        rrem[i] = i;
        crem[i] = i;
        prow_of_col[i] = -1;
    }
    nrem = n;

    while (nrem > 0 && status == GR_SUCCESS)
    {
        nmod_mat_t T;
        slong * P, * pr, * pcl, * orows, * rest, * P2, no, nrest;
        char * usedr, * usedc;
        gr_mat_t S11, LU, R, W, GO, Q;
        int unknown = 0;

        a = WORD_MAX;
        for (i = 0; i < nrem; i++)
            for (j = 0; j < nrem; j++)
            {
                const padic_radix_struct * x = CENTRY(G, rrem[i], crem[j]);
                if (x->u.size != 0)
                    a = FLINT_MIN(a, x->v);
                else if (x->N != PADIC_RADIX_EXACT)
                    unknown = 1;
            }

        if (a == WORD_MAX)
        {
            status = unknown ? GR_UNABLE : GR_DOMAIN;
            break;
        }

        /* rank profile at level a */
        nmod_mat_init(T, nrem, nrem, p);
        for (i = 0; i < nrem; i++)
            for (j = 0; j < nrem; j++)
            {
                const padic_radix_struct * x = CENTRY(G, rrem[i], crem[j]);
                ulong t = 0;
                if (x->u.size != 0 && x->v == a)
                {
                    t = nmod_set_ui(x->u.d[0], mod);
                    if (x->u.size < 0)
                        t = nmod_neg(t, mod);
                }
                nmod_mat_entry(T, i, j) = t;
            }

        P = flint_malloc(sizeof(slong) * (5 * n + nc));
        pr = P + n;
        pcl = pr + n;
        orows = pcl + n;
        P2 = orows + n;
        rest = P2 + n;
        usedr = flint_calloc(2 * n, 1);
        usedc = usedr + n;

        r = nmod_mat_lu(P, T, 0);

        for (k = 0, j = 0; k < r; k++)
        {
            while (j < nrem && nmod_mat_entry(T, k, j) == 0)
                j++;
            pr[k] = rrem[P[k]];
            pcl[k] = crem[j];
            usedr[pr[k]] = 1;
            usedc[pcl[k]] = 1;
            j++;
        }
        nmod_mat_clear(T);

        /* the other rows, and the columns still to be updated */
        for (i = 0, no = 0; i < n; i++)
            if (!usedr[i])
                orows[no++] = i;
        for (j = 0, nrest = 0; j < nrem; j++)
            if (!usedc[crem[j]])
                rest[nrest++] = crem[j];
        for (j = 0; j < c; j++)
            rest[nrest++] = n + j;

        /* W = S11^(-1) G[pr, rest] */
        gr_mat_init(S11, r, r, ctx);
        gr_mat_init(LU, r, r, ctx);
        gr_mat_init(R, r, nrest, ctx);
        gr_mat_init(W, r, nrest, ctx);
        for (i = 0; i < r; i++)
        {
            for (j = 0; j < r; j++)
                status |= padic_radix_set(ENTRY(S11, i, j), CENTRY(G, pr[i], pcl[j]), ctx);
            for (j = 0; j < nrest; j++)
                status |= padic_radix_set(ENTRY(R, i, j), CENTRY(G, pr[i], rest[j]), ctx);
        }

        status |= gr_mat_lu(&rk, P2, LU, S11, 1, ctx);
        if (status != GR_SUCCESS || rk != r)
            status |= GR_UNABLE;
        else
            status |= gr_mat_nonsingular_solve_lu_precomp(W, P2, LU, R, ctx);

        /* G[o, rest] -= G[o, pcl] W for all other rows o */
        if (status == GR_SUCCESS && no > 0 && nrest > 0)
        {
            gr_mat_init(GO, no, r, ctx);
            gr_mat_init(Q, no, nrest, ctx);
            for (i = 0; i < no; i++)
                for (j = 0; j < r; j++)
                    status |= padic_radix_set(ENTRY(GO, i, j), CENTRY(G, orows[i], pcl[j]), ctx);
            status |= gr_mat_mul(Q, GO, W, ctx);
            for (i = 0; i < no; i++)
                for (j = 0; j < nrest; j++)
                    status |= padic_radix_sub(ENTRY(G, orows[i], rest[j]),
                        CENTRY(G, orows[i], rest[j]), CENTRY(Q, i, j), ctx);
            gr_mat_clear(GO, ctx);
            gr_mat_clear(Q, ctx);
        }

        /* the pivot rows become W (identity in the pivot columns) */
        for (i = 0; i < r && status == GR_SUCCESS; i++)
        {
            for (j = 0; j < nrest; j++)
                padic_radix_swap(ENTRY(G, pr[i], rest[j]), ENTRY(W, i, j), ctx);
            prow_of_col[pcl[i]] = pr[i];
        }

        /* remaining rows and columns */
        {
            slong nr2 = 0, nc2 = 0;
            for (i = 0; i < nrem; i++)
            {
                if (!usedr[rrem[i]])
                    rrem[nr2++] = rrem[i];
                if (!usedc[crem[i]])
                    crem[nc2++] = crem[i];
            }
            nrem = nr2;
        }

        gr_mat_clear(S11, ctx);
        gr_mat_clear(LU, ctx);
        gr_mat_clear(R, ctx);
        gr_mat_clear(W, ctx);
        flint_free(P);
        flint_free(usedr);
    }

    if (status == GR_SUCCESS)
        for (j = 0; j < n; j++)
            for (k = 0; k < c; k++)
                padic_radix_swap(ENTRY(X, j, k), ENTRY(G, prow_of_col[j], n + k), ctx);

    flint_free(rrem);
    gr_mat_clear(G, ctx);
    return status;
}

int
padic_radix_mat_nonsingular_solve(gr_mat_t X, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    slong n = A->r;
    gr_mat_t Y;
    int status;

    if (n != A->c || B->r != n || X->r != n || X->c != B->c)
        return GR_DOMAIN;

    if (n == 0)
        return GR_SUCCESS;

    gr_mat_init(Y, n, B->c, ctx);
    if (PADIC_RADIX_CTX_FLAGS(ctx) & PADIC_RADIX_FAST_LINALG)
        status = gr_mat_nonsingular_solve_lu(Y, A, B, ctx);    /* speed over precision */
    else
        status = _solve_gauss_jordan_levels(Y, A, B, ctx);
    if (status == GR_SUCCESS)
        gr_mat_swap(X, Y, ctx);     /* (X may alias A or B) */
    gr_mat_clear(Y, ctx);
    return status;
}

int
padic_radix_mat_inv(gr_mat_t X, const gr_mat_t A, gr_ctx_t ctx)
{
    slong n = A->r;
    gr_mat_t I;
    int status;

    if (n != A->c || X->r != n || X->c != n)
        return GR_DOMAIN;

    gr_mat_init(I, n, n, ctx);
    status = gr_mat_one(I, ctx);
    status |= padic_radix_mat_nonsingular_solve(X, A, I, ctx);
    gr_mat_clear(I, ctx);
    return status;
}
