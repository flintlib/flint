/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "nmod_mat/impl.h"

/* columns of the triangular matrix copied to rows at a time, for the dot products */
#define TRSM_RIGHT_BLOCK 64

/*
    With the delayed reductions of solve_tri_delayed.c, unless the vectors
    would be long and few (see impl.h). Otherwise, X U = B row by row:
    x_j = (b_j - sum_{t < j} x_t U[t][j]) / U[j][j]
    for j = 0, ..., r - 1. For the dot products, the columns of U above the
    diagonal are copied to rows, TRSM_RIGHT_BLOCK of them at a time, each
    block serving all the rows of X.
*/
void
nmod_mat_solve_triu_right_classical(nmod_mat_t X, const nmod_mat_t U,
                                    const nmod_mat_t B, int unit)
{
    slong i, j, t, j0, j1, r = U->r, k = B->r;
    nmod_t mod = U->mod;
    nn_ptr inv = NULL;
    nmod_mat_t T;
    TMP_INIT;

    if (r == 0 || k == 0)
        return;

    /* with delayed reductions (impl.h), unless the vectors would be long */
    if (NMOD_MAT_SOLVE_TRI_RIGHT_ROWS(r, k, mod.n))
    {
        _nmod_mat_solve_tri_rows_delayed(X, U, B, unit, 1, 1);
        return;
    }
    if (NMOD_MAT_SOLVE_TRI_RIGHT_AXPY(r, k, mod.n))
    {
        _nmod_mat_solve_tri_right_axpy_delayed(X, U, B, unit, 1);
        return;
    }

    TMP_START;

    /* inverses of the diagonal entries: one inversion */
    if (!unit)
    {
        inv = TMP_ALLOC(r * sizeof(ulong));
        for (j = 0; j < r; j++)
            inv[j] = nmod_mat_entry(U, j, j);
        _nmod_vec_invert(inv, inv, r, mod);
    }

    nmod_mat_init(T, FLINT_MIN(r, TRSM_RIGHT_BLOCK), r, mod.n);

    const dot_params_t params = _nmod_vec_dot_params(r, mod);

    for (j0 = 0; j0 < r; j0 = j1)
    {
        j1 = FLINT_MIN(r, j0 + TRSM_RIGHT_BLOCK);

        /* T[j - j0][t] = U[t][j] for j0 <= j < j1 and t < j */
        for (j = j0; j < j1; j++)
            for (t = 0; t < j; t++)
                nmod_mat_entry(T, j - j0, t) = nmod_mat_entry(U, t, j);

        for (i = 0; i < k; i++)
        {
            nn_ptr x = nmod_mat_entry_ptr(X, i, 0);
            nn_srcptr b = nmod_mat_entry_ptr(B, i, 0);

            for (j = j0; j < j1; j++)
            {
                ulong s = _nmod_vec_dot(nmod_mat_entry_ptr(T, j - j0, 0), x,
                                        j, mod, params);
                s = nmod_sub(b[j], s, mod);
                if (!unit)
                    s = nmod_mul(s, inv[j], mod);
                x[j] = s;
            }
        }
    }

    nmod_mat_clear(T);
    TMP_END;
}

void
nmod_mat_solve_triu_right_recursive(nmod_mat_t X, const nmod_mat_t U,
                                    const nmod_mat_t B, int unit)
{
    nmod_mat_t UA, UB, UD, XX, XY, BX, BY;
    slong r, n, k;

    n = U->r;
    k = B->r;
    r = n / 2;

    if (n == 0 || k == 0)
        return;

    /*
        With U = [A B; 0 D], X = [XX XY] and B = [BX BY]:
        XX = BX A^-1, then XY = (BY - XX B) D^-1.
    */
    nmod_mat_window_init(UA, U, 0, 0, r, r);
    nmod_mat_window_init(UB, U, 0, r, r, n);
    nmod_mat_window_init(UD, U, r, r, n, n);
    nmod_mat_window_init(BX, B, 0, 0, k, r);
    nmod_mat_window_init(BY, B, 0, r, k, n);
    nmod_mat_window_init(XX, X, 0, 0, k, r);
    nmod_mat_window_init(XY, X, 0, r, k, n);

    nmod_mat_solve_triu_right(XX, UA, BX, unit);
    nmod_mat_submul(XY, BY, XX, UB);
    nmod_mat_solve_triu_right(XY, UD, XY, unit);

    nmod_mat_window_clear(UA);
    nmod_mat_window_clear(UB);
    nmod_mat_window_clear(UD);
    nmod_mat_window_clear(BX);
    nmod_mat_window_clear(BY);
    nmod_mat_window_clear(XX);
    nmod_mat_window_clear(XY);
}

void
nmod_mat_solve_triu_right(nmod_mat_t X, const nmod_mat_t U, const nmod_mat_t B, int unit)
{
    /* see impl.h */
    if (NMOD_MAT_SOLVE_TRI_RIGHT_REC(U->r, B->r, U->mod.n))
        nmod_mat_solve_triu_right_recursive(X, U, B, unit);
    else
        nmod_mat_solve_triu_right_classical(X, U, B, unit);
}
