/*
    Copyright (C) 2010, 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_mat.h"

/* the substitution, out of line: inlined in the caller, its short dot
   products were compiled to slower code (seen on x86-64) */
FLINT_STATIC_NOINLINE void
_solve_triu_classical(nmod_mat_t X, const nmod_mat_t U, const nmod_mat_t B,
        int unit, nn_srcptr inv, nn_ptr tmp)
{
    slong i, j, n = U->r, m = B->c;
    nmod_t mod = U->mod;

    const dot_params_t params = _nmod_vec_dot_params(n, mod);

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
            tmp[j] = nmod_mat_entry(X, j, i);

        for (j = n - 1; j >= 0; j--)
        {
            ulong s;
            s = _nmod_vec_dot(nmod_mat_entry_ptr(U, j, j + 1),
                              tmp + j + 1, n - j - 1, mod, params);
            s = nmod_sub(nmod_mat_entry(B, j, i), s, mod);
            if (!unit)
                s = nmod_mul(s, inv[j], mod);
            tmp[j] = s;
        }

        for (j = 0; j < n; j++)
            nmod_mat_entry(X, j, i) = tmp[j];
    }
}

void
nmod_mat_solve_triu_classical(nmod_mat_t X, const nmod_mat_t U, const nmod_mat_t B, int unit)
{
    slong i, n;
    nmod_t mod;
    nn_ptr inv, tmp;
    TMP_INIT;

    n = U->r;
    mod = U->mod;

    TMP_START;

    /* inverses of the diagonal entries: one inversion */
    if (!unit)
    {
        inv = TMP_ALLOC(n * sizeof(ulong));
        for (i = 0; i < n; i++)
            inv[i] = nmod_mat_entry(U, i, i);
        _nmod_vec_invert(inv, inv, n, mod);
    }
    else
        inv = NULL;

    tmp = TMP_ALLOC(n * sizeof(ulong));
    _solve_triu_classical(X, U, B, unit, inv, tmp);

    TMP_END;
}

void
nmod_mat_solve_triu_recursive(nmod_mat_t X, const nmod_mat_t U, const nmod_mat_t B, int unit)
{
    nmod_mat_t UA, UB, UD, XX, XY, BX, BY;
    slong r, n, m;

    n = U->r;
    m = B->c;
    r = n / 2;

    if (n == 0 || m == 0)
        return;

    /*
    Denoting inv(M) by M^, we have:

    [A B]^ [X]  ==  [A^ (X - B D^ Y)]
    [0 D]  [Y]  ==  [    D^ Y       ]
    */

    nmod_mat_window_init(UA, U, 0, 0, r, r);
    nmod_mat_window_init(UB, U, 0, r, r, n);
    nmod_mat_window_init(UD, U, r, r, n, n);
    nmod_mat_window_init(BX, B, 0, 0, r, m);
    nmod_mat_window_init(BY, B, r, 0, n, m);
    nmod_mat_window_init(XX, X, 0, 0, r, m);
    nmod_mat_window_init(XY, X, r, 0, n, m);

    nmod_mat_solve_triu(XY, UD, BY, unit);
    nmod_mat_submul(XX, BX, UB, XY);
    nmod_mat_solve_triu(XX, UA, XX, unit);

    nmod_mat_window_clear(UA);
    nmod_mat_window_clear(UB);
    nmod_mat_window_clear(UD);
    nmod_mat_window_clear(BX);
    nmod_mat_window_clear(BY);
    nmod_mat_window_clear(XX);
    nmod_mat_window_clear(XY);
}

void
nmod_mat_solve_triu(nmod_mat_t X, const nmod_mat_t U, const nmod_mat_t B, int unit)
{
    if (B->r < NMOD_MAT_SOLVE_TRI_ROWS_CUTOFF ||
        B->c < NMOD_MAT_SOLVE_TRI_COLS_CUTOFF)
    {
        nmod_mat_solve_triu_classical(X, U, B, unit);
    }
    else
    {
        nmod_mat_solve_triu_recursive(X, U, B, unit);
    }
}
