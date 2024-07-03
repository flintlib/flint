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

void
nmod_mat_solve_tril_classical(nmod_mat_t X, const nmod_mat_t L, const nmod_mat_t B, int unit)
{
    slong i, j, n, m;
    nmod_t mod;
    nn_ptr inv, tmp;

    n = L->r;
    m = B->c;
    mod = L->mod;

    if (!unit)
    {
        inv = _nmod_vec_init(n);
        for (i = 0; i < n; i++)
            inv[i] = n_invmod(nmod_mat_entry(L, i, i), mod.n);
    }
    else
        inv = NULL;

    const dot_params_t params = _nmod_vec_dot_params(n, mod);
    tmp = _nmod_vec_init(n);

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
            tmp[j] = nmod_mat_entry(X, j, i);

        for (j = 0; j < n; j++)
        {
            ulong s;
            s = _nmod_vec_dot(L->rows[j], tmp, j, mod, params);
            s = nmod_sub(nmod_mat_entry(B, j, i), s, mod);
            if (!unit)
                s = n_mulmod2_preinv(s, inv[j], mod.n, mod.ninv);
            tmp[j] = s;
        }

        for (j = 0; j < n; j++)
            nmod_mat_entry(X, j, i) = tmp[j];
    }

    _nmod_vec_clear(tmp);
    if (!unit)
        _nmod_vec_clear(inv);
}

void
nmod_mat_solve_tril_recursive(nmod_mat_t X, const nmod_mat_t L, const nmod_mat_t B, int unit)
{
    nmod_mat_t LA, LC, LD, XX, XY, BX, BY;
    slong r, n, m;

    n = L->r;
    m = B->c;
    r = n / 2;

    if (n == 0 || m == 0)
        return;

    /*
    Denoting inv(M) by M^, we have:

    [A 0]^ [X]  ==  [A^          0 ] [X]  ==  [A^ X]
    [C D]  [Y]  ==  [-D^ C A^    D^] [Y]  ==  [D^ (Y - C A^ X)]
    */

    nmod_mat_window_init(LA, L, 0, 0, r, r);
    nmod_mat_window_init(LC, L, r, 0, n, r);
    nmod_mat_window_init(LD, L, r, r, n, n);
    nmod_mat_window_init(BX, B, 0, 0, r, m);
    nmod_mat_window_init(BY, B, r, 0, n, m);
    nmod_mat_window_init(XX, X, 0, 0, r, m);
    nmod_mat_window_init(XY, X, r, 0, n, m);

    nmod_mat_solve_tril(XX, LA, BX, unit);
    nmod_mat_submul(XY, BY, LC, XX);
    nmod_mat_solve_tril(XY, LD, XY, unit);

    nmod_mat_window_clear(LA);
    nmod_mat_window_clear(LC);
    nmod_mat_window_clear(LD);
    nmod_mat_window_clear(BX);
    nmod_mat_window_clear(BY);
    nmod_mat_window_clear(XX);
    nmod_mat_window_clear(XY);
}

void
nmod_mat_solve_tril(nmod_mat_t X, const nmod_mat_t L, const nmod_mat_t B, int unit)
{
    if (B->r < NMOD_MAT_SOLVE_TRI_ROWS_CUTOFF ||
        B->c < NMOD_MAT_SOLVE_TRI_COLS_CUTOFF)
    {
        nmod_mat_solve_tril_classical(X, L, B, unit);
    }
    else
    {
        nmod_mat_solve_tril_recursive(X, L, B, unit);
    }
}
