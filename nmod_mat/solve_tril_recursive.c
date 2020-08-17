/*
    Copyright (C) 2010,2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_mat.h"
#include "nmod_vec.h"


void
nmod_mat_solve_tril_recursive(nmod_mat_t X,
                                    const nmod_mat_t L, const nmod_mat_t B,
                                                                    int unit)
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
