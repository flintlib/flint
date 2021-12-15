/*
    Copyright (C) 2010, 2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mat.h"

void fmpz_mod_mat_solve_tril_recursive(fmpz_mod_mat_t X,
                      const fmpz_mod_mat_t L, const fmpz_mod_mat_t B, int unit)
{
    fmpz_mod_mat_t LA, LC, LD, XX, XY, BX, BY;
    slong r, n, m;

    n = L->mat->r;
    m = B->mat->c;
    r = n / 2;

    if (n == 0 || m == 0)
        return;

    /*
       Denoting inv(M) by M^, we have:

       [A 0]^ [X]  ==  [A^          0 ] [X]  ==  [A^ X]
       [C D]  [Y]  ==  [-D^ C A^    D^] [Y]  ==  [D^ (Y - C A^ X)]
     */

    fmpz_mod_mat_window_init(LA, L, 0, 0, r, r);
    fmpz_mod_mat_window_init(LC, L, r, 0, n, r);
    fmpz_mod_mat_window_init(LD, L, r, r, n, n);
    fmpz_mod_mat_window_init(BX, B, 0, 0, r, m);
    fmpz_mod_mat_window_init(BY, B, r, 0, n, m);
    fmpz_mod_mat_window_init(XX, X, 0, 0, r, m);
    fmpz_mod_mat_window_init(XY, X, r, 0, n, m);

    fmpz_mod_mat_solve_tril(XX, LA, BX, unit);
    fmpz_mod_mat_submul(XY, BY, LC, XX);
    fmpz_mod_mat_solve_tril(XY, LD, XY, unit);

    fmpz_mod_mat_window_clear(LA);
    fmpz_mod_mat_window_clear(LC);
    fmpz_mod_mat_window_clear(LD);
    fmpz_mod_mat_window_clear(BX);
    fmpz_mod_mat_window_clear(BY);
    fmpz_mod_mat_window_clear(XX);
    fmpz_mod_mat_window_clear(XY);
}

