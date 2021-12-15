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

void fmpz_mod_mat_solve_triu_recursive(fmpz_mod_mat_t X,
                      const fmpz_mod_mat_t U, const fmpz_mod_mat_t B, int unit)
{
    fmpz_mod_mat_t UA, UB, UD, XX, XY, BX, BY;
    slong r, n, m;

    n = U->mat->r;
    m = B->mat->c;
    r = n / 2;

    if (n == 0 || m == 0)
        return;

    /*
       Denoting inv(M) by M^, we have:

       [A B]^ [X]  ==  [A^ (X - B D^ Y)]
       [0 D]  [Y]  ==  [    D^ Y       ]
     */

    fmpz_mod_mat_window_init(UA, U, 0, 0, r, r);
    fmpz_mod_mat_window_init(UB, U, 0, r, r, n);
    fmpz_mod_mat_window_init(UD, U, r, r, n, n);
    fmpz_mod_mat_window_init(BX, B, 0, 0, r, m);
    fmpz_mod_mat_window_init(BY, B, r, 0, n, m);
    fmpz_mod_mat_window_init(XX, X, 0, 0, r, m);
    fmpz_mod_mat_window_init(XY, X, r, 0, n, m);

    fmpz_mod_mat_solve_triu(XY, UD, BY, unit);
    fmpz_mod_mat_submul(XX, BX, UB, XY);
    fmpz_mod_mat_solve_triu(XX, UA, XX, unit);

    fmpz_mod_mat_window_clear(UA);
    fmpz_mod_mat_window_clear(UB);
    fmpz_mod_mat_window_clear(UD);
    fmpz_mod_mat_window_clear(BX);
    fmpz_mod_mat_window_clear(BY);
    fmpz_mod_mat_window_clear(XX);
    fmpz_mod_mat_window_clear(XY);
}

