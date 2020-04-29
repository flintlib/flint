/*
    Copyright (C) 2015 Tommy Hofmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_sparse_mat.h"
#include "nmod_sparse_mat.h"

slong
fmpz_sparse_mat_hnf_modular_eldiv(fmpz_sparse_mat_t M, const fmpz_t n)
{
    slong i;

    if (fmpz_sparse_mat_is_zero(M)) return 0;

    if (fmpz_abs_fits_ui(n))
    {
        /* TODO: have nmod version */
/*         nmod_init(&mod, fmpz_get_ui(n));
        nmod_sparse_mat_init(Mmod, M->r, M->c, mod);
        fmpz_sparse_mat_get_nmod_sparse_mat(Mmod, M);
        rank = nmod_sparse_mat_strong_echelon_form(Mmod);
        fmpz_sparse_mat_set_nmod_sparse_mat_unsigned(M, Mmod);
        nmod_sparse_mat_clear(Mmod);
 */    }

    fmpz_sparse_mat_strong_echelon_form_mod(M, n);
    
    for (i = 0; i < M->r; ++i)
        if (fmpz_sparse_vec_is_zero(&M->rows[i]))
            fmpz_sparse_vec_set_entry(&M->rows[i], i, n);
    return M->r;
}
    
