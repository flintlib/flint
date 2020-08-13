/*
    Copyright (C) 2015 Tommy Hofmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

void
fmpz_mat_hnf_modular_eldiv(fmpz_mat_t A, const fmpz_t D)
{
    slong i;
    mp_limb_t Dlimbt;
    nmod_mat_t AmodD;

    if (fmpz_mat_is_empty(A))
        return;

    if (fmpz_abs_fits_ui(D))
    {
        Dlimbt = fmpz_get_ui(D);
        nmod_mat_init(AmodD, A->r, A->c, Dlimbt);
        fmpz_mat_get_nmod_mat(AmodD, A);
        nmod_mat_strong_echelon_form(AmodD);
        fmpz_mat_set_nmod_mat_unsigned(A, AmodD);
        nmod_mat_clear(AmodD);
    }
    else
    {
        fmpz_mat_strong_echelon_form_mod(A, D);
    }

    for (i = 0; i < A->c; i++)
    {
        if (fmpz_is_zero(fmpz_mat_entry(A, i, i)))
        {
            fmpz_set(fmpz_mat_entry(A, i, i), D);
        }
    }
}
    
