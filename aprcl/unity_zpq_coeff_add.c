/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "aprcl.h"

void
unity_zpq_coeff_add(unity_zpq f, ulong i, ulong j, const fmpz_t x)
{
    fmpz_t coeff;
    fmpz_init(coeff);
    fmpz_mod_poly_get_coeff_fmpz(coeff, f->polys[j], i);
    if (fmpz_is_zero(coeff))
    {
        fmpz_clear(coeff);
        fmpz_mod_poly_set_coeff_fmpz(f->polys[j], i, x);
        return;
    }
    fmpz_clear(coeff);
    fmpz_add(f->polys[j]->coeffs + i, f->polys[j]->coeffs + i, x);           
    if (fmpz_cmp(f->polys[j]->coeffs + i, f->n) >= 0)
        fmpz_sub(f->polys[j]->coeffs + i, f->polys[j]->coeffs + i, f->n);
}

void
unity_zpq_coeff_add_ui(unity_zpq f, ulong i, ulong j, ulong x)
{
    fmpz_t coeff;
    fmpz_init(coeff);
    fmpz_mod_poly_get_coeff_fmpz(coeff, f->polys[j], i);
    if (fmpz_is_zero(coeff))
    {
        fmpz_clear(coeff);
        fmpz_mod_poly_set_coeff_ui(f->polys[j], i, x);
        return;
    }
    fmpz_clear(coeff);
    fmpz_add_ui(f->polys[j]->coeffs + i, f->polys[j]->coeffs + i, x);          
    if (fmpz_cmp(f->polys[j]->coeffs + i, f->n) >= 0)
        fmpz_sub(f->polys[j]->coeffs + i, f->polys[j]->coeffs + i, f->n);
}

