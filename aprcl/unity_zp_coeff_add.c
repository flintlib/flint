/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2015 Vladimir Glazachev
   
******************************************************************************/

#include "aprcl.h"

void
unity_zp_coeff_add_fmpz(unity_zp f, ulong ind, const fmpz_t x)
{
    fmpz_t coeff;
    fmpz_init(coeff);
    fmpz_mod_poly_get_coeff_fmpz(coeff, f->poly, ind);
    if (fmpz_is_zero(coeff))
    {
        fmpz_clear(coeff);
        fmpz_mod_poly_set_coeff_fmpz(f->poly, ind, x);
        return;
    }
    fmpz_clear(coeff);
    fmpz_add(f->poly->coeffs + ind, f->poly->coeffs + ind, x);           
    if (fmpz_cmp(f->poly->coeffs + ind, f->n) >= 0)
        fmpz_sub(f->poly->coeffs + ind, f->poly->coeffs + ind, f->n);
}

void
unity_zp_coeff_add_ui(unity_zp f, ulong ind, ulong x)
{
    fmpz_t coeff;
    fmpz_init(coeff);
    fmpz_mod_poly_get_coeff_fmpz(coeff, f->poly, ind);
    if (fmpz_is_zero(coeff))
    {
        fmpz_clear(coeff);
        fmpz_mod_poly_set_coeff_ui(f->poly, ind, x);
        return;
    }
    fmpz_clear(coeff);
    fmpz_add_ui(f->poly->coeffs + ind, f->poly->coeffs + ind, x);           
    if (fmpz_cmp(f->poly->coeffs + ind, f->n) >= 0)
        fmpz_sub(f->poly->coeffs + ind, f->poly->coeffs + ind, f->n);
}

