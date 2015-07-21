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
unity_zp_mul(unity_zp f, const unity_zp g, const unity_zp h)
{
    ulong i, p;
    const slong len1 = g->poly->length;
    const slong len2 = h->poly->length;

    if (len1 == 0 || len2 == 0)
    {
        fmpz_mod_poly_zero(f->poly);
        return;
    }

    fmpz_mod_poly_fit_length(f->poly, len1 + len2 - 1);
    _fmpz_poly_mul(f->poly->coeffs, g->poly->coeffs, len1, h->poly->coeffs, len2);
    _fmpz_mod_poly_set_length(f->poly, len1 + len2 - 1);
    _fmpz_mod_poly_normalise(f->poly);

    p = n_pow(f->p, f->exp);
    for (i = f->poly->length - 1; i >= p; i--)
    {
        fmpz_add(f->poly->coeffs + i - p,
                f->poly->coeffs + i - p, f->poly->coeffs + i);

        fmpz_set_ui(f->poly->coeffs + i, 0);
    }

    _unity_zp_reduce_cyclotomic(f);
}


