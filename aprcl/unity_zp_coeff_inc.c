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
unity_zp_coeff_inc(unity_zp f, ulong ind)
{
    if (ind >= f->poly->length)
    {
        fmpz_mod_poly_set_coeff_ui(f->poly, ind, 1);
        return;
    }

    fmpz_add_ui(f->poly->coeffs + ind, f->poly->coeffs + ind, 1);
    if (fmpz_equal(f->poly->coeffs + ind, f->n))
        fmpz_set_ui(f->poly->coeffs + ind, 0);
}

void
unity_zp_coeff_dec(unity_zp f, ulong ind)
{
    if (ind >= f->poly->length)
    {
        fmpz_mod_poly_set_coeff_si(f->poly, ind, -1);
        return;
    }

    fmpz_sub_ui(f->poly->coeffs + ind, f->poly->coeffs + ind, 1);
    if (fmpz_equal_si(f->poly->coeffs + ind, -1))
        fmpz_add(f->poly->coeffs + ind, f->poly->coeffs + ind, f->n);
}

