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
unity_zp_mul11(unity_zp f, const unity_zp g, const unity_zp h, fmpz_t * t)
{
    int i;

    for (i = 0; i < 10; i++)
        fmpz_mod_poly_get_coeff_fmpz(t[40 + i], g->poly, i);
    for (i = 0; i < 10; i++)
        fmpz_mod_poly_get_coeff_fmpz(t[50 + i], h->poly, i);
    for (i = 0; i < 5; i++)
    {
        fmpz_add(t[i], t[40 + i], t[45 + i]);
        fmpz_add(t[5 + i], t[50 + i], t[55 + i]);
    }

    unity_zp_ar3(t);
    for (i = 0; i < 9; i++)
        fmpz_set(t[40 + i], t[10 + i]);

    for (i = 0; i < 5; i++)
    {
        fmpz_mod_poly_get_coeff_fmpz(t[i], g->poly, i);
        fmpz_mod_poly_get_coeff_fmpz(t[5 + i], h->poly, i);
    }

    unity_zp_ar3(t);
    for (i = 0; i < 9; i++)
        fmpz_set(t[50 + i], t[10 + i]);

    for (i = 0; i < 5; i++)
    {
        fmpz_mod_poly_get_coeff_fmpz(t[i], g->poly, 5 + i);
        fmpz_mod_poly_get_coeff_fmpz(t[5 + i], h->poly, 5 + i);
    }

    unity_zp_ar3(t);
    for (i = 0; i < 9; i++)
    {
        fmpz_sub(t[40 + i], t[40 + i], t[10 + i]);
        fmpz_sub(t[40 + i], t[40 + i], t[50 + i]);
    }

    fmpz_add(t[1], t[10], t[45]);
    for (i = 0; i < 8; i++)
        fmpz_add(t[50 + i], t[50 + i], t[11 + i]);

    for (i = 0; i < 3; i++)
        fmpz_add(t[50 + i], t[50 + i], t[46 + i]);

    for (i = 5; i < 9; i++)
        fmpz_add(t[50 + i], t[50 + i], t[35 + i]);

    for (i = 0; i < 9; i++)
    {
        fmpz_sub(t[0], t[50 + i], t[1]);
        unity_zp_coeff_set_fmpz(f, i, t[0]);
    }
    fmpz_sub(t[0], t[44], t[1]);
    unity_zp_coeff_set_fmpz(f, 9, t[0]);
}

