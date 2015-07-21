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

/*
    Computes f = g * g for p = 11. 
    g must be reduced by F_11 cyclotomic polynomial.
    t is the memory for fmpz_t;
    Resulting f reduced by F_11 cyclotomic polynomial.
*/
void
unity_zp_sqr11(unity_zp f, const unity_zp g, fmpz_t * t)
{
    int i;

    for (i = 0; i < 10; i++)
        fmpz_mod_poly_get_coeff_fmpz(t[30 + i], g->poly, i);

    fmpz_set(t[0], t[30]);
    fmpz_set(t[1], t[31]);
    fmpz_set(t[2], t[32]);
    fmpz_set(t[3], t[33]);
    fmpz_set(t[4], t[34]);

    unity_zp_ar4(t);
    fmpz_set(t[50], t[5]);
    fmpz_set(t[51], t[6]);
    fmpz_set(t[52], t[7]);
    fmpz_set(t[53], t[8]);
    fmpz_set(t[54], t[9]);
    fmpz_set(t[55], t[10]);
    fmpz_set(t[56], t[11]);
    fmpz_set(t[57], t[12]);
    fmpz_set(t[58], t[13]);

    fmpz_set(t[0], t[35]);
    fmpz_set(t[1], t[36]);
    fmpz_set(t[2], t[37]);
    fmpz_set(t[3], t[38]);
    fmpz_set(t[4], t[39]);

    unity_zp_ar4(t);
    fmpz_set(t[60], t[5]);
    fmpz_set(t[61], t[6]);
    fmpz_set(t[62], t[7]);
    fmpz_set(t[63], t[8]);
    fmpz_set(t[64], t[9]);
    fmpz_set(t[65], t[10]);
    fmpz_set(t[66], t[11]);
    fmpz_set(t[67], t[12]);
    fmpz_set(t[68], t[13]);

    fmpz_set(t[0], t[35]);
    fmpz_set(t[1], t[36]);
    fmpz_set(t[2], t[37]);
    fmpz_set(t[3], t[38]);
    fmpz_set(t[4], t[39]);

    fmpz_mul_2exp(t[5], t[30], 1);
    fmpz_mul_2exp(t[6], t[31], 1);
    fmpz_mul_2exp(t[7], t[32], 1);
    fmpz_mul_2exp(t[8], t[33], 1);
    fmpz_mul_2exp(t[9], t[34], 1);

    unity_zp_ar3(t);

    fmpz_add(t[1], t[60], t[15]);
    for (i = 0; i < 8; i++)
        fmpz_add(t[50 + i], t[50 + i], t[61 + i]);

    for (i = 0; i < 3; i++)
        fmpz_add(t[50 + i], t[50 + i], t[16 + i]);

    for (i = 5; i < 9; i++)
        fmpz_add(t[50 + i], t[50 + i], t[5 + i]);

    for (i = 0; i < 9; i++)
    {
        fmpz_sub(t[0], t[50 + i], t[1]);
        unity_zp_coeff_set_fmpz(f, i, t[0]);
    }
    fmpz_sub(t[0], t[14], t[1]);
    unity_zp_coeff_set_fmpz(f, 9, t[0]);
}

