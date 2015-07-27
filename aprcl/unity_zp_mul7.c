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
    Computes f = g * h for p = 7. 
    g and h must be reduced by F_7 cyclotomic polynomial.
    t is the memory for fmpz_t; size of t must be > 68.
    Resulting f reduced by F_7 cyclotomic polynomial.
*/
void
unity_zp_mul7(unity_zp f, const unity_zp g, const unity_zp h, fmpz_t * t)
{
    fmpz_mod_poly_get_coeff_fmpz(t[30], g->poly, 0);
    fmpz_mod_poly_get_coeff_fmpz(t[31], g->poly, 1);
    fmpz_mod_poly_get_coeff_fmpz(t[32], g->poly, 2);
    fmpz_mod_poly_get_coeff_fmpz(t[33], g->poly, 3);
    fmpz_mod_poly_get_coeff_fmpz(t[34], g->poly, 4);
    fmpz_mod_poly_get_coeff_fmpz(t[35], g->poly, 5);

    fmpz_mod_poly_get_coeff_fmpz(t[40], h->poly, 0);
    fmpz_mod_poly_get_coeff_fmpz(t[41], h->poly, 1);
    fmpz_mod_poly_get_coeff_fmpz(t[42], h->poly, 2);
    fmpz_mod_poly_get_coeff_fmpz(t[43], h->poly, 3);
    fmpz_mod_poly_get_coeff_fmpz(t[44], h->poly, 4);
    fmpz_mod_poly_get_coeff_fmpz(t[45], h->poly, 5);

    fmpz_set(t[0], t[30]);
    fmpz_set(t[1], t[31]);
    fmpz_set(t[2], t[32]);
    fmpz_set(t[3], t[40]);
    fmpz_set(t[4], t[41]);
    fmpz_set(t[5], t[42]);
    unity_zp_ar1(t);

    fmpz_set(t[50], t[6]);
    fmpz_set(t[51], t[7]);
    fmpz_set(t[52], t[8]);
    fmpz_set(t[53], t[9]);
    fmpz_set(t[54], t[10]);

    fmpz_set(t[0], t[33]);
    fmpz_set(t[1], t[34]);
    fmpz_set(t[2], t[35]);
    fmpz_set(t[3], t[43]);
    fmpz_set(t[4], t[44]);
    fmpz_set(t[5], t[45]);
    unity_zp_ar1(t);

    fmpz_set(t[56], t[6]);
    fmpz_set(t[57], t[7]);
    fmpz_set(t[58], t[8]);
    fmpz_set(t[59], t[9]);
    fmpz_set(t[60], t[10]);

    fmpz_sub(t[0], t[30], t[33]);
    fmpz_sub(t[1], t[31], t[34]);
    fmpz_sub(t[2], t[32], t[35]);
    fmpz_sub(t[3], t[43], t[40]);
    fmpz_sub(t[4], t[44], t[41]);
    fmpz_sub(t[5], t[45], t[42]);
    unity_zp_ar1(t);

    fmpz_set(t[61], t[6]);
    fmpz_set(t[62], t[7]);
    fmpz_set(t[63], t[8]);
    fmpz_set(t[64], t[9]);
    fmpz_set(t[65], t[10]);

    fmpz_add(t[68], t[56], t[64]);
    fmpz_add(t[66], t[68], t[53]);
    fmpz_add(t[68], t[57], t[65]);
    fmpz_add(t[67], t[68], t[54]);
    fmpz_add(t[68], t[53], t[61]);
    fmpz_add(t[53], t[68], t[50]);
    fmpz_add(t[68], t[54], t[62]);
    fmpz_add(t[54], t[68], t[51]);

    fmpz_add(t[55], t[52], t[63]);
    fmpz_add(t[63], t[53], t[56]);
    fmpz_add(t[64], t[54], t[57]);
    fmpz_add(t[65], t[55], t[58]);
    fmpz_add(t[56], t[66], t[59]);
    fmpz_add(t[57], t[67], t[60]);

    fmpz_add(t[68], t[50], t[57]);
    fmpz_sub(t[0], t[68], t[56]);
    unity_zp_coeff_set_fmpz(f, 0, t[0]);
    fmpz_add(t[68], t[51], t[58]);
    fmpz_sub(t[0], t[68], t[56]);
    unity_zp_coeff_set_fmpz(f, 1, t[0]);
    fmpz_add(t[68], t[52], t[59]);
    fmpz_sub(t[0], t[68], t[56]);
    unity_zp_coeff_set_fmpz(f, 2, t[0]);
    fmpz_add(t[68], t[63], t[60]);
    fmpz_sub(t[0], t[68], t[56]);
    unity_zp_coeff_set_fmpz(f, 3, t[0]);
    fmpz_sub(t[0], t[64], t[56]);
    unity_zp_coeff_set_fmpz(f, 4, t[0]);
    fmpz_sub(t[0], t[65], t[56]);
    unity_zp_coeff_set_fmpz(f, 5, t[0]);
}

