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
unity_zp_mul4(unity_zp f, const unity_zp g, const unity_zp h, fmpz_t * t)
{
    fmpz_mod_poly_get_coeff_fmpz(t[0], g->poly, 0);
    fmpz_mod_poly_get_coeff_fmpz(t[1], g->poly, 1);
    fmpz_mod_poly_get_coeff_fmpz(t[2], h->poly, 0);
    fmpz_mod_poly_get_coeff_fmpz(t[3], h->poly, 1);

    fmpz_add(t[4], t[0], t[1]);
    fmpz_add(t[5], t[2], t[3]);
    fmpz_sub(t[6], t[3], t[2]);
    fmpz_mul(t[7], t[4], t[2]);
    fmpz_mul(t[8], t[5], t[1]);
    fmpz_mul(t[9], t[6], t[0]);
    fmpz_sub(t[0], t[7], t[8]);

    unity_zp_coeff_set_fmpz(f, 0, t[0]);
    fmpz_add(t[0], t[7], t[9]);
    unity_zp_coeff_set_fmpz(f, 1, t[0]);
}

void
unity_zp_mul8(unity_zp f, const unity_zp g, const unity_zp h, fmpz_t * t)
{
    fmpz_mod_poly_get_coeff_fmpz(t[0], g->poly, 0);
    fmpz_mod_poly_get_coeff_fmpz(t[1], g->poly, 1);
    fmpz_mod_poly_get_coeff_fmpz(t[2], g->poly, 2);
    fmpz_mod_poly_get_coeff_fmpz(t[3], g->poly, 3);

    fmpz_mod_poly_get_coeff_fmpz(t[4], h->poly, 0);
    fmpz_mod_poly_get_coeff_fmpz(t[5], h->poly, 1);
    fmpz_mod_poly_get_coeff_fmpz(t[6], h->poly, 2);
    fmpz_mod_poly_get_coeff_fmpz(t[7], h->poly, 3);

    fmpz_add(t[8], t[1], t[3]);
    fmpz_add(t[9], t[5], t[7]);
    fmpz_add(t[10], t[2], t[3]);
    fmpz_add(t[11], t[6], t[7]);
    fmpz_add(t[12], t[0], t[1]);
    fmpz_add(t[13], t[4], t[5]);

    fmpz_add(t[14], t[0], t[2]);
    fmpz_add(t[15], t[4], t[6]);
    fmpz_mul(t[16], t[0], t[4]);
    fmpz_mul(t[17], t[1], t[5]);
    fmpz_mul(t[18], t[2], t[6]);
    fmpz_mul(t[19], t[3], t[7]);
    fmpz_mul(t[22], t[12], t[13]);

    fmpz_mul(t[23], t[14], t[15]);
    fmpz_mul(t[24], t[8], t[9]);
    fmpz_mul(t[25], t[10], t[11]);
    fmpz_add(t[10], t[8], t[14]);
    fmpz_add(t[11], t[9], t[15]);
    fmpz_mul(t[20], t[10], t[11]);

    fmpz_add(t[26], t[16], t[17]);
    fmpz_add(t[27], t[18], t[19]);
    fmpz_add(t[28], t[26], t[19]);
    fmpz_add(t[21], t[24], t[18]);
    fmpz_sub(t[0], t[28], t[21]);
    unity_zp_coeff_set_fmpz(f, 0, t[0]);

    fmpz_add(t[28], t[22], t[27]);
    fmpz_add(t[21], t[26], t[25]);
    fmpz_sub(t[0], t[28], t[21]);
    unity_zp_coeff_set_fmpz(f, 1, t[0]);
    fmpz_add(t[28], t[17], t[23]);
    fmpz_add(t[21], t[16], t[27]);

    fmpz_sub(t[0], t[28], t[21]);
    unity_zp_coeff_set_fmpz(f, 2, t[0]);
    fmpz_add(t[28], t[23], t[22]);
    fmpz_add(t[21], t[28], t[24]);
    fmpz_add(t[28], t[21], t[25]);
    fmpz_add(t[19], t[26], t[27]);

    fmpz_add(t[21], t[19], t[20]);
    fmpz_sub(t[0], t[21], t[28]);
    unity_zp_coeff_set_fmpz(f, 3, t[0]);
}

void
unity_zp_mul16(unity_zp f, const unity_zp g, const unity_zp h, fmpz_t * t)
{
}

