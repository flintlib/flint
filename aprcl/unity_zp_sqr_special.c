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

/* computes f = g * g for p = 7 */
void
unity_zp_sqr7(unity_zp f, const unity_zp g, fmpz_t * t)
{
    fmpz_mod_poly_get_coeff_fmpz(t[0], g->poly, 0);
    fmpz_mod_poly_get_coeff_fmpz(t[1], g->poly, 1);
    fmpz_mod_poly_get_coeff_fmpz(t[2], g->poly, 2);
    fmpz_mod_poly_get_coeff_fmpz(t[3], g->poly, 3);
    fmpz_mod_poly_get_coeff_fmpz(t[4], g->poly, 4);
    fmpz_mod_poly_get_coeff_fmpz(t[5], g->poly, 5);

    fmpz_sub(t[6], t[0], t[1]);
    fmpz_sub(t[7], t[1], t[2]);
    fmpz_sub(t[8], t[2], t[3]);
    fmpz_sub(t[9], t[3], t[4]);
    fmpz_sub(t[10], t[5], t[4]);
    fmpz_add(t[11], t[6], t[7]);

    fmpz_add(t[12], t[7], t[8]);
    fmpz_add(t[13], t[8], t[9]);
    fmpz_sub(t[14], t[3], t[5]);
    fmpz_add(t[15], t[8], t[11]);
    fmpz_add(t[16], t[9], t[12]);

    fmpz_add(t[18], t[11], t[13]);
    fmpz_add(t[19], t[12], t[14]);
    fmpz_add(t[21], t[0], t[1]);
    fmpz_add(t[22], t[0], t[15]);
    fmpz_mul(t[24], t[3], t[22]);
    
    fmpz_sub(t[22], t[19], t[4]);
    fmpz_add(t[23], t[19], t[4]);
    fmpz_mul(t[25], t[22], t[23]);
    fmpz_sub(t[22], t[13], t[7]);
    fmpz_mul(t[26], t[16], t[22]);

    fmpz_add(t[22], t[19], t[14]);
    fmpz_mul(t[27], t[22], t[12]);
    fmpz_add(t[22], t[1], t[1]);
    fmpz_mul(t[28], t[22], t[11]);
    fmpz_mul(t[29], t[6], t[21]);
    fmpz_add(t[22], t[8], t[8]);

    fmpz_add(t[7], t[0], t[18]);
    fmpz_mul(t[30], t[22], t[10]);
    fmpz_add(t[31], t[24], t[25]);
    fmpz_add(t[24], t[31], t[26]);
    unity_zp_coeff_set_fmpz(f, 3, t[24]);
    fmpz_add(t[31], t[26], t[27]);

    fmpz_add(t[24], t[31], t[28]);
    unity_zp_coeff_set_fmpz(f, 1, t[24]);
    fmpz_add(t[31], t[27], t[29]);
    fmpz_add(t[24], t[31], t[30]);
    unity_zp_coeff_set_fmpz(f, 0, t[24]);
    fmpz_add(t[22], t[12], t[19]);

    fmpz_mul(t[24], t[14], t[22]);
    fmpz_sub(t[22], t[13], t[5]);
    fmpz_add(t[23], t[2], t[10]);
    fmpz_mul(t[25], t[22], t[23]);
    fmpz_mul(t[26], t[7], t[4]);
    fmpz_add(t[22], t[8], t[13]);

    fmpz_mul(t[27], t[22], t[9]);
    fmpz_add(t[22], t[6], t[6]);
    fmpz_mul(t[28], t[22], t[10]);
    fmpz_sub(t[22], t[19], t[10]);
    fmpz_mul(t[29], t[22], t[16]);

    fmpz_add(t[22], t[2], t[2]);
    fmpz_mul(t[30], t[22], t[15]);
    fmpz_add(t[31], t[24], t[25]);
    fmpz_add(t[24], t[31], t[26]);
    unity_zp_coeff_set_fmpz(f, 4, t[24]);
    fmpz_add(t[31], t[26], t[27]);

    fmpz_add(t[24], t[31], t[28]);
    unity_zp_coeff_set_fmpz(f, 5, t[24]);
    fmpz_add(t[31], t[27], t[29]);
    fmpz_add(t[24], t[31], t[30]);
    unity_zp_coeff_set_fmpz(f, 2, t[24]);
}


void
unity_zp_ar2(fmpz_t * t)
{
    /*
        a0 = t[0]; a1 = t[1]; a2 = t[2]; a3 = t[3];
        b0 = t[4]; b1 = t[5]; b2 = t[6]; b3 = t[7];

        c0 = t[8]; c1 = t[9]; c2 = t[10]; c3 = t[11];
        c4 = t[12]; c5 = t[13]; c6 = t[14];

        m1 = t[15]; m2 = t[16]; m3 = t[17]; m4 = t[18];
        m5 = t[19];

        d1 = t[20]; d2 = t[21]; d3 = t[22]; d4 = t[23];
        d5 = t[24]; d6 = t[25]; d7 = t[26]; d8 = t[27];
        d9 = t[28].
    */

    fmpz_mul(t[8], t[0], t[4]);
    fmpz_mul(t[20], t[1], t[5]);
    fmpz_mul(t[21], t[2], t[6]);
    fmpz_mul(t[14], t[3], t[7]);
    fmpz_add(t[15], t[0], t[1]);
    fmpz_add(t[16], t[4], t[5]);
    fmpz_mul(t[22], t[15], t[16]);

    fmpz_add(t[15], t[0], t[2]);
    fmpz_add(t[16], t[4], t[6]);
    fmpz_mul(t[23], t[15], t[16]);
    fmpz_add(t[17], t[2], t[3]);
    fmpz_add(t[18], t[6], t[7]);
    fmpz_mul(t[24], t[17], t[18]);

    fmpz_add(t[17], t[1], t[3]);
    fmpz_add(t[18], t[5], t[7]);
    fmpz_mul(t[25], t[17], t[18]);
    fmpz_add(t[26], t[8], t[20]);
    fmpz_sub(t[9], t[22], t[26]);
    fmpz_add(t[26], t[8], t[21]);

    fmpz_add(t[27], t[20], t[23]);
    fmpz_sub(t[10], t[27], t[26]);
    fmpz_add(t[19], t[15], t[17]);
    fmpz_add(t[17], t[16], t[18]);
    fmpz_add(t[26], t[21], t[14]);
    fmpz_sub(t[13], t[24], t[26]);

    fmpz_mul(t[26], t[17], t[19]);
    fmpz_add(t[27], t[9], t[13]);
    fmpz_add(t[28], t[27], t[25]);
    fmpz_add(t[27], t[28], t[23]);
    fmpz_sub(t[11], t[26], t[27]);
    fmpz_add(t[26], t[25], t[21]);

    fmpz_add(t[27], t[20], t[14]);
    fmpz_sub(t[12], t[26], t[27]);
}

