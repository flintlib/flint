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
unity_zp_mont_sqr7(unity_zp_mont f, const unity_zp_mont g, fmpz_t * t)
{
    fmpz_poly_get_coeff_fmpz(t[0], g->poly, 0);
    fmpz_poly_get_coeff_fmpz(t[1], g->poly, 1);
    fmpz_poly_get_coeff_fmpz(t[2], g->poly, 2);
    fmpz_poly_get_coeff_fmpz(t[3], g->poly, 3);
    fmpz_poly_get_coeff_fmpz(t[4], g->poly, 4);
    fmpz_poly_get_coeff_fmpz(t[5], g->poly, 5);

    mod_sub(t[6], t[0], t[1], f->n);
    mod_sub(t[7], t[1], t[2], f->n);
    mod_sub(t[8], t[2], t[3], f->n);
    mod_sub(t[9], t[3], t[4], f->n);
    mod_sub(t[10], t[5], t[4], f->n);
    mod_add(t[11], t[6], t[7], f->n);

    mod_add(t[12], t[7], t[8], f->n);
    mod_add(t[13], t[8], t[9], f->n);
    mod_sub(t[14], t[3], t[5], f->n);
    mod_add(t[15], t[8], t[11], f->n);
    mod_add(t[16], t[9], t[12], f->n);

    mod_add(t[18], t[11], t[13], f->n);
    mod_add(t[19], t[12], t[14], f->n);
    mod_add(t[21], t[0], t[1], f->n);
    mod_add(t[22], t[0], t[15], f->n);
    mod_mul(t[24], t[3], t[22], f->n, f->ninv, t[48], t[49], f->r);
    
    mod_sub(t[22], t[19], t[4], f->n);
    mod_add(t[23], t[19], t[4], f->n);
    mod_mul(t[25], t[22], t[23], f->n, f->ninv, t[48], t[49], f->r);
    mod_sub(t[22], t[13], t[7], f->n);
    mod_mul(t[26], t[16], t[22], f->n, f->ninv, t[48], t[49], f->r);

    mod_add(t[22], t[19], t[14], f->n);
    mod_mul(t[27], t[22], t[12], f->n, f->ninv, t[48], t[49], f->r);
    mod_add(t[22], t[1], t[1], f->n);
    mod_mul(t[28], t[22], t[11], f->n, f->ninv, t[48], t[49], f->r);
    mod_mul(t[29], t[6], t[21], f->n, f->ninv, t[48], t[49], f->r);
    mod_add(t[22], t[8], t[8], f->n);

    mod_add(t[7], t[0], t[18], f->n);
    mod_mul(t[30], t[22], t[10], f->n, f->ninv, t[48], t[49], f->r);
    mod_add(t[31], t[24], t[25], f->n);
    mod_add(t[24], t[31], t[26], f->n);
    unity_zp_coeff_set_fmpz(f, 3, t[24]);
    mod_add(t[31], t[26], t[27], f->n);

    mod_add(t[24], t[31], t[28], f->n);
    fmpz_poly_get_coeff_fmpz(f->poly, 1, t[24]);
    mod_add(t[31], t[27], t[29], f->n);
    mod_add(t[24], t[31], t[30], f->n);
    fmpz_poly_get_coeff_fmpz(f->poly, 0, t[24]);
    mod_add(t[22], t[12], t[19], f->n);

    mod_mul(t[24], t[14], t[22], f->n, f->ninv, t[48], t[49], f->r);
    mod_sub(t[22], t[13], t[5], f->n);
    mod_add(t[23], t[2], t[10], f->n);
    mod_mul(t[25], t[22], t[23], f->n, f->ninv, t[48], t[49], f->r);
    mod_mul(t[26], t[7], t[4], f->n, f->ninv, t[48], t[49], f->r);
    mod_add(t[22], t[8], t[13], f->n);

    mod_mul(t[27], t[22], t[9], f->n, f->ninv, t[48], t[49], f->r);
    mod_add(t[22], t[6], t[6], f->n);
    mod_mul(t[28], t[22], t[10], f->n, f->ninv, t[48], t[49], f->r);
    mod_sub(t[22], t[19], t[10], f->n);
    mod_mul(t[29], t[22], t[16], f->n, f->ninv, t[48], t[49], f->r);

    mod_add(t[22], t[2], t[2], f->n);
    mod_mul(t[30], t[22], t[15], f->n, f->ninv, t[48], t[49], f->r);
    mod_add(t[31], t[24], t[25], f->n);
    mod_add(t[24], t[31], t[26], f->n);
    fmpz_poly_get_coeff_fmpz(f->poly, 4, t[24]);
    mod_add(t[31], t[26], t[27], f->n);

    mod_add(t[24], t[31], t[28], f->n);
    fmpz_poly_get_coeff_fmpz(f->poly, 5, t[24]);
    mod_add(t[31], t[27], t[29], f->n);
    mod_add(t[24], t[31], t[30], f->n);
    fmpz_poly_get_coeff_fmpz(f->poly, 2, t[24]);
}

