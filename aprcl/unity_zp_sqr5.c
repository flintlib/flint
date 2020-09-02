/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "aprcl.h"

/*
    Computes f = g * g for p = 5. 
    g must be reduced by F_5 cyclotomic polynomial.
    t is the memory for fmpz_t; size of t must be > 14.
    Resulting f reduced by F_5 cyclotomic polynomial.
*/
void
unity_zp_sqr5(unity_zp f, const unity_zp g, fmpz_t * t)
{
    /*
        g = (x0, x1, x2, x3);
        f = (y0, y1, y2, y3);

        x0 = t[0]; x1 = t[1]; x2 = t[2]; x3 = t[3];
        m1 = t[4]; m2 = t[5]; m3 = t[6]; m4 = t[7];
        m5 = t[8]; m6 = t[9]; m7 = t[10]; m8 = t[11];
        d1 = t[12]; d2 = t[13]; d3 = t[14]; d4 = t[15].
    */

    /* set xi */
    fmpz_mod_poly_get_coeff_fmpz(t[0], g->poly, 0, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[1], g->poly, 1, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[2], g->poly, 2, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[3], g->poly, 3, g->ctx);

    fmpz_sub(t[4], t[0], t[2]);             /*  m1 = x0 - x2    */
    fmpz_add(t[5], t[0], t[2]);             /*  m2 = x0 + x2    */
    fmpz_sub(t[6], t[2], t[1]);             /*  m3 = x2 - x1    */
    fmpz_sub(t[7], t[0], t[3]);             /*  m4 = x0 - x3    */
    fmpz_sub(t[8], t[1], t[0]);             /*  m5 = x1 - x0    */
    fmpz_sub(t[9], t[2], t[3]);             /*  m6 = x2 - x3    */

    fmpz_sub(t[10], t[1], t[3]);            /*  m7 = x1 - x3    */
    fmpz_add(t[11], t[3], t[3]);            /*  m8 = x3 + x3    */
    fmpz_mul(t[12], t[4], t[5]);            /*  d1 = m1 * m2    */
    fmpz_mul(t[13], t[6], t[11]);           /*  d2 = m3 * m8    */
    fmpz_add(t[14], t[12], t[13]);          /*  d3 = d1 + d2    */
    unity_zp_coeff_set_fmpz(f, 0, t[14]);   /*  y0 = d3 mod n   */
    fmpz_add(t[11], t[8], t[10]);           /*  m8 = m5 + m7    */

    fmpz_mul(t[13], t[7], t[11]);           /*  d2 = m4 * m8    */
    fmpz_add(t[15], t[12], t[13]);          /*  d4 = d1 + d2    */
    unity_zp_coeff_set_fmpz(f, 1, t[15]);   /*  y1 = d4 mod n   */
    fmpz_add(t[6], t[4], t[0]);             /*  m3 = m1 + x0    */
    fmpz_mul(t[12], t[2], t[6]);            /*  d1 = x2 * m3    */
    fmpz_sub(t[5], t[10], t[3]);            /*  m2 = m7 - x3    */
    fmpz_mul(t[13], t[5], t[1]);            /*  d2 = m2 * x1    */

    fmpz_add(t[14], t[12], t[13]);          /*  d3 = d1 + d2    */
    unity_zp_coeff_set_fmpz(f, 2, t[14]);   /*  y2 = d3 mod n   */
    fmpz_add(t[10], t[9], t[9]);            /*  m7 = m6 + m6    */
    fmpz_mul(t[13], t[10], t[8]);           /*  d2 = m7 * m5    */
    fmpz_add(t[14], t[12], t[13]);          /*  d3 = d1 + d2    */
    unity_zp_coeff_set_fmpz(f, 3, t[14]);   /*  y3 = d3 mod n   */
}

