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
    Computes f = g * h for p = 5. 
    g and h must be reduced by F_5 cyclotomic polynomial.
    t is the memory for fmpz_t; size of t must be > 24.
    Resulting f reduced by F_5 cyclotomic polynomial.
*/
void
unity_zp_mul5(unity_zp f, const unity_zp g, const unity_zp h, fmpz_t * t)
{
    /*
        g = (x0, x1, x2, x3);
        h = (y0, y1, y2, y3);
        f = (z0, z1, z2, z3);

        x0 = t[0]; x1 = t[1]; x2 = t[2]; x3 = t[3];
        y0 = t[4]; y1 = t[5]; y2 = t[6]; y3 = t[7];

        m1 = t[8]; m2 = t[9]; m3 = t[10]; m4 = t[11];
        m5 = t[12]; m6 = t[13]; m7 = t[14]; m8 = t[15];

        d0 = t[16]; d1 = t[17]; d2 = t[18]; d3 = t[19];
        d4 = t[20]; d5 = t[21]; d6 = t[22]; d7 = t[23];
        d8 = t[24]. 
    */

    /* set xi */
    fmpz_mod_poly_get_coeff_fmpz(t[0], g->poly, 0, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[1], g->poly, 1, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[2], g->poly, 2, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[3], g->poly, 3, g->ctx);

    /* set yi */
    fmpz_mod_poly_get_coeff_fmpz(t[4], h->poly, 0, h->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[5], h->poly, 1, h->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[6], h->poly, 2, h->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[7], h->poly, 3, h->ctx);

    fmpz_sub(t[8], t[1], t[3]);             /* m1 = x1 - x3  */
    fmpz_sub(t[9], t[5], t[7]);             /* m2 = y1 - y3  */
    fmpz_sub(t[10], t[2], t[3]);            /* m3 = x2 - x3  */
    fmpz_sub(t[11], t[7], t[6]);            /* m4 = y3 - y2  */
    fmpz_sub(t[12], t[0], t[1]);            /* m5 = x0 - x1  */
    fmpz_sub(t[13], t[5], t[4]);            /* m6 = y1 - y0  */

    fmpz_sub(t[14], t[0], t[2]);            /* m7 = x0 - x2  */
    fmpz_sub(t[15], t[6], t[4]);            /* m8 = y2 - y0  */
    fmpz_mul(t[16], t[0], t[4]);            /* d0 = x0 * y0  */
    fmpz_mul(t[18], t[8], t[9]);            /* d2 = m1 * m2  */
    fmpz_add(t[17], t[16], t[18]);          /* d1 = d0 + d2  */
    fmpz_mul(t[18], t[10], t[11]);          /* d2 = m3 * m4  */

    fmpz_mul(t[19], t[12], t[13]);          /* d3 = m5 * m6  */
    fmpz_mul(t[20], t[14], t[15]);          /* d4 = m7 * m8  */
    fmpz_mul(t[21], t[1], t[5]);            /* d5 = x1 * y1  */
    fmpz_mul(t[22], t[2], t[6]);            /* d6 = x2 * y2  */
    fmpz_mul(t[23], t[3], t[7]);            /* d7 = x3 * y3  */
    fmpz_add(t[24], t[17], t[18]);          /* d8 = d1 + d2  */

    fmpz_sub(t[0], t[24], t[21]);           /* x0 = d8 - d5  */
    unity_zp_coeff_set_fmpz(f, 0, t[0]);    /* z0 = x0 mod n */

    fmpz_add(t[24], t[17], t[19]);          /* d8 = d1 + d3  */
    fmpz_sub(t[0], t[24], t[22]);           /* x0 = d8 - d6  */
    unity_zp_coeff_set_fmpz(f, 1, t[0]);    /* z1 = x0 mod n */

    fmpz_add(t[24], t[17], t[20]);          /* d8 = d1 + d4  */
    fmpz_sub(t[0], t[24], t[23]);           /* x0 = d8 - d7  */
    unity_zp_coeff_set_fmpz(f, 2, t[0]);    /* z2 = x0 mod n */

    fmpz_sub(t[10], t[8], t[14]);           /* m3 = m1 - m7  */
    fmpz_add(t[11], t[9], t[15]);           /* m4 = m2 + m8  */
    fmpz_mul(t[17], t[10], t[11]);          /* d1 = m3 * m4  */
    fmpz_add(t[24], t[16], t[17]);          /* d8 = d0 + d1  */
    fmpz_add(t[23], t[24], t[18]);          /* d7 = d8 + d2  */
    fmpz_add(t[24], t[23], t[19]);          /* d8 = d7 + d3  */
    fmpz_add(t[0], t[24], t[20]);           /* d7 = d8 + d4  */
    unity_zp_coeff_set_fmpz(f, 3, t[0]);    /* z3 = d7 mod n */
}

