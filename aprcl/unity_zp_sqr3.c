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
    Computes f = g * g for p = 3. 
    g must be reduced by F_3 cyclotomic polynomial.
    t is the memory for fmpz_t; size of t must be > 3.
    Resulting f reduced by F_3 cyclotomic polynomial.
*/
void
unity_zp_sqr3(unity_zp f, const unity_zp g, fmpz_t * t)
{
    /*
        g = (x0, x1);
        f = (y0, y1);

        x0 = t[0]; x1 = t[1];
        m1 = t[2]; m2 = t[3];
        d1 = t[4].
    */
    fmpz_mod_poly_get_coeff_fmpz(t[0], g->poly, 0, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[1], g->poly, 1, g->ctx);

    fmpz_sub(t[2], t[0], t[1]);             /*  m1 = x0 - x1    */
    fmpz_add(t[3], t[0], t[1]);             /*  m2 = x0 + x1    */
    fmpz_mul(t[4], t[2], t[3]);             /*  d1 = m1 * m2    */
    fmpz_add(t[3], t[2], t[0]);             /*  m2 = m1 + m0    */
    unity_zp_coeff_set_fmpz(f, 0, t[4]);    /*  y0 = d1 mod n   */
    fmpz_mul(t[4], t[1], t[3]);             /*  d1 = x1 * m2    */
    unity_zp_coeff_set_fmpz(f, 1, t[4]);    /*  y1 = d1 mod n   */
}

/*
    Computes f = g * g for p = 3. 
    g must be reduced by F_3 cyclotomic polynomial.
    t is the memory for fmpz_t; size of t must be > 30.
    Resulting f reduced by F_3 cyclotomic polynomial.
*/
void
unity_zp_sqr9(unity_zp f, const unity_zp g, fmpz_t * t)
{
    /*
        g = (x0, x1, x2, x3, x4, x5);
        f = (y0, y1, y2, y3, y4, y5);

        x0 = t[20]; x1 = t[21]; x2 = t[22];
        x3 = t[23]; x4 = t[24]; x5 = t[25];

        a0 = t[0]; a1 = t[1]; a2 = t[2];
        b0 = t[3]; b1 = t[4]; b2 = t[5];
        c0 = t[6]; c1 = t[7]; c2 = t[8];
        c3 = t[9]; c4 = t[10];

        d0 = t[26]; d1 = t[27]; d2 = t[28];
        d3 = t[29]; d4 = t[30].
    */
    fmpz_mod_poly_get_coeff_fmpz(t[20], g->poly, 0, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[21], g->poly, 1, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[22], g->poly, 2, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[23], g->poly, 3, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[24], g->poly, 4, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[25], g->poly, 5, g->ctx);

    fmpz_sub(t[0], t[20], t[23]);           /*  a0 = x0 - x3    */
    fmpz_sub(t[1], t[21], t[24]);           /*  a1 = x1 - x4    */
    fmpz_sub(t[2], t[22], t[25]);           /*  a2 = x2 - x5    */
    fmpz_add(t[3], t[20], t[23]);           /*  b0 = x0 + x3    */
    fmpz_add(t[4], t[21], t[24]);           /*  b1 = x1 + x4    */
    fmpz_add(t[5], t[22], t[25]);           /*  b2 = x2 + x5    */

    /* 
        apply auxiliary routine 2 with (a0, a1, a2) and (b0, b1, b2)
        store result in (c0, .. , c4)
    */
    unity_zp_ar1(t);

    fmpz_set(t[26], t[6]);                  /*  set d0 = c0     */
    fmpz_set(t[27], t[7]);                  /*  set d1 = c1     */
    fmpz_set(t[28], t[8]);                  /*  set d2 = c2     */
    fmpz_set(t[29], t[9]);                  /*  set d3 = c3     */
    fmpz_set(t[30], t[10]);                 /*  set d4 = c4     */

    fmpz_add(t[3], t[20], t[0]);            /*  b0 = x0 + a0    */        
    fmpz_add(t[4], t[21], t[1]);            /*  b1 = x1 + a1    */
    fmpz_add(t[5], t[22], t[2]);            /*  b2 = x2 + a2    */
    fmpz_set(t[0], t[23]);                  /*  set a0 = x3     */
    fmpz_set(t[1], t[24]);                  /*  set a1 = x4     */
    fmpz_set(t[2], t[25]);                  /*  set a2 = x5     */

    /* 
        apply auxiliary routine 2 with (a0, a1, a2) and (b0, b1, b2)
        store result in (c0, .. , c4)
    */
    unity_zp_ar1(t);

    fmpz_sub(t[0], t[26], t[9]);            /*  a0 = d0 - c3    */
    unity_zp_coeff_set_fmpz(f, 0, t[0]);    /*  y0 = a0 mod n   */
    fmpz_sub(t[0], t[27], t[10]);           /*  a0 = d1 - c4    */
    unity_zp_coeff_set_fmpz(f, 1, t[0]);    /*  y1 = a0 mod n   */
    unity_zp_coeff_set_fmpz(f, 2, t[28]);   /*  y2 = d2 mod n   */
    fmpz_add(t[0], t[29], t[6]);            /*  a0 = d3 + c0    */
    fmpz_sub(t[1], t[0], t[9]);             /*  a1 = a0 - c3    */
    unity_zp_coeff_set_fmpz(f, 3, t[1]);    /*  y3 = a1 mod n   */
    fmpz_add(t[0], t[30], t[7]);            /*  a0 = d4 + c1    */
    fmpz_sub(t[1], t[0], t[10]);            /*  a1 = a0 - c4    */
    unity_zp_coeff_set_fmpz(f, 4, t[1]);    /*  y4 = a1 mod n   */
    unity_zp_coeff_set_fmpz(f, 5, t[8]);    /*  y5 = c2 mod n   */
}

