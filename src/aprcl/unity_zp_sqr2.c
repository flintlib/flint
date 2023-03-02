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
    Computes f = g * g for p = 2^2. 
    g must be reduced by F_4 cyclotomic polynomial.
    t is the memory for fmpz_t; size of t must be > 3.
    Resulting f reduced by F_4 cyclotomic polynomial.
*/
void
unity_zp_sqr4(unity_zp f, const unity_zp g, fmpz_t * t)
{
    /*
        g = (x0, x1);
        f = (y0, y1);

        x0 = t[0]; x1 = t[1];
        m1 = t[2]; m2 = t[3];
        d1 = t[4].
    */

    /* set xi */
    fmpz_mod_poly_get_coeff_fmpz(t[0], g->poly, 0, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[1], g->poly, 1, g->ctx);

    fmpz_sub(t[2], t[0], t[1]);             /*  m1 = x0 - x1    */
    fmpz_add(t[3], t[0], t[1]);             /*  m2 = x0 + x1    */
    fmpz_mul(t[4], t[2], t[3]);             /*  d1 = m1 * m2    */
    fmpz_add(t[2], t[0], t[0]);             /*  m1 = x0 + x0    */
    unity_zp_coeff_set_fmpz(f, 0, t[4]);    /*  y0 = d1 mod n   */
    fmpz_mul(t[4], t[2], t[1]);             /*  d1 = m1 * x1    */
    unity_zp_coeff_set_fmpz(f, 1, t[4]);    /*  y1 = d1 mod n   */
}

/*
    Computes f = g * g for p = 2^3. 
    g must be reduced by F_8 cyclotomic polynomial.
    t is the memory for fmpz_t; size of t must be > 16.
    Resulting f reduced by F_8 cyclotomic polynomial.
*/
void
unity_zp_sqr8(unity_zp f, const unity_zp g, fmpz_t * t)
{
    /*
        g = (x0, x1, x2, x3);
        f = (y0, y1, y2, y3);

        x0 = t[0]; x1 = t[1]; x2 = t[2]; x3 = t[3];
        m1 = t[4]; m2 = t[5]; m3 = t[6]; m4 = t[7];
        m5 = t[8]; m6 = t[9]; m7 = t[10]; m8 = t[11];
        d1 = t[12]; d2 = t[13]; d3 = t[14]; d4 = t[15];
        d5 = t[16]; d6 = t[17].
    */

    /* set xi */
    fmpz_mod_poly_get_coeff_fmpz(t[0], g->poly, 0, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[1], g->poly, 1, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[2], g->poly, 2, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[3], g->poly, 3, g->ctx);

    fmpz_sub(t[4], t[0], t[2]);             /*  m1 = x0 - x2    */
    fmpz_add(t[5], t[0], t[2]);             /*  m2 = x0 + x2    */                     
    fmpz_sub(t[6], t[1], t[3]);             /*  m3 = x1 - x3    */
    fmpz_add(t[7], t[1], t[3]);             /*  m4 = x1 + x3    */
    fmpz_add(t[8], t[0], t[0]);             /*  m5 = x0 + x0    */
    fmpz_add(t[9], t[1], t[1]);             /*  m6 = x1 + x1    */

    fmpz_add(t[10], t[4], t[6]);            /*  m7 = m1 + m3    */
    fmpz_add(t[11], t[5], t[7]);            /*  m8 = m2 + m4    */
    fmpz_mul(t[12], t[4], t[5]);            /*  d1 = m1 * m2    */
    fmpz_mul(t[13], t[6], t[7]);            /*  d2 = m3 * m4    */
    fmpz_mul(t[14], t[9], t[3]);            /*  d3 = m6 * x3    */
    fmpz_mul(t[15], t[8], t[2]);            /*  d4 = m5 * x2    */

    fmpz_add(t[5], t[2], t[3]);             /*  m2 = x2 + x3    */
    fmpz_sub(t[16], t[12], t[14]);          /*  d5 = d1 - d3    */
    unity_zp_coeff_set_fmpz(f, 0, t[16]);   /*  y0 = d5 mod n   */
    fmpz_add(t[17], t[13], t[15]);          /*  d6 = d2 + d4    */
    unity_zp_coeff_set_fmpz(f, 2, t[17]);   /*  y2 = d6 mod n   */
    fmpz_mul(t[16], t[10], t[11]);          /*  d5 = m7 * m8    */

    fmpz_add(t[17], t[12], t[13]);          /*  d6 = d1 + d2    */
    fmpz_sub(t[13], t[16], t[17]);          /*  d2 = d5 - d6    */
    unity_zp_coeff_set_fmpz(f, 1, t[13]);   /*  y1 = d2 mod n   */
    fmpz_add(t[4], t[8], t[9]);             /*  m1 = m5 + m6    */
    fmpz_mul(t[12], t[4], t[5]);            /*  d1 = m1 * m2    */
    fmpz_add(t[17], t[14], t[15]);          /*  d6 = d3 + d4    */

    fmpz_sub(t[13], t[12], t[17]);          /*  d2 = d1 - d6    */
    unity_zp_coeff_set_fmpz(f, 3, t[13]);   /*  y3 = d2 mod n   */
}

/*
    Computes f = g * g for p = 2^4. 
    g must be reduced by F_16 cyclotomic polynomial.
    t is the memory for fmpz_t; size of t must be > 44.
    Resulting f reduced by F_16 cyclotomic polynomial.
*/
void
unity_zp_sqr16(unity_zp f, const unity_zp g, fmpz_t * t)
{
    ulong i;

    /*
        g = (x0, x1, x2, x3, x4, x5, x6, x7);
        f = (y0, y1, y2, y3, y4, y5, y6, y7);

        x0 = t[30]; x1 = t[31]; x2 = t[32]; x3 = t[33];
        x4 = t[34]; x5 = t[35]; x6 = t[37];

        a0 = t[0]; a1 = t[1]; a2 = t[2]; a3 = t[3];
        b0 = t[4]; b1 = t[5]; b2 = t[6]; b3 = t[7];
        c0 = t[8]; c1 = t[9]; c2 = t[10]; c3 = t[11];
        c4 = t[12]; c5 = t[13]; c6 = t[14];

        d0 = t[38]; d1 = t[39]; d2 = t[40]; d3 = t[41];
        d4 = t[42]; d5 = t[43]; d6 = t[44]; d7 = t[45].
    */

    /* set xi */
    for (i = 0; i < 8; i++)
        fmpz_mod_poly_get_coeff_fmpz(t[30 + i], g->poly, i, g->ctx);

    fmpz_add(t[0], t[30], t[34]);           /*  a0 = x0 + x4    */
    fmpz_add(t[1], t[31], t[35]);           /*  a1 = x1 + x5    */
    fmpz_add(t[2], t[32], t[36]);           /*  a2 = x2 + x6    */
    fmpz_add(t[3], t[33], t[37]);           /*  a3 = x3 + x7    */
    fmpz_sub(t[4], t[30], t[34]);           /*  b0 = x0 - x4    */
    fmpz_sub(t[5], t[31], t[35]);           /*  b1 = x1 - x5    */
    fmpz_sub(t[6], t[32], t[36]);           /*  b2 = x2 - x6    */
    fmpz_sub(t[7], t[33], t[37]);           /*  b3 = x3 - x7    */

    /* 
        apply auxiliary routine 2 with (a0, .. , a3) and (b0, .. , b3)
        store result in (c0, .. , c6)
    */
    unity_zp_ar2(t);

    /* set d_i = c_i */
    for (i = 8; i < 15; i++)
        fmpz_set(t[30 + i], t[i]);

    fmpz_add(t[0], t[30], t[30]);           /*  a0 = x0 + x0    */
    fmpz_add(t[1], t[31], t[31]);           /*  a1 = x1 + x1    */
    fmpz_add(t[2], t[32], t[32]);           /*  a2 = x2 + x2    */
    fmpz_add(t[3], t[33], t[33]);           /*  a3 = x3 + x3    */
    fmpz_set(t[4], t[34]);                  /*  set b0 = x4     */ 
    fmpz_set(t[5], t[35]);                  /*  set b1 = x5     */
    fmpz_set(t[6], t[36]);                  /*  set b2 = x6     */
    fmpz_set(t[7], t[37]);                  /*  set b3 = x7     */
    
    /* 
        apply auxiliary routine 2 with (a0, .. , a3) and (b0, .. , b3)
        store result in (c0, .. , c6)
    */
    unity_zp_ar2(t);

    fmpz_sub(t[16], t[38], t[12]);          /*  d7 = d0 - c4    */
    unity_zp_coeff_set_fmpz(f, 0, t[16]);   /*  y0 = d7 mod n   */
    fmpz_sub(t[16], t[39], t[13]);          /*  d7 = d1 - c5    */
    unity_zp_coeff_set_fmpz(f, 1, t[16]);   /*  y1 = d7 mod n   */
    fmpz_sub(t[16], t[40], t[14]);          /*  d7 = d2 - c6    */
    unity_zp_coeff_set_fmpz(f, 2, t[16]);   /*  y2 = d7 mod n   */
    unity_zp_coeff_set_fmpz(f, 3, t[41]);   /*  y3 = d3 mod n   */
    fmpz_add(t[16], t[42], t[8]);           /*  d7 = d4 + c0    */
    unity_zp_coeff_set_fmpz(f, 4, t[16]);   /*  y4 = d7 mod n   */
    fmpz_add(t[16], t[43], t[9]);           /*  d7 = d5 + c1    */
    unity_zp_coeff_set_fmpz(f, 5, t[16]);   /*  y5 = d7 mod n   */
    fmpz_add(t[16], t[44], t[10]);          /*  d7 = d6 + c2    */
    unity_zp_coeff_set_fmpz(f, 6, t[16]);   /*  y6 = d7 mod n   */
    unity_zp_coeff_set_fmpz(f, 7, t[11]);   /*  y7 = c3 mod n   */
}

