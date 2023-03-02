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
    Computes f = g * h for p = 2^2. 
    g and h must be reduced by F_4 cyclotomic polynomial.
    t is the memory for fmpz_t; size of t must be > 9.
    Resulting f reduced by F_4 cyclotomic polynomial.
*/
void
unity_zp_mul4(unity_zp f, const unity_zp g, const unity_zp h, fmpz_t * t)
{
    /*
        g = (x0, x1);
        h = (y0, y1);
        f = (z0, z1);

        x0 = t[0]; x1 = t[1];
        y0 = t[2]; y1 = t[3];

        m1 = t[4]; m2 = t[5]; m3 = t[6];
        d1 = t[7]; d2 = t[8]; d3 = t[9].
    */

    /* set xi */
    fmpz_mod_poly_get_coeff_fmpz(t[0], g->poly, 0, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[1], g->poly, 1, g->ctx);

    /* set yi */
    fmpz_mod_poly_get_coeff_fmpz(t[2], h->poly, 0, h->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[3], h->poly, 1, h->ctx);

    fmpz_add(t[4], t[0], t[1]);             /*  m1 = x0 + x1    */
    fmpz_add(t[5], t[2], t[3]);             /*  m2 = y0 + y1    */
    fmpz_sub(t[6], t[3], t[2]);             /*  m3 = y1 - y0    */
    fmpz_mul(t[7], t[4], t[2]);             /*  d1 = m1 * y0    */
    fmpz_mul(t[8], t[5], t[1]);             /*  d2 = m2 * x1    */
    fmpz_mul(t[9], t[6], t[0]);             /*  d3 = m3 * x0    */
    fmpz_sub(t[0], t[7], t[8]);             /*  t[0] = d1 - d2  */

    unity_zp_coeff_set_fmpz(f, 0, t[0]);    /*  z0 = t[0] mod n */
    fmpz_add(t[0], t[7], t[9]);             /*  t[0] = d1 + d3  */
    unity_zp_coeff_set_fmpz(f, 1, t[0]);    /*  z1 = t[0] mod n */
}

/*
    Computes f = g * h for p = 2^3. 
    g and h must be reduced by F_8 cyclotomic polynomial.
    t is the memory for fmpz_t; size of t must be > 27.
    Resulting f reduced by F_8 cyclotomic polynomial.
*/
void
unity_zp_mul8(unity_zp f, const unity_zp g, const unity_zp h, fmpz_t * t)
{
    /*
        g = (x0, x1, x2, x3);
        h = (y0, y1, y2, y3);
        f = (z0, z1, z2, z3);

        x0 = t[0]; x1 = t[1]; x2 = t[2]; x3 = t[3];
        y0 = t[4]; y1 = t[5]; y2 = t[6]; y3 = t[7];

        m1 = t[8]; ... ; m8 = t[15];
        d1 = t[16]; ... ; d12 = t[28].
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

    fmpz_add(t[8], t[1], t[3]);             /*  m1 = x1 + x3    */
    fmpz_add(t[9], t[5], t[7]);             /*  m2 = y1 + y3    */
    fmpz_add(t[10], t[2], t[3]);            /*  m3 = x2 + x3    */
    fmpz_add(t[11], t[6], t[7]);            /*  m4 = y2 + y3    */
    fmpz_add(t[12], t[0], t[1]);            /*  m5 = x0 + x1    */
    fmpz_add(t[13], t[4], t[5]);            /*  m6 = y0 + y1    */

    fmpz_add(t[14], t[0], t[2]);            /*  m7 = x0 + x2    */
    fmpz_add(t[15], t[4], t[6]);            /*  m8 = y0 + y2    */
    fmpz_mul(t[16], t[0], t[4]);            /*  d0 = x0 * y0    */
    fmpz_mul(t[17], t[1], t[5]);            /*  d1 = x1 * y1    */
    fmpz_mul(t[18], t[2], t[6]);            /*  d2 = x2 * y2    */
    fmpz_mul(t[19], t[3], t[7]);            /*  d3 = x3 * y3    */
    fmpz_mul(t[22], t[12], t[13]);          /*  d6 = m5 * m6    */

    fmpz_mul(t[23], t[14], t[15]);          /*  d7 = m7 * m8    */
    fmpz_mul(t[24], t[8], t[9]);            /*  d8 = m1 * m2    */
    fmpz_mul(t[25], t[10], t[11]);          /*  d9 = m3 * m4    */
    fmpz_add(t[10], t[8], t[14]);           /*  m3 = m1 + m7    */
    fmpz_add(t[11], t[9], t[15]);           /*  m4 = m2 + m8    */
    fmpz_mul(t[20], t[10], t[11]);          /*  d4 = m3 * m4    */

    fmpz_add(t[26], t[16], t[17]);          /*  d10 = d0 + d1   */
    fmpz_add(t[27], t[18], t[19]);          /*  d11 = d2 + d3   */
    fmpz_add(t[28], t[26], t[19]);          /*  d12 = d10 + d3  */
    fmpz_add(t[21], t[24], t[18]);          /*  d5 = d8 + d2    */
    fmpz_sub(t[0], t[28], t[21]);           /*  t[0] = d12 - d5 */
    unity_zp_coeff_set_fmpz(f, 0, t[0]);    /*  z0 = t[0] mod n */

    fmpz_add(t[28], t[22], t[27]);          /*  d12 = d6 + d11  */
    fmpz_add(t[21], t[26], t[25]);          /*  d5 = d10 + d9   */
    fmpz_sub(t[0], t[28], t[21]);           /*  t[0] = d12 - d5 */
    unity_zp_coeff_set_fmpz(f, 1, t[0]);    /*  z1 = t[0] mod n */
    fmpz_add(t[28], t[17], t[23]);          /*  d12 = d1 + d7   */
    fmpz_add(t[21], t[16], t[27]);          /*  d5 = d0 + d11   */

    fmpz_sub(t[0], t[28], t[21]);           /*  t[0] = d12 - d5 */
    unity_zp_coeff_set_fmpz(f, 2, t[0]);    /*  z2 = t[0] mod n */
    fmpz_add(t[28], t[23], t[22]);          /*  d12 = d7 + d6   */
    fmpz_add(t[21], t[28], t[24]);          /*  d5 = d12 + d8   */
    fmpz_add(t[28], t[21], t[25]);          /*  d12 = d5 + d8   */
    fmpz_add(t[19], t[26], t[27]);          /*  d3 = d10 + d11  */

    fmpz_add(t[21], t[19], t[20]);          /*  d5 = d3 + d4    */
    fmpz_sub(t[0], t[21], t[28]);           /*  t[0] = d5 - d12 */
    unity_zp_coeff_set_fmpz(f, 3, t[0]);    /*  z3 = t[0] mod n */
}

/*
    Computes f = g * h for p = 2^4. 
    g and h must be reduced by F_16 cyclotomic polynomial.
    t is the memory for fmpz_t; size of t must be > 63.
    Resulting f reduced by F_16 cyclotomic polynomial.
*/
void
unity_zp_mul16(unity_zp f, const unity_zp g, const unity_zp h, fmpz_t * t)
{
    int i;

    /*
        g = (x0, ... , x7);
        h = (y0, ... , y7);
        f = (z0, ... , z7);

        x0 = t[30]; x1 = t[31]; x2 = t[32]; x3 = t[33];
        x4 = t[34]; x5 = t[35]; x6 = t[36]; x7 = t[37];

        y0 = t[40]; y1 = t[41]; y2 = t[42]; y3 = t[43];
        y4 = t[44]; y5 = t[45]; y6 = t[46]; y7 = t[47];

        a0 = t[0]; ... ; a3 = t[3];
        b0 = t[4]; ... ; b3 = t[7];
        c0 = t[8]; ... ; t[14];

        d0 = t[50]; ... ; d12 = t[62];
    */

    /* set xi */
    for (i = 0; i < 8; i++)
        fmpz_mod_poly_get_coeff_fmpz(t[30 + i], g->poly, i, g->ctx);

    /* set yi */
    for (i = 0; i < 8; i++)
        fmpz_mod_poly_get_coeff_fmpz(t[40 + i], h->poly, i, h->ctx);

    fmpz_add(t[0], t[30], t[34]);           /*  a0 = x0 + x4    */
    fmpz_add(t[1], t[31], t[35]);           /*  a1 = x1 + x5    */
    fmpz_add(t[2], t[32], t[36]);           /*  a2 = x2 + x6    */
    fmpz_add(t[3], t[33], t[37]);           /*  a3 = x3 + x7    */
    fmpz_set(t[4], t[40]);                  /*  set b0 = y0     */
    fmpz_set(t[5], t[41]);                  /*  set b1 = y1     */
    fmpz_set(t[6], t[42]);                  /*  set b2 = y2     */
    fmpz_set(t[7], t[43]);                  /*  set b3 = y3     */

    /* 
        apply auxiliary routine 2 with (a0, .. , a3) and (b0, .. , b3)
        store result in (c0, .. , c6)
    */
    unity_zp_ar2(t);

    /* set d_i = c_i */
    for (i = 0; i < 7; i++)
        fmpz_set(t[50 + i], t[8 + i]);

    fmpz_add(t[0], t[40], t[44]);           /*  a0 = y0 + y4    */
    fmpz_add(t[1], t[41], t[45]);           /*  a1 = y1 + y5    */
    fmpz_add(t[2], t[42], t[46]);           /*  a2 = y2 + y6    */
    fmpz_add(t[3], t[43], t[47]);           /*  a3 = y3 + y7    */
    fmpz_set(t[4], t[34]);                  /*  set b0 = x4     */
    fmpz_set(t[5], t[35]);                  /*  set b1 = x5     */
    fmpz_set(t[6], t[36]);                  /*  set b2 = x6     */
    fmpz_set(t[7], t[37]);                  /*  set b3 = x7     */

    /* 
        apply auxiliary routine 2 with (a0, .. , a3) and (b0, .. , b3)
        store result in (c0, .. , c6)
    */
    unity_zp_ar2(t);

    /* set d_{7 + i} = c_i */
    for (i = 0; i < 7; i++)
        fmpz_set(t[57 + i], t[8 + i]);

    fmpz_sub(t[0], t[44], t[40]);           /*  a0 = y4 - y0    */
    fmpz_sub(t[1], t[45], t[41]);           /*  a1 = y5 - y1    */
    fmpz_sub(t[2], t[46], t[42]);           /*  a2 = y6 - y2    */
    fmpz_sub(t[3], t[47], t[43]);           /*  a3 = y7 - y3    */
    fmpz_set(t[4], t[30]);                  /*  set b0 = x0     */
    fmpz_set(t[5], t[31]);                  /*  set b1 = x1     */
    fmpz_set(t[6], t[32]);                  /*  set b2 = x2     */
    fmpz_set(t[7], t[33]);                  /*  set b3 = x3     */

    /* 
        apply auxiliary routine 2 with (a0, .. , a3) and (b0, .. , b3)
        store result in (c0, .. , c6)
    */
    unity_zp_ar2(t);

    fmpz_add(t[1], t[54], t[57]);           /*  a1 = d4 + d7    */
    fmpz_add(t[2], t[1], t[12]);            /*  a2 = a1 + c3    */
    fmpz_sub(t[0], t[50], t[2]);            /*  a0 = d0 - a2    */
    unity_zp_coeff_set_fmpz(f, 0, t[0]);    /*  z0 = a0 mod n   */

    fmpz_add(t[1], t[55], t[58]);           /*  a1 = d5 + d8    */
    fmpz_add(t[2], t[1], t[13]);            /*  a2 = a1 + c4    */
    fmpz_sub(t[0], t[51], t[2]);            /*  a0 = d1 - a2    */
    unity_zp_coeff_set_fmpz(f, 1, t[0]);    /*  z1 = a0 mod n   */

    fmpz_add(t[1], t[56], t[59]);           /*  a1 = d6 + d9    */
    fmpz_add(t[2], t[1], t[14]);            /*  a2 = a1 + c5    */
    fmpz_sub(t[0], t[52], t[2]);            /*  a0 = d2 - a2    */
    unity_zp_coeff_set_fmpz(f, 2, t[0]);    /*  z2 = a0 mod n   */

    fmpz_sub(t[0], t[53], t[60]);           /*  a0 = d3 - d10   */
    unity_zp_coeff_set_fmpz(f, 3, t[0]);    /*  z3 = a0 mod n   */

    fmpz_add(t[1], t[54], t[50]);           /*  a1 = d4 + d0    */
    fmpz_add(t[2], t[1], t[8]);             /*  a2 = a1 + c0    */
    fmpz_sub(t[0], t[2], t[61]);            /*  a0 = a2 - d11   */
    unity_zp_coeff_set_fmpz(f, 4, t[0]);    /*  z4 = a0 mod n   */

    fmpz_add(t[1], t[55], t[51]);           /*  a1 = d5 + d1    */
    fmpz_add(t[2], t[1], t[9]);             /*  a2 = a1 + c0    */
    fmpz_sub(t[0], t[2], t[62]);            /*  a0 = a2 - d12   */
    unity_zp_coeff_set_fmpz(f, 5, t[0]);    /*  z5 = a0 mod n   */

    fmpz_add(t[1], t[56], t[52]);           /*  a1 = d6 + d2    */
    fmpz_add(t[2], t[1], t[10]);            /*  a2 = a1 + c1    */
    fmpz_sub(t[0], t[2], t[63]);            /*  a0 = a2 - a13   */
    unity_zp_coeff_set_fmpz(f, 6, t[0]);    /*  z6 = a0 mod n   */

    fmpz_add(t[1], t[53], t[11]);           /*  a1 = d3 + c2    */
    unity_zp_coeff_set_fmpz(f, 7, t[1]);    /*  z7 = a1 mod n   */
}

