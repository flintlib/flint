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
    Computes f = g * h for p = 3. 
    g and h must be reduced by F_3 cyclotomic polynomial.
    t is the memory for fmpz_t; size of t must be > 3.
    Resulting f reduced by F_3 cyclotomic polynomial.
*/
void
unity_zp_mul3(unity_zp f, const unity_zp g, const unity_zp h, fmpz_t * t)
{
    /*
        g = (x0, x1);
        h = (y0, y1);
        f = (z0, z1);

        x0 = t[0]; x1 = t[1];
        y0 = t[2]; y1 = t[3];
        m1 = t[4]; m2 = t[5];
        d1 = t[6]; d2 = t[7]; d3 = t[8].
    */

    /* set xi */
    fmpz_mod_poly_get_coeff_fmpz(t[0], g->poly, 0, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[1], g->poly, 1, g->ctx);

    /* set yi */
    fmpz_mod_poly_get_coeff_fmpz(t[2], h->poly, 0, h->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[3], h->poly, 1, h->ctx);

    fmpz_mul(t[6], t[0], t[2]);             /*  d1 = x0 * y0    */
    fmpz_mul(t[7], t[1], t[3]);             /*  d2 = x1 * y1    */
    fmpz_sub(t[4], t[0], t[1]);             /*  m1 = x0 - x1    */
    fmpz_sub(t[5], t[3], t[2]);             /*  m2 = y1 - y0    */
    fmpz_mul(t[8], t[4], t[5]);             /*  d3 = m1 * m2    */
    fmpz_add(t[8], t[8], t[6]);             /*  d3 = d3 + d1    */
    
    unity_zp_coeff_set_fmpz(f, 1, t[8]);    /*  z1 = d3 mod n   */
    fmpz_sub(t[0], t[6], t[7]);             /*  x0 = d1 - d2    */
    unity_zp_coeff_set_fmpz(f, 0, t[0]);    /*  z0 = x0 mod n   */
}

/*
    Computes f = g * h for p = 3^2. 
    g and h must be reduced by F_9 cyclotomic polynomial.
    t is the memory for fmpz_t; size of t must be > 50.
    Resulting f reduced by F_9 cyclotomic polynomial.
*/
void
unity_zp_mul9(unity_zp f, const unity_zp g, const unity_zp h, fmpz_t * t)
{
    /*
        g = (x0, ... , x5);
        h = (y0, ... , y5);
        f = (z0, ... , z5);

        x0 = t[20]; ... ; x5 = t[25];
        y0 = t[26]; ... ; y5 = t[31];

        a0 = t[0]; a1 = t[1]; a2 = t[2];
        b0 = t[3]; b1 = t[4]; b2 = t[5];
        c0 = t[6]; ... ; c4 = t[10];

        d0 = t[32]; ... ; d19 = t[51].
    */

    /* set xi */
    fmpz_mod_poly_get_coeff_fmpz(t[20], g->poly, 0, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[21], g->poly, 1, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[22], g->poly, 2, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[23], g->poly, 3, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[24], g->poly, 4, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[25], g->poly, 5, g->ctx);

    /* set yi */
    fmpz_mod_poly_get_coeff_fmpz(t[26], h->poly, 0, h->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[27], h->poly, 1, h->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[28], h->poly, 2, h->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[29], h->poly, 3, h->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[30], h->poly, 4, h->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[31], h->poly, 5, h->ctx);

    fmpz_set(t[0], t[20]);                  /*  set a0 = x0     */
    fmpz_set(t[1], t[21]);                  /*  set a1 = x1     */
    fmpz_set(t[2], t[22]);                  /*  set a2 = x2     */
    fmpz_set(t[3], t[26]);                  /*  set b0 = y0     */
    fmpz_set(t[4], t[27]);                  /*  set b1 = y1     */
    fmpz_set(t[5], t[28]);                  /*  set b2 = y2     */

    /* 
        apply auxiliary routine 2 with (a0, a1, a2) and (b0, b1, b2)
        store result in (c0, .. , c4)
    */
    unity_zp_ar1(t);

    fmpz_set(t[32], t[6]);                  /*  set d0 = c0     */
    fmpz_set(t[33], t[7]);                  /*  set d1 = c1     */
    fmpz_set(t[34], t[8]);                  /*  set d2 = c2     */
    fmpz_set(t[35], t[9]);                  /*  set d3 = c3     */
    fmpz_set(t[36], t[10]);                 /*  set d4 = c4     */

    fmpz_set(t[0], t[23]);                  /*  set a0 = x3     */
    fmpz_set(t[1], t[24]);                  /*  set a1 = x4     */
    fmpz_set(t[2], t[25]);                  /*  set a2 = x5     */
    fmpz_set(t[3], t[29]);                  /*  set b0 = y3     */
    fmpz_set(t[4], t[30]);                  /*  set b1 = y4     */
    fmpz_set(t[5], t[31]);                  /*  set b2 = y5     */

    /* 
        apply auxiliary routine 2 with (a0, a1, a2) and (b0, b1, b2)
        store result in (c0, .. , c4)
    */
    unity_zp_ar1(t);

    fmpz_set(t[38], t[6]);                  /*  set d6 = c0     */
    fmpz_set(t[39], t[7]);                  /*  set d7 = c1     */
    fmpz_set(t[40], t[8]);                  /*  set d8 = c2     */
    fmpz_set(t[41], t[9]);                  /*  set d9 = c3     */
    fmpz_set(t[42], t[10]);                 /*  set d10 = c4    */

    fmpz_sub(t[0], t[20], t[23]);           /*  a0 = x0 - x3    */
    fmpz_sub(t[1], t[21], t[24]);           /*  a1 = x1 - x4    */
    fmpz_sub(t[2], t[22], t[25]);           /*  a2 = x2 - x5    */
    fmpz_sub(t[3], t[29], t[26]);           /*  b0 = y3 - y0    */
    fmpz_sub(t[4], t[30], t[27]);           /*  b1 = y4 - y1    */
    fmpz_sub(t[5], t[31], t[28]);           /*  b2 = y5 - y2    */

    /* 
        apply auxiliary routine 2 with (a0, a1, a2) and (b0, b1, b2)
        store result in (c0, .. , c4)
    */
    unity_zp_ar1(t);

    fmpz_set(t[43], t[6]);                  /*  set d11 = c0    */
    fmpz_set(t[44], t[7]);                  /*  set d12 = c1    */
    fmpz_set(t[45], t[8]);                  /*  set d13 = c2    */
    fmpz_set(t[46], t[9]);                  /*  set d14 = c3    */
    fmpz_set(t[47], t[10]);                 /*  set d15 = c4    */

    fmpz_add(t[50], t[38], t[46]);          /*  d18 = d6 + d14  */
    fmpz_add(t[48], t[50], t[35]);          /*  d16 = d18 + d3  */
    fmpz_add(t[50], t[39], t[47]);          /*  d18 = d7 + d15  */
    fmpz_add(t[49], t[50], t[36]);          /*  d17 = d18 + d4  */
    fmpz_add(t[50], t[35], t[43]);          /*  d18 = d3 + d11  */
    fmpz_add(t[35], t[50], t[32]);          /*  d3 = d18 + d0   */
    fmpz_add(t[50], t[36], t[44]);          /*  d18 = d4 + d12  */

    fmpz_add(t[36], t[50], t[33]);          /*  d4 = d18 + d1   */
    fmpz_add(t[37], t[34], t[45]);          /*  d5 = d2 + d13   */
    fmpz_sub(t[0], t[32], t[48]);           /*  a0 = d0 - d16   */
    unity_zp_coeff_set_fmpz(f, 0, t[0]);    /*  z0 = a0 mod n   */

    fmpz_sub(t[0], t[33], t[49]);           /*  a0 = d1 - d17   */
    unity_zp_coeff_set_fmpz(f, 1, t[0]);    /*  z1 = a0 mod n   */

    fmpz_sub(t[0], t[34], t[40]);           /*  a0 = d2 - d8    */
    unity_zp_coeff_set_fmpz(f, 2, t[0]);    /*  z2 = a0 mod n   */

    unity_zp_coeff_set_fmpz(f, 5, t[37]);   /*  z5 = d5 mod n   */

    fmpz_add(t[50], t[35], t[38]);          /*  d18 = d3 + d6   */
    fmpz_add(t[51], t[48], t[41]);          /*  d19 = d16 + d9  */
    fmpz_sub(t[0], t[50], t[51]);           /*  a0 = d18 - d19  */
    unity_zp_coeff_set_fmpz(f, 3, t[0]);    /*  z3 = a0 mod n   */

    fmpz_add(t[50], t[36], t[39]);          /*  d18 = d4 + d7   */
    fmpz_add(t[51], t[42], t[49]);          /*  d19 = d10 + d17 */
    fmpz_sub(t[0], t[50], t[51]);           /*  a0 = d18 - d19  */
    unity_zp_coeff_set_fmpz(f, 4, t[0]);    /*  z4 = a0 mod n   */
}

