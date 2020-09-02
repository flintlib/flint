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
    Computes f = g * h for p = 7. 
    g and h must be reduced by F_7 cyclotomic polynomial.
    t is the memory for fmpz_t; size of t must be > 68.
    Resulting f reduced by F_7 cyclotomic polynomial.
*/
void
unity_zp_mul7(unity_zp f, const unity_zp g, const unity_zp h, fmpz_t * t)
{
    /*
        g = (x0, ... , x5);
        h = (y0, ... , y5);
        f = (z0, ... , z5);

        x0 = t[30]; x1 = t[31]; x2 = t[32]; x3 = t[33];
        x4 = t[34]; x5 = t[35];

        y0 = t[40]; y1 = t[41]; y2 = t[42]; y3 = t[43];
        y4 = t[44]; y5 = t[45];

        a0 = t[0]; a1 = t[1] ; a2 = t[2];
        b0 = t[3]; b1 = t[4] ; b2 = t[5];
        c0 = t[6]; ... ; c4 = t[10];

        d0 = t[50]; ... ; d18 = t[68];
    */

    /* set xi */
    fmpz_mod_poly_get_coeff_fmpz(t[30], g->poly, 0, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[31], g->poly, 1, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[32], g->poly, 2, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[33], g->poly, 3, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[34], g->poly, 4, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[35], g->poly, 5, g->ctx);

    /* set yi */
    fmpz_mod_poly_get_coeff_fmpz(t[40], h->poly, 0, h->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[41], h->poly, 1, h->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[42], h->poly, 2, h->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[43], h->poly, 3, h->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[44], h->poly, 4, h->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[45], h->poly, 5, h->ctx);

    fmpz_set(t[0], t[30]);                  /*  set a0 = x0     */
    fmpz_set(t[1], t[31]);                  /*  set a1 = x1     */
    fmpz_set(t[2], t[32]);                  /*  set a2 = x2     */
    fmpz_set(t[3], t[40]);                  /*  set b0 = y0     */
    fmpz_set(t[4], t[41]);                  /*  set b1 = y1     */
    fmpz_set(t[5], t[42]);                  /*  set b2 = y2     */

    /* 
        apply auxiliary routine 2 with (a0, a1, a2) and (b0, b1, b2)
        store result in (c0, .. , c4)
    */
    unity_zp_ar1(t);

    fmpz_set(t[50], t[6]);                  /*  set d0 = c0     */
    fmpz_set(t[51], t[7]);                  /*  set d1 = c1     */
    fmpz_set(t[52], t[8]);                  /*  set d2 = c2     */
    fmpz_set(t[53], t[9]);                  /*  set d3 = c3     */
    fmpz_set(t[54], t[10]);                 /*  set d4 = c4     */

    fmpz_set(t[0], t[33]);                  /*  set a0 = x3     */
    fmpz_set(t[1], t[34]);                  /*  set a1 = x4     */
    fmpz_set(t[2], t[35]);                  /*  set a2 = x5     */
    fmpz_set(t[3], t[43]);                  /*  set b0 = y3     */
    fmpz_set(t[4], t[44]);                  /*  set b1 = y4     */
    fmpz_set(t[5], t[45]);                  /*  set b2 = y5     */

    /* 
        apply auxiliary routine 2 with (a0, a1, a2) and (b0, b1, b2)
        store result in (c0, .. , c4)
    */
    unity_zp_ar1(t);

    fmpz_set(t[56], t[6]);                  /*  set d6 = c0     */
    fmpz_set(t[57], t[7]);                  /*  set d7 = c1     */
    fmpz_set(t[58], t[8]);                  /*  set d8 = c2     */
    fmpz_set(t[59], t[9]);                  /*  set d9 = c3     */
    fmpz_set(t[60], t[10]);                 /*  set d10 = c4    */

    fmpz_sub(t[0], t[30], t[33]);           /*  a0 = x0 - x3    */
    fmpz_sub(t[1], t[31], t[34]);           /*  a1 = x1 - x4    */
    fmpz_sub(t[2], t[32], t[35]);           /*  a2 = x2 - x5    */
    fmpz_sub(t[3], t[43], t[40]);           /*  b0 = y3 - y0    */
    fmpz_sub(t[4], t[44], t[41]);           /*  b1 = y4 - y1    */
    fmpz_sub(t[5], t[45], t[42]);           /*  b2 = y5 - y2    */

    /* 
        apply auxiliary routine 2 with (a0, a1, a2) and (b0, b1, b2)
        store result in (c0, .. , c4)
    */
    unity_zp_ar1(t);

    fmpz_set(t[61], t[6]);                  /*  set d11 = c0    */
    fmpz_set(t[62], t[7]);                  /*  set d12 = c1    */
    fmpz_set(t[63], t[8]);                  /*  set d13 = c2    */
    fmpz_set(t[64], t[9]);                  /*  set d14 = c3    */
    fmpz_set(t[65], t[10]);                 /*  set d15 = c4    */

    fmpz_add(t[68], t[56], t[64]);          /*  d18 = d6 + d14  */
    fmpz_add(t[66], t[68], t[53]);          /*  d16 = d18 + d3  */
    fmpz_add(t[68], t[57], t[65]);          /*  d18 = d7 + d15  */
    fmpz_add(t[67], t[68], t[54]);          /*  d17 = d18 + d4  */
    fmpz_add(t[68], t[53], t[61]);          /*  d18 = d3 + d11  */
    fmpz_add(t[53], t[68], t[50]);          /*  d3 = d18 + d0   */
    fmpz_add(t[68], t[54], t[62]);          /*  d18 = d4 + d12  */
    fmpz_add(t[54], t[68], t[51]);          /*  d4 = d18 + d1   */

    fmpz_add(t[55], t[52], t[63]);          /*  d5 = d2 + d13   */
    fmpz_add(t[63], t[53], t[56]);          /*  d13 = d3 + d6   */
    fmpz_add(t[64], t[54], t[57]);          /*  d14 = d4 + d7   */
    fmpz_add(t[65], t[55], t[58]);          /*  d15 = d5 + d8   */
    fmpz_add(t[56], t[66], t[59]);          /*  d6 = d16 + d9   */
    fmpz_add(t[57], t[67], t[60]);          /*  d7 = d17 + d10  */

    fmpz_add(t[68], t[50], t[57]);          /*  d18 = d10 + d7  */
    fmpz_sub(t[0], t[68], t[56]);           /*  a0 = d18 - d6   */
    unity_zp_coeff_set_fmpz(f, 0, t[0]);    /*  z0 = a0 mod n   */

    fmpz_add(t[68], t[51], t[58]);          /*  d18 = d1 + d8   */
    fmpz_sub(t[0], t[68], t[56]);           /*  a0 = d18 - d6   */
    unity_zp_coeff_set_fmpz(f, 1, t[0]);    /*  z1 = a0 mod n   */

    fmpz_add(t[68], t[52], t[59]);          /*  d18 = d2 + d9   */
    fmpz_sub(t[0], t[68], t[56]);           /*  a0 = d18 - d6   */
    unity_zp_coeff_set_fmpz(f, 2, t[0]);    /*  z2 = a0 mod n   */

    fmpz_add(t[68], t[63], t[60]);          /*  d18 = d13 + d10 */
    fmpz_sub(t[0], t[68], t[56]);           /*  a0 = d18 - d6   */
    unity_zp_coeff_set_fmpz(f, 3, t[0]);    /*  z3 = a0 mod n   */

    fmpz_sub(t[0], t[64], t[56]);           /*  a0 = d14 - d6   */
    unity_zp_coeff_set_fmpz(f, 4, t[0]);    /*  z4 = a0 mod n   */

    fmpz_sub(t[0], t[65], t[56]);           /*  a0 = d15 - d6   */
    unity_zp_coeff_set_fmpz(f, 5, t[0]);    /*  z5 = a0 mod n   */
}

