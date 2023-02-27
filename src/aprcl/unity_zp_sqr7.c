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
    Computes f = g * g for p = 7. 
    g must be reduced by F_7 cyclotomic polynomial.
    t is the memory for fmpz_t; size of t must be > 31.
    Resulting f reduced by F_7 cyclotomic polynomial.
*/
void
unity_zp_sqr7(unity_zp f, const unity_zp g, fmpz_t * t)
{
    /*
        g = (x0, ... , x5);
        f = (y0, ... , y5);

        x0 = t[0]; x1 = t[1]; x2 = t[2];
        x3 = t[3]; x4 = t[4]; x5 = t[5];

        m1 = t[6]; m2 = t[7]; ... ; m18 = t[23];
        d1 = t[24]; ... ; d8 = t[31].
    */

    /* set xi */
    fmpz_mod_poly_get_coeff_fmpz(t[0], g->poly, 0, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[1], g->poly, 1, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[2], g->poly, 2, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[3], g->poly, 3, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[4], g->poly, 4, g->ctx);
    fmpz_mod_poly_get_coeff_fmpz(t[5], g->poly, 5, g->ctx);

    fmpz_sub(t[6], t[0], t[1]);             /*  m1 = x0 - x1    */
    fmpz_sub(t[7], t[1], t[2]);             /*  m2 = x1 - x2    */
    fmpz_sub(t[8], t[2], t[3]);             /*  m3 = x2 - x3    */
    fmpz_sub(t[9], t[3], t[4]);             /*  m4 = x3 - x4    */
    fmpz_sub(t[10], t[5], t[4]);            /*  m5 = x5 - x4    */
    fmpz_add(t[11], t[6], t[7]);            /*  m6 = x6 + x7    */

    fmpz_add(t[12], t[7], t[8]);            /*  m7 = m2 + m3    */
    fmpz_add(t[13], t[8], t[9]);            /*  m8 = m3 + m4    */
    fmpz_sub(t[14], t[3], t[5]);            /*  m9 = x3 - x5    */
    fmpz_add(t[15], t[8], t[11]);           /*  m10 = m3 + m6   */
    fmpz_add(t[16], t[9], t[12]);           /*  m11 = m4 + m7   */

    fmpz_add(t[18], t[11], t[13]);          /*  m13 = m6 + m8   */
    fmpz_add(t[19], t[12], t[14]);          /*  m14 = m7 + m9   */
    fmpz_add(t[21], t[0], t[1]);            /*  m16 = x0 + x1   */
    fmpz_add(t[22], t[0], t[15]);           /*  m17 = m8 + m2   */
    fmpz_mul(t[24], t[3], t[22]);           /*  d1 = x3 * m17   */
    
    fmpz_sub(t[22], t[19], t[4]);           /*  m17 = m14 - m9  */
    fmpz_add(t[23], t[19], t[4]);           /*  m18 = m14 + m9  */ 
    fmpz_mul(t[25], t[22], t[23]);          /*  d2 = m17 * m18  */
    fmpz_sub(t[22], t[13], t[7]);           /*  m17 = m8 - m2   */
    fmpz_mul(t[26], t[16], t[22]);          /*  d3 = m11 * m17  */

    fmpz_add(t[22], t[19], t[14]);          /*  m17 = m14 + m9  */
    fmpz_mul(t[27], t[22], t[12]);          /*  d4 = m17 * m7   */
    fmpz_add(t[22], t[1], t[1]);            /*  m17 = x1 + x1   */
    fmpz_mul(t[28], t[22], t[11]);          /*  d5 = m17 * m6   */
    fmpz_mul(t[29], t[6], t[21]);           /*  d6 = m1 * m16   */
    fmpz_add(t[22], t[8], t[8]);            /*  m17 = m3 + m3   */

    fmpz_add(t[7], t[0], t[18]);            /*  m2 = x0 + m13   */
    fmpz_mul(t[30], t[22], t[10]);          /*  d7 = m17 * m5   */
    fmpz_add(t[31], t[24], t[25]);          /*  d8 = d1 + d2    */
    fmpz_add(t[24], t[31], t[26]);          /*  d1 = d8 + d5    */ 
    unity_zp_coeff_set_fmpz(f, 3, t[24]);   /*  y3 = d1 mod n   */
    fmpz_add(t[31], t[26], t[27]);          /*  d8 = d3 + d4    */

    fmpz_add(t[24], t[31], t[28]);          /*  d1 = d8 + d5    */
    unity_zp_coeff_set_fmpz(f, 1, t[24]);   /*  y1 = d1 mod n   */
    fmpz_add(t[31], t[27], t[29]);          /*  d8 = d4 + d6    */
    fmpz_add(t[24], t[31], t[30]);          /*  d1 = d8 + d7    */
    unity_zp_coeff_set_fmpz(f, 0, t[24]);   /*  y0 = d1 mod n   */
    fmpz_add(t[22], t[12], t[19]);          /*  m17 = m7 + m14  */

    fmpz_mul(t[24], t[14], t[22]);          /*  d1 = m9 * m17   */
    fmpz_sub(t[22], t[13], t[5]);           /*  m17 = m8 - x5   */
    fmpz_add(t[23], t[2], t[10]);           /*  m18 = x2 + m5   */
    fmpz_mul(t[25], t[22], t[23]);          /*  d2 = m17 * m18  */
    fmpz_mul(t[26], t[7], t[4]);            /*  d3 = m2 * x4    */
    fmpz_add(t[22], t[8], t[13]);           /*  m17 = m3 + m8   */

    fmpz_mul(t[27], t[22], t[9]);           /*  d4 = m17 * m4   */
    fmpz_add(t[22], t[6], t[6]);            /*  m17 = m1 + m1   */
    fmpz_mul(t[28], t[22], t[10]);          /*  d5 = m17 * m5   */
    fmpz_sub(t[22], t[19], t[10]);          /*  m17 = m14 - m5  */
    fmpz_mul(t[29], t[22], t[16]);          /*  d6 = m17 * m11  */

    fmpz_add(t[22], t[2], t[2]);            /*  m17 = x2 + x2   */
    fmpz_mul(t[30], t[22], t[15]);          /*  d7 = m17 * m10  */
    fmpz_add(t[31], t[24], t[25]);          /*  d8 = d1 + d2    */
    fmpz_add(t[24], t[31], t[26]);          /*  d1 = d8 + d3    */
    unity_zp_coeff_set_fmpz(f, 4, t[24]);   /*  y4 = d1 mod n   */
    fmpz_add(t[31], t[26], t[27]);          /*  d8 = d3 + d4    */

    fmpz_add(t[24], t[31], t[28]);          /*  d1 = d8 + d5    */
    unity_zp_coeff_set_fmpz(f, 5, t[24]);   /*  y5 = d1 mod n   */
    fmpz_add(t[31], t[27], t[29]);          /*  d8 = d4 + d6    */
    fmpz_add(t[24], t[31], t[30]);          /*  d1 = d8 + d7    */ 
    unity_zp_coeff_set_fmpz(f, 2, t[24]);   /*  y2 = d1 mod n   */
}

