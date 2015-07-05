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
    fmpz_mod_poly_get_coeff_fmpz(t[0], g->poly, 0);
    fmpz_mod_poly_get_coeff_fmpz(t[1], g->poly, 1);

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
    t is the memory for fmpz_t; size of t must be > 3.
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
    fmpz_mod_poly_get_coeff_fmpz(t[0], g->poly, 0);
    fmpz_mod_poly_get_coeff_fmpz(t[1], g->poly, 1);
    fmpz_mod_poly_get_coeff_fmpz(t[2], g->poly, 2);
    fmpz_mod_poly_get_coeff_fmpz(t[3], g->poly, 3);

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

void
unity_zp_sqr16(unity_zp f, const unity_zp g, fmpz_t * t)
{
    fmpz_mod_poly_get_coeff_fmpz(t[30], g->poly, 0);
    fmpz_mod_poly_get_coeff_fmpz(t[31], g->poly, 1);
    fmpz_mod_poly_get_coeff_fmpz(t[32], g->poly, 2);
    fmpz_mod_poly_get_coeff_fmpz(t[33], g->poly, 3);
    fmpz_mod_poly_get_coeff_fmpz(t[34], g->poly, 4);
    fmpz_mod_poly_get_coeff_fmpz(t[35], g->poly, 5);
    fmpz_mod_poly_get_coeff_fmpz(t[36], g->poly, 6);
    fmpz_mod_poly_get_coeff_fmpz(t[37], g->poly, 7);

    fmpz_add(t[0], t[30], t[34]);
    fmpz_add(t[1], t[31], t[35]);
    fmpz_add(t[2], t[32], t[36]);
    fmpz_add(t[3], t[33], t[37]);
    fmpz_sub(t[4], t[30], t[34]);
    fmpz_sub(t[5], t[31], t[35]);
    fmpz_sub(t[6], t[32], t[36]);
    fmpz_sub(t[7], t[33], t[37]);

    unity_zp_ar2(t);

    fmpz_set(t[38], t[8]);
    fmpz_set(t[39], t[9]);
    fmpz_set(t[40], t[10]);
    fmpz_set(t[41], t[11]);
    fmpz_set(t[42], t[12]);
    fmpz_set(t[43], t[13]);
    fmpz_set(t[44], t[14]);

    fmpz_add(t[0], t[30], t[30]);
    fmpz_add(t[1], t[31], t[31]);
    fmpz_add(t[2], t[32], t[32]);
    fmpz_add(t[3], t[33], t[33]);
    
    unity_zp_ar2(t);

    fmpz_sub(t[16], t[38], t[12]);
    unity_zp_coeff_set_fmpz(f, 0, t[16]);
    fmpz_sub(t[16], t[39], t[13]);
    unity_zp_coeff_set_fmpz(f, 1, t[16]);
    fmpz_sub(t[16], t[30], t[14]);
    unity_zp_coeff_set_fmpz(f, 2, t[16]);
    unity_zp_coeff_set_fmpz(f, 3, t[41]);
    fmpz_add(t[16], t[42], t[8]);
    unity_zp_coeff_set_fmpz(f, 4, t[16]);
    fmpz_add(t[16], t[43], t[9]);
    unity_zp_coeff_set_fmpz(f, 5, t[16]);
    fmpz_add(t[16], t[44], t[10]);
    unity_zp_coeff_set_fmpz(f, 6, t[16]);
    unity_zp_coeff_set_fmpz(f, 7, t[11]);
}

