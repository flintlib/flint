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
    fmpz_mod_poly_get_coeff_fmpz(t[0], g->poly, 0); /*  set x0          */
    fmpz_mod_poly_get_coeff_fmpz(t[1], g->poly, 1); /*  set x1          */

    fmpz_mod_poly_get_coeff_fmpz(t[2], h->poly, 0); /*  set y0          */
    fmpz_mod_poly_get_coeff_fmpz(t[3], h->poly, 1); /*  set y1          */

    fmpz_mul(t[6], t[0], t[2]);                     /*  d1 = x0 * y0    */
    fmpz_mul(t[7], t[1], t[3]);                     /*  d2 = x1 * y1    */
    fmpz_sub(t[4], t[0], t[1]);                     /*  m1 = x0 - x1    */
    fmpz_sub(t[5], t[3], t[2]);                     /*  m2 = y1 - y0    */
    fmpz_mul(t[8], t[4], t[5]);                     /*  d3 = m1 * m2    */
    fmpz_add(t[8], t[8], t[6]);                     /*  d3 = d3 + d1    */
    
    unity_zp_coeff_set_fmpz(f, 1, t[8]);            /*  z1 = d3 mod n   */
    fmpz_sub(t[0], t[6], t[7]);                     /*  x0 = d1 - d2    */
    unity_zp_coeff_set_fmpz(f, 0, t[0]);            /*  z0 = x0 mod n   */
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
    fmpz_mod_poly_get_coeff_fmpz(t[20], g->poly, 0);
    fmpz_mod_poly_get_coeff_fmpz(t[21], g->poly, 1);
    fmpz_mod_poly_get_coeff_fmpz(t[22], g->poly, 2);
    fmpz_mod_poly_get_coeff_fmpz(t[23], g->poly, 3);
    fmpz_mod_poly_get_coeff_fmpz(t[24], g->poly, 4);
    fmpz_mod_poly_get_coeff_fmpz(t[25], g->poly, 5);

    fmpz_mod_poly_get_coeff_fmpz(t[26], h->poly, 0);
    fmpz_mod_poly_get_coeff_fmpz(t[27], h->poly, 1);
    fmpz_mod_poly_get_coeff_fmpz(t[28], h->poly, 2);
    fmpz_mod_poly_get_coeff_fmpz(t[29], h->poly, 3);
    fmpz_mod_poly_get_coeff_fmpz(t[30], h->poly, 4);
    fmpz_mod_poly_get_coeff_fmpz(t[31], h->poly, 5);

    fmpz_set(t[0], t[20]);                  /*  set a0 = x0 */
    fmpz_set(t[1], t[21]);                  /*  set a1 = x1 */
    fmpz_set(t[2], t[22]);                  /*  set a2 = x2 */
    fmpz_set(t[3], t[26]);                  /*  set b0 = y0 */
    fmpz_set(t[4], t[27]);                  /*  set b1 = y1 */
    fmpz_set(t[5], t[28]);                  /*  set b2 = y2 */

    /* 
        apply auxiliary routine 2 with (a0, a1, a2) and (b0, b1, b2)
        store result in (c0, .. , c4)
    */
    unity_zp_ar1(t);

    fmpz_set(t[32], t[6]);                  /*  set d0 = c0 */
    fmpz_set(t[33], t[7]);                  /*  set d1 = c1 */
    fmpz_set(t[34], t[8]);                  /*  set d2 = c2 */
    fmpz_set(t[35], t[9]);                  /*  set d3 = c3 */
    fmpz_set(t[36], t[10]);                 /*  set d4 = c4 */

    fmpz_set(t[0], t[23]);                  /*  set a0 = x3 */
    fmpz_set(t[1], t[24]);                  /*  set a1 = x4 */
    fmpz_set(t[2], t[25]);
    fmpz_set(t[3], t[29]);
    fmpz_set(t[4], t[30]);
    fmpz_set(t[5], t[31]);
    unity_zp_ar1(t);

    fmpz_set(t[38], t[6]);
    fmpz_set(t[39], t[7]);
    fmpz_set(t[40], t[8]);
    fmpz_set(t[41], t[9]);
    fmpz_set(t[42], t[10]);

    fmpz_sub(t[0], t[20], t[23]);
    fmpz_sub(t[1], t[21], t[24]);
    fmpz_sub(t[2], t[22], t[25]);
    fmpz_sub(t[3], t[29], t[26]);
    fmpz_sub(t[4], t[30], t[27]);
    fmpz_sub(t[5], t[31], t[28]);
    unity_zp_ar1(t);

    fmpz_set(t[43], t[6]);
    fmpz_set(t[44], t[7]);
    fmpz_set(t[45], t[8]);
    fmpz_set(t[46], t[9]);
    fmpz_set(t[47], t[10]);

    fmpz_add(t[50], t[38], t[46]);
    fmpz_add(t[48], t[50], t[35]);
    fmpz_add(t[50], t[39], t[47]);
    fmpz_add(t[49], t[50], t[36]);
    fmpz_add(t[50], t[35], t[43]);
    fmpz_add(t[35], t[50], t[32]);
    fmpz_add(t[50], t[36], t[44]);

    fmpz_add(t[36], t[50], t[33]);
    fmpz_add(t[37], t[34], t[45]);
    fmpz_sub(t[0], t[32], t[48]);
    unity_zp_coeff_set_fmpz(f, 0, t[0]);
    fmpz_sub(t[0], t[33], t[49]);
    unity_zp_coeff_set_fmpz(f, 1, t[0]);
    fmpz_sub(t[0], t[34], t[40]);
    unity_zp_coeff_set_fmpz(f, 2, t[0]);

    unity_zp_coeff_set_fmpz(f, 5, t[37]);
    fmpz_add(t[50], t[35], t[38]);
    fmpz_add(t[51], t[48], t[41]);
    fmpz_sub(t[0], t[50], t[51]);
    unity_zp_coeff_set_fmpz(f, 3, t[0]);
    fmpz_add(t[50], t[36], t[39]);
    fmpz_add(t[51], t[42], t[49]);
    fmpz_sub(t[0], t[50], t[51]);
    unity_zp_coeff_set_fmpz(f, 4, t[0]);
}

