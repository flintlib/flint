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
    fmpz_mod_poly_get_coeff_fmpz(t[0], g->poly, 0); /* set x0        */
    fmpz_mod_poly_get_coeff_fmpz(t[1], g->poly, 1); /* set x1        */
    fmpz_mod_poly_get_coeff_fmpz(t[2], h->poly, 0); /* set y0        */
    fmpz_mod_poly_get_coeff_fmpz(t[3], h->poly, 1); /* set y1        */

    fmpz_mul(t[6], t[0], t[2]);                     /* d1 = x0 * y0  */
    fmpz_mul(t[7], t[1], t[3]);                     /* d2 = x1 * y1  */
    fmpz_sub(t[4], t[0], t[1]);                     /* m1 = x0 - x1  */
    fmpz_sub(t[5], t[3], t[2]);                     /* m2 = y1 - y0  */
    fmpz_mul(t[8], t[4], t[5]);                     /* d3 = m1 * m2  */
    fmpz_add(t[8], t[8], t[6]);                     /* d3 = d3 + d1  */
    
    unity_zp_coeff_set_fmpz(f, 1, t[8]);            /* z1 = d3 mod n */
    fmpz_sub(t[0], t[6], t[7]);                     /* x0 = d1 - d2  */
    unity_zp_coeff_set_fmpz(f, 0, t[0]);            /* z0 = x0 mod n */
}

void
unity_zp_mul9(unity_zp f, const unity_zp g, const unity_zp h, fmpz_t * t)
{

}

