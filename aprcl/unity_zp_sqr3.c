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
    fmpz_mod_poly_get_coeff_fmpz(t[0], g->poly, 0);
    fmpz_mod_poly_get_coeff_fmpz(t[1], g->poly, 1);

    fmpz_sub(t[2], t[0], t[1]);             /*  m1 = x0 - x1    */
    fmpz_add(t[3], t[0], t[1]);             /*  m2 = x0 + x1    */
    fmpz_mul(t[4], t[2], t[3]);             /*  d1 = m1 * m2    */
    fmpz_add(t[3], t[2], t[0]);             /*  m2 = m1 + m0    */
    unity_zp_coeff_set_fmpz(f, 0, t[4]);    /*  y0 = d1 mod n   */
    fmpz_mul(t[4], t[1], t[3]);             /*  d1 = x1 * m2    */
    unity_zp_coeff_set_fmpz(f, 1, t[4]);    /*  y1 = d1 mod n   */
}

