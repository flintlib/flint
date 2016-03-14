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

void
unity_zp_mul(unity_zp f, const unity_zp g, const unity_zp h)
{
    slong glen, hlen;

    glen = g->poly->length;
    hlen = h->poly->length;

    if (glen == 0 || hlen == 0)
    {
        fmpz_mod_poly_zero(f->poly);
        return;
    }

    fmpz_mod_poly_fit_length(f->poly, glen + hlen - 1);

    if (glen >= hlen)
        _fmpz_poly_mul(f->poly->coeffs, g->poly->coeffs, glen, h->poly->coeffs, hlen);
    else
        _fmpz_poly_mul(f->poly->coeffs, h->poly->coeffs, hlen, g->poly->coeffs, glen);

    _fmpz_mod_poly_set_length(f->poly, glen + hlen - 1);

    _unity_zp_reduce_cyclotomic_divmod(f);
}

void
unity_zp_mul_inplace(unity_zp f, const unity_zp g, const unity_zp h, fmpz_t * t)
{
    /* multiplication for p^k = 4 */
    if (f->p == 2 && f->exp == 2)
    {
        unity_zp_mul4(f, g, h, t);
        return;
    }

    /* multiplication for p^k = 8 */
    if (f->p == 2 && f->exp == 3)
    {
        unity_zp_mul8(f, g, h, t);
        return;
    }

    /* multiplication for p^k = 16 */
    if (f->p == 2 && f->exp == 4)
    {
        unity_zp_mul16(f, g, h, t);
        return;
    }

    /* multiplication for p^k = 3 */
    if (f->p == 3 && f->exp == 1)
    {
        unity_zp_mul3(f, g, h, t);
        return;
    }

    /* multiplicatiom for p^k = 9 */
    if (f->p == 3 && f->exp == 2)
    {
        unity_zp_mul9(f, g, h, t);
        return;
    }

    /* multiplication for p^k = 5 */
    if (f->p == 5 && f->exp == 1)
    {
        unity_zp_mul5(f, g, h, t);
        return;
    }

    /* multiplication for p^k = 7 */
    if (f->p == 7 && f->exp == 1)
    {
        unity_zp_mul7(f, g, h, t);
        return;
    }

    /* multiplication for p^k = 11 */
    if (f->p == 11 && f->exp == 1)
    {
        unity_zp_mul11(f, g, h, t);
        return;
    }

    /* traditional multiplication */
    unity_zp_mul(f, g, h);      
}


