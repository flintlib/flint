/*============================================================================

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

    Copyright (C) 2011 Sebastian Pancratz
 
******************************************************************************/

#include "padic_poly.h"

/*
    Brings f in the form p^v * poly, no reduction.
 */
void _padic_poly_canonicalise(fmpz *poly, long *v, long len, const fmpz_t p)
{
    const long min = _fmpz_vec_ord_p(poly, len, p);

    if (min == 0)
    {
        if (_fmpz_vec_is_zero(poly, len))
            *v = 0;
    }
    else  /* min > 0 */
    {
        fmpz_t pow;

        fmpz_init(pow);
        fmpz_pow_ui(pow, p, min);
        _fmpz_vec_scalar_divexact_fmpz(poly, poly, len, pow);
        fmpz_clear(pow);

        *v += min;
    }
}

void padic_poly_canonicalise(padic_poly_t poly, const fmpz_t p)
{
    _padic_poly_canonicalise(poly->coeffs, &(poly->val), poly->length, p);
}

