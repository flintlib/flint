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

    Copyright (C) 2008-2011 William Hart
    
******************************************************************************/

#include "fmpz_poly.h"
#include "fft.h"
#include "fft_tuning.h"

void _fmpz_poly_mul_SS(fmpz *output, const fmpz *input1, len_t len1, 
                       const fmpz *input2, len_t len2)
{
    const len_t rlen = len1 + len2 - 1;

    _fmpz_poly_mullow_SS(output, input1, len1, input2, len2, rlen);
}

void
fmpz_poly_mul_SS(fmpz_poly_t res,
                 const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
    const len_t len1 = poly1->length, len2 = poly2->length;
    len_t rlen;

    if (len1 == 0 || len2 == 0)
    {
        fmpz_poly_zero(res);
        return;
    }
    if (len1 == 1 || len2 == 1)
    {
        fmpz_poly_mul_classical(res, poly1, poly2);
        return;
    }

    rlen = len1 + len2 - 1;

    fmpz_poly_fit_length(res, rlen);
    if (len1 >= len2)
        _fmpz_poly_mullow_SS(res->coeffs, poly1->coeffs, len1,
                          poly2->coeffs, len2, rlen);
    else
        _fmpz_poly_mullow_SS(res->coeffs, poly2->coeffs, len2,
                          poly1->coeffs, len1, rlen);
    _fmpz_poly_set_length(res, rlen);
}
