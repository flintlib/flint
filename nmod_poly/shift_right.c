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

    Copyright (C) 2010 William Hart

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

void _nmod_poly_shift_right(mp_ptr res, mp_srcptr poly, long len, long k)
{
    flint_mpn_copyi(res, poly + k, len);
}

void nmod_poly_shift_right(nmod_poly_t res, const nmod_poly_t poly, long k)
{
    if (k >= poly->length) /* shift all coeffs out */
        res->length = 0;
    else
    {
        const long len = poly->length - k;
        nmod_poly_fit_length(res, len);

        _nmod_poly_shift_right(res->coeffs, poly->coeffs, len, k);

        res->length = len;
    }
}
