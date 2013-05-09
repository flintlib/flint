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

void _nmod_poly_shift_left(mp_ptr res, mp_srcptr poly, len_t len, len_t k)
{
    flint_mpn_copyd(res + k, poly, len);
    flint_mpn_zero(res, k);
}

void nmod_poly_shift_left(nmod_poly_t res, const nmod_poly_t poly, len_t k)
{
    nmod_poly_fit_length(res, poly->length + k);
   
    _nmod_poly_shift_left(res->coeffs, poly->coeffs, poly->length, k);
   
    res->length = poly->length + k;
}

