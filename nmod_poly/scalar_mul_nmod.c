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

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

void
nmod_poly_scalar_mul_nmod(nmod_poly_t res, const nmod_poly_t poly1, mp_limb_t c)
{
    if ((poly1->length == 0) || (c == 0))
    {
        nmod_poly_zero(res);
        return;
    }

    nmod_poly_fit_length(res, poly1->length);

    _nmod_vec_scalar_mul_nmod(res->coeffs, poly1->coeffs, poly1->length,
                         c, poly1->mod);

    res->length = poly1->length;
    _nmod_poly_normalise(res);  /* there may have been cancellation */
}
