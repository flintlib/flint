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
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

void
nmod_poly_randtest(nmod_poly_t poly, flint_rand_t state, len_t len)
{
    nmod_poly_fit_length(poly, len);
    _nmod_vec_randtest(poly->coeffs, state, len, poly->mod);
    poly->length = len;
    _nmod_poly_normalise(poly);
}

void
nmod_poly_randtest_irreducible(nmod_poly_t poly, flint_rand_t state, len_t len)
{
    do {
        nmod_poly_randtest(poly, state, len);
    } while (nmod_poly_is_zero(poly) || !(nmod_poly_is_irreducible(poly)));
}
