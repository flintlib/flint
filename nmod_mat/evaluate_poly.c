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
#include "flint.h"
#include "nmod_mat.h"

/*Function to evaluate a polynomial with Matrix as an argument. It assumes that the dest and c matrix has same
    dimensions. dest and c cannot be aliased*/
void
_nmod_mat_evaluate_poly(nmod_mat_t dest,const mp_srcptr poly, slong len, const nmod_mat_t c, nmod_t mod)
{
    if (len == 0)
    {
        nmod_mat_zero(dest);
        return;
    }

    slong m = len-1;
    nmod_mat_t IDENTITY;
    nmod_mat_init_set(IDENTITY, c);
    nmod_mat_one(IDENTITY);
    nmod_mat_scalar_mul(dest, IDENTITY, poly[m]);

    if (len == 1 || nmod_mat_is_zero(c))
    {
        nmod_mat_clear(IDENTITY);
        return;
    }

    nmod_mat_print_pretty(dest);
    nmod_mat_t temp;
    nmod_mat_init_set(temp, c);

    for( m-- ; m >= 0 ; m--)
    {
        nmod_mat_scalar_mul(temp, IDENTITY, poly[m]);
        nmod_mat_addmul(temp, temp, dest, c);
        nmod_mat_set(dest, temp);
        nmod_mat_print_pretty(dest);
    }

    nmod_mat_clear(IDENTITY);
    nmod_mat_clear(temp);
}

void
nmod_mat_evaluate_poly(nmod_mat_t dest, const nmod_poly_t poly, const nmod_mat_t c)
{
    return _nmod_mat_evaluate_poly(dest, poly->coeffs, poly->length, c, poly->mod);
}
