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

    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include "nmod_poly.h"

void
nmod_poly_factor_insert(nmod_poly_factor_t fac,
                        const nmod_poly_t poly, len_t exp)
{
    len_t i;

    if (poly->length <= 1)
        return;

    for (i = 0; i < fac->num; i++)
    {
        if (nmod_poly_equal(poly, fac->p + i))
        {
            fac->exp[i] += exp;
            return;
        }
    }

    if (fac->alloc == fac->num)
    {
        len_t new_size = 2 * fac->alloc;

        fac->p   = flint_realloc(fac->p, sizeof(nmod_poly_struct) * new_size);
        fac->exp = flint_realloc(fac->exp, sizeof(len_t) * new_size);

        for (i = fac->alloc; i < new_size; i++)
            nmod_poly_init_preinv(fac->p + i, 0, 0);

        fac->alloc = new_size;
    }

    nmod_poly_set(fac->p + (fac->num), poly);
    (fac->p + (fac->num))->mod = poly->mod;
    fac->exp[fac->num] = exp;
    fac->num++;
}
