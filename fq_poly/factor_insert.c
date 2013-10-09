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

    Copyright (C) 2013 Mike Hansen

******************************************************************************/
#include "fq_poly.h"

void
fq_poly_factor_insert(fq_poly_factor_t fac, const fq_poly_t poly, slong exp,
                      const fq_ctx_t ctx)
{
    slong i;

    if (poly->length <= 1)
        return;

    for (i = 0; i < fac->num; i++)
    {
        if (fq_poly_equal(poly, fac->poly + i))
        {
            fac->exp[i] += exp;
            return;
        }
    }

    if (fac->alloc == fac->num)
    {
        slong new_size = 2 * fac->alloc;

        fac->poly =
            flint_realloc(fac->poly, sizeof(fq_poly_struct) * new_size);
        fac->exp = flint_realloc(fac->exp, sizeof(slong) * new_size);

        for (i = fac->alloc; i < new_size; i++)
            fq_poly_init(fac->poly + i);

        fac->alloc = new_size;
    }

    fq_poly_set(fac->poly + fac->num, poly, ctx);
    fac->exp[fac->num] = exp;
    fac->num++;
}
