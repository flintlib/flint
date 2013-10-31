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

#include "math.h"
#include "fq_poly.h"

int
fq_poly_is_irreducible_ddf(const fq_poly_t f, const fq_ctx_t ctx)
{
    slong i, n;
    slong *degs;
    fq_poly_factor_t dist_deg;

    n = fq_poly_degree(f, ctx);

    if (n < 2)
        return 1;

    if (!fq_poly_is_squarefree(f, ctx))
        return 0;

    if (!(degs = (slong *) flint_malloc(n * sizeof(slong))))
    {
        flint_printf("Exception (fq_poly_is_irreducible_ddf): \n");
        flint_printf("Not enough memory.\n");
        abort();
    }

    fq_poly_factor_init(dist_deg, ctx);
    fq_poly_factor_distinct_deg(dist_deg, f, &degs, ctx);
    for (i = 0; i < dist_deg->num; i++)
    {
        if (degs[i] == n)
        {
            flint_free(degs);
            fq_poly_factor_clear(dist_deg, ctx);
            return 1;
        }
        if (degs[i] > 0)
        {
            flint_free(degs);
            fq_poly_factor_clear(dist_deg, ctx);
            return 0;
        }
    }

    flint_free(degs);
    fq_poly_factor_clear(dist_deg, ctx);

    return 1;
}
