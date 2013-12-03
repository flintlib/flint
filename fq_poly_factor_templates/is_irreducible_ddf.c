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

    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Martin Lee
    Copyright (C) 2013 Mike Hansen

******************************************************************************/


#ifdef T

#include "templates.h"

#include "math.h"
int
TEMPLATE(T, poly_is_irreducible_ddf) (const TEMPLATE(T, poly_t) f,
                                      const TEMPLATE(T, ctx_t) ctx)
{
    slong i, n;
    slong *degs;
    TEMPLATE(T, poly_factor_t) dist_deg;

    n = TEMPLATE(T, poly_degree) (f, ctx);

    if (n < 2)
        return 1;

    if (!TEMPLATE(T, poly_is_squarefree) (f, ctx))
        return 0;

    if (!(degs = (slong *) flint_malloc(n * sizeof(slong))))
    {
        TEMPLATE_PRINTF("Exception (%s_poly_is_irreducible_ddf): \n", T);
        flint_printf("Not enough memory.\n");
        abort();
    }

    TEMPLATE(T, poly_factor_init) (dist_deg, ctx);
    TEMPLATE(T, poly_factor_distinct_deg) (dist_deg, f, &degs, ctx);
    for (i = 0; i < dist_deg->num; i++)
    {
        if (degs[i] == n)
        {
            flint_free(degs);
            TEMPLATE(T, poly_factor_clear) (dist_deg, ctx);
            return 1;
        }
        if (degs[i] > 0)
        {
            flint_free(degs);
            TEMPLATE(T, poly_factor_clear) (dist_deg, ctx);
            return 0;
        }
    }

    flint_free(degs);
    TEMPLATE(T, poly_factor_clear) (dist_deg, ctx);

    return 1;
}


#endif
