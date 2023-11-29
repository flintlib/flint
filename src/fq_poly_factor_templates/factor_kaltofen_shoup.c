/*
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

#include <math.h>

void
TEMPLATE(T, poly_factor_kaltofen_shoup) (TEMPLATE(T, poly_factor_t) res,
                                         const TEMPLATE(T, poly_t) poly,
                                         const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, poly_t) v;
    TEMPLATE(T, poly_factor_t) sq_free, dist_deg;
    slong i, j, k, l, res_num, dist_deg_num;
    slong *degs;

    degs = flint_malloc(TEMPLATE(T, poly_degree) (poly, ctx) * sizeof(slong));

    TEMPLATE(T, poly_init) (v, ctx);

    TEMPLATE(T, poly_make_monic) (v, poly, ctx);

    /* compute squarefree factorisation */
    TEMPLATE(T, poly_factor_init) (sq_free, ctx);
    TEMPLATE(T, poly_factor_squarefree) (sq_free, v, ctx);

    /* compute distinct-degree factorisation */
    TEMPLATE(T, poly_factor_init) (dist_deg, ctx);
    for (i = 0; i < sq_free->num; i++)
    {
        dist_deg_num = dist_deg->num;

        TEMPLATE(T, poly_factor_distinct_deg) (dist_deg, sq_free->poly + i,
                                               &degs, ctx);

        /* compute equal-degree factorisation */
        for (j = dist_deg_num, l = 0; j < dist_deg->num; j++, l++)
        {
            res_num = res->num;

            TEMPLATE(T, poly_factor_equal_deg) (res, dist_deg->poly + j,
                                                degs[l], ctx);
            for (k = res_num; k < res->num; k++)
                res->exp[k] = TEMPLATE(T, poly_remove) (v, res->poly + k, ctx);
        }
    }

    flint_free(degs);
    TEMPLATE(T, poly_clear) (v, ctx);
    TEMPLATE(T, poly_factor_clear) (dist_deg, ctx);
    TEMPLATE(T, poly_factor_clear) (sq_free, ctx);
}


#endif
