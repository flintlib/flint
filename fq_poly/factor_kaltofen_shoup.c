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
#include <math.h>

void
fq_poly_factor_kaltofen_shoup(fq_poly_factor_t res, const fq_poly_t poly,
                              const fq_ctx_t ctx)
{
    fq_poly_t v;
    fq_poly_factor_t sq_free, dist_deg;
    slong i, j, k, l, res_num, dist_deg_num;
    slong *degs;

    if (!(degs = flint_malloc(fq_poly_degree(poly) * sizeof(slong))))
    {
        flint_printf("Exception (fq_poly_factor_kaltofen_shoup): \n");
        flint_printf("Not enough memory.\n");
        abort();
    }

    fq_poly_init(v);

    fq_poly_make_monic(v, poly, ctx);

    /* compute squarefree factorisation */
    fq_poly_factor_init(sq_free, ctx);
    fq_poly_factor_squarefree(sq_free, v, ctx);

    /* compute distinct-degree factorisation */
    fq_poly_factor_init(dist_deg, ctx);
    for (i = 0; i < sq_free->num; i++)
    {
        dist_deg_num = dist_deg->num;

        fq_poly_factor_distinct_deg(dist_deg, sq_free->poly + i, &degs, ctx);

        /* compute equal-degree factorisation */
        for (j = dist_deg_num, l = 0; j < dist_deg->num; j++, l++)
        {
            res_num = res->num;

            fq_poly_factor_equal_deg(res, dist_deg->poly + j, degs[l], ctx);
            for (k = res_num; k < res->num; k++)
                res->exp[k] = fq_poly_remove(v, res->poly + k, ctx);
        }
    }

    flint_free(degs);
    fq_poly_clear(v);
    fq_poly_factor_clear(dist_deg);
    fq_poly_factor_clear(sq_free);
}
