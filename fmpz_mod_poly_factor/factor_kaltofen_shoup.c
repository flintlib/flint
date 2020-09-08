/*
    Copyright (C) 2012 Lina Kulakova

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "fmpz_mod_poly.h"
#include "thread_support.h"

void
fmpz_mod_poly_factor_kaltofen_shoup(fmpz_mod_poly_factor_t res,
                          const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_t v;
    fmpz_mod_poly_factor_t sq_free, dist_deg;
    slong i, j, k, l, res_num, dist_deg_num;
    slong *degs;

    fmpz_mod_poly_init(v, ctx);

    fmpz_mod_poly_make_monic(v, poly, ctx);
    
    if (poly->length <= 2)
    {
        fmpz_mod_poly_factor_insert(res, v, 1, ctx);
        fmpz_mod_poly_clear(v, ctx);
        return;        
    }
    
    if (!(degs = flint_malloc(fmpz_mod_poly_degree(poly, ctx) * sizeof(slong))))
    {
        flint_printf("Exception (fmpz_mod_poly_factor_kaltofen_shoup): \n");
        flint_printf("Not enough memory.\n");
        flint_abort();
    }

    /* compute squarefree factorisation */
    fmpz_mod_poly_factor_init(sq_free, ctx);
    fmpz_mod_poly_factor_squarefree(sq_free, v, ctx);

    /* compute distinct-degree factorisation */
    fmpz_mod_poly_factor_init(dist_deg, ctx);
    for (i = 0; i < sq_free->num; i++)
    {
        dist_deg_num = dist_deg->num;

        if ((flint_get_num_threads() > 1) &&
            ((sq_free->poly + i)->length > (1024*flint_get_num_threads())/4))
            fmpz_mod_poly_factor_distinct_deg_threaded(dist_deg,
                                                sq_free->poly + i, &degs, ctx);
        else
            fmpz_mod_poly_factor_distinct_deg(dist_deg, sq_free->poly + i,
                                              &degs, ctx);

        /* compute equal-degree factorisation */
        for (j = dist_deg_num, l = 0; j < dist_deg->num; j++, l++)
        {
            res_num = res->num;

            fmpz_mod_poly_factor_equal_deg(res, dist_deg->poly + j, degs[l], ctx);
            for (k = res_num; k < res->num; k++)
                res->exp[k] = fmpz_mod_poly_remove(v, res->poly + k, ctx);
        }
    }

    flint_free(degs);
    fmpz_mod_poly_clear(v, ctx);
    fmpz_mod_poly_factor_clear(dist_deg, ctx);
    fmpz_mod_poly_factor_clear(sq_free, ctx);
}
