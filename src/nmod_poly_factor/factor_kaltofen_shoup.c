/*
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2022 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "nmod_poly.h"
#include "nmod_poly_factor.h"

void nmod_poly_factor_kaltofen_shoup(nmod_poly_factor_t res,
                                     const nmod_poly_t poly)
{
    nmod_poly_t v;
    nmod_poly_factor_t sq_free, dist_deg;
    slong i, j, k, l, res_num, dist_deg_num;
    slong * degs;

    nmod_poly_init_mod(v, poly->mod);

    nmod_poly_make_monic(v, poly);

    if (poly->length <= 2)
    {
        nmod_poly_factor_insert (res, v, 1);
        nmod_poly_clear (v);

        return;
    }

    degs = flint_malloc(nmod_poly_degree(poly) * sizeof(slong));

    /* compute squarefree factorisation */
    nmod_poly_factor_init(sq_free);
    nmod_poly_factor_squarefree(sq_free, v);

    /* compute distinct-degree factorisation */
    nmod_poly_factor_init(dist_deg);
    for (i = 0; i < sq_free->num; i++)
    {
        dist_deg_num = dist_deg->num;

        if ((flint_get_num_threads() > 1) &&
            ((sq_free->p + i)->length > (1024*flint_get_num_threads())/4))
            nmod_poly_factor_distinct_deg_threaded(dist_deg, sq_free->p + i,
                                                                        &degs);
        else
            nmod_poly_factor_distinct_deg(dist_deg, sq_free->p + i, &degs);

        /* compute equal-degree factorisation */
        for (j = dist_deg_num, l = 0; j < dist_deg->num; j++, l++)
        {
            res_num = res->num;

            nmod_poly_factor_equal_deg(res, dist_deg->p + j, degs[l]);
            for (k = res_num; k < res->num; k++)
                res->exp[k] = nmod_poly_remove(v, res->p + k);
        }
    }

    flint_free(degs);
    nmod_poly_clear(v);
    nmod_poly_factor_clear(dist_deg);
    nmod_poly_factor_clear(sq_free);
}
