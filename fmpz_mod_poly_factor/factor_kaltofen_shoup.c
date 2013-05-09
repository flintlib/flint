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

******************************************************************************/

#include <math.h>
#include "fmpz_mod_poly_factor.h"

void
fmpz_mod_poly_factor_kaltofen_shoup(fmpz_mod_poly_factor_t res,
                                    const fmpz_mod_poly_t poly)
{
    fmpz_mod_poly_t v;
    fmpz_mod_poly_factor_t sq_free, dist_deg;
    len_t i, j, k, l, res_num, dist_deg_num;
    len_t *degs;

    if (!(degs = flint_malloc(fmpz_mod_poly_degree(poly) * sizeof(len_t))))
    {
        printf("Exception (fmpz_mod_poly_factor_kaltofen_shoup): \n");
        printf("Not enough memory.\n");
        abort();
    }

    fmpz_mod_poly_init(v, &poly->p);

    fmpz_mod_poly_make_monic(v, poly);

    /* compute squarefree factorisation */
    fmpz_mod_poly_factor_init(sq_free);
    fmpz_mod_poly_factor_squarefree(sq_free, v);

    /* compute distinct-degree factorisation */
    fmpz_mod_poly_factor_init(dist_deg);
    for (i = 0; i < sq_free->num; i++)
    {
        dist_deg_num = dist_deg->num;

        fmpz_mod_poly_factor_distinct_deg(dist_deg, sq_free->poly + i, &degs);

        /* compute equal-degree factorisation */
        for (j = dist_deg_num, l = 0; j < dist_deg->num; j++, l++)
        {
            res_num = res->num;

            fmpz_mod_poly_factor_equal_deg(res, dist_deg->poly + j, degs[l]);
            for (k = res_num; k < res->num; k++)
                res->exp[k] = fmpz_mod_poly_remove(v, res->poly + k);
        }
    }

    flint_free(degs);
    fmpz_mod_poly_clear(v);
    fmpz_mod_poly_factor_clear(dist_deg);
    fmpz_mod_poly_factor_clear(sq_free);
}
