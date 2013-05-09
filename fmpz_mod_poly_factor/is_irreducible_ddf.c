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

int fmpz_mod_poly_is_irreducible_ddf(const fmpz_mod_poly_t f)
{
    len_t i, n;
    len_t *degs;
    fmpz_mod_poly_factor_t dist_deg;

    n = fmpz_mod_poly_degree(f);

    if (n < 2)
        return 1;

    if (!fmpz_mod_poly_is_squarefree(f))
        return 0;

    if (!(degs = (len_t *)flint_malloc(n * sizeof(len_t))))
    {
        printf("Exception (fmpz_mod_poly_is_irreducible_ddf): \n");
        printf("Not enough memory.\n");
        abort();
    }

    fmpz_mod_poly_factor_init(dist_deg);
    fmpz_mod_poly_factor_distinct_deg(dist_deg, f, &degs);
    for (i = 0; i < dist_deg->num; i++)
    {
        if (degs[i] == n)
        {
            flint_free(degs);
            fmpz_mod_poly_factor_clear(dist_deg);
            return 1;
        }
        if (degs[i] > 0)
        {
            flint_free(degs);
            fmpz_mod_poly_factor_clear(dist_deg);
            return 0;
        }
    }

    flint_free(degs);
    fmpz_mod_poly_factor_clear(dist_deg);

    return 1;
}
