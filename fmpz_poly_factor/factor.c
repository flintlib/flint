/*
    Copyright (C) 2011 Andy Novocin
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "fmpz_poly.h"

void fmpz_poly_factor(fmpz_poly_factor_t fac, const fmpz_poly_t G)
{
    const slong lenG = G->length;
    fmpz_poly_t g;

    if (lenG == 0)
    {
        fmpz_set_ui(&fac->c, 0);
        return;
    }
    if (lenG == 1)
    {
        fmpz_set(&fac->c, G->coeffs);
        return;
    }

    fmpz_poly_init(g);

    if (lenG == 2)
    {
        fmpz_poly_content(&fac->c, G);
        if (fmpz_sgn(fmpz_poly_lead(G)) < 0)
            fmpz_neg(&fac->c, &fac->c);
        fmpz_poly_scalar_divexact_fmpz(g, G, &fac->c);
        fmpz_poly_factor_insert(fac, g, 1);
    }
    else
    {
        slong j, k;
        fmpz_poly_factor_t sq_fr_fac;

        /* Does a presearch for a factor of form x^k */
        for (k = 0; fmpz_is_zero(G->coeffs + k); k++) ;

        if (k != 0)
        {
            fmpz_poly_t t;

            fmpz_poly_init(t);
            fmpz_poly_set_coeff_ui(t, 1, 1);
            fmpz_poly_factor_insert(fac, t, k);
            fmpz_poly_clear(t);
        }

        fmpz_poly_shift_right(g, G, k);

        /* Could make other tests for x-1 or simple things 
           maybe take advantage of the composition algorithm */
        fmpz_poly_factor_init(sq_fr_fac);
        fmpz_poly_factor_squarefree(sq_fr_fac, g);

        fmpz_set(&fac->c, &sq_fr_fac->c);

        /* Factor each square-free part */
        for (j = 0; j < sq_fr_fac->num; j++)
            _fmpz_poly_factor_zassenhaus(fac, sq_fr_fac->exp[j], sq_fr_fac->p + j, 10, 1);

        fmpz_poly_factor_clear(sq_fr_fac);
    }
    fmpz_poly_clear(g);
}


