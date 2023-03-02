/*
    Copyright (C) 2011 Andy Novocin
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "fmpz_poly.h"

void fmpz_poly_factor_deflation(fmpz_poly_factor_t fac, const fmpz_poly_t G, int deflation)
{
    const slong lenG = G->length;
    fmpz_poly_t g;

    fac->num = 0;

    if (lenG <= 1)
    {
        if (lenG < 1)
            fmpz_zero(&fac->c);
        else
            fmpz_set(&fac->c, G->coeffs + 0);
        return;
    }

    fmpz_poly_init(g);

    if (lenG < 5)
    {
        fmpz_poly_content(&fac->c, G);
        if (fmpz_sgn(fmpz_poly_lead(G)) < 0)
            fmpz_neg(&fac->c, &fac->c);
        fmpz_poly_scalar_divexact_fmpz(g, G, &fac->c);

        if (lenG < 3)
            fmpz_poly_factor_insert(fac, g, 1);
        else if (lenG == 3)
            _fmpz_poly_factor_quadratic(fac, g, 1);
        else
            _fmpz_poly_factor_cubic(fac, g, 1);
    }
    else
    {
        slong i, j, k, d;
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

        if (deflation && (d = fmpz_poly_deflation(G)) > 1)
        {
            fmpz_poly_factor_t gfac;
            fmpz_poly_factor_init(gfac);

            fmpz_poly_deflate(g, g, d);
            fmpz_poly_factor(gfac, g);
            fmpz_set(&fac->c, &gfac->c);

            for (i = 0; i < gfac->num; i++)
            {
                fmpz_poly_factor_t hfac;
                fmpz_poly_factor_init(hfac);
                fmpz_poly_inflate(gfac->p + i, gfac->p + i, d);
                fmpz_poly_factor_deflation(hfac, gfac->p + i, 0);

                for (j = 0; j < hfac->num; j++)
                    fmpz_poly_factor_insert(fac, hfac->p + j, gfac->exp[i] * hfac->exp[j]);

                fmpz_poly_factor_clear(hfac);
            }

            fmpz_poly_factor_clear(gfac);
        }
        else
        {
            /* Could make other tests for x-1 or simple things 
               maybe take advantage of the composition algorithm */
            fmpz_poly_factor_init(sq_fr_fac);
            fmpz_poly_factor_squarefree(sq_fr_fac, g);

            fmpz_set(&fac->c, &sq_fr_fac->c);

            /* Factor each square-free part */
            for (j = 0; j < sq_fr_fac->num; j++)
            {
                _fmpz_poly_factor_zassenhaus(fac, sq_fr_fac->exp[j],
                                                           sq_fr_fac->p + j, 8, 1);
            }

            fmpz_poly_factor_clear(sq_fr_fac);
        }
    }
    fmpz_poly_clear(g);
}

void fmpz_poly_factor(fmpz_poly_factor_t fac, const fmpz_poly_t G)
{
    fmpz_poly_factor_deflation(fac, G, 1);
}
