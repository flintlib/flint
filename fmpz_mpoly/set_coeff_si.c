/*
    Copyright (C) 2008, 2009, 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"

void
fmpz_mpoly_set_coeff_si(fmpz_mpoly_t poly,
                                  slong n, slong x, const fmpz_mpoly_ctx_t ctx)
{
    if (x == 0)
    {
       fmpz ptr;
       slong i, m;

       if (n >= poly->length)
          return;

       fmpz_zero(poly->coeffs + n);
       ptr = poly->coeffs[n];

       for (i = n; i < poly->length - 1; i++)
          poly->coeffs[i] = poly->coeffs[i + 1];

       poly->coeffs[i] = ptr;

       m = (poly->bits*ctx->n - 1)/FLINT_BITS + 1;

       for (i = n*m; i < (poly->length - 1)*m; i++)
          poly->exps[i] = poly->exps[i + m];

       poly->length--;
    }
    else
    {
        fmpz_mpoly_fit_length(poly, n + 1, ctx);

        if (n + 1 > poly->length)
        {
           slong i;
           
           for (i = poly->length; i < n; i++)
               fmpz_zero(poly->coeffs + i);
           
           poly->length = n + 1;
        }

        fmpz_set_si(poly->coeffs + n, x);
    }
}
