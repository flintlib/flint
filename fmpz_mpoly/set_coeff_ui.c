/*
    Copyright (C) 2008, 2009 William Hart

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
fmpz_mpoly_set_coeff_ui(fmpz_mpoly_t poly,
                                  slong n, ulong x, const fmpz_mpoly_ctx_t ctx)
{
    if (x == 0)
    {
       fmpz ptr;
       slong i, N;

       if (n >= poly->length)
          return;

       fmpz_zero(poly->coeffs + n);
       ptr = poly->coeffs[n];

       for (i = n; i < poly->length - 1; i++)
          poly->coeffs[i] = poly->coeffs[i + 1];

       poly->coeffs[i] = ptr;

       N = words_per_exp(ctx->n, poly->bits);

       for (i = n*N; i < (poly->length - 1)*N; i++)
          poly->exps[i] = poly->exps[i + N];

       poly->length--;
    }
    else
    {
        fmpz_mpoly_fit_length(poly, n + 1, ctx);

        if (n == poly->length)
           poly->length++;
        else if (n > poly->length)
           flint_throw(FLINT_ERROR,
                                   "Invalid index in fmpz_mpoly_set_coeff_ui");

        fmpz_set_ui(poly->coeffs + n, x);
    }
}
