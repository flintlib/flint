/*
    Copyright (C) 2016 William Hart

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

void _fmpz_mpoly_gen1(fmpz * poly, ulong * exps, slong i,
                                         slong bits, slong n, int deg, int rev)
{
    slong k = FLINT_BITS/bits;

    fmpz_set_ui(poly + 0, 1);
    if (rev)
       exps[0] = (UWORD(1) << ((k - n + i)*bits));
    else
       exps[0] = (UWORD(1) << ((k - i - deg - 1)*bits));
    
    if (deg)
       exps[0] |= (UWORD(1) << ((k - 1)*bits));
}

void _fmpz_mpoly_gen(fmpz * poly, ulong * exps, slong i,
                                slong bits, slong n, int deg, int rev, slong N)
{
    slong j;
    ulong * mon;

    TMP_INIT;

    if (N == 1)
    {
       _fmpz_mpoly_gen1(poly, exps, i, bits, n, deg, rev);

       return;
    }
    
    fmpz_set_ui(poly + 0, 1);
    
    TMP_START;

    mon = (ulong *) TMP_ALLOC((n - deg)*sizeof(ulong));
    
    for (j = 0; j < n - deg; j++)
       mon[j] = 0;

    mon[i] = 1;

    _fmpz_mpoly_set_monomial(exps, mon, bits, n, deg, rev);

    TMP_END;
}

void fmpz_mpoly_gen(fmpz_mpoly_t poly, slong i, const fmpz_mpoly_ctx_t ctx)
{
   int deg, rev;

   slong N = (poly->bits*ctx->n - 1)/FLINT_BITS + 1;

   fmpz_mpoly_fit_length(poly, 1, ctx);

   degrev_from_ord(deg, rev, ctx->ord);

   _fmpz_mpoly_gen(poly->coeffs, poly->exps, i, poly->bits,
                                                          ctx->n, deg, rev, N);
      
   _fmpz_mpoly_set_length(poly, 1, ctx);
}
