/*
    Copyright (C) 2017 Daniel Schultz

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

slong _fmpz_mpoly_derivative(fmpz * coeff1, ulong * exp1,
                         const fmpz * coeff2, const ulong * exp2, slong len,
                slong var, int deg, int rev, slong nfields, slong bits, slong N)
{
    slong off, shift, fpw, i, k = 0;
    ulong c, mask;
    ulong*one;
    TMP_INIT;

    TMP_START;

    fpw = FLINT_BITS/bits;
    mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    mpoly_off_shift(&off, &shift, var, deg, rev, fpw, nfields, bits);
    one = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_univar_exp(one, var, deg, N, off, shift, fpw, bits);

    /* x^c -> c*x^(c-1) */
    for (i = 0; i < len; i++)
    {
        c = (exp2[N*i + off] >> shift) & mask;
        if (c != 0)
        {
            mpoly_monomial_sub(exp1 + N*k, exp2 + N*i, one, N);
            fmpz_mul_ui(coeff1 + k, coeff2 + i, c);
            k++;
        }
    }
    TMP_END;
    return k;
}

void fmpz_mpoly_derivative(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                                         slong var, const fmpz_mpoly_ctx_t ctx)
{
    int deg, rev;
    slong N, len;

    N = words_per_exp(ctx->n, poly2->bits);
    degrev_from_ord(deg, rev, ctx->ord);

    fmpz_mpoly_fit_length(poly1, poly2->length, ctx);
    fmpz_mpoly_fit_bits(poly1, poly2->bits, ctx);
    poly1->bits = poly2->bits;

    len = _fmpz_mpoly_derivative(poly1->coeffs, poly1->exps,
                   poly2->coeffs, poly2->exps, poly2->length,
                        var, deg, rev, ctx->n, poly2->bits, N);

    _fmpz_mpoly_set_length(poly1, len, ctx);
}
