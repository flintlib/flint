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

slong _fmpz_mpoly_integral(fmpz_t s, fmpz * coeff1, ulong * exp1,
                         const fmpz * coeff2, const ulong * exp2, slong len,
                slong idx, int deg, int rev, slong nfields, slong bits, slong N)
{
    slong off, shift, fpw, i;
    ulong c, mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    ulong * one;
    fmpz_t d, g;
    TMP_INIT;

    TMP_START;

    fmpz_init(d);
    fmpz_init(g);
    fmpz_set_si(s, WORD(1));

    fpw = FLINT_BITS/bits;
    if (rev)
    {
        off   = (nfields - 1 - idx)/fpw;
        shift = (nfields - 1 - idx)%fpw;
    } else
    {
        off   = (deg + idx)/fpw;
        shift = (deg + idx)%fpw;
    }
    shift = (fpw - 1 - shift) * bits;

    /* get exponent to add */
    one = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    for (i = 0; i < N; i++)
    {
        one[i] = 0;
    }
    one[off] = WORD(1) << shift;
    if (deg)
    {
        one[0] |= WORD(1) << ((fpw - 1)*bits);
    }

    /* scan once to find required denominator */
    for (i = 0; i < len; i++)
    {
        c = (exp2[N*i + off] >> shift) & mask;
        fmpz_set_ui(d, c + 1);
        fmpz_gcd(g, d, coeff2 + i);
        fmpz_divexact(d, d, g);
        fmpz_lcm(s, s, d);
    }

    /* then scan again to compute the terms */
    /* x^n -> x^(n+1)/(n+1) */
    for (i = 0; i < len; i++)
    {
        c = (exp2[N*i + off] >> shift) & mask;
        mpoly_monomial_add(exp1 + N*i, exp2 + N*i, one, N);
        fmpz_set_ui(d, c + 1);
        fmpz_mul(g, coeff2 + i, s);
        fmpz_mul(coeff1 + i, coeff2 + i, g);
        fmpz_divexact(coeff1 + i, g, d);
    }

    fmpz_clear(g);
    fmpz_clear(d);
    TMP_END;
    return len;
}

void fmpz_mpoly_integral(fmpz_mpoly_t poly1, fmpz_t scale,
               const fmpz_mpoly_t poly2, slong idx, const fmpz_mpoly_ctx_t ctx)
{
    int deg, rev;
    slong len, N, exp_bits;
    ulong max_exp, * max_deg, * exp2 = poly2->exps;
    int free2 = 0;
    TMP_INIT;

    TMP_START;

    degrev_from_ord(deg, rev, ctx->ord);

    max_deg = (ulong *) TMP_ALLOC(ctx->n*sizeof(ulong));

    /* compute bits required to represent result */
    mpoly_max_degrees(max_deg, poly2->exps, poly2->length, poly2->bits, ctx->n);
    max_exp = max_deg[rev ? idx : ctx->n - 1 - deg - idx];
    if (deg)
    {
        max_exp = FLINT_MAX(max_exp, max_deg[ctx->n - 1]);
    }
    exp_bits = FLINT_MAX(WORD(8), 1 + FLINT_BIT_COUNT(max_exp + 1));
    exp_bits = FLINT_MAX(exp_bits, poly2->bits);
    if (exp_bits > FLINT_BITS)
    {
        flint_throw(FLINT_EXPOF, "Exponent overflow in fmpz_mpoly_integrate");
    }
    exp_bits = mpoly_optimize_bits(exp_bits, ctx->n);
    N = words_per_exp(ctx->n, exp_bits);

    /* ensure input exponents are packed into same sized fields as output */
    if (exp_bits > poly2->bits)
    {
        free2 = 1;
        exp2 = (ulong *) flint_malloc(N*poly2->length*sizeof(ulong));
        mpoly_unpack_monomials(exp2, exp_bits, poly2->exps, poly2->bits,
                                                        poly2->length, ctx->n);
    }

    /* deal with aliasing and do integration */
    if (poly1 == poly2)
    {
        fmpz_mpoly_t temp;

        fmpz_mpoly_init2(temp, poly2->length, ctx);
        fmpz_mpoly_fit_bits(temp, exp_bits, ctx);
        temp->bits = exp_bits;

        len = _fmpz_mpoly_integral(scale, temp->coeffs, temp->exps,
                       poly2->coeffs, exp2, poly2->length,
                            idx, deg, rev, ctx->n, exp_bits, N);
        _fmpz_mpoly_set_length(temp, len, ctx);

        fmpz_mpoly_swap(temp, poly1, ctx);
        fmpz_mpoly_clear(temp, ctx);
    } else
    {
        fmpz_mpoly_fit_length(poly1, poly2->length, ctx);
        fmpz_mpoly_fit_bits(poly1, exp_bits, ctx);
        poly1->bits = exp_bits;

        len = _fmpz_mpoly_integral(scale, poly1->coeffs, poly1->exps,
                       poly2->coeffs, exp2, poly2->length,
                            idx, deg, rev, ctx->n, exp_bits, N);
        _fmpz_mpoly_set_length(poly1, len, ctx);
    }

    if (free2)
    {
        flint_free(exp2);
    }

    TMP_END
}
