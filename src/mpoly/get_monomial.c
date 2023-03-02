/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "mpoly.h"

/* !!! this file DOES need to change with new orderings */

/*
    poly_exps (when unpacked) must contain the exponent of each variable in
    some field. These functions should place the exponent of each variable
    in the corresponding entry of user_exps.
*/

void mpoly_get_monomial_ui_unpacked_ffmpz(ulong * user_exps,
                                const fmpz * poly_exps, const mpoly_ctx_t mctx)
{
    slong i;

    /* poly_exps is already unpacked, just read of the correct fields */
    for (i = 0; i < mctx->nvars; i++)
    {
        slong off = mctx->rev ? i : mctx->nvars - 1 - i;
        FLINT_ASSERT(fmpz_abs_fits_ui(poly_exps + off));
        user_exps[i] = fmpz_get_ui(poly_exps + off);
    }
}

void mpoly_get_monomial_ffmpz_unpacked_ffmpz(fmpz * user_exps,
                                const fmpz * poly_exps, const mpoly_ctx_t mctx)
{
    slong i;

    /* poly_exps is already unpacked, just read of the correct fields */
    for (i = 0; i < mctx->nvars; i++)
    {
        slong off = mctx->rev ? i : mctx->nvars - 1 - i;
        fmpz_set(user_exps + i, poly_exps + off);
    }
}

void mpoly_get_monomial_pfmpz_unpacked_ffmpz(fmpz ** user_exps,
                                const fmpz * poly_exps, const mpoly_ctx_t mctx)
{
    slong i;

    /* poly_exps is already unpacked, just read of the correct fields */
    for (i = 0; i < mctx->nvars; i++)
    {
        slong off = mctx->rev ? i : mctx->nvars - 1 - i;
        fmpz_set(user_exps[i], poly_exps + off);
    }
}

void mpoly_get_monomial_ui_unpacked_ui(ulong * user_exps,
                               const ulong * poly_exps, const mpoly_ctx_t mctx)
{
    slong i;

    /* poly_exps is already unpacked, just read of the correct fields */
    for (i = 0; i < mctx->nvars; i++)
    {
        slong off = mctx->rev ? i : mctx->nvars - 1 - i;
        user_exps[i] = poly_exps[off];
    }
}

void mpoly_get_monomial_ui_sp(ulong * user_exps, const ulong * poly_exps,
                                   flint_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    slong nvars = mctx->nvars;
    slong i, shift;
    ulong u, mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    ulong * exp1;
    const ulong * exp2;
    slong dir;

    FLINT_ASSERT(bits <= FLINT_BITS);

    exp2 = poly_exps;
    exp1 = user_exps + nvars - 1;
    dir = -WORD(1);
    if (mctx->rev)
    {
        exp1 = user_exps;
        dir = UWORD(1);
    }

    if (nvars < 1)
        return;

    i = 0;
    u = *exp2++;
    shift = 0;
    *exp1 = u & mask;
    exp1 += dir;
    u = u >> bits;      /* number of bits to encode 0th field */
    shift += bits;      /* number of bits to encode 0th field */
    while (++i < nvars)
    {
        if (shift + bits > FLINT_BITS)
        {
            u = *exp2++;
            shift = 0;
        }
        *exp1 = u & mask;
        exp1 += dir;
        u = u >> bits;      /* number of bits to encode ith field */
        shift += bits;      /* number of bits to encode ith field */
    }
}

void mpoly_get_monomial_ui_mp(ulong * user_exps, const ulong * poly_exps,
                                   flint_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    slong nvars = mctx->nvars;
    slong i, j;
    ulong * exp1;
    const ulong * exp2;
    ulong words_per_field = bits/FLINT_BITS;
    ulong check_mask;
    slong dir;

    FLINT_ASSERT(bits%FLINT_BITS == 0);
    FLINT_ASSERT(bits > FLINT_BITS);

    exp2 = poly_exps;
    exp1 = user_exps + nvars - 1;
    dir = -WORD(1);
    if (mctx->rev)
    {
        exp1 = user_exps;
        dir = UWORD(1);
    }

    check_mask = 0;
    for (i = 0; i < nvars; i++)
    {
        *exp1 = *exp2;
        exp1 += dir;

        for (j = 1; j < words_per_field; j++)
            check_mask |= exp2[j];

        exp2 += words_per_field;
    }

    if (check_mask != 0)
        flint_throw(FLINT_ERROR, "Exponent vector does not fit a ulong.");
}

void mpoly_get_monomial_si_mp(slong * user_exps, const ulong * poly_exps,
                                   flint_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    slong nvars = mctx->nvars;
    slong i, j;
    slong * exp1;
    const ulong * exp2;
    ulong words_per_field = bits/FLINT_BITS;
    ulong check_mask;
    slong dir;

    FLINT_ASSERT(bits%FLINT_BITS == 0);
    FLINT_ASSERT(bits > FLINT_BITS);

    exp2 = poly_exps;
    exp1 = user_exps + nvars - 1;
    dir = -WORD(1);
    if (mctx->rev)
    {
        exp1 = user_exps;
        dir = UWORD(1);
    }

    check_mask = 0;
    for (i = 0; i < nvars; i++)
    {
        *exp1 = (slong) *exp2;
        exp1 += dir;

        check_mask |= FLINT_SIGN_EXT(exp2[0]);
        for (j = 1; j < words_per_field; j++)
            check_mask |= exp2[j];

        exp2 += words_per_field;
    }

    if (check_mask != 0)
        flint_throw(FLINT_ERROR, "Exponent vector does not fit an slong.");
}

void mpoly_get_monomial_ffmpz(fmpz * user_exps, const ulong * poly_exps,
                                   flint_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    slong i, nvars = mctx->nvars;
    fmpz * tmp_exps;
    TMP_INIT;

    TMP_START;
    tmp_exps = (fmpz *) TMP_ALLOC(mctx->nfields*sizeof(fmpz));
    for (i = 0; i < mctx->nfields; i++)
        fmpz_init(tmp_exps + i);

    mpoly_unpack_vec_fmpz(tmp_exps, poly_exps, bits, mctx->nfields, 1);

    for (i = 0; i < nvars; i++)
        fmpz_swap(user_exps + i, tmp_exps + (mctx->rev ? i : nvars - 1 - i));

    for (i = 0; i < mctx->nfields; i++)
        fmpz_clear(tmp_exps + i);

    TMP_END;
}

void mpoly_get_monomial_pfmpz(fmpz ** user_exps, const ulong * poly_exps,
                                   flint_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    slong i, nvars = mctx->nvars;
    fmpz * tmp_exps;
    TMP_INIT;

    TMP_START;
    tmp_exps = (fmpz *) TMP_ALLOC(mctx->nfields*sizeof(fmpz));
    for (i = 0; i < mctx->nfields; i++)
        fmpz_init(tmp_exps + i);

    mpoly_unpack_vec_fmpz(tmp_exps, poly_exps, bits, mctx->nfields, 1);

    for (i = 0; i < nvars; i++)
        fmpz_swap(user_exps[i], tmp_exps + (mctx->rev ? i : nvars - 1 - i));

    for (i = 0; i < mctx->nfields; i++)
        fmpz_clear(tmp_exps + i);

    TMP_END;
}
