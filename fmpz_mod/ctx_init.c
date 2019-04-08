/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mod.h"

void fmpz_mod_ctx_init(fmpz_mod_ctx_t ctx, const fmpz_t n)
{
    mp_bitcnt_t bits;

    if  (fmpz_sgn(n) <= 0)
    {
        flint_throw(FLINT_ERROR, "Nonpositive modulus in fmpz_mod_ctx_init");
    }

    fmpz_init_set(ctx->n, n);
    fmpz_preinvn_init(ctx->ninv, n);

    bits = fmpz_bits(n);

    if (bits < FLINT_BITS)
    {
        ctx->sub_fxn = fmpz_mod_sub1;
        ctx->n_limbs[0] = fmpz_get_ui(n);
        ctx->n_limbs[1] = 0;
        ctx->n_limbs[2] = 0;
    }
    else if (bits < 2*FLINT_BITS)
    {
        ctx->sub_fxn = fmpz_mod_sub2;
        fmpz_get_uiui(ctx->n_limbs + 1, ctx->n_limbs + 0, n);
        ctx->n_limbs[2] = 0;
    }
    else if (bits < 3*FLINT_BITS)
    {
        ctx->sub_fxn = fmpz_mod_sub3;
        fmpz_get_uiuiui(ctx->n_limbs + 2, ctx->n_limbs + 1, ctx->n_limbs + 0, n);
    }
    else
    {
        ctx->sub_fxn = fmpz_mod_subN;
        ctx->n_limbs[0] = 0;
        ctx->n_limbs[1] = 0;
        ctx->n_limbs[2] = 0;
    }
}
