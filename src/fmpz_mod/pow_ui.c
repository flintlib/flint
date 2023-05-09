/*
    Copyright (C) 2019 Daniel Schultz
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mod.h"
#include "gr.h"

#ifdef FLINT_HAVE_FFT_SMALL

void
_fmpz_mod_pow_ui(fmpz_t res, const fmpz_t x, ulong e, const fmpz_mod_ctx_t ctx)
{
    if (e <= 2)
    {
        if (e == 0)
            fmpz_mod_set_ui(res, 1, ctx);
        else if (e == 1)
            fmpz_set(res, x);
        else
            fmpz_mod_mul(res, x, x, ctx);
    }
    else if (fmpz_is_zero(x) || fmpz_is_one(x))
    {
        fmpz_set(res, x);
    }
    else if (fmpz_bits(ctx->n) < 70000)
    {
        fmpz_powm_ui(res, x, e, ctx->n);
    }
    else
    {
        gr_ctx_t gctx;
        _gr_ctx_init_fmpz_mod_from_ref(gctx, ctx);

        if (!COEFF_IS_MPZ(*x) || FLINT_BIT_COUNT(e) < 20)
            GR_MUST_SUCCEED(gr_generic_pow_ui_binexp(res, x, e, gctx));
        else
            GR_MUST_SUCCEED(gr_generic_pow_ui_sliding(res, x, e, gctx));
    }
}

#else

void
_fmpz_mod_pow_ui(fmpz_t res, const fmpz_t x, ulong e, const fmpz_mod_ctx_t ctx)
{
    fmpz_powm_ui(res, x, e, ctx->n);
}

#endif

void fmpz_mod_pow_ui(fmpz_t a, const fmpz_t b, ulong pow,
                                                     const fmpz_mod_ctx_t ctx)
{
    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));
    _fmpz_mod_pow_ui(a, b, pow, ctx);
    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}
