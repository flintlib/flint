/*
    Copyright (C) 2019 Daniel Schultz
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "longlong.h"
#include "fmpz.h"
#include "fmpz_mod.h"

#include "mpn_extras.h"

/* canonical x in [0, n) as mn limbs */
static void
_fmpz_mod_to_mpn(nn_ptr bp, const fmpz_t x, mp_size_t mn)
{
    if (!COEFF_IS_MPZ(*x))
    {
        bp[0] = *x;
        flint_mpn_zero(bp + 1, mn - 1);
    }
    else
    {
        mpz_ptr zx = COEFF_TO_PTR(*x);
        mp_size_t c = zx->_mp_size;
        flint_mpn_copyi(bp, zx->_mp_d, c);
        flint_mpn_zero(bp + c, mn - c);
    }
}

/* x canonical, exponent (ep, en) positive; uses the context's
   precomputed inverse so that small exponents need no setup at all */
static void
_fmpz_mod_pow_mpn(fmpz_t res, const fmpz_t x, nn_srcptr ep, mp_size_t en,
                  const fmpz_mod_ctx_t ctx)
{
    mpz_ptr zn = COEFF_TO_PTR(*ctx->n);
    mp_size_t mn = zn->_mp_size;
    nn_ptr bp, rp;
    TMP_INIT;

    TMP_START;
    bp = TMP_ALLOC(2 * mn * sizeof(ulong));
    rp = bp + mn;
    _fmpz_mod_to_mpn(bp, x, mn);
    flint_mpn_powm_preinvn(rp, bp, ep, en, zn->_mp_d, mn,
                           ctx->ninv_huge->dinv, ctx->ninv_huge->norm);
    fmpz_set_ui_array(res, rp, mn);
    TMP_END;
}

static void
_fmpz_mod_pow_fmpz(fmpz_t res, const fmpz_t x, const fmpz_t e, const fmpz_mod_ctx_t ctx)
{
    if (*e <= 2)
    {
        if (*e == 0)
            fmpz_mod_set_ui(res, 1, ctx);
        else if (*e == 1)
            fmpz_set(res, x);
        else
            fmpz_mod_mul(res, x, x, ctx);
    }
    else if (fmpz_is_zero(x) || fmpz_is_one(x))
    {
        fmpz_set(res, x);
    }
    else if (COEFF_IS_MPZ(*ctx->n) && ctx->ninv_huge != NULL)
    {
        if (!COEFF_IS_MPZ(*e))
        {
            ulong ee = *e;
            _fmpz_mod_pow_mpn(res, x, &ee, 1, ctx);
        }
        else
        {
            mpz_ptr ze = COEFF_TO_PTR(*e);
            _fmpz_mod_pow_mpn(res, x, ze->_mp_d, ze->_mp_size, ctx);
        }
    }
    else
    {
        fmpz_powm(res, x, e, ctx->n);
    }
}

static void
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
    else if (COEFF_IS_MPZ(*ctx->n) && ctx->ninv_huge != NULL)
    {
        _fmpz_mod_pow_mpn(res, x, &e, 1, ctx);
    }
    else
    {
        fmpz_powm_ui(res, x, e, ctx->n);
    }
}

int fmpz_mod_pow_fmpz(fmpz_t a, const fmpz_t b, const fmpz_t pow,
                                                      const fmpz_mod_ctx_t ctx)
{
    int success = 1;
    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));

    if (fmpz_sgn(pow) >= 0)
    {
        _fmpz_mod_pow_fmpz(a, b, pow, ctx);
    }
    else
    {
        fmpz_t d;
        fmpz_init(d);
        fmpz_gcdinv(d, a, b, ctx->n);
        if (fmpz_is_one(d))
        {
            fmpz_neg(d, pow);
            _fmpz_mod_pow_fmpz(a, a, d, ctx);
        }
        else
        {
            success = 0;
        }
        fmpz_clear(d);
    }

    FLINT_ASSERT(!success || fmpz_mod_is_canonical(a, ctx));
    return success;
}

void fmpz_mod_pow_ui(fmpz_t a, const fmpz_t b, ulong pow,
                                                     const fmpz_mod_ctx_t ctx)
{
    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));
    _fmpz_mod_pow_ui(a, b, pow, ctx);
    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}
