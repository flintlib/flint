/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "mpn_extras.h"
#include "mpn_mod.h"

/* powering via flint_mpn_powm_preinvn on the context's precomputed
   inverse: no per-call setup for small exponents */

int
mpn_mod_pow_ui(nn_ptr res, nn_srcptr x, ulong e, gr_ctx_t ctx)
{
    mp_size_t n = MPN_MOD_CTX_NLIMBS(ctx);
    ulong rp[MPN_MOD_MAX_LIMBS];

    flint_mpn_powm_preinvn(rp, x, &e, 1, MPN_MOD_CTX_MODULUS(ctx), n,
                           MPN_MOD_CTX_MODULUS_PREINV(ctx),
                           MPN_MOD_CTX_NORM(ctx));
    flint_mpn_copyi(res, rp, n);
    return GR_SUCCESS;
}

int
mpn_mod_pow_si(nn_ptr res, nn_srcptr x, slong e, gr_ctx_t ctx)
{
    if (e >= 0)
        return mpn_mod_pow_ui(res, x, (ulong) e, ctx);
    else
    {
        int status = mpn_mod_inv(res, x, ctx);
        if (status != GR_SUCCESS)
            return status;
        return mpn_mod_pow_ui(res, res, -(ulong) e, ctx);
    }
}

int
mpn_mod_pow_fmpz(nn_ptr res, nn_srcptr x, const fmpz_t e, gr_ctx_t ctx)
{
    if (!COEFF_IS_MPZ(*e))
        return mpn_mod_pow_si(res, x, *e, ctx);
    else
    {
        mpz_ptr ze = COEFF_TO_PTR(*e);
        mp_size_t n = MPN_MOD_CTX_NLIMBS(ctx);
        mp_size_t en = ze->_mp_size;
        ulong rp[MPN_MOD_MAX_LIMBS];
        nn_srcptr bp = x;
        int status = GR_SUCCESS;

        if (en < 0)
        {
            status = mpn_mod_inv(res, x, ctx);
            if (status != GR_SUCCESS)
                return status;
            bp = res;
            en = -en;
        }

        flint_mpn_powm_preinvn(rp, bp, ze->_mp_d, en,
                               MPN_MOD_CTX_MODULUS(ctx), n,
                               MPN_MOD_CTX_MODULUS_PREINV(ctx),
                               MPN_MOD_CTX_NORM(ctx));
        flint_mpn_copyi(res, rp, n);
        return GR_SUCCESS;
    }
}
