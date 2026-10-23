/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "mpn_mod.h"

/*
    The context already holds everything flint_mpn_sqrtmod_preinv wants:
    the modulus, the precomputed inverse of its normalisation, and the
    normalisation itself.
*/
truth_t
mpn_mod_is_square(nn_srcptr x, gr_ctx_t ctx)
{
    if (mpn_mod_is_zero(x, ctx) == T_TRUE || mpn_mod_is_one(x, ctx) == T_TRUE)
        return T_TRUE;

    if (MPN_MOD_CTX_MODULUS(ctx)[0] & 1)
    {
        /* a Jacobi symbol of -1 rules out a square, prime modulus or not */
        if (!flint_mpn_is_square_mod(x, MPN_MOD_CTX_MODULUS(ctx),
                    MPN_MOD_CTX_NLIMBS(ctx)))
            return T_FALSE;

        if (MPN_MOD_CTX_IS_PRIME(ctx) == T_TRUE)
            return T_TRUE;
    }

    return T_UNKNOWN;
}

int
mpn_mod_sqrt(nn_ptr res, nn_srcptr x, gr_ctx_t ctx)
{
    int found;

    if (mpn_mod_is_zero(x, ctx) == T_TRUE || mpn_mod_is_one(x, ctx) == T_TRUE)
        return mpn_mod_set(res, x, ctx);

    /* an even modulus is not prime here, whatever the context claims */
    if (!(MPN_MOD_CTX_MODULUS(ctx)[0] & 1))
        return GR_UNABLE;

    if (MPN_MOD_CTX_IS_PRIME(ctx) != T_TRUE)
    {
        /* the Jacobi symbol can still rule a root out */
        if (mpn_mod_is_square(x, ctx) == T_FALSE)
        {
            mpn_mod_zero(res, ctx);
            return GR_DOMAIN;
        }

        return GR_UNABLE;
    }

    /* the symbol is computed once, inside, and res is zeroed on failure */
    found = flint_mpn_sqrtmod_preinv(res, x, MPN_MOD_CTX_MODULUS(ctx),
                MPN_MOD_CTX_NLIMBS(ctx), MPN_MOD_CTX_MODULUS_PREINV(ctx),
                MPN_MOD_CTX_NORM(ctx));

    if (found == 1)
        return GR_SUCCESS;

    /* a failure of the algorithm means the modulus is not prime after all */
    return (found == 0) ? GR_DOMAIN : GR_UNABLE;
}
