/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T
#ifdef B

#include "templates.h"

void TEMPLATE(T, embed_gens)(TEMPLATE(T, t) gen_sub,
                             TEMPLATE(T, t) gen_sup,
                             TEMPLATE(B, poly_t) minpoly,
                             const TEMPLATE(T, ctx_t) sub_ctx,
                             const TEMPLATE(T, ctx_t) sup_ctx)
{
    if (TEMPLATE(T, ctx_degree)(sub_ctx) == 1) 
    {
        TEMPLATE(T, gen)(gen_sub, sub_ctx);
        TEMPLATE(T, set)(gen_sup, gen_sub, sup_ctx);
    }
    else 
    {
        _TEMPLATE(T, embed_gens_naive)(gen_sub, gen_sup, minpoly, sub_ctx, sup_ctx);
    }
}

void _TEMPLATE(T, embed_gens_naive)(TEMPLATE(T, t) gen_sub,
                                    TEMPLATE(T, t) gen_sup,
                                    TEMPLATE(B, poly_t) minpoly,
                                    const TEMPLATE(T, ctx_t) sub_ctx,
                                    const TEMPLATE(T, ctx_t) sup_ctx)
{
    TEMPLATE(T, poly_t) modulus, fact;
    flint_rand_t state;

    TEMPLATE(T, poly_init) (modulus, sup_ctx);
    TEMPLATE(T, poly_init) (fact, sup_ctx);
    TEMPLATE4(T, poly_set, B, poly)(modulus,
				    TEMPLATE(T, ctx_modulus)(sub_ctx),
				    sup_ctx);
    
    flint_randinit(state);

    /* Get one linear factor of sub_ctx->modulus in sup_ctx */
    while (TEMPLATE(T, poly_degree)(modulus, sup_ctx) != 1) 
    {
        while (!TEMPLATE(T, poly_factor_equal_deg_prob)
               (fact, state, modulus, 1, sup_ctx))
        {
        };
        TEMPLATE(T, poly_set)(modulus, fact, sup_ctx);
    }

    TEMPLATE(T, gen)(gen_sub, sub_ctx);
    TEMPLATE(T, set)(gen_sup, modulus->coeffs, sup_ctx);
    TEMPLATE(T, neg)(gen_sup, gen_sup, sup_ctx);

    TEMPLATE(B, poly_set)(minpoly, TEMPLATE(T, ctx_modulus)(sub_ctx));

    TEMPLATE(T, poly_clear) (modulus, sup_ctx);
    TEMPLATE(T, poly_clear) (fact, sup_ctx);
}


#endif
#endif
