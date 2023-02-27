/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

#include "math.h"

int
TEMPLATE(T, poly_is_irreducible_ben_or) (const TEMPLATE(T, poly_t) f,
                                         const TEMPLATE(T, ctx_t) ctx)
{
    int result;
    slong i, n;
    fmpz_t q;
    TEMPLATE(T, poly_t) g, x, xq, xqimx;
    TEMPLATE(T, poly_t) v, vinv;

    n = TEMPLATE(T, poly_degree) (f, ctx);

    if (n < 2)
        return 1;

    if (!TEMPLATE(T, poly_is_squarefree) (f, ctx))
        return 0;

    TEMPLATE(T, poly_init) (v, ctx);
    TEMPLATE(T, poly_init) (vinv, ctx);
    TEMPLATE(T, poly_make_monic) (v, f, ctx);
    TEMPLATE(T, poly_reverse) (vinv, v, v->length, ctx);
    TEMPLATE(T, poly_inv_series_newton) (vinv, vinv, v->length, ctx);

    TEMPLATE(T, poly_init) (x, ctx);
    TEMPLATE(T, poly_init) (xq, ctx);
    TEMPLATE(T, poly_init) (xqimx, ctx);

    /* Compute xq = x^q mod f */
    fmpz_init(q);
    fmpz_pow_ui(q, TEMPLATE(T, ctx_prime) (ctx),
                TEMPLATE(T, ctx_degree) (ctx));
    TEMPLATE(T, poly_gen) (x, ctx);
    TEMPLATE(T, poly_powmod_fmpz_sliding_preinv) (xq, x, q, 0, v, vinv, ctx);
    TEMPLATE(T, poly_set) (xqimx, xq, ctx);

    result = 1;
    TEMPLATE(T, poly_init) (g, ctx);
    for (i = 1; i <= n / 2; i++)
    {
        TEMPLATE(T, poly_sub) (xqimx, xqimx, x, ctx);
        TEMPLATE(T, poly_gcd) (g, f, xqimx, ctx);
        if (!TEMPLATE(T, poly_is_one) (g, ctx))
        {
            result = 0;
            break;
        }
        else if (i == n / 2)
        {
            /* We don't need to compute the last step */
            break;
        }

        TEMPLATE(T, poly_add) (xqimx, xqimx, x, ctx);

        if (TEMPLATE(CAP_T, POLY_ITERATED_FROBENIUS_CUTOFF) (ctx, v->length))
        {
            TEMPLATE(T, poly_powmod_fmpz_sliding_preinv) (xqimx, xqimx, q, 0,
                                                          v, vinv, ctx);
        }
        else
        {
            TEMPLATE(T, poly_compose_mod_preinv) (xqimx, xqimx, xq, v, vinv,
                                                  ctx);
        }
    }

    TEMPLATE(T, poly_clear) (g, ctx);
    TEMPLATE(T, poly_clear) (x, ctx);
    TEMPLATE(T, poly_clear) (xq, ctx);
    TEMPLATE(T, poly_clear) (xqimx, ctx);
    TEMPLATE(T, poly_clear) (v, ctx);
    TEMPLATE(T, poly_clear) (vinv, ctx);
    fmpz_clear(q);

    return result;
}


#endif
