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

void
TEMPLATE(T, poly_randtest_irreducible) (TEMPLATE(T, poly_t) f,
                                        flint_rand_t state, slong len,
                                        const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, poly_t) xq, xqi, x, g_i, finv;
    fmpz_t q;
    slong i, restart;

    /* Compute q */
    fmpz_init_set(q, TEMPLATE(T, ctx_prime) (ctx));
    fmpz_pow_ui(q, q, TEMPLATE(T, ctx_degree) (ctx));

    TEMPLATE(T, poly_init) (x, ctx);
    TEMPLATE(T, poly_gen) (x, ctx);
    TEMPLATE(T, poly_init) (xq, ctx);
    TEMPLATE(T, poly_init) (xqi, ctx);
    TEMPLATE(T, poly_init) (g_i, ctx);
    TEMPLATE(T, poly_init) (finv, ctx);

    while (1)
    {
        restart = 0;

        /* Generate random monic polynomial of length len */
        TEMPLATE(T, poly_randtest_monic) (f, state, len, ctx);

        TEMPLATE(T, poly_reverse) (finv, f, f->length, ctx);
        TEMPLATE(T, poly_inv_series_newton) (finv, finv, f->length, ctx);

        /* Compute xq = x^q mod f */
        TEMPLATE(T, poly_powmod_fmpz_binexp_preinv) (xq, x, q, f, finv, ctx);
        TEMPLATE(T, poly_set) (xqi, xq, ctx);

        for (i = 1; i <= (len - 1) / 2; i++)
        {
            TEMPLATE(T, poly_sub) (xqi, xqi, x, ctx);
            TEMPLATE(T, poly_gcd) (g_i, xqi, f, ctx);
            TEMPLATE(T, poly_add) (xqi, xqi, x, ctx);
            if (!TEMPLATE(T, poly_is_one) (g_i, ctx))
            {
                restart = 1;
                break;
            }
            TEMPLATE(T, poly_compose_mod_brent_kung_preinv) (xqi, xqi, xq, f,
                                                             finv, ctx);

        }
        if (!restart)
        {
            break;
        }
    }

    TEMPLATE(T, poly_clear) (x, ctx);
    TEMPLATE(T, poly_clear) (xq, ctx);
    TEMPLATE(T, poly_clear) (xqi, ctx);
    TEMPLATE(T, poly_clear) (g_i, ctx);
    TEMPLATE(T, poly_clear) (finv, ctx);
    fmpz_clear(q);
}


#endif
