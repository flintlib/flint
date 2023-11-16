/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"

int
ca_poly_factor_squarefree(ca_t c, ca_poly_vec_t fac, ulong * exp, const ca_poly_t F, ca_ctx_t ctx)
{
    ca_poly_t f, d, t1;
    ca_poly_t v, w, s;
    slong i;
    int success;

    if (F->length == 0)
    {
        ca_zero(c, ctx);
        ca_poly_vec_set_length(fac, 0, ctx);
        return 1;
    }

    if (!ca_poly_is_proper(F, ctx))
    {
        return 0;
    }

    ca_set(c, F->coeffs + F->length - 1, ctx);

    if (F->length == 1)
    {
        ca_poly_vec_set_length(fac, 0, ctx);
        return 1;
    }

    ca_poly_init(f, ctx);
    ca_poly_init(d, ctx);
    ca_poly_init(t1, ctx);
    ca_poly_init(v, ctx);
    ca_poly_init(w, ctx);
    ca_poly_init(s, ctx);

    success = 0;

    ca_poly_make_monic(f, F, ctx);
    ca_poly_derivative(t1, f, ctx);
    if (!ca_poly_gcd(d, f, t1, ctx))
        goto cleanup;

    ca_poly_vec_set_length(fac, 0, ctx);

    if (d->length == 1)
    {
        ca_poly_vec_append(fac, f, ctx);
        exp[fac->length - 1] = 1;
    }
    else
    {
        ca_poly_div(v, f, d, ctx);
        ca_poly_div(w, t1, d, ctx);

        /* invariant: v is monic, so we don't need to check if it is proper */
        for (i = 1; ; i++)
        {
            ca_poly_derivative(t1, v, ctx);
            ca_poly_sub(s, w, t1, ctx);

            if (!ca_poly_is_proper(s, ctx))
                goto cleanup;

            if (s->length == 0)
            {
                if (v->length > 1)
                {
                    ca_poly_vec_append(fac, v, ctx);
                    exp[fac->length - 1] = i;
                }
                break;
            }

            if (!ca_poly_gcd(d, v, s, ctx))
                goto cleanup;

            ca_poly_div(v, v, d, ctx);
            ca_poly_div(w, s, d, ctx);

            if (d->length > 1)
            {
                ca_poly_vec_append(fac, d, ctx);
                exp[fac->length - 1] = i;
            }
        }

    }

    success = 1;

cleanup:
    ca_poly_clear(f, ctx);
    ca_poly_clear(d, ctx);
    ca_poly_clear(t1, ctx);
    ca_poly_clear(v, ctx);
    ca_poly_clear(w, ctx);
    ca_poly_clear(s, ctx);

    return success;
}
