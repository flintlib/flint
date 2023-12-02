/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

#include "ulong_extras.h"
int
TEMPLATE(T, poly_factor_equal_deg_prob) (TEMPLATE(T, poly_t) factor,
                                         flint_rand_t state,
                                         const TEMPLATE(T, poly_t) pol,
                                         slong d,
                                         const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, poly_t) a, b, c, polinv;
    TEMPLATE(T, t) t;
    fmpz_t exp, q;
    int res = 1;
    slong i, k;

    if (pol->length <= 1)
    {
        flint_throw(FLINT_ERROR, "Exception (" TEMPLATE_STR(T) "_poly_factor_equal_deg_prob): Input polynomial is linear.\n");
    }

    fmpz_init(q);
    TEMPLATE(T, ctx_order) (q, ctx);

    TEMPLATE(T, poly_init) (a, ctx);

    do
    {
        TEMPLATE(T, poly_randtest) (a, state, pol->length - 1, ctx);
    } while (a->length <= 1);

    TEMPLATE(T, poly_gcd) (factor, a, pol, ctx);

    if (factor->length != 1)
    {
        TEMPLATE(T, poly_clear) (a, ctx);
        fmpz_clear(q);
        return 1;
    }

    TEMPLATE(T, poly_init) (b, ctx);
    TEMPLATE(T, poly_init) (polinv, ctx);

    TEMPLATE(T, poly_reverse) (polinv, pol, pol->length, ctx);
    TEMPLATE(T, poly_inv_series_newton) (polinv, polinv, polinv->length, ctx);

    fmpz_init(exp);
    if (fmpz_cmp_ui(TEMPLATE(T, ctx_prime) (ctx), 2) > 0)
    {
        /* compute a^{(q^d-1)/2} rem pol */
        fmpz_pow_ui(exp, q, d);
        fmpz_sub_ui(exp, exp, 1);
        fmpz_fdiv_q_2exp(exp, exp, 1);

        TEMPLATE(T, poly_powmod_fmpz_sliding_preinv) (b, a, exp, 0, pol,
                                                      polinv, ctx);
    }
    else
    {
        /* compute b = (a^{2^{k*d-1}}+a^{2^{k*d-2}}+...+a^4+a^2+a) rem pol */
        k = d * TEMPLATE(T, ctx_degree) (ctx);  /* TODO: Handle overflow? */
        TEMPLATE(T, poly_rem) (b, a, pol, ctx);
        TEMPLATE(T, poly_init) (c, ctx);
        TEMPLATE(T, poly_set) (c, b, ctx);
        for (i = 1; i < k; i++)
        {
            /* c = a^{2^i} = (a^{2^{i-1}})^2 */
            TEMPLATE(T, poly_powmod_ui_binexp_preinv) (c, c, 2, pol, polinv,
                                                       ctx);
            TEMPLATE(T, poly_add) (b, b, c, ctx);
        }
        TEMPLATE(T, poly_rem) (b, b, pol, ctx);
        TEMPLATE(T, poly_clear) (c, ctx);
    }
    fmpz_clear(exp);

    TEMPLATE(T, init) (t, ctx);
    TEMPLATE(T, sub_one) (t, b->coeffs + 0, ctx);
    TEMPLATE(T, poly_set_coeff) (b, 0, t, ctx);
    TEMPLATE(T, clear) (t, ctx);

    TEMPLATE(T, poly_gcd) (factor, b, pol, ctx);

    if ((factor->length <= 1) || (factor->length == pol->length))
        res = 0;

    TEMPLATE(T, poly_clear) (a, ctx);
    TEMPLATE(T, poly_clear) (b, ctx);
    TEMPLATE(T, poly_clear) (polinv, ctx);
    fmpz_clear(q);

    return res;
}


#endif
