/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"
#include "ulong_extras.h"

int
fmpz_mod_poly_factor_equal_deg_prob(fmpz_mod_poly_t factor, flint_rand_t state,
                  const fmpz_mod_poly_t pol, slong d, const fmpz_mod_ctx_t ctx)
{
    const fmpz * p = fmpz_mod_ctx_modulus(ctx);
    fmpz_mod_poly_t a, b, c, polinv;
    fmpz_t exp, t;
    int res = 1;
    slong i;

    if (pol->length <= 1)
    {
        flint_printf("Exception (fmpz_mod_poly_factor_equal_deg_prob): \n");
        flint_printf("Input polynomial is linear.\n");
        flint_abort();
    }

    fmpz_mod_poly_init(a, ctx);

    do
    {
        fmpz_mod_poly_randtest(a, state, pol->length - 1, ctx);
    } while (a->length <= 1);

    fmpz_mod_poly_gcd(factor, a, pol, ctx);

    if (factor->length != 1)
    {
        fmpz_mod_poly_clear(a, ctx);
        return 1;
    }

    fmpz_mod_poly_init(b, ctx);
    fmpz_mod_poly_init(polinv, ctx);

    fmpz_mod_poly_reverse(polinv, pol, pol->length, ctx);
    fmpz_mod_poly_inv_series_newton(polinv, polinv, polinv->length, ctx);

    fmpz_init(exp);
    if (fmpz_cmp_ui(p, 2) > 0)
    {
        /* compute a^{(p^d-1)/2} rem pol */
        fmpz_pow_ui(exp, p, d);
        fmpz_sub_ui(exp, exp, 1);
        fmpz_fdiv_q_2exp(exp, exp, 1);

        fmpz_mod_poly_powmod_fmpz_binexp_preinv(b, a, exp, pol, polinv, ctx);
    }
    else
    {
        /* compute b = (a^{2^{d-1}}+a^{2^{d-2}}+...+a^4+a^2+a) rem pol */
        fmpz_mod_poly_rem(b, a, pol, ctx);
        fmpz_mod_poly_init(c, ctx);
        fmpz_mod_poly_set(c, b, ctx);
        for (i = 1; i < d; i++)
        {
            /* c = a^{2^i} = (a^{2^{i-1}})^2 */
            fmpz_mod_poly_powmod_ui_binexp_preinv(c, c, 2, pol, polinv, ctx);
            fmpz_mod_poly_add(b, b, c, ctx);
        }
        fmpz_mod_poly_rem(b, b, pol, ctx);
        fmpz_mod_poly_clear(c, ctx);
    }
    fmpz_clear(exp);

    fmpz_init(t);
    fmpz_sub_ui(t, &(b->coeffs[0]), 1);
    fmpz_mod(t, t, p);
    fmpz_mod_poly_set_coeff_fmpz(b, 0, t, ctx);
    fmpz_clear(t);

    fmpz_mod_poly_gcd(factor, b, pol, ctx);

    if ((factor->length <= 1) || (factor->length == pol->length))
        res = 0;

    fmpz_mod_poly_clear(a, ctx);
    fmpz_mod_poly_clear(b, ctx);
    fmpz_mod_poly_clear(polinv, ctx);

    return res;
}
