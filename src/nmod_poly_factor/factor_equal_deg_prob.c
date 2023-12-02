/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2022 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "gmpcompat.h"
#include "ulong_extras.h"
#include "nmod_poly.h"
#include "nmod_poly_factor.h"

int
nmod_poly_factor_equal_deg_prob(nmod_poly_t factor,
    flint_rand_t state, const nmod_poly_t pol, slong d)
{
    nmod_poly_t a, b, c, polinv;
    fmpz_t exp;
    int res = 1;
    slong i;

    if (pol->length <= 1)
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_factor_equal_deg_prob): Input polynomial is linear.\n");
    }

    nmod_poly_init_mod(a, pol->mod);

    do
    {
        nmod_poly_randtest(a, state, pol->length - 1);
    } while (a->length <= 1);

    nmod_poly_gcd(factor, a, pol);

    if (factor->length != 1)
    {
        nmod_poly_clear(a);
        return 1;
    }

    nmod_poly_init_mod(b, pol->mod);
    nmod_poly_init_mod(polinv, pol->mod);

    nmod_poly_reverse(polinv, pol, pol->length);
    nmod_poly_inv_series(polinv, polinv, polinv->length);

    fmpz_init(exp);

    if (pol->mod.n > 2)
    {
        /* compute a^{(p^d-1)/2} rem pol */
        fmpz_ui_pow_ui(exp, pol->mod.n, d);
        fmpz_sub_ui(exp, exp, 1);
        fmpz_tdiv_q_2exp(exp, exp, 1);

        nmod_poly_powmod_fmpz_binexp_preinv(b, a, exp, pol, polinv);
    }
    else
    {
        /* compute b = (a^{2^{d-1}}+a^{2^{d-2}}+...+a^4+a^2+a) rem pol */
        nmod_poly_rem(b, a, pol);
        nmod_poly_init_mod(c, pol->mod);
        nmod_poly_set(c, b);

        for (i = 1; i < d; i++)
        {
            /* c = a^{2^i} = (a^{2^{i-1}})^2 */
            nmod_poly_powmod_ui_binexp_preinv(c, c, 2, pol, polinv);
            nmod_poly_add(b, b, c);
        }

        nmod_poly_rem(b, b, pol);
        nmod_poly_clear(c);
    }

    fmpz_clear(exp);

    nmod_poly_set_coeff_ui(b, 0, n_submod(b->coeffs[0], 1, pol->mod.n));

    nmod_poly_gcd(factor, b, pol);

    if (factor->length <= 1 || factor->length == pol->length)
       res = 0;

    nmod_poly_clear(polinv);
    nmod_poly_clear(a);
    nmod_poly_clear(b);

    return res;
}
