/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "nmod_poly_factor.h"
#include "ulong_extras.h"
#include "nmod_poly_factor_gr.h"

/* Note: unlike gr_poly_is_irreducible, this module considers constants
   (and the zero polynomial) irreducible. */
static int
_nmod_poly_irreducible_gr(const nmod_poly_t f, int ddf)
{
    gr_ctx_t ctx;
    gr_poly_t P;
    nmod_t mod = f->mod;
    truth_t res;

    if (f->length <= 2)
        return 1;

    _gr_ctx_init_nmod(ctx, &mod);
    GR_MUST_SUCCEED(gr_ctx_set_is_field(ctx, T_TRUE));
    NMOD_POLY_AS_GR(P, f);

    res = ddf ? gr_poly_is_irreducible_ddf(P, ctx) : gr_poly_is_irreducible(P, ctx);

    gr_ctx_clear(ctx);

    if (res == T_UNKNOWN)
        flint_throw(FLINT_ERROR, "nmod_poly_is_irreducible: unable to decide\n");

    return (res == T_TRUE);
}

int nmod_poly_is_irreducible_ddf(const nmod_poly_t poly)
{
    return _nmod_poly_irreducible_gr(poly, 1);
}

int
nmod_poly_is_irreducible(const nmod_poly_t f)
{
    return _nmod_poly_irreducible_gr(f, 0);
}

static void
nmod_poly_powpowmod(nmod_poly_t res, const nmod_poly_t pol,
                                    ulong exp, ulong exp2, const nmod_poly_t f)
{
    nmod_poly_t pow;
    ulong i;

    nmod_poly_init_mod(pow, f->mod);

    nmod_poly_powmod_ui_binexp(pow, pol, exp, f);

    nmod_poly_set(res, pow);

    if (!nmod_poly_equal(pow, pol))
        for (i = 1; i < exp2; i++)
            nmod_poly_powmod_ui_binexp(res, res, exp, f);

    nmod_poly_clear(pow);
}

int
nmod_poly_is_irreducible_rabin(const nmod_poly_t f)
{
    if (nmod_poly_length(f) > 2)
    {
        const ulong p = nmod_poly_modulus(f);
        const slong n     = nmod_poly_degree(f);
        nmod_poly_t a, x, x_p;

        nmod_poly_init(a, p);
        nmod_poly_init(x, p);
        nmod_poly_init(x_p, p);

	    nmod_poly_set_coeff_ui(x, 1, 1);

        /* Compute x^q mod f */
        nmod_poly_powpowmod(x_p, x, p, n, f);

	    if (!nmod_poly_is_zero(x_p))
            nmod_poly_make_monic(x_p, x_p);

        /* Now do the irreducibility test */
        if (!nmod_poly_equal(x_p, x))
        {
            nmod_poly_clear(a);
            nmod_poly_clear(x);
            nmod_poly_clear(x_p);

	    return 0;
        } else
        {
            n_factor_t factors;
            slong i;

            n_factor_init(&factors);

	        n_factor(&factors, n, 1);

            for (i = 0; i < factors.num; i++)
            {
                nmod_poly_powpowmod(a, x, p, n/factors.p[i], f);
                nmod_poly_sub(a, a, x);

                if (!nmod_poly_is_zero(a))
                    nmod_poly_make_monic(a, a);

                nmod_poly_gcd(a, a, f);

                if (a->length != 1)
                {
                    nmod_poly_clear(a);
                    nmod_poly_clear(x);
                    nmod_poly_clear(x_p);

		    return 0;
                }
            }
        }

        nmod_poly_clear(a);
        nmod_poly_clear(x);
        nmod_poly_clear(x_p);
    }

    return 1;
}
