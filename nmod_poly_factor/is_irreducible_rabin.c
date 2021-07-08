/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "ulong_extras.h"

void
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
        const mp_limb_t p = nmod_poly_modulus(f);
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
