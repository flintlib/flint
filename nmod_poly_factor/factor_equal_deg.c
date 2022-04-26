/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010, 2022 William Hart
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
nmod_poly_factor_equal_deg(nmod_poly_factor_t factors,
                           const nmod_poly_t pol, slong d)
{
    if (pol->length == d + 1)
    {
        nmod_poly_factor_insert(factors, pol, 1);
    }
    else
    {
        nmod_poly_t f, g;
        flint_rand_t state;

        nmod_poly_init_mod(f, pol->mod);

        flint_randinit(state);

        while (!nmod_poly_factor_equal_deg_prob(f, state, pol, d))
           ;

        flint_randclear(state);

        nmod_poly_init_mod(g, pol->mod);
        nmod_poly_div(g, pol, f);

        nmod_poly_factor_equal_deg(factors, f, d);
        nmod_poly_clear(f);
        nmod_poly_factor_equal_deg(factors, g, d);
        nmod_poly_clear(g);
    }
}
