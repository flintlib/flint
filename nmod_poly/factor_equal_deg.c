/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include "nmod_poly.h"
#include "ulong_extras.h"

void
nmod_poly_factor_equal_deg(nmod_poly_factor_t factors,
                           const nmod_poly_t pol, len_t d)
{
    if (pol->length == d + 1)
    {
        nmod_poly_factor_insert(factors, pol, 1);
    }
    else
    {
        nmod_poly_t f, g;
        flint_rand_t state;

        nmod_poly_init_preinv(f, pol->mod.n, pol->mod.ninv);

        flint_randinit(state);

        while (!nmod_poly_factor_equal_deg_prob(f, state, pol, d)) {};

        flint_randclear(state);

        nmod_poly_init_preinv(g, pol->mod.n, pol->mod.ninv);
        nmod_poly_div(g, pol, f);

        nmod_poly_factor_equal_deg(factors, f, d);
        nmod_poly_clear(f);
        nmod_poly_factor_equal_deg(factors, g, d);
        nmod_poly_clear(g);
    }
}
