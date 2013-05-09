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
    Copyright (C) 2012 Lina Kulakova

******************************************************************************/

#include "fmpz_mod_poly_factor.h"
#include "ulong_extras.h"

void
fmpz_mod_poly_factor_equal_deg(fmpz_mod_poly_factor_t factors,
                               const fmpz_mod_poly_t pol, len_t d)
{
    if (pol->length == d + 1)
    {
        fmpz_mod_poly_factor_insert(factors, pol, 1);
    }
    else
    {
        fmpz_mod_poly_t f, g, r;
        flint_rand_t state;

        fmpz_mod_poly_init(f, &pol->p);

        flint_randinit(state);

        while (!fmpz_mod_poly_factor_equal_deg_prob(f, state, pol, d))
        {
        };

        flint_randclear(state);

        fmpz_mod_poly_init(g, &pol->p);
        fmpz_mod_poly_init(r, &pol->p);
        fmpz_mod_poly_divrem(g, r, pol, f);
        fmpz_mod_poly_clear(r);

        fmpz_mod_poly_factor_equal_deg(factors, f, d);
        fmpz_mod_poly_clear(f);
        fmpz_mod_poly_factor_equal_deg(factors, g, d);
        fmpz_mod_poly_clear(g);
    }
}
