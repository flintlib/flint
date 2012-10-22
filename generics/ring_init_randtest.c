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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "generics.h"

int
ring_init_randtest(ring_t * R, flint_rand_t state, int maxdepth)
{
    int depth, frac_pos;

    maxdepth = 1 + n_randint(state, maxdepth);

    if (n_randint(state, 2))
    {
        ring_init_fmpz(R[0]);
        frac_pos = 1;
    }
    else
    {
        ring_init_limb(R[0]);
        frac_pos = 0;
    }

    if (maxdepth < 2)
        return 1;

    depth = 1;

    /* mod p in the base ring */
    if (n_randint(state, 2))
    {
        gen_t prime;
        gen_init(prime, R[0]);
        gen_set_si(prime, 7);   /* TODO: set_ui, random prime, and TODO: provide a way to free this... */
        ring_init_mod(R[1], R[0], prime->elem);
        frac_pos = 0;
        depth++;
    }

    /* else consider inserting fractions */
    if (frac_pos != 0)
        frac_pos = n_randint(state, maxdepth);

    while (depth < maxdepth)
    {
        if (depth == frac_pos)
        {
            ring_init_frac(R[depth], R[depth-1], R[0]);
        }
        else
        {
            ring_init_poly(R[depth], R[depth-1]);
        }

        depth++;
    }

    return depth;
}

