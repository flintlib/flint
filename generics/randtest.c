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

void
elem_randtest(elem_t res, flint_rand_t state, const long * size, const ring_t ring)
{
    switch (ring->type)
    {
        case TYPE_FMPZ:
            fmpz_randtest(&res->z, state, size[0]);
            break;

        case TYPE_LIMB:
            res->n = n_randtest(state);
            break;

        case TYPE_POLY:
            {
                long i, len;
                elem_ptr ptr;

                len = n_randint(state, size[0]);
                _elem_poly_fit_length(res->poly, len, ring);
                ptr = res->poly->coeffs;

                for (i = 0; i < len; i++)
                    elem_randtest(ptr + i, state, size + 1, ring->parent);

                _elem_poly_set_length(res->poly, len, ring);
                _elem_poly_normalise(res->poly, ring);
            }
            break;

        case TYPE_MOD:
            {
                switch (RING_PARENT(ring)->type)
                {
                    case TYPE_FMPZ:
                        fmpz_randtest_mod(&res->z, state, &RING_MODULUS(ring)->z);
                        break;

                    case TYPE_LIMB:
                        res->n = n_randint(state, ring->nmod.n);
                        break;

                    default:
                        NOT_IMPLEMENTED("randtest (mod)", ring);
                }
            }
            break;

        default:
            NOT_IMPLEMENTED("randtest", ring);
    }

}

void
gen_randtest(gen_t res, flint_rand_t state, const long * size)
{
    elem_randtest(&res->elem, state, size, res->ring);
}
