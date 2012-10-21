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


static __inline__ void
elem_poly_randtest(elem_poly_struct * res, flint_rand_t state, const long * size, const ring_t ring)
{
    long i, len, _size = RING_PARENT(ring)->size;
    elem_ptr ptr;

    len = n_randint(state, size[0]);
    _elem_poly_fit_length(res, len, ring);
    ptr = res->coeffs;

    for (i = 0; i < len; i++)
        elem_randtest(INDEX(ptr, i, _size), state, size + 1, ring->parent);

    _elem_poly_set_length(res, len, ring);
    _elem_poly_normalise(res, ring);
}


void
elem_randtest(elem_ptr res, flint_rand_t state, const long * size, const ring_t ring)
{
    switch (ring->type)
    {
        case TYPE_FMPZ:
            fmpz_randtest(res, state, size[0]);
            break;

        case TYPE_LIMB:
            *((mp_ptr) res) = n_randtest(state);
            break;

        case TYPE_POLY:
            elem_poly_randtest(res, state, size, ring);
            break;

        case TYPE_MOD:
            {
                switch (RING_PARENT(ring)->type)
                {
                    case TYPE_FMPZ:
                        fmpz_randtest_mod(res, state, RING_MODULUS(ring));
                        break;

                    case TYPE_LIMB:
                        *((mp_ptr) res) = n_randint(state, ring->nmod.n);
                        break;

                    default:
                        NOT_IMPLEMENTED("randtest (mod)", ring);
                }
            }
            break;

        case TYPE_FRAC:
            elem_randtest(NUMER(res, ring), state, size, ring->numer);
            elem_randtest(DENOM(res, ring), state, size, ring->denom);
            if (elem_is_zero(DENOM(res, ring), ring->denom))
                elem_one(DENOM(res, ring), ring->denom);
            elem_frac_canonicalise(res, ring);
            break;

        default:
            NOT_IMPLEMENTED("randtest", ring);
    }

}

void
gen_randtest(gen_t res, flint_rand_t state, const long * size)
{
    elem_randtest(res->elem, state, size, res->ring);
}
