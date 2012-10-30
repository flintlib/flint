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
elem_neg(elem_ptr res, elem_srcptr src, const ring_t ring)
{
    switch (ring->type)
    {
        case TYPE_FMPZ:
            fmpz_neg(res, src);
            break;

        case TYPE_LIMB:
            *((mp_ptr) res) = -(*((mp_srcptr) src));
            break;

        case TYPE_MOD:
            switch (RING_PARENT(ring)->type)
            {
                case TYPE_LIMB:
                    *((mp_ptr) res) = n_negmod(*((mp_srcptr) src), ring->nmod.n);
                    break;

                case TYPE_FMPZ:
                    if (!fmpz_is_zero(src))
                        fmpz_sub(res, RING_MODULUS(ring), src);
                    else
                        fmpz_zero(res);
                    break;

                default:
                    NOT_IMPLEMENTED("neg (mod)", ring);
            }
            break;

        case TYPE_POLY:
            elem_poly_neg(res, src, ring);
            break;

        case TYPE_FRAC:
            elem_neg(NUMER(res, ring), NUMER(src, ring), RING_NUMER(ring));
            elem_set(DENOM(res, ring), DENOM(src, ring), RING_DENOM(ring));
            break;

        case TYPE_COMPLEX:
            elem_neg(REALPART(res, ring), REALPART(src, ring), ring->parent);
            elem_neg(IMAGPART(res, ring), IMAGPART(src, ring), ring->parent);
            break;

        default:
            NOT_IMPLEMENTED("neg", ring);
    }
}

void
gen_neg(gen_t y, const gen_t x)
{
    if (y->ring == x->ring)
    {
        elem_neg(y->elem, x->elem, y->ring);
    }
    else
    {
        NOT_IMPLEMENTED("gen_neg coercing into ", y->ring);
    }
}
