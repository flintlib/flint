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
elem_add(elem_t res, const elem_t op1, const elem_t op2, const ring_t ring)
{
    switch (ring->type)
    {
        case TYPE_FMPZ:
            fmpz_add(&res->z, &op1->z, &op2->z);
            break;

        case TYPE_LIMB:
            res->n = op1->n + op2->n;
            break;

        case TYPE_POLY:
            {
                long max = FLINT_MAX(op1->poly->length, op2->poly->length);

                _elem_poly_fit_length(res->poly, max, ring);

                _elem_poly_add(res->poly->coeffs,
                    op1->poly->coeffs, op1->poly->length,
                    op2->poly->coeffs, op2->poly->length, ring->parent);

                _elem_poly_set_length(res->poly, max, ring);
                _elem_poly_normalise(res->poly, ring);
            }
            break;

        case TYPE_MOD:
            {
                switch (RING_PARENT(ring)->type)
                {
                    case TYPE_LIMB:
                        res->n = n_addmod(op1->n, op2->n, ring->nmod.n);
                        break;

                    case TYPE_FMPZ:
                        fmpz_add(&res->z, &op1->z, &op2->z);
                        if (fmpz_cmpabs(&res->z, &RING_MODULUS(ring)->z) >= 0)
                            fmpz_sub(&res->z, &res->z, &RING_MODULUS(ring)->z);
                        break;

                    default:
                        NOT_IMPLEMENTED("add (mod)", ring);
                }
            }
            break;

        default:
            NOT_IMPLEMENTED("add", ring);
    }
}

void
gen_add(gen_t z, const gen_t x, const gen_t y)
{
    if (x->ring == y->ring && x->ring == z->ring)
    {
        elem_add(&z->elem, &x->elem, &y->elem, z->ring);
    }
    else
    {
        NOT_IMPLEMENTED("gen_add coercing into ", z->ring);
    }
}
