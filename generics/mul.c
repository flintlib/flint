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
elem_mul(elem_t res, const elem_t op1, const elem_t op2, const ring_t ring)
{
    switch (ring->type)
    {
        case TYPE_FMPZ:
            fmpz_mul(&res->z, &op1->z, &op2->z);
            break;

        case TYPE_LIMB:
            res->n = op1->n * op2->n;
            break;

        case TYPE_POLY:
            {
                long rlen, len1, len2;

                len1 = op1->poly->length;
                len2 = op2->poly->length;

                if (len1 == 0 || len2 == 0)
                {
                    elem_zero(res, ring);
                    return;
                }

                rlen = len1 + len2 - 1;

                if (res == op1 || res == op2)
                {
                    elem_t t;
                    elem_init(t, ring);
                    _elem_poly_fit_length(t->poly, rlen, ring);

                    if (len1 >= len2)
                    {
                        _elem_poly_mul(t->poly->coeffs,
                            op1->poly->coeffs, len1,
                            op2->poly->coeffs, len2, ring->parent);
                    }
                    else
                    {
                        _elem_poly_mul(t->poly->coeffs,
                            op2->poly->coeffs, len2,
                            op1->poly->coeffs, len1, ring->parent);
                    }

                    elem_swap(res, t);
                    elem_clear(t, ring);
                }
                else
                {
                    _elem_poly_fit_length(res->poly, rlen, ring);

                    if (len1 >= len2)
                    {
                        _elem_poly_mul(res->poly->coeffs,
                            op1->poly->coeffs, len1,
                            op2->poly->coeffs, len2, ring->parent);
                    }
                    else
                    {
                        _elem_poly_mul(res->poly->coeffs,
                            op2->poly->coeffs, len2,
                            op1->poly->coeffs, len1, ring->parent);
                    }
                }

                _elem_poly_set_length(res->poly, rlen, ring);
                _elem_poly_normalise(res->poly, ring);
            }
            break;

        case TYPE_MOD:
            {
                switch (RING_PARENT(ring)->type)
                {
                    case TYPE_LIMB:
                        res->n = n_mulmod2_preinv(op1->n, op2->n,
                            ring->nmod.n, ring->nmod.ninv);
                        break;

                    case TYPE_FMPZ:
                        fmpz_mul(&res->z, &op1->z, &op2->z);
                        fmpz_mod(&res->z, &res->z, &RING_MODULUS(ring)->z);
                        break;

                    default:
                        NOT_IMPLEMENTED("mul (mod)", ring);
                }
            }
            break;

        default:
            NOT_IMPLEMENTED("mul", ring);
    }
}

void
gen_mul(gen_t z, const gen_t x, const gen_t y)
{
    if (x->ring == y->ring && x->ring == z->ring)
    {
        elem_mul(&z->elem, &x->elem, &y->elem, z->ring);
    }
    else
    {
        NOT_IMPLEMENTED("gen_mul coercing into ", z->ring);
    }
}
