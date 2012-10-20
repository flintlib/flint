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
elem_poly_mul(elem_poly_struct * res, const elem_poly_struct * op1, const elem_poly_struct * op2, const ring_t ring)
{
    long rlen, len1, len2;

    len1 = op1->length;
    len2 = op2->length;

    if (len1 == 0 || len2 == 0)
    {
        elem_zero(res, ring);
        return;
    }

    rlen = len1 + len2 - 1;

    if (res == op1 || res == op2)
    {
        elem_poly_struct t[1];
        elem_init(t, ring);
        _elem_poly_fit_length(t, rlen, ring);

        if (len1 >= len2)
            _elem_poly_mul(t->coeffs, op1->coeffs, len1, op2->coeffs, len2, ring->parent);
        else
            _elem_poly_mul(t->coeffs, op2->coeffs, len2, op1->coeffs, len1, ring->parent);

        elem_poly_swap(res, t);
        elem_clear(t, ring);
    }
    else
    {
        _elem_poly_fit_length(res, rlen, ring);

        if (len1 >= len2)
            _elem_poly_mul(res->coeffs, op1->coeffs, len1, op2->coeffs, len2, ring->parent);
        else
            _elem_poly_mul(res->coeffs, op2->coeffs, len2, op1->coeffs, len1, ring->parent);
    }

    _elem_poly_set_length(res, rlen, ring);
    _elem_poly_normalise(res, ring);
}

void
elem_mul(elem_ptr res, elem_srcptr op1, elem_srcptr op2, const ring_t ring)
{
    switch (ring->type)
    {
        case TYPE_FMPZ:
            fmpz_mul(res, op1, op2);
            break;

        case TYPE_LIMB:
            *((mp_ptr) res) = (*((mp_srcptr) op1)) * (*((mp_srcptr) op2));
            break;

        case TYPE_POLY:
            elem_poly_mul(res, op1, op2, ring);
            break;

        case TYPE_MOD:
            {
                switch (RING_PARENT(ring)->type)
                {
                    case TYPE_LIMB:
                        *((mp_ptr) res) = n_mulmod2_preinv(*((mp_srcptr) op1),
                            *((mp_srcptr) op2),
                            ring->nmod.n, ring->nmod.ninv);
                        break;

                    case TYPE_FMPZ:
                        fmpz_mul(res, op1, op2);
                        fmpz_mod(res, res, RING_MODULUS(ring));
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
        elem_mul(z->elem, x->elem, y->elem, z->ring);
    }
    else
    {
        NOT_IMPLEMENTED("gen_mul coercing into ", z->ring);
    }
}
