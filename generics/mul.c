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

        case TYPE_FRAC:
            elem_frac_mul(res, op1, op2, ring);
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
