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
    Copyright (C) 2013 William Hart

******************************************************************************/

#include "generics.h"

void
elem_add(elem_ptr res, elem_srcptr op1, elem_srcptr op2, const ring_t ring)
{
    switch (ring->type)
    {
        case TYPE_FMPZ:
            fmpz_add(res, op1, op2);
            break;

        case TYPE_MPZ:
            mpz_add(res, op1, op2);
            break;

        case TYPE_LIMB:
            *((mp_ptr) res) = *((mp_srcptr) op1) + *((mp_srcptr) op2);
            break;

        case TYPE_POLY:
            elem_poly_add(res, op1, op2, ring);
            break;

        case TYPE_MOD:
            {
                switch (RING_PARENT(ring)->type)
                {
                    case TYPE_LIMB:
                        *((mp_ptr) res) = n_addmod(*((mp_srcptr) op1), *((mp_srcptr) op2), ring->nmod.n);
                        break;

                    case TYPE_FMPZ:
                        fmpz_add(res, op1, op2);
                        if (fmpz_cmpabs(res, RING_MODULUS(ring)) >= 0)
                            fmpz_sub(res, res, RING_MODULUS(ring));
                        break;

                    default:
                        NOT_IMPLEMENTED("add (mod)", ring);
                }
            }
            break;

        case TYPE_FRAC:
            elem_frac_add(res, op1, op2, ring);
            break;

        case TYPE_COMPLEX:
            elem_add(REALPART(res, ring), REALPART(op1, ring), REALPART(op2, ring), ring->parent);
            elem_add(IMAGPART(res, ring), IMAGPART(op1, ring), IMAGPART(op2, ring), ring->parent);
            break;

        default:
            NOT_IMPLEMENTED("add", ring);
    }
}

