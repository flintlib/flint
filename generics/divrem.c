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

    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "generics.h"

void
elem_divrem(elem_ptr Q, elem_ptr R, elem_srcptr A, elem_srcptr B, const ring_t ring)
{
    switch (ring->type)
    {
        case TYPE_FMPZ:
            fmpz_fdiv_qr(Q, R, A, B);
            break;

        case TYPE_LIMB:
            {
                mp_limb_t q, r;
                q = (*((mp_srcptr) A)) / (*((mp_srcptr) B));
                r = (*((mp_srcptr) A)) % (*((mp_srcptr) B));
                *((mp_ptr) Q) = q;
                *((mp_ptr) R) = r;
            }
            break;

        case TYPE_POLY:
            elem_poly_divrem(Q, R, A, B, ring);
            break;

        case TYPE_MOD:
            elem_divexact(Q, A, B, ring);
            elem_zero(R, ring);
            break;

        case TYPE_FRAC:
            /* assume that it's a field */
            if (RING_NUMER(ring) == RING_DENOM(ring))
            {
                elem_divexact(Q, A, B, ring);
                elem_zero(R, ring);
                break;
            }

        default:
            NOT_IMPLEMENTED("divrem", ring);
    }
}

void
gen_divrem(gen_t Q, gen_t R, const gen_t A, const gen_t B)
{
    if (Q->ring == A->ring && A->ring == B->ring && R->ring == A->ring)
    {
        elem_divrem(Q->elem, R->elem, A->elem, B->elem, Q->ring);
    }
    else
    {
        NOT_IMPLEMENTED("gen_divrem coercing into ", Q->ring);
    }
}
