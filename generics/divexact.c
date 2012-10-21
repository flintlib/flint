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

static void
elem_poly_divexact(elem_poly_struct * Q,
                      const elem_poly_struct * A, const elem_poly_struct * B, const ring_t ring)
{
    elem_poly_struct R[1];
    elem_init(R, ring);
    elem_poly_divrem(Q, R, A, B, ring);
    elem_clear(R, ring);
}

void
elem_divexact(elem_ptr Q, elem_srcptr A, elem_srcptr B, const ring_t ring)
{
    switch (ring->type)
    {
        case TYPE_FMPZ:
            fmpz_divexact(Q, A, B);
            break;

        case TYPE_POLY:
            elem_poly_divexact(Q, A, B, ring);
            break;

        case TYPE_MOD:
            {
                switch (RING_PARENT(ring)->type)
                {
                    case TYPE_LIMB:
                        {
                            mp_limb_t t;
                            t = n_invmod((*((mp_srcptr) B)), ring->nmod.n);
                            t = n_mulmod2_preinv(t, (*((mp_srcptr) A)), ring->nmod.n, ring->nmod.ninv);
                            *((mp_ptr) Q) = t;
                        }
                        break;

                    case TYPE_FMPZ:
                        {
                            fmpz_t t;
                            fmpz_init(t);
                            fmpz_invmod(t, B, RING_MODULUS(ring));
                            fmpz_mul(Q, t, A);
                            fmpz_mod(Q, Q, RING_MODULUS(ring));
                            fmpz_clear(t);
                        }
                        break;

                    default:
                        NOT_IMPLEMENTED("divexact", ring);
                }
            }
            break;

        default:
            NOT_IMPLEMENTED("divexact", ring);
    }
}

void
gen_divexact(gen_t Q, const gen_t A, const gen_t B)
{
    if (Q->ring == A->ring && A->ring == B->ring)
    {
        elem_divexact(Q->elem, A->elem, B->elem, Q->ring);
    }
    else
    {
        NOT_IMPLEMENTED("divexact coercing into ", Q->ring);
    }
}
