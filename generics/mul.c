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
elem_complex_mul(elem_ptr z, elem_srcptr x, elem_srcptr y, const ring_t ring)
{
    elem_srcptr a, b, c, d;
    elem_ptr e, f;
    elem_ptr t, u, v;
    const ring_struct * rring = RING_PARENT(ring);

    a = REALPART(x, ring);
    b = IMAGPART(x, ring);

    c = REALPART(y, ring);
    d = IMAGPART(y, ring);

    e = REALPART(z, ring);
    f = IMAGPART(z, ring);

    ELEM_TMP_INIT(t, rring);
    ELEM_TMP_INIT(u, rring);
    ELEM_TMP_INIT(v, rring);

    elem_add(t, a, b, rring);
    elem_add(u, c, d, rring);
    elem_mul(v, t, u, rring);

    elem_mul(t, a, c, rring);
    elem_mul(u, b, d, rring);

    elem_sub(e, t, u, rring);
    elem_sub(f, v, t, rring);
    elem_sub(f, f, u, rring);

    ELEM_TMP_CLEAR(t, rring);
    ELEM_TMP_CLEAR(u, rring);
    ELEM_TMP_CLEAR(v, rring);
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

        case TYPE_FRAC:
            elem_frac_mul(res, op1, op2, ring);
            break;

        case TYPE_COMPLEX:
            elem_complex_mul(res, op1, op2, ring);
            break;

        default:
            NOT_IMPLEMENTED("mul", ring);
    }
}

