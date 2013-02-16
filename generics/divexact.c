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
    Copyright (C) 2013 William Hart

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
elem_complex_divexact(elem_ptr q, elem_srcptr x, elem_srcptr y, const ring_t ring)
{
    elem_srcptr a, b, c, d;
    elem_ptr e, f, t, u, v, w;
    const ring_struct * rring = RING_PARENT(ring);

    a = REALPART(x, ring);
    b = IMAGPART(x, ring);
    c = REALPART(y, ring);
    d = IMAGPART(y, ring);
    e = REALPART(q, ring);
    f = IMAGPART(q, ring);

    ELEM_TMP_INIT(t, rring);
    ELEM_TMP_INIT(u, rring);
    ELEM_TMP_INIT(v, rring);
    ELEM_TMP_INIT(w, rring);

    elem_mul(t, c, c, rring);
    elem_mul(u, d, d, rring);
    elem_add(v, t, u, rring);

    /* re = (ac + bd) / (c^2 + d^2) */
    elem_mul(t, a, c, rring);
    elem_mul(u, b, d, rring);
    elem_add(t, t, u, rring);

    /* im = (bc - ad) / (c^2 + d^2) */
    elem_mul(u, b, c, rring);
    elem_mul(w, a, d, rring);
    elem_sub(u, u, w, rring);

    elem_divexact(e, t, v, rring);
    elem_divexact(f, u, v, rring);

    ELEM_TMP_CLEAR(t, rring);
    ELEM_TMP_CLEAR(u, rring);
    ELEM_TMP_CLEAR(v, rring);
    ELEM_TMP_CLEAR(w, rring);
}

void
elem_divexact(elem_ptr Q, elem_srcptr A, elem_srcptr B, const ring_t ring)
{
    switch (ring->type)
    {
        case TYPE_FMPZ:
            fmpz_divexact(Q, A, B);
            break;

        case TYPE_MPZ:
            mpz_divexact(Q, A, B);
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

        case TYPE_COMPLEX:
            elem_complex_divexact(Q, A, B, ring);
            break;

        case TYPE_FRAC:
            if (RING_NUMER(ring) == RING_DENOM(ring))
            {
                elem_frac_div(Q, A, B, ring);
                break;
            }
            /* fall through */

        default:
            {
                elem_ptr tmp;
                ELEM_TMP_INIT(tmp, ring);
                elem_divrem(Q, tmp, A, B, ring);
                if (!elem_is_zero(tmp, ring))
                {
                    printf("divexact: divrem gave nonzero remainder\n");
                    abort();
                }
                ELEM_TMP_CLEAR(tmp, ring);
                break;
            }
    }
}

