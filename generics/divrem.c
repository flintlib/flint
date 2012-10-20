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

static int
euclidean_test(elem_srcptr X, elem_srcptr L, const ring_t ring)
{
    switch (ring->type)
    {
        case TYPE_FMPZ:
            return fmpz_cmpabs(X, L) < 0;

        case TYPE_MOD:
            return elem_is_zero(X, ring);

        case TYPE_POLY:
            return ((const elem_poly_struct *) X)->length < ((const elem_poly_struct *) L)->length;

        case TYPE_LIMB:
        default:
            NOT_IMPLEMENTED("euclidean test", ring); 
    }
}



void _elem_poly_divrem(elem_ptr Q, elem_ptr R, elem_srcptr A,
    long lenA, elem_srcptr B, long lenB, const ring_t ring)
{
    elem_ptr tmp;
    elem_srcptr leadB;
    long i, iQ, iR, size;

    size = ring->size;
    leadB = SRC_INDEX(B, lenB - 1, size);
    tmp = flint_malloc(size);       /* should be tmp init */
    elem_init(tmp, ring);

    if (R != A)
    {
        _elem_vec_set(R, A, lenA, ring);
    }

    for (iQ = lenA - lenB, iR = lenA - 1; iQ >= 0; iQ--, iR--)
    {
        /* TODO: euclidean test */
        if (euclidean_test(INDEX(R, iR, size), leadB, ring))
        {
            elem_zero(INDEX(Q, iQ, size), ring);
        }
        else
        {
            /* TODO: this could use div, not writing out the remainder */
            elem_divrem(INDEX(Q, iQ, size), tmp, SRC_INDEX(R, iR, size), leadB, ring);

            for (i = 0; i < lenB; i++)
            {
                /* R[iQ + i] -= B[i] * Q[iQ] */
                elem_mul(tmp, INDEX(B, i, size), SRC_INDEX(Q, iQ, size), ring);
                elem_sub(INDEX(R, iQ + i, size), SRC_INDEX(R, iQ + i, size), tmp, ring);
            }
        }
    }

    elem_clear(tmp, ring);
    flint_free(tmp);
}

void
elem_poly_divrem(elem_poly_struct * Q, elem_poly_struct * R,
                      const elem_poly_struct * A, const elem_poly_struct * B, const ring_t ring)
{
    long lenA = A->length;
    long lenB = B->length;

    if (lenB == 0)
    {
        printf("Exception: division by zero in elem_poly_divrem\n");
        abort();
    }

    if (lenA < lenB)
    {
        elem_set(R, A, ring);
        elem_zero(Q, ring);
        return;
    }

    if (Q == A || Q == B || R == A || R == B)
    {
        printf("not implemented: aliasing in divrem\n");
        abort();
    }

    /* TODO: lenB - 1, if a field... */
    _elem_poly_fit_length(R, lenA, ring);
    _elem_poly_fit_length(Q, lenA - lenB + 1, ring);

    _elem_poly_divrem(Q->coeffs, R->coeffs, A->coeffs, lenA, B->coeffs, lenB, ring->parent);

    _elem_poly_set_length(Q, lenA - lenB + 1, ring);
    _elem_poly_set_length(R, lenA, ring);
    _elem_poly_normalise(Q, ring);
    _elem_poly_normalise(R, ring);
}

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
            {
                switch (RING_PARENT(ring)->type)
                {
                    case TYPE_LIMB:
                        {
                            mp_limb_t t;
                            t = n_invmod((*((mp_srcptr) B)), ring->nmod.n);
                            t = n_mulmod2_preinv(t, (*((mp_srcptr) A)), ring->nmod.n, ring->nmod.ninv);
                            *((mp_ptr) Q) = t;
                            *((mp_ptr) R) = 0;
                        }
                        break;

                    case TYPE_FMPZ:
                        {
                            fmpz_t t;
                            fmpz_init(t);
                            fmpz_invmod(t, B, RING_MODULUS(ring));
                            fmpz_mul(Q, t, A);
                            fmpz_mod(Q, Q, RING_MODULUS(ring));
                            fmpz_zero(R);
                            fmpz_clear(t);
                        }
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
