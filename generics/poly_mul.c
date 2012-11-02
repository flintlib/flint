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
_elem_poly_mul(elem_ptr res, elem_srcptr poly1, long len1,
    elem_srcptr poly2, long len2, const ring_t ring)
{
    long size = ring->size;

    if (ring->type == TYPE_MOD && RING_PARENT(ring)->type == TYPE_LIMB)
    {
        _nmod_poly_mul(res, poly1, len1, poly2, len2, ring->nmod);
        return;
    }

    if (ring->type == TYPE_FMPZ && 0)
    {
        _fmpz_poly_mul(res, poly1, len1, poly2, len2);
        return;
    }

/*
    if (ring->type == TYPE_MOD && RING_PARENT(ring)->type == TYPE_FMPZ)
    {
        _fmpz_mod_poly_mul(res, poly1, len1, poly2, len2, RING_MODULUS(ring));
        return;
    }
*/

    if (len1 == 1 && len2 == 1)
    {
        elem_mul(res, poly1, poly2, ring);
    }
    else
    {
        long i;

        _elem_vec_scalar_mul(res, poly1, len1, poly2, ring);

        _elem_vec_scalar_mul(INDEX(res, len1, size),
                            SRC_INDEX(poly2, 1, size), len2 - 1,
                            SRC_INDEX(poly1, len1 - 1, size), ring);

        for (i = 0; i < len1 - 1; i++)
            _elem_vec_scalar_addmul(INDEX(res, i + 1, size),
                    SRC_INDEX(poly2, 1, size), len2 - 1,
                    SRC_INDEX(poly1, i, size), ring);
    }
}

void
elem_poly_mul(elem_poly_t res, const elem_poly_t op1, const elem_poly_t op2, const ring_t ring)
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
        elem_poly_t t;
        elem_init(t, ring);
        elem_poly_fit_length(t, rlen, ring);

        if (len1 >= len2)
            _elem_poly_mul(t->coeffs, op1->coeffs, len1, op2->coeffs, len2, ring->parent);
        else
            _elem_poly_mul(t->coeffs, op2->coeffs, len2, op1->coeffs, len1, ring->parent);

        elem_poly_swap(res, t);
        elem_clear(t, ring);
    }
    else
    {
        elem_poly_fit_length(res, rlen, ring);

        if (len1 >= len2)
            _elem_poly_mul(res->coeffs, op1->coeffs, len1, op2->coeffs, len2, ring->parent);
        else
            _elem_poly_mul(res->coeffs, op2->coeffs, len2, op1->coeffs, len1, ring->parent);
    }

    elem_poly_set_length(res, rlen, ring);
    elem_poly_normalise(res, ring);
}

