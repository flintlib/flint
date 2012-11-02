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

static int is_nmod_poly(const ring_t ring)
{
    if (ring->type == TYPE_POLY)
        return is_nmod_poly(RING_PARENT(ring));
    return ring->type == TYPE_MOD && RING_PARENT(ring)->type == TYPE_LIMB;
}

void
_elem_poly_mul(elem_ptr res, elem_srcptr poly1, long len1,
    elem_srcptr poly2, long len2, const ring_t ring)
{
    if (is_nmod_poly(ring))
    {
        _elem_poly_nmod_mul(res, poly1, len1, poly2, len2, ring);
        return;
    }

    if (ring->type == TYPE_FMPZ && 0)
    {
        _fmpz_poly_mul(res, poly1, len1, poly2, len2);
        return;
    }

    _elem_poly_mul_classical(res, poly1, len1, poly2, len2, ring);
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

