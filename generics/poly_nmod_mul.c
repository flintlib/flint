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

long _elem_poly_rect_hull(long * bounds, const elem_poly_struct * coeffs, long len, const ring_t ring);

static long
ring_get_vars(const ring_t ring)
{
    if (ring->type == TYPE_POLY)
        return ring_get_vars(RING_PARENT(ring)) + 1;
    else
        return 0;
}


#define FLAT_MAX_LENGTH 25
#define KS_MIN_SIZE 50
#define KS_MIN_DENSITY 0.1

void
_elem_poly_nmod_mul(elem_ptr z, elem_srcptr x, long xlen, elem_srcptr y, long ylen, const ring_t ring)
{
    long i, vars, xterms, yterms, * xbounds, * ybounds, xmax, ymax;
    double xrect, yrect;

    if (ring->type != TYPE_POLY)
    {
        _nmod_poly_mul(z, x, xlen, y, ylen, ring->nmod);
        return;
    }

    vars = ring_get_vars(ring) + 1;

    xbounds = flint_calloc(vars, sizeof(long));
    ybounds = flint_calloc(vars, sizeof(long));

    xterms = _elem_poly_rect_hull(xbounds, x, xlen, ring);
    yterms = _elem_poly_rect_hull(ybounds, y, ylen, ring);

    xrect = xmax = xbounds[0];
    for (i = 1; i < vars; i++)
    {
        xrect *= xbounds[i];
        xmax = FLINT_MAX(xmax, xbounds[i]);
    }

    yrect = ymax = ybounds[0];
    for (i = 1; i < vars; i++)
    {
        yrect *= ybounds[i];
        ymax = FLINT_MAX(ymax, ybounds[i]);
    }

    if ((xterms > xrect * KS_MIN_DENSITY) &&
        (yterms > yrect * KS_MIN_DENSITY) &&
        (xmax > KS_MIN_SIZE) &&
        (ymax > KS_MIN_SIZE))
    {
        _elem_poly_nmod_mul_KS(z, x, xlen, y, ylen, ring);
    }
    else if (ymax < FLAT_MAX_LENGTH && xmax < FLAT_MAX_LENGTH && vars == 2)
    {
        _elem_poly_nmod_mul_flat(z, x, xlen, y, ylen, ring);
    }
    else
    {
        _elem_poly_mul_classical(z, x, xlen, y, ylen, ring);
    }

    flint_free(xbounds);
    flint_free(ybounds);
}

void
elem_poly_nmod_mul(elem_poly_t res, const elem_poly_t op1, const elem_poly_t op2, const ring_t ring)
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
            _elem_poly_nmod_mul(t->coeffs, op1->coeffs, len1, op2->coeffs, len2, ring->parent);
        else
            _elem_poly_nmod_mul(t->coeffs, op2->coeffs, len2, op1->coeffs, len1, ring->parent);

        elem_poly_swap(res, t);
        elem_clear(t, ring);
    }
    else
    {
        elem_poly_fit_length(res, rlen, ring);

        if (len1 >= len2)
            _elem_poly_nmod_mul(res->coeffs, op1->coeffs, len1, op2->coeffs, len2, ring->parent);
        else
            _elem_poly_nmod_mul(res->coeffs, op2->coeffs, len2, op1->coeffs, len1, ring->parent);
    }

    elem_poly_set_length(res, rlen, ring);
    elem_poly_normalise(res, ring);
}

