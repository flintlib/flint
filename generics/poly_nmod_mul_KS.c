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

long ring_get_vars(const ring_t ring)
{
    if (ring->type == TYPE_POLY)
        return ring_get_vars(RING_PARENT(ring)) + 1;
    else
        return 0;
}

mp_limb_t ring_get_modulus_ui(const ring_t ring)
{
    if (ring->type == TYPE_POLY)
        return ring_get_modulus_ui(RING_PARENT(ring));
    else
        return ring->nmod.n;
}

void
_elem_poly_nmod_bit_pack(mp_ptr xp, elem_srcptr x, long xlen, const long * stride, long bits, const ring_t ring)
{
    if (ring->type != TYPE_POLY)
    {
        _nmod_poly_bit_pack(xp, (mp_srcptr) x, xlen, bits);
    }
    else
    {
        long i;
        const elem_poly_struct * polys = x;

        for (i = 0; i < xlen; i++)
        {
            _elem_poly_nmod_bit_pack(xp + stride[0] * i,
                (polys + i)->coeffs,
                (polys + i)->length, stride + 1, bits, RING_PARENT(ring));
        }
    }
}

/* TODO: this currently allocates space for leading zeros, which is madness */
void
_elem_poly_nmod_bit_unpack(elem_ptr z, mp_srcptr zp,
    const long * zbounds, const long * stride, long bits, const ring_t ring)
{
    if (ring->type != TYPE_POLY)
    {
        _nmod_poly_bit_unpack((mp_ptr) z, zbounds[0], zp, bits, ring->nmod);
    }
    else
    {
        long i;
        elem_poly_struct * x = z;

        for (i = 0; i < zbounds[0]; i++)
        {
            elem_poly_fit_length(x + i, zbounds[1], ring);
            _elem_poly_nmod_bit_unpack((x + i)->coeffs,
                zp + stride[0] * i, zbounds + 1, stride + 1, bits, RING_PARENT(ring));
            (x + i)->length = zbounds[1];
            elem_poly_normalise(x + i, ring);
        }
    }
}

void
_elem_poly_nmod_mul_KS(elem_ptr z, elem_srcptr x, long xlen, elem_srcptr y, long ylen, const ring_t ring)
{
    mp_ptr zp, xp, yp;
    mp_limb_t m;
    long i, maxterms, bits, vars;
    long *stride, *xbounds, *ybounds, *zbounds;
    long xlimbs, ylimbs;

    if (ring->type != TYPE_POLY)
    {
        _nmod_poly_mul_KS(z, x, xlen, y, ylen, 0, ring->nmod);
        return;
    }

    /* variables including the top level polynomials */
    vars = ring_get_vars(ring) + 1;
    m = ring_get_modulus_ui(ring);

    xbounds = flint_calloc(vars, sizeof(long));
    ybounds = flint_calloc(vars, sizeof(long));

    _elem_poly_rect_hull(xbounds, x, xlen, ring);
    _elem_poly_rect_hull(ybounds, y, ylen, ring);

    stride = flint_calloc(vars, sizeof(long));
    zbounds = flint_calloc(vars, sizeof(long));

    for (i = 0; i < vars; i++)
        zbounds[i] = (xbounds[i] + ybounds[i] - 1);

    /* bound coefficients.
    TODO: compute max bits instead of using modulus.
    TODO: can overflow if polynomial is extremely sparse (but then
    we shouldn't be using KS...) */
    maxterms = 1;
    for (i = 0; i < vars; i++)
        maxterms *= FLINT_MIN(xbounds[i], ybounds[i]);
    bits = FLINT_BIT_COUNT(maxterms) + 2*FLINT_BIT_COUNT(m - 1);

    /* number of limbs at the bottom level */
    stride[vars-1] = (zbounds[vars-1] * bits + FLINT_BITS - 1) / FLINT_BITS;
    for (i = vars - 2; i >= 0; i--)
        stride[i] = stride[i + 1] * zbounds[i];

    xlimbs = stride[1] * xbounds[0];
    ylimbs = stride[1] * ybounds[0];

    xp = flint_calloc(xlimbs + 1, sizeof(mp_limb_t));
    yp = flint_calloc(ylimbs + 1, sizeof(mp_limb_t));
    zp = flint_calloc(xlimbs + ylimbs + 1, sizeof(mp_limb_t));

    _elem_poly_nmod_bit_pack(xp, x, xbounds[0], stride + 1, bits, ring);
    _elem_poly_nmod_bit_pack(yp, y, ybounds[0], stride + 1, bits, ring);

    if (xlimbs < ylimbs)
        mpn_mul(zp, yp, ylimbs, xp, xlimbs);
    else
        mpn_mul(zp, xp, xlimbs, yp, ylimbs);

    _elem_poly_nmod_bit_unpack(z, zp, zbounds, stride + 1, bits, ring);

    flint_free(xp);
    flint_free(yp);
    flint_free(zp);
    flint_free(stride);
    flint_free(xbounds);
    flint_free(ybounds);
    flint_free(zbounds);
}

void
elem_poly_nmod_mul_KS(elem_poly_t res, const elem_poly_t op1, const elem_poly_t op2, const ring_t ring)
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

    elem_poly_fit_length(res, rlen, ring);

    if (len1 >= len2)
        _elem_poly_nmod_mul_KS(res->coeffs, op1->coeffs, len1, op2->coeffs, len2, ring->parent);
    else
        _elem_poly_nmod_mul_KS(res->coeffs, op2->coeffs, len2, op1->coeffs, len1, ring->parent);

    elem_poly_set_length(res, rlen, ring);
    elem_poly_normalise(res, ring);
}

