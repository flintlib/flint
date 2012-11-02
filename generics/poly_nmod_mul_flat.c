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

/*
Computes the maximum length in each variable, assuming that
the lengths vector has been pre-initialised to zero.
Returns the total length of all polynomials at the bottom level.
*/

long
_elem_poly_rect_hull(long * bounds, const elem_poly_struct * coeffs, long len, const ring_t ring)
{
    bounds[0] = FLINT_MAX(bounds[0], len);

    if (ring->type != TYPE_POLY)
    {
        return len;
    }
    else if (RING_PARENT(ring)->type != TYPE_POLY)
    {
        long total, clen, i;

        for (i = total = 0; i < len; i++)
        {
            clen = (coeffs + i)->length;
            bounds[1] = FLINT_MAX(bounds[1], clen);
            total += clen;
        }

        return total;
    }
    else
    {
        long total, i;

        for (i = total = 0; i < len; i++)
            total += _elem_poly_rect_hull(bounds + 1, (coeffs + i)->coeffs, (coeffs + i)->length, RING_PARENT(ring));

        return total;
    }
}

static void
fit_length(elem_poly_t poly, long len)
{
    long alloc = poly->alloc;

    if (len > alloc)
    {
        if (len < 2 * alloc)
            len = 2 * alloc;

        poly->coeffs = flint_realloc(poly->coeffs, len * sizeof(mp_limb_t));
        poly->alloc = len;
    }
}

static void
normalise(elem_poly_t poly)
{
    long i;
    mp_ptr ptr = poly->coeffs;
    for (i = poly->length - 1; (i >= 0) && ptr[i] == 0; i--);
    poly->length = i + 1;
}


void
_elem_poly_nmod_mul_flat_2d(elem_poly_struct * z,
        const elem_poly_struct * x, long xlen,
        const elem_poly_struct * y, long ylen,
        long xybound, long bits, const ring_t ring)
{
    long xi, yi, i, j, zlen, alen, blen, nlimbs;
    mp_limb_t a, b, c, d;
    mp_ptr PZ;
    mp_srcptr PX, PY;
    mp_ptr tmp, t;
    nmod_t mod = RING_PARENT(ring)->nmod;

    zlen = xlen + ylen - 1;

    if (bits <= FLINT_BITS)
        nlimbs = 1;
    else if (bits <= 2 * FLINT_BITS)
        nlimbs = 2;
    else
        nlimbs = 3;

    tmp = calloc(nlimbs * zlen * xybound, sizeof(mp_limb_t));

    for (i = 0; i < zlen; i++)
        z[i].length = 0;

    for (xi = 0; xi < xlen; xi++)
    {
        alen = x[xi].length;

        for (yi = 0; yi < ylen; yi++)
        {
            blen = y[yi].length;

            PX = (x + xi)->coeffs;
            PY = (y + yi)->coeffs;
            PZ = (z + xi + yi)->coeffs;

            if (nlimbs == 1)
                for (i = 0; i < alen; i++)
                    for (j = 0; j < blen; j++)
                        tmp[xybound*(xi+yi) + i+j] += PX[i] * PY[j];
            else if (nlimbs == 2)
                for (i = 0; i < alen; i++)
                    for (j = 0; j < blen; j++)
                    {
                        t = tmp + 2 * (xybound * (xi+yi) + i+j);
                        umul_ppmm(b, a, PX[i], PY[j]);
                        add_ssaaaa(t[1], t[0], t[1], t[0], b, a);
                    }
            else
                for (i = 0; i < alen; i++)
                    for (j = 0; j < blen; j++)
                    {
                        t = tmp + 3 * (xybound * (xi+yi) + i+j);
                        umul_ppmm(b, a, PX[i], PY[j]);
                        add_sssaaaaaa(t[2], t[1], t[0],
                                      t[2], t[1], t[0], 0, b, a);
                    }

            z[xi+yi].length = FLINT_MAX(z[xi+yi].length, alen + blen - 1);
        }
    }

    for (i = 0; i < zlen; i++)
    {
        fit_length(z + i, z[i].length);

        if (nlimbs == 1)
            _nmod_vec_reduce((z + i)->coeffs,
                (mp_srcptr) tmp + xybound * i, z[i].length, mod);
        else if (nlimbs == 2)
            for (j = 0; j < z[i].length; j++)
            {
                a = tmp[2 * (xybound * i + j)];
                b = tmp[2 * (xybound * i + j) + 1];
                c = n_ll_mod_preinv(b, a, mod.n, mod.ninv);
                ((mp_ptr) (z + i)->coeffs)[j] = c;
            }
        else
            for (j = 0; j < z[i].length; j++)
            {
                a = tmp[3 * (xybound * i + j)];
                b = tmp[3 * (xybound * i + j) + 1];
                c = tmp[3 * (xybound * i + j) + 2];
                d = n_lll_mod_preinv(c, b, a, mod.n, mod.ninv);
                ((mp_ptr) (z + i)->coeffs)[j] = d;
            }

        normalise(z + i);
    }

    free(tmp);
}

void
_elem_poly_nmod_mul_flat(elem_ptr z, elem_srcptr x, long xlen, elem_srcptr y, long ylen, const ring_t ring)
{
    long zlen;
    long xbounds[2];
    long ybounds[2];
    long xybound, bits;

    xbounds[0] = 0;
    xbounds[1] = 0;
    _elem_poly_rect_hull(xbounds, x, xlen, ring);

    ybounds[0] = 0;
    ybounds[1] = 0;
    _elem_poly_rect_hull(ybounds, y, ylen, ring);

    zlen = xlen + ylen - 1;

    /* bound inner length */
    xybound = xbounds[1] + ybounds[1] - 1;

    /* bound coefficients */
    bits = FLINT_BIT_COUNT(FLINT_MIN(xbounds[0], ybounds[0]) *
           FLINT_MIN(xbounds[1], ybounds[1])) +
        2*FLINT_BIT_COUNT(RING_PARENT(ring)->nmod.n - 1);

    _elem_poly_nmod_mul_flat_2d(z, x, xlen, y, ylen, xybound, bits, ring);
}

void
elem_poly_nmod_mul_flat(elem_poly_t res, const elem_poly_t op1, const elem_poly_t op2, const ring_t ring)
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
            _elem_poly_nmod_mul_flat(t->coeffs, op1->coeffs, len1, op2->coeffs, len2, ring->parent);
        else
            _elem_poly_nmod_mul_flat(t->coeffs, op2->coeffs, len2, op1->coeffs, len1, ring->parent);

        elem_poly_swap(res, t);
        elem_clear(t, ring);
    }
    else
    {
        elem_poly_fit_length(res, rlen, ring);

        if (len1 >= len2)
            _elem_poly_nmod_mul_flat(res->coeffs, op1->coeffs, len1, op2->coeffs, len2, ring->parent);
        else
            _elem_poly_nmod_mul_flat(res->coeffs, op2->coeffs, len2, op1->coeffs, len1, ring->parent);
    }

    elem_poly_set_length(res, rlen, ring);
    elem_poly_normalise(res, ring);
}

