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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

mp_ptr * _nmod_poly_tree_alloc(long len)
{
    mp_ptr * tree = NULL;

    if (len)
    {
        long i, height = FLINT_CLOG2(len);

        tree = flint_malloc(sizeof(mp_ptr) * (height + 1));
        for (i = 0; i <= height; i++)
            tree[i] = _nmod_vec_init(len + (len >> i) + 1);
    }

    return tree;
}

void _nmod_poly_tree_free(mp_ptr * tree, long len)
{
    if (len)
    {
        long i, height = FLINT_CLOG2(len);

        for (i = 0; i <= height; i++)
            flint_free(tree[i]);

        flint_free(tree);
    }
}

void
_nmod_poly_tree_build(mp_ptr * tree, mp_srcptr roots, long len, nmod_t mod)
{
    long height, pow, left, i;
    mp_ptr pa, pb;

    if (len == 0)
        return;

    height = FLINT_CLOG2(len);

    /* zeroth level, (x-a) */
    for (i = 0; i < len; i++)
    {
        tree[0][2 * i + 1] = 1;
        tree[0][2 * i] = nmod_neg(roots[i], mod);
    }

    /* first level, (x-a)(x-b) = x^2 + (-a-b)*x + a*b */
    if (height > 1)
    {
        pa = tree[1];

        for (i = 0; i < len / 2; i++)
        {
            mp_limb_t a, b;

            a = roots[2 * i];
            b = roots[2 * i + 1];

            pa[3 * i] = nmod_mul(a, b, mod);
            pa[3 * i + 1] = nmod_neg(nmod_add(a, b, mod), mod);
            pa[3 * i + 2] = 1;
        }

        if (len & 1)
        {
            pa[3 * (len / 2)] = nmod_neg(roots[len-1], mod);
            pa[3 * (len / 2) + 1] = 1;
        }
    }

    for (i = 1; i < height - 1; i++)
    {
        left = len;
        pow = 1L << i;
        pa = tree[i];
        pb = tree[i + 1];

        while (left >= 2 * pow)
        {
            _nmod_poly_mul(pb, pa, pow + 1, pa + pow + 1, pow + 1, mod);
            left -= 2 * pow;
            pa += 2 * pow + 2;
            pb += 2 * pow + 1;
        }

        if (left > pow)
            _nmod_poly_mul(pb, pa, pow + 1, pa + pow + 1, left - pow + 1, mod);
        else if (left > 0)
            _nmod_vec_set(pb, pa, left + 1);
    }
}
