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
    Copyright (C) 2012 William Hart

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"

void
_fmpz_mod_poly_evaluate_fmpz_vec_fast_precomp(fmpz * vs, const fmpz * poly,
    len_t plen, fmpz_poly_struct ** tree, len_t len, const fmpz_t mod)
{
    len_t height, i, j, pow, left;
    len_t tree_height;
    fmpz_t temp, inv;
    fmpz * t, * u, * pb, * pc, * swap;
    fmpz_poly_struct * pa;

    fmpz_init(temp);
    fmpz_init(inv);

    /* avoid worrying about some degenerate cases */
    if (len < 2 || plen < 2)
    {
        if (len == 1)
        {
            fmpz_negmod(temp, tree[0]->coeffs, mod);
            _fmpz_mod_poly_evaluate_fmpz(vs, poly, plen, temp, mod);
        } else if (len != 0 && plen == 0)
            _fmpz_vec_zero(vs, len);
        else if (len != 0 && plen == 1)
            for (i = 0; i < len; i++)
                fmpz_set(vs + i, poly);
        
        fmpz_clear(temp);
        return;
    }

    t = _fmpz_vec_init(2*len);
    u = _fmpz_vec_init(2*len);

    left = len;

    /* Initial reduction. We allow the polynomial to be larger
       or smaller than the number of points. */
    height = FLINT_BIT_COUNT(plen - 1) - 1;
    tree_height = FLINT_CLOG2(len);
    while (height >= tree_height)
        height--;
    pow = 1L << height;

    for (i = j = 0; i < len; i += pow, j++)
    {
        pa = tree[height] + j;
        fmpz_invmod(inv, pa->coeffs + pa->length - 1, mod);
        _fmpz_mod_poly_rem(t + i, poly, plen, pa->coeffs, pa->length, inv, mod);
    }

    for (i = height - 1; i >= 0; i--)
    {
        pow = 1L << i;
        left = len;
        pa = tree[i];
        pb = t;
        pc = u;

        left = len;
        while (left >= 2 * pow)
        {
            fmpz_invmod(inv, pa->coeffs + pa->length - 1, mod);
            _fmpz_mod_poly_rem(pc, pb, 2 * pow, pa->coeffs, pa->length, inv, mod);
            
            pa++;
            fmpz_invmod(inv, pa->coeffs + pa->length - 1, mod);
            _fmpz_mod_poly_rem(pc + pow, pb, 2 * pow, pa->coeffs, pa->length, inv, mod);
            
            pa++;
            pb += 2 * pow;
            pc += 2 * pow;
            left -= 2 * pow;
        }
        
        if (left > pow)
        {
            fmpz_invmod(inv, pa->coeffs + pa->length - 1, mod);
            _fmpz_mod_poly_rem(pc, pb, left, pa->coeffs, pa->length, inv, mod);
            
            pa ++;
            fmpz_invmod(inv, pa->coeffs + pa->length - 1, mod);
            _fmpz_mod_poly_rem(pc + pow, pb, left, pa->coeffs, pa->length, inv, mod);
        }
        else if (left > 0)
           _fmpz_vec_set(pc, pb, left);

        swap = t;
        t = u;
        u = swap;
    }

    fmpz_clear(temp);
    fmpz_clear(inv);

    _fmpz_vec_set(vs, t, len);

    _fmpz_vec_clear(t, 2*len);
    _fmpz_vec_clear(u, 2*len);
}

void _fmpz_mod_poly_evaluate_fmpz_vec_fast(fmpz * ys, const fmpz * poly, len_t plen,
    const fmpz * xs, len_t n, const fmpz_t mod)
{
    fmpz_poly_struct ** tree;

    tree = _fmpz_mod_poly_tree_alloc(n);
    _fmpz_mod_poly_tree_build(tree, xs, n, mod);
    _fmpz_mod_poly_evaluate_fmpz_vec_fast_precomp(ys, poly, plen, tree, n, mod);
    _fmpz_mod_poly_tree_free(tree, n);
}

void
fmpz_mod_poly_evaluate_fmpz_vec_fast(fmpz * ys,
        const fmpz_mod_poly_t poly, const fmpz * xs, len_t n)
{
    _fmpz_mod_poly_evaluate_fmpz_vec_fast(ys, poly->coeffs,
                                        poly->length, xs, n, &(poly->p));
}
