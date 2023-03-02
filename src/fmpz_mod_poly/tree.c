/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"

fmpz_poly_struct ** _fmpz_mod_poly_tree_alloc(slong len)
{
    fmpz_poly_struct ** tree = NULL;

    if (len)
    {
        slong i, j, height = FLINT_CLOG2(len);

        tree = flint_malloc(sizeof(fmpz_poly_struct *) * (height + 1));
        for (i = 0; i <= height; i++, len = (len + 1)/2)
        {
           tree[i] = flint_malloc(sizeof(fmpz_poly_struct) * len);
           for (j = 0; j < len; j++)
              fmpz_poly_init(tree[i] + j);
        }
            
    }

    return tree;
}

void _fmpz_mod_poly_tree_free(fmpz_poly_struct ** tree, slong len)
{
    if (len)
    {
        slong i, j, height = FLINT_CLOG2(len);

        for (i = 0; i <= height; i++, len = (len + 1)/2)
        {
           for (j = 0; j < len; j++)
              fmpz_poly_clear(tree[i] + j);
           flint_free(tree[i]);
        }

        flint_free(tree);
    }
}

void
_fmpz_mod_poly_tree_build(fmpz_poly_struct ** tree, const fmpz * roots, slong len, const fmpz_t mod)
{
    slong height, pow, left, i;
    fmpz_poly_struct * pa, * pb;

    if (len == 0)
        return;

    height = FLINT_CLOG2(len);

    /* zeroth level, (x-a) */
    for (i = 0; i < len; i++)
    {
        fmpz_poly_set_coeff_ui(tree[0] + i, 1, 1);
        fmpz_negmod((tree[0] + i)->coeffs, roots + i, mod);
    }

    for (i = 0; i < height - 1; i++)
    {
        left = len;
        pow = WORD(1) << i;
        pa = tree[i];
        pb = tree[i + 1];

        while (left >= 2 * pow)
        {
            fmpz_poly_fit_length(pb, pa->length + (pa + 1)->length - 1);
            _fmpz_mod_poly_mul(pb->coeffs, pa->coeffs, pa->length, (pa + 1)->coeffs, (pa + 1)->length, mod);
            _fmpz_poly_set_length(pb, pa->length + (pa + 1)->length - 1);
            left -= 2 * pow;
            pa += 2;
            pb += 1;
        }

        if (left > pow)
        {
            fmpz_poly_fit_length(pb, pa->length + (pa + 1)->length - 1);
            _fmpz_mod_poly_mul(pb->coeffs, pa->coeffs, pa->length, (pa + 1)->coeffs, (pa + 1)->length, mod);
            _fmpz_poly_set_length(pb, pa->length + (pa + 1)->length - 1);
        } else if (left > 0)
            fmpz_poly_set(pb, pa);
    }
}
