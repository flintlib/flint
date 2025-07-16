/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2025 Rémi Prébet

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "longlong.h"
#include "fmpz.h"
#include "fmpq.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"


fmpz ** _fmpz_poly_tree_alloc(slong len)
{
    fmpz ** tree = NULL;

    if (len)
    {
        slong i, height = FLINT_CLOG2(len);

        tree = flint_malloc(sizeof(fmpz_t) * (height + 1));
        for (i = 0; i <= height; i++)
            tree[i] = _fmpz_vec_init(len + (len >> i) + 1);
    }

    return tree;
}

void _fmpz_poly_tree_free(fmpz ** tree, slong len)
{
    if (len)
    {
        slong i, height = FLINT_CLOG2(len);

        for (i = 0; i <= height; i++)
            flint_free(tree[i]);

        flint_free(tree);
    }
}

/* Given (a1/b1,...,an/bn), construct the subproduct tree of (b1*x-a1)...(bn*x-an) */
void
_fmpz_poly_tree_build_fmpq_vec(fmpz ** tree, const fmpq * roots, slong len)
{
    slong height, pow, left, i;
    fmpz * pa, * pb;
    fmpq_t ab, cd;

    if (len == 0)
        return;

    height = FLINT_CLOG2(len);

    /* zeroth level, (b*x-a) */
    for (i = 0; i < len; i++)
    {
        fmpz_neg(tree[0] + 2 * i,     fmpq_numref(roots + i));
        fmpz_set(tree[0] + 2 * i + 1, fmpq_denref(roots + i));
    }

    /* first level, (b*x-a)(d*x-c) = bd*x^2 + (-ad-bc)*x + ac */
    if (height > 1)
    {
        fmpq_init(ab); fmpq_init(cd);
        pa = tree[1];

        for (i = 0; i < len / 2; i++)
        {
            fmpq_set(ab, roots + 2 * i);
            fmpq_set(cd, roots + 2 * i + 1);

            // deg 0
            fmpz_mul(pa + 3 * i, fmpq_numref(ab), fmpq_numref(cd));
            // deg 1
            fmpz_mul(   pa + 3 * i + 1, fmpq_numref(ab), fmpq_denref(cd));
            fmpz_addmul(pa + 3 * i + 1, fmpq_denref(ab), fmpq_numref(cd));
            fmpz_neg(pa + 3 * i + 1, pa + 3 * i + 1);
            // deg 2
            fmpz_mul(pa + 3 * i + 2, fmpq_denref(ab), fmpq_denref(cd));
        }

        if (len & 1)
        {
            fmpz_neg(pa + 3 * (len / 2),     fmpq_numref(roots + len - 1));
            fmpz_set(pa + 3 * (len / 2) + 1, fmpq_denref(roots + len - 1));
        }
    }

    for (i = 1; i < height - 1; i++)
    {
        left = len;
        pow = WORD(1) << i;
        pa = tree[i];
        pb = tree[i + 1];

        while (left >= 2 * pow)
        {
            _fmpz_poly_mul(pb, pa, pow + 1, pa + pow + 1, pow + 1);
            left -= 2 * pow;
            pa = pa + 2 * pow + 2;
            pb = pb + 2 * pow + 1;
        }

        if (left > pow)
        {
            _fmpz_poly_mul(pb, pa, pow + 1, pa + pow + 1, left - pow + 1);
        }
        else if (left > 0)
        {
            _fmpz_vec_set(pb, pa, left + 1);
        }
    }
}
