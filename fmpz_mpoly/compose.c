/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"
#include "ulong_extras.h"
#include "profiler.h"
#include "assert.h"


void _mpoly_rbnode_clear_comp(mpoly_rbtree_t tree, mpoly_rbnode_t node, slong s,
                         fmpz_mpoly_t v, fmpz_mpoly_t x, fmpz_mpoly_ctx_t ctx2)
{
    fmpz_mpoly_t l, r, xp;
    slong e = node->key;
    assert(e >= s);

    fmpz_mpoly_init(r, ctx2);
    fmpz_mpoly_zero(r, ctx2);
    if (node->right != tree->null)
        _mpoly_rbnode_clear_comp(tree, node->right, e, r, x, ctx2);

    fmpz_mpoly_init(l, ctx2);
    fmpz_mpoly_zero(l, ctx2);
    if (node->left != tree->null)
        _mpoly_rbnode_clear_comp(tree, node->left, s, l, x, ctx2);

    fmpz_mpoly_init(xp, ctx2);
    fmpz_mpoly_pow_fps(xp, x, e - s, ctx2);
    fmpz_mpoly_add(v, r, node->data, ctx2);
    fmpz_mpoly_mul_johnson(r, xp, v, ctx2);
    fmpz_mpoly_add(v, l, r, ctx2);

    fmpz_mpoly_clear(r, ctx2);
    fmpz_mpoly_clear(l, ctx2);
    fmpz_mpoly_clear(xp, ctx2);
    fmpz_mpoly_clear(node->data, ctx2);
    flint_free(node->data);
    flint_free(node);
}


/*
    evaluate a f(xbar) at xbar = polys2,
*/
void fmpz_mpoly_compose(fmpz_mpoly_t res, fmpz_mpoly_t poly1,
     fmpz_mpoly_struct ** polys2, fmpz_mpoly_ctx_t ctx1, fmpz_mpoly_ctx_t ctx2)
{
    int deg, rev, new;
    slong i, j, k, N, nvars, bits;
    slong main_exp, main_var, main_shift, main_off, shift, off, fpw;
    ulong mask;
    slong entries, k_len;
    slong p_len;
    fmpz * p_coeff;
    ulong * p_exp;
    slong * degrees;
    slong * offs;
    ulong * masks;
    fmpz_mpoly_struct * powers;
    fmpz_mpoly_t t;
    mpoly_rbtree_t tree;
    mpoly_rbnode_struct * node;
    TMP_INIT;

    p_len = poly1->length;
    p_coeff = poly1->coeffs;
    p_exp = poly1->exps;
    bits = poly1->bits;

    assert(res != poly1);
    if (p_len == 0)
    {
        fmpz_mpoly_zero(res, ctx1);
        return;
    }

    TMP_START;

    N = words_per_exp(ctx1->n, bits);
    degrev_from_ord(deg, rev, ctx1->ord);
    nvars = ctx1->n - deg;

    degrees = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    fmpz_mpoly_degrees(degrees, poly1, ctx1);

    /* pick main variable with highest degree */
    main_var = 0;
    for (i = 1; i < nvars; i++)
    {
        if (degrees[i] > degrees[main_var])
            main_var = i;
    }

    fpw = FLINT_BITS/bits;

    /* compute how many masks are needed */
    entries = 0;
    for (i = 0; i < nvars; i++)
    {
        if (i == main_var)
            continue;
        entries += FLINT_BIT_COUNT(degrees[i]);
    }
    offs = (slong *) TMP_ALLOC(entries*sizeof(slong));
    masks = (ulong *) TMP_ALLOC(entries*sizeof(slong));
    powers = (fmpz_mpoly_struct *) TMP_ALLOC(entries*sizeof(fmpz_mpoly_struct));

    /* store bit masks for each power of two of the non-main variables */
    k = 0;
    for (i = 0; i < nvars; i++)
    {
        if (i == main_var)
            continue;

        mpoly_off_shift(&off, &shift, i, deg, rev, fpw, ctx1->n, bits);
        for (j = 1; j <= degrees[i]; j *= 2)
        {
            offs[k] = off;
            masks[k] = j << shift;
            fmpz_mpoly_init(powers + k, ctx2);
            if (j == 1)
                fmpz_mpoly_set(powers + k, polys2[i], ctx2);
            else
                fmpz_mpoly_mul_johnson(powers + k, powers + k - 1, powers + k - 1, ctx2);
            k++;
        }
    }
    k_len = k;
    assert(k_len == entries);

    /* accumulate coefficients of the main variable */
    mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    mpoly_off_shift(&main_off, &main_shift, main_var, deg, rev, fpw, ctx1->n, bits);
    mpoly_rbtree_init(tree);
    fmpz_mpoly_init(t, ctx2);
    for (i = 0; i < p_len; i++)
    {
        main_exp = (p_exp[N*i + main_off] >> main_shift) & mask;
        node = mpoly_rbtree_get(&new, tree, main_exp);
        if (new)
        {
            node->data = flint_malloc(sizeof(fmpz_mpoly_struct));
            fmpz_mpoly_init(node->data, ctx2);
            fmpz_mpoly_zero(node->data, ctx2);
        }

        fmpz_mpoly_set_ui(t, WORD(1), ctx2);
        for (k = 0; k < k_len; k++)
        {
            if ((p_exp[N*i + offs[k]] & masks[k]) != WORD(0))
                fmpz_mpoly_mul_johnson(t, t, powers + k, ctx2);
        }
        fmpz_mpoly_scalar_mul_fmpz(t, t, p_coeff + i, ctx2);
        fmpz_mpoly_add(node->data, node->data, t, ctx2);
    }
    fmpz_mpoly_clear(t, ctx2);

    /* use tree method to evaluate in the main variable */
    _mpoly_rbnode_clear_comp(tree, tree->head->left, WORD(0), res, polys2[main_var], ctx2);

    for (k = 0; k < k_len; k++)
        fmpz_mpoly_clear(powers + k, ctx2);

    TMP_END;
}

