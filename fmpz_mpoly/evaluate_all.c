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
#include "assert.h"


/*
    Given a polynomial tree with exponents stored in the keys and
    coefficients stored in the data member,
    the function _mpoly_rbnode_clear_evalall_tree_fmpz clears the tree
    and stores the evaluation of the polynomial in "l".

    poly = a0 x^0 + a1 x^1 + a2 x^2 + a3 x^3 + a4 x^4 + a5 x^5 + a6 x^6

    tree =                 a3 x^3
                 a1 x^1              a5 x^5

            a0 x^0    a2 x^2    a4 x^4    a6 x^6

    l = a0*x^0 + x^1*(a1 + a2*x^1) + x^3*(a3 + a4*x^1 + x^2*(a5 + a6*x^1))

*/
void _mpoly_rbnode_clear_evalall_tree_fmpz(mpoly_rbtree_t tree,
                              mpoly_rbnode_t node, slong s, fmpz_t l, fmpz_t x)
{
    fmpz_t r, xp;
    slong e = node->key;
    assert(e >= s);

    fmpz_init(r);
    fmpz_zero(r);
    if (node->right != tree->null)
        _mpoly_rbnode_clear_evalall_tree_fmpz(tree, node->right, e, r, x);

    fmpz_zero(l);
    if (node->left != tree->null)
        _mpoly_rbnode_clear_evalall_tree_fmpz(tree, node->left, s, l, x);

    fmpz_init(xp);
    fmpz_pow_ui(xp, x, e - s);
    fmpz_add(r, r, (fmpz*)(&node->data));
    fmpz_addmul(l, xp, r);

    fmpz_clear(r);
    fmpz_clear(xp);
    fmpz_clear((fmpz*)(&node->data));
    flint_free(node);
}


/*
    evaluate a f(xbar) at xbar = val,
*/
void fmpz_mpoly_evaluate_all_tree_fmpz(fmpz_t ev, fmpz_mpoly_t poly,
                                             fmpz ** val, fmpz_mpoly_ctx_t ctx)
{
    int new;
    slong i, j, k, N, bits, nvars = ctx->minfo->nvars;
    slong main_exp, main_var, main_shift, main_off, shift, off;
    ulong mask;
    slong entries, k_len;
    slong p_len;
    fmpz * p_coeff;
    ulong * p_exp;
    slong * degrees;
    slong * offs;
    ulong * masks;
    fmpz * powers;
    fmpz_t t;
    mpoly_rbtree_t tree;
    mpoly_rbnode_struct * node;
    TMP_INIT;

    p_len = poly->length;
    p_coeff = poly->coeffs;
    p_exp = poly->exps;
    bits = poly->bits;

    if (p_len == 0)
    {
        fmpz_zero(ev);
        return;
    }

    TMP_START;

    degrees = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    fmpz_mpoly_degrees_si(degrees, poly, ctx);

    /* pick main variable with highest degree */
    main_var = 0;
    for (i = 1; i < nvars; i++)
    {
        if (degrees[i] > degrees[main_var])
            main_var = i;
    }

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
    powers = (fmpz *) TMP_ALLOC(entries*sizeof(fmpz));

    N = mpoly_words_per_exp(bits, ctx->minfo);

    /* store bit masks for each power of two of the non-main variables */
    k = 0;
    for (i = 0; i < nvars; i++)
    {
        if (i == main_var)
            continue;

        mpoly_gen_offset_shift(&off, &shift, i, N, bits, ctx->minfo);
        for (j = 1; j <= degrees[i]; j *= 2)
        {
            offs[k] = off;
            masks[k] = j << shift;
            fmpz_init(powers + k);
            if (j == 1)
                fmpz_set(powers + k, val[i]);
            else
                fmpz_mul(powers + k, powers + k - 1, powers + k - 1);
            k++;
        }
    }
    k_len = k;
    assert(k_len == entries);

    /* accumulate coefficients of the main variable */
    mpoly_gen_offset_shift(&main_off, &main_shift, main_var, N, bits, ctx->minfo);
    mpoly_rbtree_init(tree);
    fmpz_init(t);
    mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    for (i = 0; i < p_len; i++)
    {
        main_exp = (p_exp[N*i + main_off] >> main_shift) & mask;
        node = mpoly_rbtree_get(&new, tree, main_exp);
        if (new)
        {
            fmpz_init((fmpz*)(&node->data));
            fmpz_zero((fmpz*)(&node->data));
        }

        fmpz_set_ui(t, WORD(1));
        for (k = 0; k < k_len; k++)
        {
            if ((p_exp[N*i + offs[k]] & masks[k]) != WORD(0))
                fmpz_mul(t, t, powers + k);
        }
        fmpz_addmul((fmpz*)(&node->data), p_coeff + i, t);
    }
    fmpz_clear(t);

    /* use tree method to evaluate in the main variable */
    _mpoly_rbnode_clear_evalall_tree_fmpz(tree, tree->head->left, WORD(0), ev, val[main_var]);

    for (k = 0; k < k_len; k++)
        fmpz_clear(powers + k);

    TMP_END;
}
