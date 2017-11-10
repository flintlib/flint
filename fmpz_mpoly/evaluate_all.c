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


/*
    evaluate a f(xbar) at xbar = val,
*/
void fmpz_mpoly_evaluate_all_fmpz_straight(fmpz_t ev, fmpz_mpoly_t poly,
                                             fmpz ** val, fmpz_mpoly_ctx_t ctx)
{
    int deg, rev;
    slong i, j, k, N, nvars, bits;
    slong shift, off, fpw;
    slong entries, k_len;
    slong p_len;
    fmpz * p_coeff;
    ulong * p_exp;
    slong * degrees;
    slong * offs;
    ulong * masks;
    fmpz * powers;
    fmpz_t t;
    TMP_INIT;

    p_len = poly->length;
    p_coeff = poly->coeffs;
    p_exp = poly->exps;
    bits = poly->bits;

    if (p_len == 0) {
        fmpz_zero(ev);
        return;
    }

    TMP_START;

    N = words_per_exp(ctx->n, bits);
    degrev_from_ord(deg, rev, ctx->ord);
    nvars = ctx->n - deg;

    degrees = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    fmpz_mpoly_degrees(degrees, poly, ctx);
    
    fpw = FLINT_BITS/bits;

    entries = 1;
    for (i = 0; i < nvars; i++)
        entries += FLINT_BIT_COUNT(degrees[i]);

    offs = (slong *) TMP_ALLOC(entries*sizeof(slong));
    masks = (ulong *) TMP_ALLOC(entries*sizeof(slong));
    powers = (fmpz *) TMP_ALLOC(entries*sizeof(fmpz));

    k = 0;
    for (i = 0; i < nvars; i++)
    {
        mpoly_off_shift(&off, &shift, i, deg, rev, fpw, ctx->n, bits);
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
    assert(k_len + 1 == entries);

    fmpz_zero(ev);
    fmpz_init(t);
    for (i = 0; i < p_len; i++)
    {
        fmpz_set_ui(t, WORD(1));
        for (k = 0; k < k_len; k++)
        {
            if ((p_exp[N*i + offs[k]] & masks[k]) != WORD(0))
                fmpz_mul(t, t, powers + k);
        }
        fmpz_addmul(ev, p_coeff + i, t);
    }

    fmpz_clear(t);
    for (k = 0; k < k_len; k++)
        fmpz_clear(powers + k);

    TMP_END;
}



/*
    Given a polynomial tree with exponents stored in the keys and
    coefficients stored in the data member,
    the function mpoly_rbtree_clear_eval clears the tree
    and stores the evaluation of the polynomial in ev.

    poly = a0 x^0 + a1 x^1 + a2 x^2 + a3 x^3 + a4 x^4 + a5 x^5 + a6 x^6

    tree =                 a3 x^3
                 a1 x^1              a5 x^5

            a0 x^0    a2 x^2    a4 x^4    a6 x^6

    ev = a0*x^0 + x^1*(a1 + a2*x^1) + x^3*(a3 + a4*x^1 + x^2*(a5 + a6*x^1))

*/
void _mpoly_rbnode_clear_evalall_fmpz_tree(mpoly_rbtree_t tree,
                              mpoly_rbnode_t node, slong s, fmpz_t l, fmpz_t x)
{
    fmpz_t r, xp;
    slong e = node->key;
    assert(e >= s);

    fmpz_init(r);
    fmpz_zero(r);
    if (node->right != tree->null)
        _mpoly_rbnode_clear_evalall_fmpz_tree(tree, node->right, e, r, x);

    fmpz_zero(l);
    if (node->left != tree->null)
        _mpoly_rbnode_clear_evalall_fmpz_tree(tree, node->left, s, l, x);

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
void fmpz_mpoly_evaluate_all_fmpz_tree(fmpz_t ev, fmpz_mpoly_t poly,
                                             fmpz ** val, fmpz_mpoly_ctx_t ctx)
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

    N = words_per_exp(ctx->n, bits);
    degrev_from_ord(deg, rev, ctx->ord);
    nvars = ctx->n - deg;

    degrees = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    fmpz_mpoly_degrees(degrees, poly, ctx);

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
    powers = (fmpz *) TMP_ALLOC(entries*sizeof(fmpz));

    /* store bit masks for each power of two of the non-main variables */
    k = 0;
    for (i = 0; i < nvars; i++)
    {
        if (i == main_var)
            continue;

        mpoly_off_shift(&off, &shift, i, deg, rev, fpw, ctx->n, bits);
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
    mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    mpoly_off_shift(&main_off, &main_shift, main_var, deg, rev, fpw, ctx->n, bits);
    mpoly_rbtree_init(tree);
    fmpz_init(t);
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
    _mpoly_rbnode_clear_evalall_fmpz_tree(tree, tree->head->left, WORD(0), ev, val[main_var]);

    for (k = 0; k < k_len; k++)
        fmpz_clear(powers + k);

    TMP_END;
}
