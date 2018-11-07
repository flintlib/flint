/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

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
static void _mpoly_rbnode_clear_sp(mpoly_rbtree_t tree, mpoly_rbnode_t node,
                                             slong s, fmpz_t l, const fmpz_t x)
{
    fmpz_t r, xp;
    slong e = node->key;
    FLINT_ASSERT(e >= s);

    fmpz_init(r);
    fmpz_zero(r);
    if (node->right != tree->null)
        _mpoly_rbnode_clear_sp(tree, node->right, e, r, x);

    fmpz_zero(l);
    if (node->left != tree->null)
        _mpoly_rbnode_clear_sp(tree, node->left, s, l, x);

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
void _fmpz_mpoly_evaluate_all_fmpz_sp(fmpz_t ev, const fmpz_mpoly_t A,
                                fmpz * const * val, const fmpz_mpoly_ctx_t ctx)
{
    int new;
    slong i, j, k, N, bits, nvars = ctx->minfo->nvars;
    slong main_exp, main_var, main_shift, main_off, shift, off;
    ulong mask;
    slong entries, k_len;
    slong Alen;
    fmpz * Acoeff;
    ulong * Aexp;
    slong * degrees;
    slong * offs;
    ulong * masks;
    fmpz * powers;
    fmpz_t t;
    mpoly_rbtree_t tree;
    mpoly_rbnode_struct * node;
    TMP_INIT;

    Alen = A->length;
    Acoeff = A->coeffs;
    Aexp = A->exps;
    bits = A->bits;

    FLINT_ASSERT(Alen > 0);

    TMP_START;

    degrees = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    mpoly_degrees_si(degrees, Aexp, Alen, bits, ctx->minfo);

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
        mp_bitcnt_t varibits;

        if (i == main_var)
            continue;

        mpoly_gen_offset_shift(&off, &shift, i, N, bits, ctx->minfo);
        varibits = FLINT_BIT_COUNT(degrees[i]);
        for (j = 0; j < varibits; j++)
        {
            offs[k] = off;
            masks[k] = UWORD(1) << (shift + j);
            fmpz_init(powers + k);
            if (j == 0)
                fmpz_set(powers + k, val[i]);
            else
                fmpz_mul(powers + k, powers + k - 1, powers + k - 1);
            k++;
        }
    }
    k_len = k;
    FLINT_ASSERT(k_len == entries);

    /* accumulate coefficients of the main variable */
    mpoly_gen_offset_shift(&main_off, &main_shift, main_var, N, bits, ctx->minfo);
    mpoly_rbtree_init(tree);
    fmpz_init(t);
    mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    for (i = 0; i < Alen; i++)
    {
        main_exp = (Aexp[N*i + main_off] >> main_shift) & mask;
        node = mpoly_rbtree_get(&new, tree, main_exp);
        if (new)
        {
            fmpz_init((fmpz*)(&node->data));
        }

        fmpz_set(t, Acoeff + i);
        for (k = 0; k < k_len; k++)
        {
            if ((Aexp[N*i + offs[k]] & masks[k]) != WORD(0))
                fmpz_mul(t, t, powers + k);
        }
        fmpz_add((fmpz*)(&node->data), (fmpz*)(&node->data), t);
    }
    fmpz_clear(t);

    /* use tree method to evaluate in the main variable */
    _mpoly_rbnode_clear_sp(tree, tree->head->left, WORD(0), ev, val[main_var]);

    for (k = 0; k < k_len; k++)
        fmpz_clear(powers + k);

    TMP_END;
}


static void _mpoly_rbnode_clear_mp(mpoly_rbtree_t tree, mpoly_rbnode_t node,
                                      const fmpz_t s, fmpz_t l, const fmpz_t x)
{
    fmpz_t r, xp;
    FLINT_ASSERT(fmpz_cmp(&node->key, s) >= 0);

    fmpz_init(r);
    if (node->right != tree->null)
        _mpoly_rbnode_clear_mp(tree, node->right, &node->key, r, x);

    fmpz_zero(l);
    if (node->left != tree->null)
        _mpoly_rbnode_clear_mp(tree, node->left, s, l, x);

    fmpz_init(xp);

    fmpz_sub((fmpz*)(&node->key), (fmpz*)(&node->key), s);
    fmpz_pow_fmpz(xp, x, (fmpz*)(&node->key));
    fmpz_add(r, r, (fmpz*)(&node->data));
    fmpz_addmul(l, xp, r);

    fmpz_clear(r);
    fmpz_clear(xp);
    fmpz_clear((fmpz*)(&node->data));
    fmpz_clear((fmpz*)(&node->key));
    flint_free(node);
}


/*
    evaluate a f(xbar) at xbar = val,
*/
void _fmpz_mpoly_evaluate_all_fmpz_mp(fmpz_t ev, const fmpz_mpoly_t A,
                               fmpz * const * vals, const fmpz_mpoly_ctx_t ctx)
{
    int new;
    slong i, j, k, N, bits, nvars = ctx->minfo->nvars;
    slong main_var, main_off, off;
    fmpz_t main_exp;
    slong entries, k_len;
    slong Alen;
    fmpz * Acoeff;
    ulong * Aexp;
    slong * degrees;
    slong * offs;
    ulong * masks;
    fmpz * powers;
    fmpz_t t;
    fmpz_t s;
    mpoly_rbtree_t tree;
    mpoly_rbnode_struct * node;
    TMP_INIT;

    Alen = A->length;
    Acoeff = A->coeffs;
    Aexp = A->exps;
    bits = A->bits;

    FLINT_ASSERT(Alen > 0);

    TMP_START;

    degrees = _fmpz_vec_init(nvars);
    mpoly_degrees_ffmpz(degrees, Aexp, Alen, bits, ctx->minfo);

    /* pick main variable with highest degree */
    main_var = 0;
    for (i = 1; i < nvars; i++)
    {
        if (fmpz_cmp(degrees + i, degrees + main_var) > 0)
            main_var = i;
    }

    /* compute how many masks are needed */
    entries = 0;
    for (i = 0; i < nvars; i++)
    {
        if (i == main_var)
            continue;
        entries += fmpz_bits(degrees + i);
    }
    offs = (slong *) TMP_ALLOC(entries*sizeof(slong));
    masks = (ulong *) TMP_ALLOC(entries*sizeof(slong));
    powers = (fmpz *) TMP_ALLOC(entries*sizeof(fmpq));

    N = mpoly_words_per_exp(bits, ctx->minfo);

    /* store bit masks for each power of two of the non-main variables */
    k = 0;
    for (i = 0; i < nvars; i++)
    {
        mp_bitcnt_t varibits;

        if (i == main_var)
            continue;

        off = mpoly_gen_offset_mp(i, N, bits, ctx->minfo);
        varibits = fmpz_bits(degrees + i);
        if (varibits >= FLINT_BITS && !fmpz_is_zero(vals[i])
                                  && !fmpz_is_pm1(vals[i]))
        {
            flint_throw(FLINT_ERROR, "Exponent too large in fmpz_mpoly_evaluate_all_fmpz");
        }
        for (j = 0; j < varibits; j++)
        {
            offs[k] = off + (j / FLINT_BITS);
            masks[k] = UWORD(1) << (j % FLINT_BITS);
            fmpz_init(powers + k);
            if (j == 0)
                fmpz_set(powers + k, vals[i]);
            else
                fmpz_mul(powers + k, powers + k - 1, powers + k - 1);
            k++;
        }
    }
    k_len = k;
    FLINT_ASSERT(k_len == entries);

    /* accumulate coefficients of the main variable */
    main_off = mpoly_gen_offset_mp(main_var, N, bits, ctx->minfo);
    mpoly_rbtree_init(tree);
    fmpz_init(t);
    fmpz_init(main_exp);
    for (i = 0; i < Alen; i++)
    {
        fmpz_set_ui_array(main_exp, Aexp + N*i + main_off, bits/FLINT_BITS);
        node = mpoly_rbtree_get_fmpz(&new, tree, main_exp);
        if (new)
        {
            fmpz_init((fmpz*)(&node->data));
        }

        fmpz_set(t, Acoeff + i);
        for (k = 0; k < k_len; k++)
        {
            if ((Aexp[N*i + offs[k]] & masks[k]) != WORD(0))
                fmpz_mul(t, t, powers + k);
        }
        fmpz_add((fmpz*)(&node->data), (fmpz*)(&node->data), t);
    }
    fmpz_clear(main_exp);
    fmpz_clear(t);

    /* use tree method to evaluate in the main variable */
    fmpz_init(s);
    _mpoly_rbnode_clear_mp(tree, tree->head->left, s, ev, vals[main_var]);
    fmpz_clear(s);

    _fmpz_vec_clear(degrees, nvars);
    for (k = 0; k < k_len; k++)
        fmpz_clear(powers + k);

    TMP_END;
}


void fmpz_mpoly_evaluate_all_fmpz(fmpz_t ev, const fmpz_mpoly_t A,
                               fmpz * const * vals, const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_mpoly_is_zero(A, ctx))
    {
        fmpz_zero(ev);
        return;
    }

    if (A->bits <= FLINT_BITS)
    {
        _fmpz_mpoly_evaluate_all_fmpz_sp(ev, A, vals, ctx);
    }
    else
    {
        _fmpz_mpoly_evaluate_all_fmpz_mp(ev, A, vals, ctx);
    }
}
