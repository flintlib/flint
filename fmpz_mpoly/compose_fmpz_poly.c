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
static void _rbnode_clear_sp(mpoly_rbtree_t tree, mpoly_rbnode_t node,
                             slong s, fmpz_poly_t l, const fmpz_poly_t x,
                                                   const fmpz_mpoly_ctx_t ctx)
{
    fmpz_poly_t r, xp;
    slong e = node->key;
    FLINT_ASSERT(e >= s);

    fmpz_poly_init(r);
    fmpz_poly_zero(r);
    if (node->right != tree->null)
        _rbnode_clear_sp(tree, node->right, e, r, x, ctx);

    fmpz_poly_zero(l);
    if (node->left != tree->null)
        _rbnode_clear_sp(tree, node->left, s, l, x, ctx);

    fmpz_poly_init(xp);
    fmpz_poly_pow(xp, x, e - s);

    fmpz_poly_add(r, r, node->data);
    fmpz_poly_mul(r, xp, r);
    fmpz_poly_add(l, l, r);

    fmpz_poly_clear(r);
    fmpz_poly_clear(xp);
    fmpz_poly_clear(node->data);
    flint_free(node->data);
    flint_free(node);
}



void _fmpz_mpoly_compose_fmpz_poly_sp(fmpz_poly_t A, const fmpz_mpoly_t B,
                      fmpz_poly_struct * const * C, const fmpz_mpoly_ctx_t ctx)
{
    int new;
    slong i, j, k, N, bits, nvars = ctx->minfo->nvars;
    slong main_exp, main_var, main_shift, main_off, shift, off;
    ulong mask;
    slong entries, k_len;
    slong Blen;
    fmpz * Bcoeff;
    ulong * Bexp;
    fmpz * degrees;
    slong * offs;
    ulong * masks;
    fmpz_poly_struct * powers;
    fmpz_poly_t t, t2;
    mpoly_rbtree_t tree;
    mpoly_rbnode_struct * node;
    TMP_INIT;

    Blen = B->length;
    Bcoeff = B->coeffs;
    Bexp = B->exps;
    bits = B->bits;

    FLINT_ASSERT(Blen != 0);

    TMP_START;

    degrees = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    mpoly_degrees_si(degrees, Bexp, Blen, bits, ctx->minfo);

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
    powers = (fmpz_poly_struct *) TMP_ALLOC(entries*sizeof(fmpz_poly_struct));

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
            fmpz_poly_init(powers + k);
            if (j == 0)
                fmpz_poly_set(powers + k, C[i]);
            else
                fmpz_poly_mul(powers + k, powers + k - 1, powers + k - 1);
            k++;
        }
    }
    k_len = k;
    FLINT_ASSERT(k_len == entries);

    /* accumulate coefficients of the main variable */
    mpoly_gen_offset_shift(&main_off, &main_shift, main_var, N, bits, ctx->minfo);
    mpoly_rbtree_init(tree);
    fmpz_poly_init(t);
    fmpz_poly_init(t2);
    mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    for (i = 0; i < Blen; i++)
    {
        main_exp = (Bexp[N*i + main_off] >> main_shift) & mask;
        node = mpoly_rbtree_get(&new, tree, main_exp);
        if (new)
        {
            node->data = flint_malloc(sizeof(nmod_poly_struct));
            fmpz_poly_init(node->data);
        }

        fmpz_poly_set_fmpz(t, Bcoeff + i);
        for (k = 0; k < k_len; k++)
        {
            if ((Bexp[N*i + offs[k]] & masks[k]) != WORD(0))
            {
                fmpz_poly_mul(t2, t, powers + k);
                fmpz_poly_swap(t, t2);
            }
        }
        fmpz_poly_add(t2, t, node->data);
        fmpz_poly_swap(t2, node->data);
    }
    fmpz_poly_clear(t);
    fmpz_poly_clear(t2);

    for (k = 0; k < k_len; k++)
        fmpz_poly_clear(powers + k);

    /* use tree method to evaluate in the main variable */
    _rbnode_clear_sp(tree, tree->head->left, WORD(0), A, C[main_var], ctx);

    TMP_END;
}


static void _rbnode_clear_mp(mpoly_rbtree_t tree, mpoly_rbnode_t node,
                       const fmpz_t s, fmpz_poly_t l, const fmpz_poly_t x,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    fmpz_poly_t r, xp;
    FLINT_ASSERT(fmpz_cmp(&node->key, s) >= 0);

    fmpz_poly_init(r);
    if (node->right != tree->null)
        _rbnode_clear_mp(tree, node->right, &node->key, r, x, ctx);

    fmpz_poly_zero(l);
    if (node->left != tree->null)
        _rbnode_clear_mp(tree, node->left, s, l, x, ctx);

    fmpz_poly_init(xp);
    fmpz_sub(&node->key, &node->key, s);
    FLINT_ASSERT(fmpz_sgn(&node->key) >= 0);
    if (fmpz_fits_si(&node->key))
    {
        fmpz_poly_pow(xp, x, fmpz_get_si(&node->key));
    }
    else
    {
        slong degree = fmpz_poly_degree(x);
        fmpz_poly_zero(xp);
        if (degree == WORD(0))
        {
            fmpz_t t;
            fmpz_init(t);
            fmpz_poly_get_coeff_fmpz(t, x, 0);
            fmpz_pow_fmpz(t ,t, &node->key);
            fmpz_poly_set_fmpz(xp, t);
            fmpz_clear(t);
        }
        else if (degree > WORD(0))
        {
            /* lets not try to power a multinomial */
            flint_throw(FLINT_EXPOF,
                      "Exponent overflow in nmod_mpoly_evaluate_nmod_poly");
        }
    
    }
    fmpz_poly_add(r, r, node->data);
    fmpz_poly_mul(r, xp, r);
    fmpz_poly_add(l, l, r);

    fmpz_clear(&node->key);
    fmpz_poly_clear(r);
    fmpz_poly_clear(xp);
    fmpz_poly_clear(node->data);
    flint_free(node->data);
    flint_free(node);
}



void _fmpz_mpoly_compose_fmpz_poly_mp(fmpz_poly_t A, const fmpz_mpoly_t B,
                      fmpz_poly_struct * const * C, const fmpz_mpoly_ctx_t ctx)
{
    int new;
    ulong l;
    slong i, k, N, bits, nvars = ctx->minfo->nvars;
    fmpz_t main_exp;
    slong main_var, main_off, off;
    slong entries, k_len;
    slong Blen;
    fmpz * Bcoeff;
    ulong * Bexp;
    fmpz * degrees;
    slong * offs;
    ulong * masks;
    mp_bitcnt_t * bitcounts;
    fmpz_t s;
    fmpz_poly_struct * powers;
    fmpz_poly_t t, t2;
    mpoly_rbtree_t tree;
    mpoly_rbnode_struct * node;
    TMP_INIT;

    Blen = B->length;
    Bcoeff = B->coeffs;
    Bexp = B->exps;
    bits = B->bits;

    FLINT_ASSERT(Blen != 0);

    TMP_START;

    bitcounts = (mp_bitcnt_t * ) TMP_ALLOC(nvars*sizeof(mp_bitcnt_t));
    degrees = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
    for (i = 0; i < nvars; i++)
    {
        fmpz_init(degrees + i);
    }
    mpoly_degrees_ffmpz(degrees, Bexp, Blen, bits, ctx->minfo);

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
        bitcounts[i] = fmpz_bits(degrees + i);
        if (i == main_var)
            continue;
        entries += bitcounts[i];
    }
    offs = (slong *) TMP_ALLOC(entries*sizeof(slong));
    masks = (ulong *) TMP_ALLOC(entries*sizeof(slong));
    powers = (fmpz_poly_struct *) TMP_ALLOC(entries*sizeof(fmpz_poly_struct));

    N = mpoly_words_per_exp(bits, ctx->minfo);

    /* store bit masks for each power of two of the non-main variables */
    k = 0;
    for (i = 0; i < nvars; i++)
    {
        if (i == main_var)
            continue;

        off = mpoly_gen_offset_mp(i, N, bits, ctx->minfo);

        for (l = 0; l < bitcounts[i]; l++)
        {
            offs[k] = off + (l / FLINT_BITS);
            masks[k] = UWORD(1) << (l % FLINT_BITS);
            fmpz_poly_init(powers + k);
            if (l == 0)
                fmpz_poly_set(powers + k, C[i]);
            else
                fmpz_poly_mul(powers + k, powers + k - 1, powers + k - 1);
            k++;
        }
    }
    k_len = k;
    FLINT_ASSERT(k_len == entries);

    /* accumulate coefficients of the main variable */
    main_off = mpoly_gen_offset_mp(main_var, N, bits, ctx->minfo);
    mpoly_rbtree_init(tree);
    fmpz_poly_init(t);
    fmpz_poly_init(t2);
    fmpz_init(main_exp);
    for (i = 0; i < Blen; i++)
    {
        fmpz_set_ui_array(main_exp, Bexp + N*i + main_off, bits/FLINT_BITS);
        node = mpoly_rbtree_get_fmpz(&new, tree, main_exp);
        if (new)
        {
            node->data = flint_malloc(sizeof(fmpz_poly_struct));
            fmpz_poly_init(node->data);
        }

        fmpz_poly_set_fmpz(t, Bcoeff + i);
        for (k = 0; k < k_len; k++)
        {
            if ((Bexp[N*i + offs[k]] & masks[k]) != WORD(0))
            {
                fmpz_poly_mul(t2, t, powers + k);
                fmpz_poly_swap(t, t2);
            }
        }
        fmpz_poly_add(t2, t, node->data);
        fmpz_poly_swap(t2, node->data);
    }
    fmpz_clear(main_exp);
    fmpz_poly_clear(t);
    fmpz_poly_clear(t2);

    for (k = 0; k < k_len; k++)
        fmpz_poly_clear(powers + k);

    for (i = 0; i < nvars; i++)
        fmpz_clear(degrees + i);

    /* use tree method to evaluate in the main variable */
    fmpz_init(s);
    _rbnode_clear_mp(tree, tree->head->left, s, A, C[main_var], ctx);
    fmpz_clear(s);

    TMP_END;
}


void fmpz_mpoly_compose_fmpz_poly(fmpz_poly_t A, const fmpz_mpoly_t B,
                      fmpz_poly_struct * const * C, const fmpz_mpoly_ctx_t ctx)
{
    if (B->length == 0)
    {
        fmpz_poly_zero(A);
        return;
    }
    else if (B->bits <= FLINT_BITS)
    {
        _fmpz_mpoly_compose_fmpz_poly_sp(A, B, C, ctx);
        return;
    }
    else
    {
        _fmpz_mpoly_compose_fmpz_poly_mp(A, B, C, ctx);
        return;
    }
}

