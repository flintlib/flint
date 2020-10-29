/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"


static void _rbnode_clear_sp(mpoly_rbtree_t tree, mpoly_rbnode_t node,
                   slong s, fq_nmod_poly_t l, const fq_nmod_poly_t x,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    fq_nmod_poly_t r, xp;
    slong e = node->key;
    FLINT_ASSERT(e >= s);

    fq_nmod_poly_init(r, ctx->fqctx);
    fq_nmod_poly_zero(r, ctx->fqctx);
    if (node->right != tree->null)
        _rbnode_clear_sp(tree, node->right, e, r, x, ctx);

    fq_nmod_poly_zero(l, ctx->fqctx);
    if (node->left != tree->null)
        _rbnode_clear_sp(tree, node->left, s, l, x, ctx);

    fq_nmod_poly_init(xp, ctx->fqctx);
    fq_nmod_poly_pow(xp, x, e - s, ctx->fqctx);

    fq_nmod_poly_add(r, r, node->data, ctx->fqctx);
    fq_nmod_poly_mul(r, xp, r, ctx->fqctx);
    fq_nmod_poly_add(l, l, r, ctx->fqctx);

    fq_nmod_poly_clear(r, ctx->fqctx);
    fq_nmod_poly_clear(xp, ctx->fqctx);
    fq_nmod_poly_clear(node->data, ctx->fqctx);
    flint_free(node->data);
    flint_free(node);
}


int _fq_nmod_mpoly_compose_fq_nmod_poly_sp(fq_nmod_poly_t A, const fq_nmod_mpoly_t B,
                fq_nmod_poly_struct * const * C, const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    int success = 1;
    int new;
    slong i, j, k, N, bits, nvars = ctx->minfo->nvars;
    slong main_exp, main_var, main_shift, main_off, shift, off;
    ulong mask;
    slong entries, k_len;
    slong Blen;
    const mp_limb_t * Bcoeff;
    ulong * Bexp;
    slong * degrees;
    slong * offs;
    ulong * masks;
    fq_nmod_poly_struct * powers;
    fq_nmod_poly_t t, t2;
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
        if (_ff_poly_pow_ui_is_not_feasible(C[i]->length, degrees[i]))
        {
            success = 0;
            goto cleanup_degrees;
        }

        if (i == main_var)
            continue;

        entries += FLINT_BIT_COUNT(degrees[i]);
    }
    offs = (slong *) TMP_ALLOC(entries*sizeof(slong));
    masks = (ulong *) TMP_ALLOC(entries*sizeof(ulong));
    powers = (fq_nmod_poly_struct *) TMP_ALLOC(entries*
                                                  sizeof(fq_nmod_poly_struct));

    N = mpoly_words_per_exp(bits, ctx->minfo);

    /* store bit masks for each power of two of the non-main variables */
    k = 0;
    for (i = 0; i < nvars; i++)
    {
        flint_bitcnt_t varibits;

        if (i == main_var)
            continue;

        mpoly_gen_offset_shift_sp(&off, &shift, i, bits, ctx->minfo);
        varibits = FLINT_BIT_COUNT(degrees[i]);
        for (j = 0; j < varibits; j++)
        {
            offs[k] = off;
            masks[k] = UWORD(1) << (shift + j);
            fq_nmod_poly_init(powers + k, ctx->fqctx);
            if (j == 0)
                fq_nmod_poly_set(powers + k, C[i], ctx->fqctx);
            else
                fq_nmod_poly_mul(powers + k, powers + k - 1, powers + k - 1,
                                                                   ctx->fqctx);
            k++;
        }
    }
    k_len = k;
    FLINT_ASSERT(k_len == entries);

    /* accumulate coefficients of the main variable */
    mpoly_gen_offset_shift_sp(&main_off, &main_shift, main_var, bits, ctx->minfo);
    mpoly_rbtree_init(tree);
    fq_nmod_poly_init(t, ctx->fqctx);
    fq_nmod_poly_init(t2, ctx->fqctx);
    mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    for (i = 0; i < Blen; i++)
    {
        main_exp = (Bexp[N*i + main_off] >> main_shift) & mask;
        node = mpoly_rbtree_get(&new, tree, main_exp);
        if (new)
        {
            node->data = flint_malloc(sizeof(fq_nmod_poly_struct));
            fq_nmod_poly_init(node->data, ctx->fqctx);
            fq_nmod_poly_zero(node->data, ctx->fqctx);
        }

        fq_nmod_poly_fit_length(t, 1, ctx->fqctx);
        n_fq_get_fq_nmod(t->coeffs + 0, Bcoeff + d*i, ctx->fqctx);
        t->length = 1;
        for (k = 0; k < k_len; k++)
        {
            if ((Bexp[N*i + offs[k]] & masks[k]) != WORD(0))
            {
                fq_nmod_poly_mul(t2, t, powers + k, ctx->fqctx);
                fq_nmod_poly_swap(t, t2, ctx->fqctx);
            }
        }
        fq_nmod_poly_add(t2, t, node->data, ctx->fqctx);
        fq_nmod_poly_swap(t2, node->data, ctx->fqctx);
    }
    fq_nmod_poly_clear(t, ctx->fqctx);
    fq_nmod_poly_clear(t2, ctx->fqctx);

    for (k = 0; k < k_len; k++)
        fq_nmod_poly_clear(powers + k, ctx->fqctx);

    /* use tree method to evaluate in the main variable */
    _rbnode_clear_sp(tree, tree->head->left, WORD(0), A, C[main_var], ctx);

cleanup_degrees:

    TMP_END;

    return success;
}


static int _rbnode_clear_mp(mpoly_rbtree_t tree, mpoly_rbnode_t node,
                     const fmpz_t s, fq_nmod_poly_t l, const fq_nmod_poly_t x,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    int success = 1;
    fq_nmod_poly_t r, xp;
    FLINT_ASSERT(fmpz_cmp(&node->key, s) >= 0);

    fq_nmod_poly_init(r, ctx->fqctx);
    fq_nmod_poly_zero(r, ctx->fqctx);
    if (node->right != tree->null)
    {
        if (!_rbnode_clear_mp(tree, node->right, &node->key, r, x, ctx))
            success = 0;
    }

    fq_nmod_poly_zero(l, ctx->fqctx);
    if (node->left != tree->null)
    {
        if (!_rbnode_clear_mp(tree, node->left, s, l, x, ctx))
            success = 0;
    }

    fq_nmod_poly_init(xp, ctx->fqctx);
    fmpz_sub(&node->key, &node->key, s);
    FLINT_ASSERT(fmpz_sgn(&node->key) >= 0);
    if (fmpz_fits_si(&node->key))
    {
        fq_nmod_poly_pow(xp, x, fmpz_get_si(&node->key), ctx->fqctx);
    }
    else
    {
        slong degree = fq_nmod_poly_degree(x, ctx->fqctx);
        fq_nmod_poly_zero(xp, ctx->fqctx);
        if (degree == 0)
        {
            fq_nmod_t t;
            fq_nmod_init(t, ctx->fqctx);
            fq_nmod_pow(t, x->coeffs + 0, &node->key, ctx->fqctx);
            fq_nmod_poly_set_coeff(xp, 0, t, ctx->fqctx);
            fq_nmod_clear(t, ctx->fqctx);
        }
        else if (degree > 0)
        {
            success = 0;
        }
    }
    fq_nmod_poly_add(r, r, node->data, ctx->fqctx);
    fq_nmod_poly_mul(r, xp, r, ctx->fqctx);
    fq_nmod_poly_add(l, l, r, ctx->fqctx);

    fmpz_clear(&node->key);
    fq_nmod_poly_clear(r, ctx->fqctx);
    fq_nmod_poly_clear(xp, ctx->fqctx);
    fq_nmod_poly_clear(node->data, ctx->fqctx);
    flint_free(node->data);
    flint_free(node);

    return success;
}


int _fq_nmod_mpoly_compose_fq_nmod_poly_mp(fq_nmod_poly_t A, const fq_nmod_mpoly_t B,
                fq_nmod_poly_struct * const * C, const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    int success = 1;
    int new;
    ulong l;
    slong i, k, N, bits, nvars = ctx->minfo->nvars;
    fmpz_t main_exp;
    slong main_var, main_off, off;
    slong entries, k_len;
    slong Blen;
    const mp_limb_t * Bcoeff;
    ulong * Bexp;
    fmpz * degrees;
    slong * offs;
    ulong * masks;
    flint_bitcnt_t * bitcounts;
    fmpz_t s;
    fq_nmod_poly_struct * powers;
    fq_nmod_poly_t t, t2;
    mpoly_rbtree_t tree;
    mpoly_rbnode_struct * node;
    TMP_INIT;

    Blen = B->length;
    Bcoeff = B->coeffs;
    Bexp = B->exps;
    bits = B->bits;

    FLINT_ASSERT(Blen != 0);

    TMP_START;

    bitcounts = (flint_bitcnt_t * ) TMP_ALLOC(nvars*sizeof(flint_bitcnt_t));
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
        if (_ff_poly_pow_fmpz_is_not_feasible(C[i]->length, degrees + i))
        {
            success = 0;
            goto cleanup_degrees;
        }

        bitcounts[i] = fmpz_bits(degrees + i);

        if (i == main_var)
            continue;

        entries += bitcounts[i];
    }
    offs = (slong *) TMP_ALLOC(entries*sizeof(slong));
    masks = (ulong *) TMP_ALLOC(entries*sizeof(ulong));
    powers = (fq_nmod_poly_struct *) TMP_ALLOC(entries*
                                                  sizeof(fq_nmod_poly_struct));

    N = mpoly_words_per_exp(bits, ctx->minfo);

    /* store bit masks for each power of two of the non-main variables */
    k = 0;
    for (i = 0; i < nvars; i++)
    {
        if (i == main_var)
            continue;

        off = mpoly_gen_offset_mp(i, bits, ctx->minfo);

        for (l = 0; l < bitcounts[i]; l++)
        {
            offs[k] = off + (l/FLINT_BITS);
            masks[k] = UWORD(1) << (l%FLINT_BITS);
            fq_nmod_poly_init(powers + k, ctx->fqctx);
            if (l == 0)
                fq_nmod_poly_set(powers + k, C[i], ctx->fqctx);
            else
                fq_nmod_poly_mul(powers + k, powers + k - 1, powers + k - 1,
                                                                   ctx->fqctx);
            k++;
        }
    }
    k_len = k;
    FLINT_ASSERT(k_len == entries);

    /* accumulate coefficients of the main variable */
    main_off = mpoly_gen_offset_mp(main_var, bits, ctx->minfo);
    mpoly_rbtree_init(tree);
    fq_nmod_poly_init(t, ctx->fqctx);
    fq_nmod_poly_init(t2, ctx->fqctx);
    fmpz_init(main_exp);
    for (i = 0; i < Blen; i++)
    {
        fmpz_set_ui_array(main_exp, Bexp + N*i + main_off, bits/FLINT_BITS);
        node = mpoly_rbtree_get_fmpz(&new, tree, main_exp);
        if (new)
        {
            node->data = flint_malloc(sizeof(fq_nmod_poly_struct));
            fq_nmod_poly_init(node->data, ctx->fqctx);
            fq_nmod_poly_zero(node->data, ctx->fqctx);
        }

        fq_nmod_poly_fit_length(t, 1, ctx->fqctx);
        n_fq_get_fq_nmod(t->coeffs + 0, Bcoeff + d*i, ctx->fqctx);
        t->length = 1;
        for (k = 0; k < k_len; k++)
        {
            if ((Bexp[N*i + offs[k]] & masks[k]) != WORD(0))
            {
                fq_nmod_poly_mul(t2, t, powers + k, ctx->fqctx);
                fq_nmod_poly_swap(t, t2, ctx->fqctx);
            }
        }
        fq_nmod_poly_add(t2, t, node->data, ctx->fqctx);
        fq_nmod_poly_swap(t2, node->data, ctx->fqctx);
    }
    fmpz_clear(main_exp);
    fq_nmod_poly_clear(t, ctx->fqctx);
    fq_nmod_poly_clear(t2, ctx->fqctx);

    for (k = 0; k < k_len; k++)
        fq_nmod_poly_clear(powers + k, ctx->fqctx);

    /* use tree method to evaluate in the main variable */
    fmpz_init(s);
    if (!_rbnode_clear_mp(tree, tree->head->left, s, A, C[main_var], ctx))
        success = 0;
    fmpz_clear(s);

cleanup_degrees:

    for (i = 0; i < nvars; i++)
        fmpz_clear(degrees + i);

    TMP_END;

    return success;
}


int fq_nmod_mpoly_compose_fq_nmod_poly(fq_nmod_poly_t A,
                    const fq_nmod_mpoly_t B, fq_nmod_poly_struct * const * C,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    if (B->length == 0)
    {
        fq_nmod_poly_zero(A, ctx->fqctx);
        return 1;
    }

    if (B->bits <= FLINT_BITS)
    {
        return _fq_nmod_mpoly_compose_fq_nmod_poly_sp(A, B, C, ctx);
    }
    else
    {
        return _fq_nmod_mpoly_compose_fq_nmod_poly_mp(A, B, C, ctx);
    }
}
