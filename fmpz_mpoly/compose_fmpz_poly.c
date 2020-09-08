/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

static int _fmpz_poly_pow_fmpz_is_not_feasible(const fmpz_poly_t b, const fmpz_t e)
{
    if (b->length < 2)
    {
        if (b->length == 1)
        {
            return _fmpz_pow_fmpz_is_not_feasible(fmpz_bits(b->coeffs + 0), e);
        }
        else
        {
            return 0;
        }
    }
    else
    {
        ulong limit = (ulong)(WORD_MAX)/(ulong)(2*sizeof(fmpz));
        return fmpz_cmp_ui(e, limit/(ulong)(b->length)) >= 0;
    }
}

static int _fmpz_poly_pow_ui_is_not_feasible(const fmpz_poly_t b, ulong e)
{
    if (b->length < 2)
    {
        if (b->length == 1)
        {
            return _fmpz_pow_ui_is_not_feasible(fmpz_bits(b->coeffs + 0), e);
        }
        else
        {
            return 0;
        }
    }
    else
    {
        ulong limit = (ulong)(WORD_MAX)/(ulong)(2*sizeof(fmpz));
        return e >= limit/(ulong)(b->length);
    }
}


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



int _fmpz_mpoly_compose_fmpz_poly_sp(fmpz_poly_t A, const fmpz_mpoly_t B,
                      fmpz_poly_struct * const * C, const fmpz_mpoly_ctx_t ctx)
{
    int success = 1;
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
        if (_fmpz_poly_pow_ui_is_not_feasible(C[i], degrees[i]))
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
    powers = (fmpz_poly_struct *) TMP_ALLOC(entries*sizeof(fmpz_poly_struct));

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
    mpoly_gen_offset_shift_sp(&main_off, &main_shift, main_var, bits, ctx->minfo);
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

cleanup_degrees:

    TMP_END;

    return success;
}


static int _rbnode_clear_mp(mpoly_rbtree_t tree, mpoly_rbnode_t node,
                       const fmpz_t s, fmpz_poly_t l, const fmpz_poly_t x,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    int success = 1;
    fmpz_poly_t r, xp;
    FLINT_ASSERT(fmpz_cmp(&node->key, s) >= 0);

    fmpz_poly_init(r);
    if (node->right != tree->null)
    {
        if (!_rbnode_clear_mp(tree, node->right, &node->key, r, x, ctx))
            success = 0;
    }

    fmpz_poly_zero(l);
    if (node->left != tree->null)
    {
        if (!_rbnode_clear_mp(tree, node->left, s, l, x, ctx))
            success = 0;
    }

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
        if (degree == 0)
        {
            fmpz_t t;
            fmpz_init(t);
            fmpz_poly_get_coeff_fmpz(t, x, 0);
            if (!fmpz_pow_fmpz(t ,t, &node->key))
                success = 0;
            fmpz_poly_set_fmpz(xp, t);
            fmpz_clear(t);
        }
        else if (degree > 0)
        {
            success = 0;
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

    return success;
}



int _fmpz_mpoly_compose_fmpz_poly_mp(fmpz_poly_t A, const fmpz_mpoly_t B,
                      fmpz_poly_struct * const * C, const fmpz_mpoly_ctx_t ctx)
{
    int success = 1;
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
    flint_bitcnt_t * bitcounts;
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
        if (_fmpz_poly_pow_fmpz_is_not_feasible(C[i], degrees + i))
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
    powers = (fmpz_poly_struct *) TMP_ALLOC(entries*sizeof(fmpz_poly_struct));

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
    main_off = mpoly_gen_offset_mp(main_var, bits, ctx->minfo);
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


int fmpz_mpoly_compose_fmpz_poly(fmpz_poly_t A, const fmpz_mpoly_t B,
                      fmpz_poly_struct * const * C, const fmpz_mpoly_ctx_t ctx)
{
    if (B->length == 0)
    {
        fmpz_poly_zero(A);
        return 1;
    }

    if (B->bits <= FLINT_BITS)
    {
        return _fmpz_mpoly_compose_fmpz_poly_sp(A, B, C, ctx);
    }
    else
    {
        return _fmpz_mpoly_compose_fmpz_poly_mp(A, B, C, ctx);
    }
}

