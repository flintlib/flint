/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_univar_init(
    fq_nmod_mpoly_univar_t A,
    const fq_nmod_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
}


void fq_nmod_mpoly_univar_clear(fq_nmod_mpoly_univar_t A,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = A->alloc - 1; i >= 0; i--)
    {
        fq_nmod_mpoly_clear(A->coeffs + i, ctx);
        fmpz_clear(A->exps + i);
    }

    if (A->coeffs)
        flint_free(A->coeffs);

    if (A->exps)
        flint_free(A->exps);
}


void fq_nmod_mpoly_univar_fit_length(fq_nmod_mpoly_univar_t A,
                                   slong length, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->exps = (fmpz *) flint_malloc(new_alloc*sizeof(fmpz));
            A->coeffs = (fq_nmod_mpoly_struct *) flint_malloc(
                                       new_alloc*sizeof(fq_nmod_mpoly_struct));
        }
        else
        {
            A->exps = (fmpz *) flint_realloc(A->exps, new_alloc*sizeof(fmpz));
            A->coeffs = (fq_nmod_mpoly_struct *) flint_realloc(A->coeffs,
                                       new_alloc*sizeof(fq_nmod_mpoly_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fmpz_init(A->exps + i);
            fq_nmod_mpoly_init(A->coeffs + i, ctx);
        }
        A->alloc = new_alloc;
    }
}

void fq_nmod_mpoly_univar_assert_canonical(fq_nmod_mpoly_univar_t A,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;

    for (i = 0; i + 1 < A->length; i++)
    {
        if (fmpz_cmp(A->exps + i, A->exps + i + 1) <= 0
            || fmpz_sgn(A->exps + i) < 0
            || fmpz_sgn(A->exps + i + 1) < 0)
        {
            flint_throw(FLINT_ERROR, "Univariate polynomial exponents out of order");
        }
    }

    for (i = 0; i < A->length; i++)
        fq_nmod_mpoly_assert_canonical(A->coeffs + i, ctx);
}

void fq_nmod_mpoly_univar_print_pretty(const fq_nmod_mpoly_univar_t A,
                                const char ** x, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    if (A->length == 0)
        flint_printf("0");
    for (i = 0; i < A->length; i++)
    {
        if (i != 0)
            flint_printf("+");
        flint_printf("(");
        fq_nmod_mpoly_print_pretty(A->coeffs + i,x,ctx);
        flint_printf(")*X^");
        fmpz_print(A->exps + i);
    }
}

static void _mpoly_rbnode_clear_sp(
    fq_nmod_mpoly_univar_t A,
    mpoly_rbtree_t tree,
    mpoly_rbnode_t node)
{
    mpoly_rbnode_struct * left = node->left;

    if (node->right != tree->null)
        _mpoly_rbnode_clear_sp(A, tree, node->right);

    FLINT_ASSERT(A->length < A->alloc);

    fmpz_set_si(A->exps + A->length, node->key);
    fq_nmod_mpoly_swap(A->coeffs + A->length, (fq_nmod_mpoly_struct *)(node->data), NULL);
    A->length++;

    fq_nmod_mpoly_clear(node->data, NULL);
    flint_free(node->data);
    flint_free(node);

    if (left != tree->null)
        _mpoly_rbnode_clear_sp(A, tree, left);
}

static void _mpoly_rbnode_clear_mp(
    fq_nmod_mpoly_univar_t A,
    mpoly_rbtree_t tree,
    mpoly_rbnode_t node)
{
    mpoly_rbnode_struct * left = node->left;

    if (node->right != tree->null)
        _mpoly_rbnode_clear_mp(A, tree, node->right);

    FLINT_ASSERT(A->length < A->alloc);

    fmpz_swap(A->exps + A->length, (fmpz*)(&node->key));
    fq_nmod_mpoly_swap(A->coeffs + A->length, (fq_nmod_mpoly_struct *)(node->data), NULL);
    A->length++;

    fmpz_clear((fmpz*)(&node->key));
    fq_nmod_mpoly_clear(node->data, NULL);
    flint_free(node->data);
    flint_free(node);

    if (left != tree->null)
        _mpoly_rbnode_clear_mp(A, tree, left);
}

/* the coefficients of A should be constructed with the same bits as B */
void fq_nmod_mpoly_to_univar(
    fq_nmod_mpoly_univar_t A,
    const fq_nmod_mpoly_t B,
    slong var,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    flint_bitcnt_t bits = B->bits;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    slong shift, off;
    slong Blen = B->length;
    const mp_limb_t * Bcoeffs = B->coeffs;
    const ulong * Bexp = B->exps;
    slong i;
    int new;
    ulong * one;
    mpoly_rbtree_t tree;
    mpoly_rbnode_struct * node;
    fq_nmod_mpoly_struct * t;
#define LUT_limit (48)
    fq_nmod_mpoly_struct LUT[LUT_limit];
    TMP_INIT;

    if (B->length == 0)
    {
        A->length = 0;
        return;
    }

    TMP_START;

    one = (ulong*) TMP_ALLOC(N*sizeof(ulong));

    mpoly_rbtree_init(tree);

    if (bits <= FLINT_BITS)
    {
        ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
        mpoly_gen_monomial_offset_shift_sp(one, &off, &shift,
                                                        var, bits, ctx->minfo);
        for (i = 0; i < LUT_limit; i++)
            fq_nmod_mpoly_init3(LUT + i, 4, bits, ctx);

        /* fill in tree/LUT from B */
        for (i = 0; i < Blen; i++)
        {
            ulong k = (Bexp[N*i + off] >> shift) & mask;
            if (k < LUT_limit)
            {
                t = LUT + k;
            }
            else
            {
                node = mpoly_rbtree_get(&new, tree, k);
                if (new)
                {
                    t = flint_malloc(sizeof(fq_nmod_mpoly_struct));
                    fq_nmod_mpoly_init3(t, 4, bits, ctx);
                    node->data = t;
                }
                else
                {
                    t = node->data;
                }
            }
            fq_nmod_mpoly_fit_length(t, t->length + 1, ctx);
            _n_fq_set(t->coeffs + d*t->length, Bcoeffs + d*i, d);
            mpoly_monomial_msub(t->exps + N*t->length, Bexp + N*i, k, one, N);
            t->length++;
        }

        /* clear out tree to A */
        fq_nmod_mpoly_univar_fit_length(A, tree->size + LUT_limit, ctx);
        A->length = 0;
        if (tree->size > 0)
            _mpoly_rbnode_clear_sp(A, tree, tree->head->left);

        for (i = LUT_limit - 1; i >= 0; i--)
        {
            t = LUT + i;
            if (t->length > 0)
            {
                FLINT_ASSERT(A->length < A->alloc);
                fmpz_set_si(A->exps + A->length, i);
                fq_nmod_mpoly_swap(A->coeffs + A->length, t, ctx);
                A->length++;
            }
            fq_nmod_mpoly_clear(t, ctx);
        }
    }
    else
    {
        fmpz_t k;
        fmpz_init(k);

        off = mpoly_gen_monomial_offset_mp(one, var, bits, ctx->minfo);

        /* fill in tree from B */
        for (i = 0; i < Blen; i++)
        {
            fmpz_set_ui_array(k, Bexp + N*i + off, bits/FLINT_BITS);

            node = mpoly_rbtree_get_fmpz(&new, tree, k);
            if (new)
            {
                t = flint_malloc(sizeof(fq_nmod_mpoly_struct));
                fq_nmod_mpoly_init3(t, 4, bits, ctx);
                node->data = t;
            }
            else
            {
                t = node->data;
            }
            fq_nmod_mpoly_fit_length(t, t->length + 1, ctx);
            _n_fq_set(t->coeffs + d*t->length, Bcoeffs + d*i, d);
            mpoly_monomial_msub_ui_array(t->exps + N*t->length, Bexp + N*i,
                                    Bexp + N*i + off, bits/FLINT_BITS, one, N);
            t->length++;
        }

        /* clear out tree to A */
        fq_nmod_mpoly_univar_fit_length(A, tree->size, ctx);
        A->length = 0;
        FLINT_ASSERT(tree->size > 0);
        _mpoly_rbnode_clear_mp(A, tree, tree->head->left);

        fmpz_clear(k);
    }

    TMP_END;
}


/*
    Currently this function does not work if the coefficients depend on "var".
    The assertion x->next == NULL would need to be replaced by a loop.
    Other asserts would need to be removed as well.
*/
void fq_nmod_mpoly_from_univar_bits(
    fq_nmod_mpoly_t A,
    flint_bitcnt_t Abits,
    const fq_nmod_mpoly_univar_t B,
    slong var,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong N = mpoly_words_per_exp(Abits, ctx->minfo);
    slong i;
    slong next_loc, heap_len = 1;
    ulong * cmpmask;
    slong total_len;
    mpoly_heap_s * heap;
    slong Alen;
    ulong ** Btexp;
    ulong * exp;
    ulong * one;
    mpoly_heap_t * chain, * x;
    TMP_INIT;

    if (B->length == 0)
    {
        fq_nmod_mpoly_fit_length_reset_bits(A, 0, Abits, ctx);
        A->length = 0;
        return;
    }

    TMP_START;

    /* pack everything into Abits */

    one = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    Btexp = (ulong **) TMP_ALLOC(B->length*sizeof(ulong*));

    total_len = 0;
    for (i = 0; i < B->length; i++)
    {
        fq_nmod_mpoly_struct * Bi = B->coeffs + i;
        total_len += Bi->length;
        Btexp[i] = Bi->exps;
        if (Abits != Bi->bits)
        {
            Btexp[i] = (ulong *) flint_malloc(N*Bi->length*sizeof(ulong));
            if (!mpoly_repack_monomials(Btexp[i], Abits,
                                   Bi->exps, Bi->bits, Bi->length, ctx->minfo))
            {
                FLINT_ASSERT(0 && "repack does not fit");
            }
        }
    }

    fq_nmod_mpoly_fit_length_reset_bits(A, total_len, Abits, ctx);

    next_loc = B->length + 2;
    heap = (mpoly_heap_s *) TMP_ALLOC((B->length + 1)*sizeof(mpoly_heap_s));
    exp = (ulong *) TMP_ALLOC(B->length*N*sizeof(ulong));
    chain = (mpoly_heap_t *) TMP_ALLOC(B->length*sizeof(mpoly_heap_t));

    mpoly_get_cmpmask(cmpmask, N, Abits, ctx->minfo);

    Alen = 0;
    if (Abits <= FLINT_BITS)
    {
        mpoly_gen_monomial_sp(one, var, Abits, ctx->minfo);

        for (i = 0; i < B->length; i++)
        {
            FLINT_ASSERT(fmpz_fits_si(B->exps + i));
            x = chain + i;
            x->i = i;
            x->j = 0;
            x->next = NULL;
            mpoly_monomial_madd(exp + N*x->i, Btexp[x->i] + N*x->j,
                                             fmpz_get_si(B->exps + i), one, N);
            _mpoly_heap_insert(heap, exp + N*i, x, &next_loc, &heap_len,
                                                                   N, cmpmask);
        }

        while (heap_len > 1)
        {
            mpoly_monomial_set(A->exps + N*Alen, heap[1].exp, N);
            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
            _n_fq_set(A->coeffs + d*Alen, (B->coeffs + x->i)->coeffs + d*x->j, d);
            Alen++;

            FLINT_ASSERT(x->next == NULL);

            if (x->j + 1 < (B->coeffs + x->i)->length)
            {
                FLINT_ASSERT(fmpz_fits_si(B->exps + x->i));
                x->j = x->j + 1;
                x->next = NULL;
                mpoly_monomial_madd(exp + N*x->i, Btexp[x->i] + N*x->j,
                                          fmpz_get_si(B->exps + x->i), one, N);
                _mpoly_heap_insert(heap, exp + N*x->i, x, &next_loc, &heap_len,
                                                                   N, cmpmask);
            }
        }
    }
    else
    {
        mpoly_gen_monomial_offset_mp(one, var, Abits, ctx->minfo);

        for (i = 0; i < B->length; i++)
        {
            x = chain + i;
            x->i = i;
            x->j = 0;
            x->next = NULL;
            mpoly_monomial_madd_fmpz(exp + N*x->i, Btexp[x->i] + N*x->j,
                                                          B->exps + i, one, N);
            _mpoly_heap_insert(heap, exp + N*i, x, &next_loc, &heap_len,
                                                                   N, cmpmask);
        }

        while (heap_len > 1)
        {
            mpoly_monomial_set(A->exps + N*Alen, heap[1].exp, N);
            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
            _n_fq_set(A->coeffs + d*Alen, (B->coeffs + x->i)->coeffs + d*x->j, d);
            Alen++;

            FLINT_ASSERT(x->next == NULL);

            if (x->j + 1 < (B->coeffs + x->i)->length)
            {
                x->j = x->j + 1;
                x->next = NULL;
                mpoly_monomial_madd_fmpz(exp + N*x->i, Btexp[x->i] + N*x->j,
                                                       B->exps + x->i, one, N);
                _mpoly_heap_insert(heap, exp + N*x->i, x, &next_loc, &heap_len,
                                                                   N, cmpmask);
            }
        }
    }

    A->length = Alen;

    for (i = 0; i < B->length; i++)
    {
        if (Btexp[i] != (B->coeffs + i)->exps)
            flint_free(Btexp[i]);
    }

    TMP_END;
}

void fq_nmod_mpoly_from_univar(fq_nmod_mpoly_t A, const fq_nmod_mpoly_univar_t B,
                                      slong var, const fq_nmod_mpoly_ctx_t ctx)
{
    slong n = ctx->minfo->nfields;
    flint_bitcnt_t bits;
    slong i;
    fmpz * gen_fields, * tmp_fields, * max_fields;
    TMP_INIT;

    if (B->length == 0)
    {
        fq_nmod_mpoly_zero(A, ctx);
        return;
    }

    TMP_START;

    /* find bits required to represent result */
    gen_fields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    tmp_fields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    max_fields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_init(gen_fields + i);
        fmpz_init(tmp_fields + i);
        fmpz_init(max_fields + i);
    }

    mpoly_gen_fields_fmpz(gen_fields, var, ctx->minfo);

    for (i = 0; i < B->length; i++)
    {
        fq_nmod_mpoly_struct * Bi = B->coeffs + i;
        mpoly_max_fields_fmpz(tmp_fields, Bi->exps, Bi->length, Bi->bits, ctx->minfo);
        _fmpz_vec_scalar_addmul_fmpz(tmp_fields, gen_fields, n, B->exps + i);
        _fmpz_vec_max_inplace(max_fields, tmp_fields, n);
    }
    bits = _fmpz_vec_max_bits(max_fields, n);
    bits = FLINT_MAX(MPOLY_MIN_BITS, bits + 1);
    bits = mpoly_fix_bits(bits, ctx->minfo);

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_clear(gen_fields + i);
        fmpz_clear(tmp_fields + i);
        fmpz_clear(max_fields + i);
    }
    TMP_END;

    fq_nmod_mpoly_from_univar_bits(A, bits, B, var, ctx);
}
