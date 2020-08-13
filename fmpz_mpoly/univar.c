/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mpoly.h"


void fmpz_mpoly_univar_init(fmpz_mpoly_univar_t A, const fmpz_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
}


void fmpz_mpoly_univar_clear(fmpz_mpoly_univar_t A, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
    {
        fmpz_mpoly_clear(A->coeffs + i, ctx);
        fmpz_clear(A->exps + i);
    }

    if (A->coeffs)
        flint_free(A->coeffs);

    if (A->exps)
        flint_free(A->exps);
}


void fmpz_mpoly_univar_fit_length(fmpz_mpoly_univar_t A,
                                      slong length, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->exps = (fmpz *) flint_malloc(new_alloc*sizeof(fmpz));
            A->coeffs = (fmpz_mpoly_struct *) flint_malloc(
                                          new_alloc*sizeof(fmpz_mpoly_struct));
        }
        else
        {
            A->exps = (fmpz *) flint_realloc(A->exps, new_alloc*sizeof(fmpz));
            A->coeffs = (fmpz_mpoly_struct *) flint_realloc(A->coeffs,
                                          new_alloc*sizeof(fmpz_mpoly_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fmpz_init(A->exps + i);
            fmpz_mpoly_init(A->coeffs + i, ctx);
        }
        A->alloc = new_alloc;
    }
}

void fmpz_mpoly_univar_assert_canonical(fmpz_mpoly_univar_t A, const fmpz_mpoly_ctx_t ctx)
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
        fmpz_mpoly_assert_canonical(A->coeffs + i, ctx);
}

void fmpz_mpoly_univar_print_pretty(const fmpz_mpoly_univar_t A,
                                   const char ** x, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    if (A->length == 0)
        flint_printf("0");
    for (i = 0; i < A->length; i++)
    {
        if (i != 0)
            flint_printf("+");
        flint_printf("(");
        fmpz_mpoly_print_pretty(A->coeffs + i,x,ctx);
        flint_printf(")*X^");
        fmpz_print(A->exps + i);
    }
}

static void _mpoly_rbnode_clear_sp(
    fmpz_mpoly_univar_t A,
    mpoly_rbtree_t tree,
    mpoly_rbnode_t node)
{
    mpoly_rbnode_struct * left = node->left;

    if (node->right != tree->null)
        _mpoly_rbnode_clear_sp(A, tree, node->right);

    FLINT_ASSERT(A->length < A->alloc);

    fmpz_set_si(A->exps + A->length, node->key);
    fmpz_mpoly_swap(A->coeffs + A->length, (fmpz_mpoly_struct *)(node->data), NULL);
    A->length++;

    fmpz_mpoly_clear(node->data, NULL);
    flint_free(node->data);
    flint_free(node);

    if (left != tree->null)
        _mpoly_rbnode_clear_sp(A, tree, left);
}

static void _mpoly_rbnode_clear_mp(
    fmpz_mpoly_univar_t A,
    mpoly_rbtree_t tree,
    mpoly_rbnode_t node)
{
    mpoly_rbnode_struct * left = node->left;

    if (node->right != tree->null)
        _mpoly_rbnode_clear_mp(A, tree, node->right);

    FLINT_ASSERT(A->length < A->alloc);

    fmpz_swap(A->exps + A->length, (fmpz*)(&node->key));
    fmpz_mpoly_swap(A->coeffs + A->length, (fmpz_mpoly_struct *)(node->data), NULL);
    A->length++;

    fmpz_clear((fmpz*)(&node->key));
    fmpz_mpoly_clear(node->data, NULL);
    flint_free(node->data);
    flint_free(node);

    if (left != tree->null)
        _mpoly_rbnode_clear_mp(A, tree, left);
}


/* the coefficients of A should be constructed with the same bits as B */
void fmpz_mpoly_to_univar(fmpz_mpoly_univar_t A, const fmpz_mpoly_t B,
                                         slong var, const fmpz_mpoly_ctx_t ctx)
{
    flint_bitcnt_t bits = B->bits;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    slong shift, off;
    slong Blen = B->length;
    const fmpz * Bcoeff = B->coeffs;
    const ulong * Bexp = B->exps;
    slong i;
    int new;
    ulong * one;
    mpoly_rbtree_t tree;
    mpoly_rbnode_struct * node;
    fmpz_mpoly_struct * d;
#define LUT_limit (48)
    fmpz_mpoly_struct LUT[LUT_limit];
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
            fmpz_mpoly_init3(LUT + i, 4, bits, ctx);

        /* fill in tree/LUT from B */
        for (i = 0; i < Blen; i++)
        {
            ulong k = (Bexp[N*i + off] >> shift) & mask;
            if (k < LUT_limit)
            {
                d = LUT + k;
            }
            else
            {
                node = mpoly_rbtree_get(&new, tree, k);
                if (new)
                {
                    d = flint_malloc(sizeof(fmpz_mpoly_struct));
                    fmpz_mpoly_init3(d, 4, bits, ctx);
                    node->data = d;
                }
                else
                {
                    d = node->data;
                }
            }
            fmpz_mpoly_fit_length(d, d->length + 1, ctx);
            fmpz_set(d->coeffs + d->length, Bcoeff + i);
            mpoly_monomial_msub(d->exps + N*d->length, Bexp + N*i, k, one, N);
            d->length++;
        }

        /* clear out tree to A */
        fmpz_mpoly_univar_fit_length(A, tree->size + LUT_limit, ctx);
        A->length = 0;
        if (tree->size > 0)
            _mpoly_rbnode_clear_sp(A, tree, tree->head->left);

        for (i = LUT_limit - 1; i >= 0; i--)
        {
            d = LUT + i;
            if (d->length > 0)
            {
                FLINT_ASSERT(A->length < A->alloc);
                fmpz_set_si(A->exps + A->length, i);
                fmpz_mpoly_swap(A->coeffs + A->length, d, ctx);
                A->length++;
            }
            fmpz_mpoly_clear(d, ctx);
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
                d = flint_malloc(sizeof(fmpz_mpoly_struct));
                fmpz_mpoly_init3(d, 4, bits, ctx);
                node->data = d;
            }
            else
            {
                d = node->data;
            }
            fmpz_mpoly_fit_length(d, d->length + 1, ctx);
            fmpz_set(d->coeffs + d->length, Bcoeff + i);
            mpoly_monomial_msub_ui_array(d->exps + N*d->length, Bexp + N*i,
                                    Bexp + N*i + off, bits/FLINT_BITS, one, N);
            d->length++;
        }

        /* clear out tree to A */
        fmpz_mpoly_univar_fit_length(A, tree->size, ctx);
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
void fmpz_mpoly_from_univar_bits(fmpz_mpoly_t A, flint_bitcnt_t Abits,
            const fmpz_mpoly_univar_t B, slong var, const fmpz_mpoly_ctx_t ctx)
{
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
        fmpz_mpoly_fit_bits(A, Abits, ctx);
        A->bits = Abits;
        _fmpz_mpoly_set_length(A, 0, ctx);
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
        fmpz_mpoly_struct * Bi = B->coeffs + i;
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

    fmpz_mpoly_fit_length(A, total_len, ctx);
    fmpz_mpoly_fit_bits(A, Abits, ctx);
    A->bits = Abits;

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
            FLINT_ASSERT(Alen < A->alloc);
            mpoly_monomial_set(A->exps + N*Alen, heap[1].exp, N);
            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
            fmpz_set(A->coeffs + Alen, (B->coeffs + x->i)->coeffs + x->j);
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
            FLINT_ASSERT(Alen < A->alloc);
            mpoly_monomial_set(A->exps + N*Alen, heap[1].exp, N);
            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
            fmpz_set(A->coeffs + Alen, (B->coeffs + x->i)->coeffs + x->j);
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

    _fmpz_mpoly_set_length(A, Alen, ctx);

    for (i = 0; i < B->length; i++)
    {
        if (Btexp[i] != (B->coeffs + i)->exps)
            flint_free(Btexp[i]);
    }

    TMP_END;
}

void fmpz_mpoly_from_univar(fmpz_mpoly_t A, const fmpz_mpoly_univar_t B,
                                         slong var, const fmpz_mpoly_ctx_t ctx)
{
    slong n = ctx->minfo->nfields;
    flint_bitcnt_t bits;
    slong i;
    fmpz * gen_fields, * tmp_fields, * max_fields;
    TMP_INIT;

    if (B->length == 0)
    {
        fmpz_mpoly_zero(A, ctx);
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
        fmpz_mpoly_struct * Bi = B->coeffs + i;
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

    fmpz_mpoly_from_univar_bits(A, bits, B, var, ctx);
}

