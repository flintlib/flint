/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
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

void fq_nmod_mpoly_univar_set_coeff_ui(
    fq_nmod_mpoly_univar_t A,
    ulong e,
    const fq_nmod_mpoly_t c,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, j;

    for (i = A->length; i >= 0; i--)
    {
        int cmp = i > 0 ? fmpz_cmp_ui(A->exps + i - 1, e) : 1;

        if (cmp > 0)
        {
            if (fq_nmod_mpoly_is_zero(c, ctx))
                return;

            fq_nmod_mpoly_univar_fit_length(A, A->length + 1, ctx);

            for (j = A->length; j > i; j--)
            {
                fq_nmod_mpoly_swap(A->coeffs + j, A->coeffs + j + 1, ctx);
                fmpz_swap(A->exps + j, A->exps + j + 1);
            }

            A->length++;

            fmpz_set_ui(A->exps + i, e);
            fq_nmod_mpoly_set(A->coeffs + i, c, ctx);
            return;
        }
        else if (cmp == 0)
        {
            fq_nmod_mpoly_set(A->coeffs + i, c, ctx);

            if (!fq_nmod_mpoly_is_zero(A->coeffs + i, ctx))
                return;

            A->length--;

            for (j = i; j < A->length; j++)
            {
                fq_nmod_mpoly_swap(A->coeffs + j, A->coeffs + j + 1, ctx);
                fmpz_swap(A->exps + j, A->exps + j + 1);
            }
        }
    }

    FLINT_ASSERT(0 && "unreachable");
    return;
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

static void _tree_data_clear_sp(
    fq_nmod_mpoly_univar_t A,
    mpoly_rbtree_ui_t tree,
    slong idx,
    const fq_nmod_mpoly_ctx_t ctx)
{
    mpoly_rbnode_ui_struct * nodes = tree->nodes + 2;
    fq_nmod_mpoly_struct * data = (fq_nmod_mpoly_struct *) tree->data;

    if (idx < 0)
        return;

    _tree_data_clear_sp(A, tree, nodes[idx].right, ctx);

    FLINT_ASSERT(A->length < A->alloc);

    fmpz_set_ui(A->exps + A->length, nodes[idx].key);
    fq_nmod_mpoly_swap(A->coeffs + A->length, data + idx, ctx);
    A->length++;

    fq_nmod_mpoly_clear(data + idx, ctx);

    _tree_data_clear_sp(A, tree, nodes[idx].left, ctx);
}

static void _tree_data_clear_mp(
    fq_nmod_mpoly_univar_t A,
    mpoly_rbtree_fmpz_t tree,
    slong idx,
    const fq_nmod_mpoly_ctx_t ctx)
{
    mpoly_rbnode_fmpz_struct * nodes = tree->nodes + 2;
    fq_nmod_mpoly_struct * data = (fq_nmod_mpoly_struct *) tree->data;

    if (idx < 0)
        return;

    _tree_data_clear_mp(A, tree, nodes[idx].right, ctx);

    FLINT_ASSERT(A->length < A->alloc);

    fmpz_set(A->exps + A->length, nodes[idx].key);
    fq_nmod_mpoly_swap(A->coeffs + A->length, data + idx, ctx);
    A->length++;

    fq_nmod_mpoly_clear(data + idx, ctx);

    _tree_data_clear_mp(A, tree, nodes[idx].left, ctx);
}


/* the coefficients of A should be constructed with the same bits as B */
void fq_nmod_mpoly_to_univar(fq_nmod_mpoly_univar_t A, const fq_nmod_mpoly_t B,
                                         slong var, const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    flint_bitcnt_t bits = B->bits;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    slong shift, off;
    slong Blen = B->length;
    const mp_limb_t * Bcoeff = B->coeffs;
    const ulong * Bexp = B->exps;
    slong i;
    int its_new;
    ulong * one;
#define LUT_limit (48)
    fq_nmod_mpoly_struct LUT[LUT_limit];

    if (B->length == 0)
    {
        A->length = 0;
        return;
    }

    one = FLINT_ARRAY_ALLOC(N, ulong);

    if (bits <= FLINT_BITS)
    {
        slong Alen;
        mpoly_rbtree_ui_t tree;
        fq_nmod_mpoly_struct * t;
        ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);

        mpoly_rbtree_ui_init(tree, sizeof(fq_nmod_mpoly_struct));

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
                t = mpoly_rbtree_ui_lookup(tree, &its_new, k);
                if (its_new)
                    fq_nmod_mpoly_init3(t, 4, bits, ctx);
            }
            fq_nmod_mpoly_fit_length(t, t->length + 1, ctx);
            _n_fq_set(t->coeffs + d*t->length, Bcoeff + d*i, d);
            mpoly_monomial_msub(t->exps + N*t->length, Bexp + N*i, k, one, N);
            t->length++;
        }

        /* clear out tree to A */

        Alen = tree->length;
        for (i = LUT_limit - 1; i >= 0; i--)
            Alen += (LUT[i].length > 0);

        fq_nmod_mpoly_univar_fit_length(A, Alen, ctx);
        A->length = 0;

        _tree_data_clear_sp(A, tree, mpoly_rbtree_ui_head(tree), ctx);

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

        mpoly_rbtree_ui_clear(tree);
    }
    else
    {
        mpoly_rbtree_fmpz_t tree;
        fq_nmod_mpoly_struct * t;
        fmpz_t k;

        fmpz_init(k);
        mpoly_rbtree_fmpz_init(tree, sizeof(fq_nmod_mpoly_struct));

        off = mpoly_gen_monomial_offset_mp(one, var, bits, ctx->minfo);

        /* fill in tree from B */
        for (i = 0; i < Blen; i++)
        {
            fmpz_set_ui_array(k, Bexp + N*i + off, bits/FLINT_BITS);

            t = mpoly_rbtree_fmpz_lookup(tree, &its_new, k);
            if (its_new)
                fq_nmod_mpoly_init3(t, 4, bits, ctx);

            fq_nmod_mpoly_fit_length(t, t->length + 1, ctx);
            _n_fq_set(t->coeffs + d*t->length, Bcoeff + d*i, d);
            mpoly_monomial_msub_ui_array(t->exps + N*t->length, Bexp + N*i,
                                    Bexp + N*i + off, bits/FLINT_BITS, one, N);
            t->length++;
        }

        /* clear out tree to A */
        fq_nmod_mpoly_univar_fit_length(A, tree->length, ctx);
        A->length = 0;

        _tree_data_clear_mp(A, tree, mpoly_rbtree_fmpz_head(tree), ctx);

        fmpz_clear(k);
        mpoly_rbtree_fmpz_clear(tree);
    }

    flint_free(one);
}


/*
    Currently this function does not work if the coefficients depend on "var".
    The assertion x->next == NULL would need to be replaced by a loop.
    Other asserts would need to be removed as well.
*/
void _fq_nmod_mpoly_from_univar(
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

    _fq_nmod_mpoly_from_univar(A, bits, B, var, ctx);
}

#define COEFF(A, i) ((void*)(A->coeffs + (i)*R->elem_size))

static void mpoly_univar_set_fq_nmod_mpoly_univar(
    mpoly_univar_t A,
    mpoly_void_ring_t R,
    const fq_nmod_mpoly_univar_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;

    mpoly_univar_fit_length(A, B->length, R);
    A->length = B->length;

    for (i = B->length - 1; i >= 0; i--)
    {
        fmpz_set(A->exps + i, B->exps + i);
        fq_nmod_mpoly_set(COEFF(A, i), B->coeffs + i, ctx);
    }
}

static void mpoly_univar_swap_fq_nmod_mpoly_univar(
    mpoly_univar_t A,
    mpoly_void_ring_t R,
    fq_nmod_mpoly_univar_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;

    mpoly_univar_fit_length(A, B->length, R);
    fq_nmod_mpoly_univar_fit_length(B, A->length, ctx);

    for (i = FLINT_MAX(A->length, B->length) - 1; i >= 0; i--)
    {
        fmpz_swap(A->exps + i, B->exps + i);
        fq_nmod_mpoly_swap(COEFF(A, i), B->coeffs + i, ctx);
    }

    FLINT_SWAP(slong, A->length, B->length);
}

int fq_nmod_mpoly_univar_pseudo_gcd(
    fq_nmod_mpoly_univar_t gx,
    const fq_nmod_mpoly_univar_t ax,
    const fq_nmod_mpoly_univar_t bx,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    mpoly_void_ring_t R;
    mpoly_univar_t Ax, Bx, Gx;

    mpoly_void_ring_init_fq_nmod_mpoly_ctx(R, ctx);
    mpoly_univar_init(Ax, R);
    mpoly_univar_init(Bx, R);
    mpoly_univar_init(Gx, R);
    mpoly_univar_set_fq_nmod_mpoly_univar(Ax, R, ax, ctx);
    mpoly_univar_set_fq_nmod_mpoly_univar(Bx, R, bx, ctx);

    success = mpoly_univar_pseudo_gcd_ducos(Gx, Ax, Bx, R);

    if (success)
        mpoly_univar_swap_fq_nmod_mpoly_univar(Gx, R, gx, ctx);

    mpoly_univar_clear(Ax, R);
    mpoly_univar_clear(Bx, R);
    mpoly_univar_clear(Gx, R);

    return success;
}

int fq_nmod_mpoly_univar_resultant(
    fq_nmod_mpoly_t d,
    const fq_nmod_mpoly_univar_t ax,
    const fq_nmod_mpoly_univar_t bx,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    mpoly_void_ring_t R;
    mpoly_univar_t Ax, Bx;

    mpoly_void_ring_init_fq_nmod_mpoly_ctx(R, ctx);
    mpoly_univar_init(Ax, R);
    mpoly_univar_init(Bx, R);
    mpoly_univar_set_fq_nmod_mpoly_univar(Ax, R, ax, ctx);
    mpoly_univar_set_fq_nmod_mpoly_univar(Bx, R, bx, ctx);

    success = mpoly_univar_resultant(d, Ax, Bx, R);

    mpoly_univar_clear(Ax, R);
    mpoly_univar_clear(Bx, R);

    return success;
}

int fq_nmod_mpoly_univar_discriminant(
    fq_nmod_mpoly_t d,
    const fq_nmod_mpoly_univar_t fx,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    mpoly_void_ring_t R;
    mpoly_univar_t Fx;

    mpoly_void_ring_init_fq_nmod_mpoly_ctx(R, ctx);
    mpoly_univar_init(Fx, R);
    mpoly_univar_set_fq_nmod_mpoly_univar(Fx, R, fx, ctx);

    success = mpoly_univar_discriminant(d, Fx, R);

    mpoly_univar_clear(Fx, R);

    return success;
}

