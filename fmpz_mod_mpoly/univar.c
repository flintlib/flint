/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"


void fmpz_mod_mpoly_univar_init(
    fmpz_mod_mpoly_univar_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
}


void fmpz_mod_mpoly_univar_clear(
    fmpz_mod_mpoly_univar_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
    {
        fmpz_mod_mpoly_clear(A->coeffs + i, ctx);
        fmpz_clear(A->exps + i);
    }

    if (A->coeffs)
        flint_free(A->coeffs);

    if (A->exps)
        flint_free(A->exps);
}


void fmpz_mod_mpoly_univar_fit_length(
    fmpz_mod_mpoly_univar_t A,
    slong length,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length <= old_alloc)
        return;

    A->exps = FLINT_ARRAY_REALLOC(A->exps, new_alloc, fmpz);
    A->coeffs = FLINT_ARRAY_REALLOC(A->coeffs, new_alloc, fmpz_mod_mpoly_struct);

    for (i = old_alloc; i < new_alloc; i++)
    {
        fmpz_init(A->exps + i);
        fmpz_mod_mpoly_init(A->coeffs + i, ctx);
    }

    A->alloc = new_alloc;
}

void fmpz_mod_mpoly_univar_set_coeff_ui(
    fmpz_mod_mpoly_univar_t A,
    ulong e,
    const fmpz_mod_mpoly_t c,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, j;

    for (i = A->length; i >= 0; i--)
    {
        int cmp = i > 0 ? fmpz_cmp_ui(A->exps + i - 1, e) : 1;

        if (cmp > 0)
        {
            if (fmpz_mod_mpoly_is_zero(c, ctx))
                return;

            fmpz_mod_mpoly_univar_fit_length(A, A->length + 1, ctx);

            for (j = A->length; j > i; j--)
            {
                fmpz_mod_mpoly_swap(A->coeffs + j, A->coeffs + j + 1, ctx);
                fmpz_swap(A->exps + j, A->exps + j + 1);
            }

            A->length++;

            fmpz_set_ui(A->exps + i, e);
            fmpz_mod_mpoly_set(A->coeffs + i, c, ctx);
            return;
        }
        else if (cmp == 0)
        {
            fmpz_mod_mpoly_set(A->coeffs + i, c, ctx);

            if (!fmpz_mod_mpoly_is_zero(A->coeffs + i, ctx))
                return;

            A->length--;

            for (j = i; j < A->length; j++)
            {
                fmpz_mod_mpoly_swap(A->coeffs + j, A->coeffs + j + 1, ctx);
                fmpz_swap(A->exps + j, A->exps + j + 1);
            }            
        }
    }

    FLINT_ASSERT(0 && "unreachable");
    return;
}



void fmpz_mod_mpoly_univar_assert_canonical(
    fmpz_mod_mpoly_univar_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;

    for (i = 0; i + 1 < A->length; i++)
    {
        if (fmpz_cmp(A->exps + i, A->exps + i + 1) <= 0 ||
            fmpz_sgn(A->exps + i) < 0 ||
            fmpz_sgn(A->exps + i + 1) < 0)
        {
            flint_throw(FLINT_ERROR,
                               "Univariate polynomial exponents out of order");
        }
    }

    for (i = 0; i < A->length; i++)
        fmpz_mod_mpoly_assert_canonical(A->coeffs + i, ctx);
}

void fmpz_mod_mpoly_univar_print_pretty(
    const fmpz_mod_mpoly_univar_t A,
    const char ** x,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;
    if (A->length == 0)
        flint_printf("0");
    for (i = 0; i < A->length; i++)
    {
        if (i != 0)
            flint_printf("+");
        flint_printf("(");
        fmpz_mod_mpoly_print_pretty(A->coeffs + i,x,ctx);
        flint_printf(")*X^");
        fmpz_print(A->exps + i);
    }
}

static void _mpoly_rbnode_clear_sp(
    fmpz_mod_mpoly_univar_t A,
    mpoly_rbtree_t tree,
    mpoly_rbnode_t node)
{
    mpoly_rbnode_struct * left = node->left;

    if (node->right != tree->null)
        _mpoly_rbnode_clear_sp(A, tree, node->right);

    FLINT_ASSERT(A->length < A->alloc);

    fmpz_set_si(A->exps + A->length, node->key);
    fmpz_mod_mpoly_swap(A->coeffs + A->length,
                                  (fmpz_mod_mpoly_struct *)(node->data), NULL);
    A->length++;

    fmpz_mod_mpoly_clear(node->data, NULL);
    flint_free(node->data);
    flint_free(node);

    if (left != tree->null)
        _mpoly_rbnode_clear_sp(A, tree, left);
}

static void _mpoly_rbnode_clear_mp(
    fmpz_mod_mpoly_univar_t A,
    mpoly_rbtree_t tree,
    mpoly_rbnode_t node)
{
    mpoly_rbnode_struct * left = node->left;

    if (node->right != tree->null)
        _mpoly_rbnode_clear_mp(A, tree, node->right);

    FLINT_ASSERT(A->length < A->alloc);

    fmpz_swap(A->exps + A->length, (fmpz*)(&node->key));
    fmpz_mod_mpoly_swap(A->coeffs + A->length,
                                  (fmpz_mod_mpoly_struct *)(node->data), NULL);
    A->length++;

    fmpz_clear((fmpz*)(&node->key));
    fmpz_mod_mpoly_clear(node->data, NULL);
    flint_free(node->data);
    flint_free(node);

    if (left != tree->null)
        _mpoly_rbnode_clear_mp(A, tree, left);
}


/* the coefficients of A should be constructed with the same bits as B */
void fmpz_mod_mpoly_to_univar(
    fmpz_mod_mpoly_univar_t A,
    const fmpz_mod_mpoly_t B,
    slong var,
    const fmpz_mod_mpoly_ctx_t ctx)
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
    fmpz_mod_mpoly_struct * d;
#define LUT_limit (48)
    fmpz_mod_mpoly_struct LUT[LUT_limit];
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
            fmpz_mod_mpoly_init3(LUT + i, 4, bits, ctx);

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
                    d = flint_malloc(sizeof(fmpz_mod_mpoly_struct));
                    fmpz_mod_mpoly_init3(d, 4, bits, ctx);
                    node->data = d;
                }
                else
                {
                    d = node->data;
                }
            }
            fmpz_mod_mpoly_fit_length(d, d->length + 1, ctx);
            fmpz_set(d->coeffs + d->length, Bcoeff + i);
            mpoly_monomial_msub(d->exps + N*d->length, Bexp + N*i, k, one, N);
            d->length++;
        }

        /* clear out tree to A */
        fmpz_mod_mpoly_univar_fit_length(A, tree->size + LUT_limit, ctx);
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
                fmpz_mod_mpoly_swap(A->coeffs + A->length, d, ctx);
                A->length++;
            }
            fmpz_mod_mpoly_clear(d, ctx);
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
                d = flint_malloc(sizeof(fmpz_mod_mpoly_struct));
                fmpz_mod_mpoly_init3(d, 4, bits, ctx);
                node->data = d;
            }
            else
            {
                d = node->data;
            }
            fmpz_mod_mpoly_fit_length(d, d->length + 1, ctx);
            fmpz_set(d->coeffs + d->length, Bcoeff + i);
            mpoly_monomial_msub_ui_array(d->exps + N*d->length, Bexp + N*i,
                                    Bexp + N*i + off, bits/FLINT_BITS, one, N);
            d->length++;
        }

        /* clear out tree to A */
        fmpz_mod_mpoly_univar_fit_length(A, tree->size, ctx);
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
void _fmpz_mod_mpoly_from_univar(
    fmpz_mod_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpz_mod_mpoly_univar_t B,
    slong var,
    const fmpz_mod_mpoly_ctx_t ctx)
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
        fmpz_mod_mpoly_fit_length_reset_bits(A, 0, Abits, ctx);
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
        fmpz_mod_mpoly_struct * Bi = B->coeffs + i;
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

    fmpz_mod_mpoly_fit_length_reset_bits(A, total_len, Abits, ctx);

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
            FLINT_ASSERT(Alen < A->coeffs_alloc);
            mpoly_monomial_set(A->exps + N*Alen, heap[1].exp, N);
            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
            fmpz_set(A->coeffs + Alen, B->coeffs[x->i].coeffs + x->j);
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
            FLINT_ASSERT(Alen < A->coeffs_alloc);
            mpoly_monomial_set(A->exps + N*Alen, heap[1].exp, N);
            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
            fmpz_set(A->coeffs + Alen, B->coeffs[x->i].coeffs + x->j);
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

void fmpz_mod_mpoly_from_univar(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_univar_t B,
    slong var,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong n = ctx->minfo->nfields;
    flint_bitcnt_t bits;
    slong i;
    fmpz * gen_fields, * tmp_fields, * max_fields;
    TMP_INIT;

    if (B->length == 0)
    {
        fmpz_mod_mpoly_zero(A, ctx);
        return;
    }

    TMP_START;

    /* find bits required to represent result */
    gen_fields = TMP_ARRAY_ALLOC(3*ctx->minfo->nfields, fmpz);
    tmp_fields = gen_fields + ctx->minfo->nfields;
    max_fields = tmp_fields + ctx->minfo->nfields;
    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_init(gen_fields + i);
        fmpz_init(tmp_fields + i);
        fmpz_init(max_fields + i);
    }

    mpoly_gen_fields_fmpz(gen_fields, var, ctx->minfo);

    for (i = 0; i < B->length; i++)
    {
        fmpz_mod_mpoly_struct * Bi = B->coeffs + i;
        mpoly_max_fields_fmpz(tmp_fields, Bi->exps, Bi->length, Bi->bits,
                                                                   ctx->minfo);
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

    _fmpz_mod_mpoly_from_univar(A, bits, B, var, ctx);
}

void fmpz_mod_mpoly_univar_set(
    fmpz_mod_mpoly_univar_t A,
    const fmpz_mod_mpoly_univar_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;

    fmpz_mod_mpoly_univar_fit_length(A, B->length, ctx);

    for (i = 0; i < B->length; i++)
    {
        fmpz_mod_mpoly_set(A->coeffs + i, B->coeffs+ i, ctx);
        fmpz_set(A->exps + i, B->exps + i);
    }

    A->length = B->length;
}

static void _fmpz_max(fmpz_t a, const fmpz_t b, const fmpz_t c)
{
    fmpz_set(a, fmpz_cmp(b, c) > 0 ? b : c);
}


/*
    A = prem(A, -B)
    C is used for working space
*/
void _fmpz_mod_mpoly_univar_prem(
    fmpz_mod_mpoly_univar_t A,
    const fmpz_mod_mpoly_univar_t B,
    fmpz_mod_mpoly_univar_t C,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, j;
    fmpz_t z1, delta, delta_org;
    fmpz_mod_mpoly_t u, v;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(B != C);
    FLINT_ASSERT(C != A);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(fmpz_cmp(A->exps + 0, B->exps + 0) >= 0);

    fmpz_init(z1);
    fmpz_init(delta);
    fmpz_init(delta_org);
    fmpz_mod_mpoly_init(u, ctx);
    fmpz_mod_mpoly_init(v, ctx);

    fmpz_sub(delta_org, A->exps + 0, B->exps + 0);
    fmpz_add_ui(delta_org, delta_org, 1);

looper:

    if (A->length < 1)
        goto done;

    fmpz_sub(delta, A->exps + 0, B->exps + 0);
    if (fmpz_sgn(delta) < 0)
        goto done;

    i = 1;
    j = 1;
    C->length = 0;
    while (i < A->length || j < B->length)
    {
        fmpz_mod_mpoly_univar_fit_length(C, C->length + 1, ctx);

        if (j < B->length)
            fmpz_add(z1, B->exps + j, delta);

        if (i < A->length && j < B->length && fmpz_equal(A->exps + i, z1))
        {
            fmpz_mod_mpoly_mul(u, A->coeffs + i, B->coeffs + 0, ctx);
            fmpz_mod_mpoly_mul(v, A->coeffs + 0, B->coeffs + j, ctx);
            fmpz_mod_mpoly_sub(C->coeffs + C->length, v, u, ctx);
            fmpz_set(C->exps + C->length, A->exps + i);
            i++;
            j++;
        }
        else if (i < A->length && (j >= B->length ||
                                                fmpz_cmp(A->exps + i, z1) > 0))
        {
            fmpz_mod_mpoly_mul(C->coeffs + C->length, A->coeffs + i,
                                                           B->coeffs + 0, ctx);
            fmpz_mod_mpoly_neg(C->coeffs + C->length,
                                                   C->coeffs + C->length, ctx);
            fmpz_set(C->exps + C->length, A->exps + i);
            i++;
        }
        else
        {
            FLINT_ASSERT(j < B->length && (i >= A->length ||
                                               fmpz_cmp(A->exps + i, z1) < 0));
            fmpz_mod_mpoly_mul(C->coeffs + C->length, A->coeffs + 0,
                                                           B->coeffs + j, ctx);
            fmpz_set(C->exps + C->length, z1);
            j++;
        }

        C->length += !fmpz_mod_mpoly_is_zero(C->coeffs + C->length, ctx);
    }

    fmpz_mod_mpoly_univar_swap(A, C, ctx);
    fmpz_sub_ui(delta_org, delta_org, 1);
    goto looper;

done:

    FLINT_ASSERT(fmpz_sgn(delta_org) >= 0);

    if (!fmpz_is_zero(delta_org))
    {
        fmpz_mod_mpoly_neg(v, B->coeffs + 0, ctx);
        fmpz_mod_mpoly_pow_fmpz(u, v, delta_org, ctx);
        for (i = 0; i < A->length; i++)
            fmpz_mod_mpoly_mul(A->coeffs + i, A->coeffs + i, u, ctx);
    }

    fmpz_clear(delta);
    fmpz_clear(delta_org);
    fmpz_mod_mpoly_clear(u, ctx);
    fmpz_mod_mpoly_clear(v, ctx);
}

/*
    poly1 = pseudo-gcd of P and Q
          = last nonzero subresultant polynomial starting with P and Q
*/
int _fmpz_mod_mpoly_univar_pgcd_ducos(
    fmpz_mod_mpoly_univar_t poly1,
    const fmpz_mod_mpoly_univar_t polyP,
    const fmpz_mod_mpoly_univar_t polyQ,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, j, k, aJ, ae;
    fmpz_t n, d, e, J, z1, alpha;
    int iexists, jexists, kexists;
    fmpz_mod_mpoly_t u, v, w, s;
    fmpz_mod_mpoly_univar_t A, B, C, D, H, T;
    fmpz_mod_mpoly_univar_struct * last;

    FLINT_ASSERT(polyP->length > 0);
    FLINT_ASSERT(polyQ->length > 0);
    FLINT_ASSERT(fmpz_cmp(polyP->exps + 0, polyQ->exps + 0) >= 0);
    FLINT_ASSERT(fmpz_sgn(polyQ->exps + 0) >= 0);

    if (fmpz_is_zero(polyQ->exps + 0))
    {
        fmpz_mod_mpoly_univar_fit_length(poly1, 1, ctx);
        poly1->length = 1;
        fmpz_zero(poly1->exps + 0);
        return fmpz_mod_mpoly_pow_fmpz(poly1->coeffs + 0, polyQ->coeffs + 0,
                                                         polyQ->exps + 0, ctx);
    }

    fmpz_init(n);
    fmpz_init(d);
    fmpz_init(e);
    fmpz_init(J);
    fmpz_init(z1);
    fmpz_init(alpha);
    fmpz_mod_mpoly_init(u,ctx);
    fmpz_mod_mpoly_init(v,ctx);
    fmpz_mod_mpoly_init(w,ctx);
    fmpz_mod_mpoly_init(s,ctx);
    fmpz_mod_mpoly_univar_init(A,ctx);
    fmpz_mod_mpoly_univar_init(B,ctx);
    fmpz_mod_mpoly_univar_init(C,ctx);
    fmpz_mod_mpoly_univar_init(D,ctx);
    fmpz_mod_mpoly_univar_init(H,ctx);
    fmpz_mod_mpoly_univar_init(T,ctx);

    i = FLINT_MAX(polyP->length, polyQ->length);
    fmpz_mod_mpoly_univar_fit_length(A, 1 + i, ctx);
    fmpz_mod_mpoly_univar_fit_length(B, 1 + i, ctx);
    fmpz_mod_mpoly_univar_fit_length(C, 1 + i, ctx);
    fmpz_mod_mpoly_univar_fit_length(D, 1 + i, ctx);
    fmpz_mod_mpoly_univar_fit_length(H, 1 + i, ctx);
    fmpz_mod_mpoly_univar_fit_length(T, 1 + i, ctx);

    fmpz_mod_mpoly_univar_set(B, polyP, ctx);
    fmpz_mod_mpoly_univar_set(A, polyQ, ctx);

    last = A;

    fmpz_sub(z1, polyP->exps + 0, polyQ->exps + 0);
    fmpz_mod_mpoly_pow_fmpz(s, polyQ->coeffs + 0, z1, ctx);

    _fmpz_mod_mpoly_univar_prem(B, A, D, ctx);

looper:

    if (B->length < 1)
        goto done;

    fmpz_set(d, A->exps + 0);
    fmpz_set(e, B->exps + 0);

    last = B;

    fmpz_sub(z1, d, e);
    if (fmpz_is_one(z1))
    {
        if (fmpz_is_zero(e))
            goto done;

        /* D = (B[e]*A - A[e]*B)/A[d] */
        /*           i        j       */
        i = 1;
        j = 1;
        if (A->length > 1 && fmpz_equal(A->exps + 1, e))
            i++;
        else
            j = B->length;
        D->length = 0;
        while (i < A->length || j < B->length)
        {
            fmpz_mod_mpoly_univar_fit_length(D, D->length + 1, ctx);

            if (i < A->length && j < B->length &&
                                          fmpz_equal(A->exps + i, B->exps + j))
            {
                fmpz_mod_mpoly_mul(u, A->coeffs + i, B->coeffs + 0, ctx);
                fmpz_mod_mpoly_mul(v, A->coeffs + 1, B->coeffs + j, ctx);
                fmpz_mod_mpoly_sub(w, u, v, ctx);
                fmpz_mod_mpoly_divexact(D->coeffs + D->length, w,
                                                           A->coeffs + 0, ctx);
                fmpz_set(D->exps + D->length, A->exps + i);
                i++;
                j++;                
            }
            else if (i < A->length && (j >= B->length ||
                                       fmpz_cmp(A->exps + i, B->exps + j) > 0))
            {
                fmpz_mod_mpoly_mul(u, A->coeffs + i, B->coeffs + 0, ctx);
                fmpz_mod_mpoly_divexact(D->coeffs + D->length, u,
                                                           A->coeffs + 0, ctx);
                fmpz_set(D->exps + D->length, A->exps + i);
                i++;
            }
            else
            {
                FLINT_ASSERT((j < B->length && (i >= A->length ||
                                     fmpz_cmp(B->exps + j, A->exps + i) > 0)));
                fmpz_mod_mpoly_mul(v, A->coeffs + 1, B->coeffs + j, ctx);
                fmpz_mod_mpoly_divexact(D->coeffs + D->length, v,
                                                           A->coeffs + 0, ctx);
                fmpz_mod_mpoly_neg(D->coeffs + D->length,
                                                   D->coeffs + D->length, ctx);
                fmpz_set(D->exps + D->length, B->exps + j);
                j++;
            }

            D->length += !fmpz_mod_mpoly_is_zero(D->coeffs + D->length, ctx);
        }

        /* A = (B[e]*(D - B*x) + B[e-1]*B)/s */
        /*            i    j            k    */
        i = 0;
        fmpz_sub_ui(z1, e, 1);
        if (B->length > 1 && fmpz_equal(B->exps + 1, z1))
        {
            j = 2;            
            k = 1;
        }
        else
        {
            j = 1;
            k = B->length;
        }

        A->length = 0;
        while (i < D->length || j < B->length || k < B->length)
        {
            fmpz * exp;

            fmpz_mod_mpoly_univar_fit_length(A, A->length + 1, ctx);

            exp = A->exps + A->length;

            fmpz_zero(exp);

            if (i < D->length)
                _fmpz_max(exp, exp, D->exps + i);

            if (j < B->length)
            {
                fmpz_add_ui(z1, B->exps + j, 1);
                _fmpz_max(exp, exp, z1);
            }

            if (k < B->length)
                _fmpz_max(exp, exp, B->exps + k);

            iexists = (i < D->length) && fmpz_equal(exp, D->exps + i);
            jexists = (j < B->length) && fmpz_equal(exp, z1);
            kexists = (k < B->length) && fmpz_equal(exp, B->exps + k);

            FLINT_ASSERT(iexists || jexists || kexists);

            if (iexists)
            {
                if (jexists)
                {
                    fmpz_mod_mpoly_sub(w, D->coeffs + i, B->coeffs + j, ctx);
                    fmpz_mod_mpoly_mul(u, B->coeffs + 0, w, ctx);
                }
                else
                {
                    fmpz_mod_mpoly_mul(u, B->coeffs + 0, D->coeffs + i, ctx);
                }

                if (kexists)
                {
                    fmpz_mod_mpoly_mul(v, B->coeffs + 1, B->coeffs + k, ctx);
                    fmpz_mod_mpoly_add(w, u, v, ctx);
                    fmpz_mod_mpoly_divexact(A->coeffs + A->length, w, s, ctx);
                }
                else
                {
                    fmpz_mod_mpoly_divexact(A->coeffs + A->length, u, s, ctx);
                }
            }
            else
            {
                if (kexists)
                {
                    fmpz_mod_mpoly_mul(u, B->coeffs + 1, B->coeffs + k, ctx);
                    if (jexists)
                    {
                        fmpz_mod_mpoly_mul(v, B->coeffs + 0, B->coeffs + j, ctx);
                        fmpz_mod_mpoly_sub(w, u, v, ctx);
                        fmpz_mod_mpoly_divexact(A->coeffs + A->length, w, s, ctx);
                    }
                    else
                    {
                        fmpz_mod_mpoly_divexact(A->coeffs + A->length, u, s, ctx);
                    }
                }
                else
                {
                    fmpz_mod_mpoly_mul(u, B->coeffs + 0, B->coeffs + j, ctx);
                    fmpz_mod_mpoly_divexact(A->coeffs + A->length, u, s, ctx);
                    fmpz_mod_mpoly_neg(A->coeffs + A->length,
                                                   A->coeffs + A->length, ctx);
                }
            }

            A->length += !fmpz_mod_mpoly_is_zero(A->coeffs + A->length, ctx);

            i += iexists;
            j += jexists;
            k += kexists;
        }

        fmpz_mod_mpoly_univar_swap(A, B, ctx);
        fmpz_mod_mpoly_set(s, A->coeffs + 0, ctx);
        last = A;
    }
    else
    {
        fmpz_sub(n, d, e);
        fmpz_sub_ui(n, n, 1);
        fmpz_one(alpha);
        while (fmpz_add(z1, alpha, alpha), fmpz_cmp(z1, n) <= 0)
            fmpz_set(alpha, z1);

        fmpz_mod_mpoly_set(u, B->coeffs + 0, ctx);
        fmpz_sub(n, n, alpha);
        while (fmpz_cmp_ui(alpha, 1) > 0)
        {
            fmpz_tdiv_q_2exp(alpha, alpha, 1);
            fmpz_mod_mpoly_mul(v, u, u, ctx);
            fmpz_mod_mpoly_divexact(u, v, s, ctx);
            if (fmpz_cmp(n, alpha) >= 0)
            {
                fmpz_mod_mpoly_mul(v, u, B->coeffs + 0, ctx);
                fmpz_mod_mpoly_divexact(u, v, s, ctx);
                fmpz_sub(n, n, alpha);
            }
        }

        fmpz_mod_mpoly_univar_fit_length(C, B->length, ctx);
        for (i = 0; i < B->length; i++)
        {
            fmpz_mod_mpoly_mul(v, u, B->coeffs + i, ctx);
            fmpz_mod_mpoly_divexact(C->coeffs + i, v, s, ctx);
            fmpz_set(C->exps + i, B->exps + i);
        }
        C->length = B->length;

        last = C;

        if (fmpz_is_zero(e))
            goto done;

        /* H = C - C[e]*x^e */
        fmpz_mod_mpoly_univar_fit_length(H, C->length, ctx);
        for (i = 1; i < C->length; i++)
        {
            fmpz_mod_mpoly_set(H->coeffs + i - 1, C->coeffs + i, ctx);
            fmpz_set(H->exps + i - 1, C->exps + i);
        }
        H->length = C->length - 1;

        /* D = C[e]*A - A[e]*H  (truncated to powers of x < e) */
        i = 0;
        j = H->length;
        ae = A->length;
        while (i < A->length && fmpz_cmp(A->exps + i, e) >= 0)
        {
            if (fmpz_equal(A->exps + i, e))
            {
                j = 0;
                ae = i;
            }
            i++;
        }
        D->length = 0;
        while (i < A->length || j < H->length)
        {
            fmpz_mod_mpoly_univar_fit_length(D, D->length + 1, ctx);

            if (i < A->length && j < H->length &&
                                          fmpz_equal(A->exps + i, H->exps + j))
            {
                fmpz_mod_mpoly_mul(u, A->coeffs + i, C->coeffs + 0, ctx);
                fmpz_mod_mpoly_mul(v, A->coeffs + ae, H->coeffs + j, ctx);
                fmpz_mod_mpoly_sub(D->coeffs + D->length, u, v, ctx);
                fmpz_set(D->exps + D->length, A->exps + i);
                i++;
                j++;                
            }
            else if (i < A->length && (j >= H->length ||
                                       fmpz_cmp(A->exps + i, H->exps + j) > 0))
            {
                fmpz_mod_mpoly_mul(D->coeffs + D->length,
                                            A->coeffs + i, C->coeffs + 0, ctx);
                fmpz_set(D->exps + D->length, A->exps + i);
                i++;
            }
            else
            {
                FLINT_ASSERT(j < H->length && (i >= A->length ||
                                      fmpz_cmp(H->exps + j, A->exps + i) > 0));
                fmpz_mod_mpoly_mul(D->coeffs + D->length,
                                           A->coeffs + ae, H->coeffs + j, ctx);
                fmpz_mod_mpoly_neg(D->coeffs + D->length,
                                                   D->coeffs + D->length, ctx);
                fmpz_set(D->exps + D->length, H->exps + j);
                j++;
            }

            D->length += !fmpz_mod_mpoly_is_zero(D->coeffs + D->length, ctx);
        }

        for (fmpz_add_ui(J, e, 1); fmpz_cmp(J, d) < 0; fmpz_add_ui(J, J, 1))
        {
            if (H->length < 1)
                break;

            /* H = H*x - H[e-1]*B/B[e] */
            fmpz_sub_ui(z1, e, 1);
            if (fmpz_equal(H->exps + 0, z1))
            {
                i = 1;
                j = 1;
                T->length = 0;
                while (i < H->length || j < B->length)
                {
                    fmpz_mod_mpoly_univar_fit_length(T, T->length + 1, ctx);

                    if (i < H->length)
                        fmpz_add_ui(z1, H->exps + i, 1);

                    if (i < H->length && j < B->length &&
                                                   fmpz_equal(z1, B->exps + j))
                    {
                        fmpz_mod_mpoly_mul(u, H->coeffs + 0, B->coeffs + j, ctx);
                        fmpz_mod_mpoly_divexact(v, u, B->coeffs + 0, ctx);
                        fmpz_mod_mpoly_sub(T->coeffs + T->length,
                                                        H->coeffs + i, v, ctx);
                        fmpz_set(T->exps + T->length, B->exps + j);
                        i++;
                        j++;                
                    }
                    else if (i < H->length && (j >= B->length ||
                                                fmpz_cmp(z1, B->exps + j) > 0))
                    {
                        fmpz_mod_mpoly_set(T->coeffs + T->length, H->coeffs + i, ctx);
                        fmpz_set(T->exps + T->length, z1);
                        i++;
                    }
                    else
                    {
                        FLINT_ASSERT(j < B->length && (i >= H->length ||
                                               fmpz_cmp(z1, B->exps + j) < 0));
                        fmpz_mod_mpoly_mul(u, H->coeffs + 0, B->coeffs + j, ctx);
                        fmpz_mod_mpoly_divexact(T->coeffs + T->length, u,
                                                           B->coeffs + 0, ctx);
                        fmpz_mod_mpoly_neg(T->coeffs + T->length,
                                                   T->coeffs + T->length, ctx);
                        fmpz_set(T->exps + T->length, B->exps + j);
                        j++;                
                    }

                    T->length += !fmpz_mod_mpoly_is_zero(T->coeffs + T->length, ctx);
                }

                fmpz_mod_mpoly_univar_swap(H, T, ctx);
            }
            else
            {
                FLINT_ASSERT(fmpz_cmp(H->exps + 0, z1) < 0);
                for (i = 0; i < H->length; i++)
                    fmpz_add_ui(H->exps + i, H->exps + i, 1);
            }

            /* find coefficient of x^J in A */
            aJ = 0;
            while (aJ < A->length && !fmpz_equal(A->exps + aJ, J))
                aJ++;
            if (aJ >= A->length)
                continue;

            /* D = D - A[J]*H */
            i = 0;
            j = 0;
            T->length = 0;
            while (i < D->length || j < H->length)
            {
                fmpz_mod_mpoly_univar_fit_length(T, T->length + 1, ctx);

                if (i < D->length && j < H->length &&
                                          fmpz_equal(D->exps + i, H->exps + j))
                {
                    fmpz_mod_mpoly_mul(u, H->coeffs + j, A->coeffs + aJ, ctx);
                    fmpz_mod_mpoly_sub(T->coeffs + T->length, D->coeffs + i,
                                                                       u, ctx);
                    fmpz_set(T->exps + T->length, D->exps + i);
                    i++;
                    j++;                
                }
                else if (i < D->length && (j >= H->length ||
                                       fmpz_cmp(D->exps + i, H->exps + j) > 0))
                {
                    fmpz_mod_mpoly_set(T->coeffs + T->length, D->coeffs + i, ctx);
                    fmpz_set(T->exps + T->length, D->exps + i);
                    i++;
                }
                else
                {
                    FLINT_ASSERT(j < H->length && (i >= D->length ||
                                      fmpz_cmp(D->exps + i, H->exps + j) < 0));
                    fmpz_mod_mpoly_mul(T->coeffs + T->length, H->coeffs + j,
                                                          A->coeffs + aJ, ctx);
                    fmpz_mod_mpoly_neg(T->coeffs + T->length,
                                                   T->coeffs + T->length, ctx);
                    fmpz_set(T->exps + T->length, H->exps + j);
                    j++;
                }

                T->length += !fmpz_mod_mpoly_is_zero(T->coeffs + T->length, ctx);
            }
            fmpz_mod_mpoly_univar_swap(D, T, ctx);
        }

        /* B = (-1)^(d-e+1) * (B[e]*(D/A[d] - H*x) +  H[e-1]*B)/s */
        i = 0;
        fmpz_sub_ui(z1, e, 1);
        if (H->length > 0 && fmpz_equal(H->exps + 0, z1))
        {
            j = 1;
            k = 1;
        }
        else
        {
            j = 0;
            k = B->length;
        }
        T->length = 0;
        while (i < D->length || j < H->length || k < B->length)
        {
            fmpz * exp;

            fmpz_mod_mpoly_univar_fit_length(T, T->length + 1, ctx);

            exp = T->exps + T->length;

            fmpz_zero(exp);

            if (i < D->length)
                _fmpz_max(exp, exp, D->exps + i);

            if (j < H->length)
            {
                fmpz_add_ui(z1, H->exps + j, 1);
                _fmpz_max(exp, exp, z1);
            }

            if (k < B->length)
                _fmpz_max(exp, exp, B->exps + k);

            iexists = (i < D->length && fmpz_equal(exp, D->exps + i));
            jexists = (j < H->length && fmpz_equal(exp, z1));
            kexists = (k < B->length && fmpz_equal(exp, B->exps + k));

            FLINT_ASSERT(iexists || jexists || kexists);

            if (iexists)
            {
                if (jexists)
                {
                    fmpz_mod_mpoly_divexact(u, D->coeffs + i, A->coeffs + 0, ctx);
                    fmpz_mod_mpoly_sub(w, u, H->coeffs + j, ctx);
                    fmpz_mod_mpoly_mul(u, B->coeffs + 0, w, ctx);
                }
                else
                {
                    fmpz_mod_mpoly_divexact(u, D->coeffs + i, A->coeffs + 0, ctx);
                    fmpz_mod_mpoly_mul(u, B->coeffs + 0, u, ctx);
                }
                if (kexists)
                {
                    fmpz_mod_mpoly_mul(v, H->coeffs + 0, B->coeffs + k, ctx);
                    fmpz_mod_mpoly_add(w, u, v, ctx);
                    fmpz_mod_mpoly_divexact(T->coeffs + T->length, w, s, ctx);
                }
                else
                {
                    fmpz_mod_mpoly_divexact(T->coeffs + T->length, u, s, ctx);
                }
            }
            else
            {
                if (kexists)
                {
                    fmpz_mod_mpoly_mul(u, H->coeffs + 0, B->coeffs + k, ctx);
                    if (jexists)
                    {
                        fmpz_mod_mpoly_mul(v, B->coeffs + 0, H->coeffs + j, ctx);
                        fmpz_mod_mpoly_sub(w, u, v, ctx);
                        fmpz_mod_mpoly_divexact(T->coeffs + T->length, w, s, ctx);
                    }
                    else
                    {
                        fmpz_mod_mpoly_divexact(T->coeffs + T->length, u, s, ctx);
                    }
                }
                else
                {
                    fmpz_mod_mpoly_mul(u, B->coeffs + 0, H->coeffs + j, ctx);
                    fmpz_mod_mpoly_divexact(T->coeffs + T->length, u, s, ctx);
                    fmpz_mod_mpoly_neg(T->coeffs + T->length,
                                                    T->coeffs + T->length, ctx);
                }
            }

            if (((fmpz_get_ui(d) - fmpz_get_ui(e)) & 1) == 0)
                fmpz_mod_mpoly_neg(T->coeffs + T->length,
                                                   T->coeffs + T->length, ctx);            

            T->length += !fmpz_mod_mpoly_is_zero(T->coeffs + T->length, ctx);

            i += iexists;
            j += jexists;
            k += kexists;
        }

        fmpz_mod_mpoly_univar_swap(B, T, ctx);
        fmpz_mod_mpoly_univar_swap(A, C, ctx);
        fmpz_mod_mpoly_set(s, A->coeffs + 0, ctx);
        last = A;
    }

    goto looper;

done:

    fmpz_mod_mpoly_univar_swap(poly1, last, ctx);

    fmpz_clear(n);
    fmpz_clear(d);
    fmpz_clear(e);
    fmpz_clear(J);
    fmpz_clear(z1);
    fmpz_clear(alpha);
    fmpz_mod_mpoly_clear(u, ctx);
    fmpz_mod_mpoly_clear(v, ctx);
    fmpz_mod_mpoly_clear(w, ctx);
    fmpz_mod_mpoly_clear(s, ctx);
    fmpz_mod_mpoly_univar_clear(A, ctx);
    fmpz_mod_mpoly_univar_clear(B, ctx);
    fmpz_mod_mpoly_univar_clear(C, ctx);
    fmpz_mod_mpoly_univar_clear(D, ctx);
    fmpz_mod_mpoly_univar_clear(H, ctx);
    fmpz_mod_mpoly_univar_clear(T, ctx);
    return 1;
}



void fmpz_mod_mpoly_univar_derivative(
    fmpz_mod_mpoly_univar_t A,
    const fmpz_mod_mpoly_univar_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong Ai, Bi;

    fmpz_mod_mpoly_univar_fit_length(A, B->length, ctx);

    Ai = 0;
    for (Bi = 0; Bi < B->length; Bi++)
    {
        if (fmpz_sgn(B->exps + Bi) <= 0)
            continue;

        fmpz_mod_mpoly_scalar_mul_fmpz(A->coeffs + Ai, B->coeffs + Bi,
                                                            B->exps + Bi, ctx);
        fmpz_sub_ui(A->exps + Ai, B->exps + Bi, 1);

        Ai += !fmpz_mod_mpoly_is_zero(A->coeffs + Ai, ctx);
    }

    A->length = Ai;
}


int fmpz_mod_mpoly_univar_resultant(
    fmpz_mod_mpoly_t r,
    const fmpz_mod_mpoly_univar_t fx,
    const fmpz_mod_mpoly_univar_t gx,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int change_sign;
    fmpz_mod_mpoly_univar_t rx;
    const fmpz_mod_mpoly_univar_struct * F, * G;

    if (fx->length < 1 || gx->length < 1)
    {
        fmpz_mod_mpoly_zero(r, ctx);
        return 1;
    }

    fmpz_mod_mpoly_univar_init(rx, ctx);

    if (fmpz_cmp(fx->exps + 0, gx->exps + 0) < 0)
    {
        change_sign = 1 & fmpz_get_ui(fx->exps + 0) & fmpz_get_ui(gx->exps + 0);
        F = gx;
        G = fx;
    }
    else
    {
        change_sign = 0;
        F = fx;
        G = gx;
    }

    if (fmpz_is_zero(G->exps + 0))
    {
        fmpz_mod_mpoly_pow_fmpz(r, G->coeffs + 0, F->exps + 0, ctx);
    }
    else
    {
        _fmpz_mod_mpoly_univar_pgcd_ducos(rx, F, G, ctx);

        if (rx->length == 1 && fmpz_is_zero(rx->exps + 0))
            fmpz_mod_mpoly_swap(r, rx->coeffs + 0, ctx);
        else
            fmpz_mod_mpoly_zero(r, ctx);
    }

    if (change_sign)
        fmpz_mod_mpoly_neg(r, r, ctx);

    fmpz_mod_mpoly_univar_clear(rx, ctx);

    return 1;
}


int fmpz_mod_mpoly_univar_discriminant(
    fmpz_mod_mpoly_t d,
    const fmpz_mod_mpoly_univar_t fx,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_mod_mpoly_univar_t rx, fxp;

    if (fx->length < 1 || fmpz_cmp_ui(fx->exps + fx->length - 1, 1) > 0)
    {
        /* the discriminant of the zero polynomial should be zero */
        fmpz_mod_mpoly_zero(d, ctx);
        return 1;
    }

    if (fmpz_is_zero(fx->exps + 0))
    {
        /* the discriminant of the constant polynomial a should be 1/a^2 */
        fmpz_mod_mpoly_one(d, ctx);
        if (!fmpz_mod_mpoly_divides(d, d, fx->coeffs + 0, ctx))
            flint_throw(FLINT_IMPINV,
                  "fmpz_mod_mpoly_discriminant: non-unit constant polynomial");
        fmpz_mod_mpoly_mul(d, d, d, ctx);
        return 1;
    }

    if (fmpz_is_one(fx->exps + 0))
    {
        /* the discriminant of a linear polynomial should be 1 */
        fmpz_mod_mpoly_one(d, ctx);
        return 1;
    }

    /* the discriminant is (-1)^(n*(n-1)/2) res(f, f')/a_n */
    fmpz_mod_mpoly_univar_init(rx, ctx);
    fmpz_mod_mpoly_univar_init(fxp, ctx);
    fmpz_mod_mpoly_univar_derivative(fxp, fx, ctx);
    _fmpz_mod_mpoly_univar_pgcd_ducos(rx, fx, fxp, ctx);

    if (rx->length == 1 && fmpz_is_zero(rx->exps + 0))
    {
        fmpz_mod_mpoly_divexact(d, rx->coeffs + 0, fx->coeffs + 0, ctx);
        if (fmpz_get_ui(fx->exps + 0) & 2)
            fmpz_mod_mpoly_neg(d, d, ctx);
    }
    else
    {
        fmpz_mod_mpoly_zero(d, ctx);
    }

    fmpz_mod_mpoly_univar_clear(rx, ctx);
    fmpz_mod_mpoly_univar_clear(fxp, ctx);

    return 1;
}

