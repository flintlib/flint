/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

/* exponents of B are not multiprecision */
static void _fq_nmod_mpoly_evaluate_one_fq_nmod_sp(
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    slong var,
    const fq_nmod_t val,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int newer;
    fq_nmod_t pp;
    slong i, j, N;
    flint_bitcnt_t bits;
    slong main_exp, main_shift, main_off;
    ulong * cmpmask, * one;
    slong Aalloc, Alen, Blen;
    fq_nmod_struct * Acoeff, * Bcoeff;
    ulong * Aexp, * Bexp;
    ulong * main_exps;
    fq_nmod_struct * powers;
    slong * inds;
    ulong * exp_array, * exp;
    slong next_loc;
    slong heap_len;
    mpoly_heap_s * heap;
    mpoly_heap_t * x;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_rbnode_struct ** stack;
    slong stack_size;
    mpoly_rbtree_t tree;
    mpoly_rbnode_struct * node;
    mpoly_rbnode_struct * root;
    TMP_INIT;

    Blen = B->length;
    Bcoeff = B->coeffs;
    Bexp = B->exps;
    bits = B->bits;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(Blen > 0);

    TMP_START;

    fq_nmod_init(pp, ctx->fqctx);

    N = mpoly_words_per_exp(bits, ctx->minfo);
    one = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_gen_monomial_offset_shift_sp(one, &main_off, &main_shift,
                                                        var, bits, ctx->minfo);
    mpoly_get_cmpmask(cmpmask, N, bits, ctx->minfo);

    /* scan poly2 and put powers of var into tree */
    /*
        suppose x > y with lex order and poly2 is
        x^3*y + x^2*y^3 + x*y + y + 1
        and we are evaluating at y = 2

        The data member of nodes in the tree of powers of y appearing in
        the polynomial contain the term number of the first such appearance
             data
        y^0  4    (1       term)
        y^1  0    (x^3*y   term)
        y^3  1    (x^2*y^3 term)

        the inds array conains the term number of the next term with
        the same power of y

        inds[0] = 2  (x*y term)
        inds[1] = -1 
        inds[2] = 3  (y   term)
        inds[3] = -1
        inds[4] = -1
    */
    inds = (slong *) TMP_ALLOC(Blen*sizeof(slong));
    mpoly_rbtree_init(tree);
    for (i = 0; i < Blen; i++)
    {
        ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
        main_exp = (Bexp[N*i + main_off] >> main_shift) & mask;
        node = mpoly_rbtree_get(&newer, tree, main_exp);
        if (newer)
        {
            node->data = (void *) i;
        }
        else
        {
            inds[(slong) node->data2] = i;
        }
        node->data2 = (void *) i;
        inds[i] = -WORD(1);
    }

    /* manually traverse tree and add node data to heap */
    heap_len = 1;
    next_loc = tree->size + 4;
    heap = (mpoly_heap_s *) TMP_ALLOC((tree->size + 1)*sizeof(mpoly_heap_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(tree->size*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*tree->size*sizeof(slong));
    exp_array = (ulong *) TMP_ALLOC(tree->size*N*sizeof(ulong));

    i = 0;
    stack_size = 0;
    main_exps = (ulong *) TMP_ALLOC(tree->size*N*sizeof(ulong));
    powers = (fq_nmod_struct *) TMP_ALLOC(tree->size*sizeof(fq_nmod_struct));
    stack = (mpoly_rbnode_struct **) TMP_ALLOC(tree->size
                                              * sizeof(mpoly_rbnode_struct *));
    root = tree->head->left;

looper:

    while (root != tree->null)
    {
        stack[stack_size++] = root;
        root = root->left;
    }

    if (stack_size == 0)
        goto done;

    root = stack[--stack_size];

    mpoly_monomial_mul_ui(main_exps + N*i, one, N, root->key);

    fq_nmod_init(powers + i, ctx->fqctx);
    fq_nmod_pow_ui(powers + i, val, root->key, ctx->fqctx);

    x = chain + i;
    x->i = i;
    x->j = (slong) root->data;

    x->next = NULL;
    mpoly_monomial_sub(exp_array + N*i, Bexp + N*x->j, main_exps + N*i, N);
    _mpoly_heap_insert(heap, exp_array + N*i, x,
                                      &next_loc, &heap_len, N, cmpmask);

    i++;    
    node = root->right;
    flint_free(root);
    root = node;

    goto looper;

done:

    FLINT_ASSERT(i == tree->size);

    /* take from heap and put into A */
    fq_nmod_mpoly_fit_bits(A, bits, ctx);
    A->bits = bits;
    Acoeff = A->coeffs;
    Aexp = A->exps;
    Aalloc = A->alloc;
    Alen = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        _fq_nmod_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, N, ctx->fqctx);

        mpoly_monomial_set(Aexp + Alen*N, exp, N);

        fq_nmod_zero(Acoeff + Alen, ctx->fqctx);
        do {
            x = (mpoly_heap_t *) _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
            do {
                *store++ = x->i;
                *store++ = x->j;
                fq_nmod_mul(pp, powers + x->i, Bcoeff + x->j, ctx->fqctx);
                fq_nmod_add(Acoeff + Alen, Acoeff + Alen, pp, ctx->fqctx);
            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

        Alen += !fq_nmod_is_zero(Acoeff + Alen, ctx->fqctx);

        while (store > store_base)
        {
            j = *--store;
            i = *--store;
            if (inds[j] >= 0)
            {
                x = chain + i;
                x->i = i;
                x->j = inds[j];
                x->next = NULL;
                mpoly_monomial_sub(exp_array + N*i, Bexp + N*x->j,
                                                       main_exps + N*i, N);
                _mpoly_heap_insert(heap, exp_array + i*N, x,
                                      &next_loc, &heap_len, N, cmpmask);
            }
        }
    }

    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->alloc = Aalloc;
    A->length = Alen;

    for (j = 0; j < tree->size; j++)
        fq_nmod_clear(powers + j, ctx->fqctx);

    fq_nmod_clear(pp, ctx->fqctx);

    TMP_END;
}


/* exponents of B are multiprecision */
static void _fq_nmod_mpoly_evaluate_one_fq_nmod_mp(
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    slong var,
    const fq_nmod_t val,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int newer;
    fq_nmod_t pp;
    slong i, j, N, bits;
    slong main_off;
    fmpz_t main_exp;
    ulong * cmpmask, * main_one;
    slong Aalloc, Alen, Blen;
    fq_nmod_struct * Acoeff, * Bcoeff;
    ulong * Aexp, * Bexp;
    ulong * main_exps;
    fq_nmod_struct * powers;
    slong * inds;
    ulong * exp_array, * exp;
    slong next_loc;
    slong heap_len;
    mpoly_heap_s * heap;
    mpoly_heap_t * x;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_rbnode_struct ** stack;
    slong stack_size;
    mpoly_rbtree_t tree;
    mpoly_rbnode_struct * node;
    mpoly_rbnode_struct * root;
    TMP_INIT;

    fq_nmod_init(pp, ctx->fqctx);

    Blen = B->length;
    Bcoeff = B->coeffs;
    Bexp = B->exps;
    bits = B->bits;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(bits > FLINT_BITS);
    FLINT_ASSERT(Blen > 0);

    TMP_START;

    N = mpoly_words_per_exp(bits, ctx->minfo);
    main_one = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    main_off = mpoly_gen_monomial_offset_mp(main_one, var, bits, ctx->minfo);
    mpoly_get_cmpmask(cmpmask, N, bits, ctx->minfo);

    /* scan poly2 and put powers of var into tree */
    /*
        suppose x > y with lex order and poly2 is
        x^3*y + x^2*y^3 + x*y + y + 1
        and we are evaluating at y = 2

        The data member of nodes in the tree of powers of y appearing in
        the polynomial contain the term number of the first such appearance
             data
        y^0  4    (1       term)
        y^1  0    (x^3*y   term)
        y^3  1    (x^2*y^3 term)

        the inds array conains the term number of the next term with
        the same power of y

        inds[0] = 2  (x*y term)
        inds[1] = -1 
        inds[2] = 3  (y   term)
        inds[3] = -1
        inds[4] = -1
    */
    inds = (slong *) TMP_ALLOC(Blen*sizeof(slong));
    mpoly_rbtree_init(tree);
    fmpz_init(main_exp);
    for (i = 0; i < Blen; i++)
    {
        fmpz_set_ui_array(main_exp, Bexp + N*i + main_off, bits/FLINT_BITS);
        node = mpoly_rbtree_get_fmpz(&newer, tree, main_exp);
        if (newer)
        {
            node->data = (void *) i;
        }
        else
        {
            inds[(slong) node->data2] = i;
        }
        node->data2 = (void *) i;
        inds[i] = -WORD(1);
    }
    fmpz_clear(main_exp);

    /* manually traverse tree and add node data to heap */
    heap_len = 1;
    next_loc = tree->size + 4;
    heap = (mpoly_heap_s *) TMP_ALLOC((tree->size + 1)*sizeof(mpoly_heap_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(tree->size*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*tree->size*sizeof(slong));
    exp_array = (ulong *) TMP_ALLOC(tree->size*N*sizeof(ulong));

    i = 0;
    stack_size = 0;
    main_exps = (ulong *) TMP_ALLOC(tree->size*N*sizeof(ulong));
    powers = (fq_nmod_struct *) TMP_ALLOC(tree->size*sizeof(fq_nmod_struct));
    stack = (mpoly_rbnode_struct **) TMP_ALLOC(tree->size
                                              * sizeof(mpoly_rbnode_struct *));
    root = tree->head->left;

looper:

    while (root != tree->null)
    {
        stack[stack_size++] = root;
        root = root->left;
    }
    if (stack_size == 0)
        goto done;

    root = stack[--stack_size];

    mpoly_monomial_mul_fmpz(main_exps + N*i, main_one, N, &root->key);
    fq_nmod_init(powers + i, ctx->fqctx);
    fq_nmod_pow(powers + i, val, (fmpz*)(&root->key), ctx->fqctx);

    x = chain + i;
    x->i = i;
    x->j = (slong) root->data;

    x->next = NULL;
    mpoly_monomial_sub_mp(exp_array + N*i, Bexp + N*x->j, main_exps + N*i, N);
    _mpoly_heap_insert(heap, exp_array + i*N, x,
                                      &next_loc, &heap_len, N, cmpmask);

    i++;    
    node = root->right;
    fmpz_clear((fmpz*)(&root->key));
    flint_free(root);
    root = node;

    goto looper;

done:

    FLINT_ASSERT(i == tree->size);

    /* take from heap and put into A */
    fq_nmod_mpoly_fit_bits(A, bits, ctx);
    A->bits = bits;
    Acoeff = A->coeffs;
    Aexp = A->exps;
    Aalloc = A->alloc;
    Alen = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        _fq_nmod_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, N, ctx->fqctx);

        mpoly_monomial_set(Aexp + Alen*N, exp, N);

        fq_nmod_zero(Acoeff + Alen, ctx->fqctx);
        do {
            x = (mpoly_heap_t *) _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
            do {
                *store++ = x->i;
                *store++ = x->j;
                fq_nmod_mul(pp, powers + x->i, Bcoeff + x->j, ctx->fqctx);
                fq_nmod_add(Acoeff + Alen, Acoeff + Alen, pp, ctx->fqctx);
            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

        Alen += !fq_nmod_is_zero(Acoeff + Alen, ctx->fqctx);

        while (store > store_base)
        {
            j = *--store;
            i = *--store;
            if (inds[j] >= 0)
            {
                x = chain + i;
                x->i = i;
                x->j = inds[j];
                x->next = NULL;
                mpoly_monomial_sub_mp(exp_array + N*i, Bexp + N*x->j,
                                                           main_exps + N*i, N);
                _mpoly_heap_insert(heap, exp_array + N*i, x,
                                      &next_loc, &heap_len, N, cmpmask);
            }
        }
    }

    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->alloc = Aalloc;
    A->length = Alen;

    for (j = 0; j < tree->size; j++)
        fq_nmod_clear(powers + j, ctx->fqctx);

    fq_nmod_clear(pp, ctx->fqctx);

    TMP_END;
}

void fq_nmod_mpoly_evaluate_one_fq_nmod(
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    slong var,
    const fq_nmod_t val,
    const fq_nmod_mpoly_ctx_t ctx)
{
    if (B->length == 0)
    {
        fq_nmod_mpoly_zero(A, ctx);
        return;
    }

    if (A == B)
    {
        fq_nmod_mpoly_t T;
        fq_nmod_mpoly_init(T, ctx);
        fq_nmod_mpoly_evaluate_one_fq_nmod(T, B, var, val, ctx);
        fq_nmod_mpoly_swap(A, T, ctx);
        fq_nmod_mpoly_clear(T, ctx);
        return;
    }

    if (B->bits <= FLINT_BITS)
    {
        _fq_nmod_mpoly_evaluate_one_fq_nmod_sp(A, B, var, val, ctx);
    }
    else
    {
        _fq_nmod_mpoly_evaluate_one_fq_nmod_mp(A, B, var, val, ctx);
    }

    fq_nmod_mpoly_assert_canonical(A, ctx);
}

