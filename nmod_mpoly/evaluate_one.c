/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

/* exponents of B are not multiprecision */
void _nmod_mpoly_evaluate_one_ui_sp(nmod_mpoly_t A, const nmod_mpoly_t B,
                              slong var, ulong val, const nmod_mpoly_ctx_t ctx)
{
    int new;
    ulong acc0, acc1, acc2, pp0, pp1;
    slong i, j, N;
    mp_bitcnt_t bits;
    slong main_exp, main_shift, main_off;
    ulong * cmpmask, * one;
    slong Aalloc, Alen, Blen;
    mp_limb_t * Acoeff, * Bcoeff;
    ulong * Aexp, * Bexp;
    ulong * main_exps;
    mp_limb_t * powers;
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

    N = mpoly_words_per_exp(bits, ctx->minfo);
    one = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_gen_oneexp_offset_shift(one, &main_off, &main_shift,
                                                     var, N, bits, ctx->minfo);
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
        node = mpoly_rbtree_get(&new, tree, main_exp);
        if (new)
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
    powers = (mp_limb_t *) TMP_ALLOC(tree->size*sizeof(mp_limb_t));
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
    powers[i] = nmod_pow_ui(val, root->key, ctx->ffinfo->mod);

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
    nmod_mpoly_fit_bits(A, bits, ctx);
    A->bits = bits;
    Acoeff = A->coeffs;
    Aexp = A->exps;
    Aalloc = A->alloc;
    Alen = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        _nmod_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, N);

        mpoly_monomial_set(Aexp + Alen*N, exp, N);
        
        acc0 = acc1 = acc2 = 0;
        do {
            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
            do {
                *store++ = x->i;
                *store++ = x->j;
                umul_ppmm(pp1, pp0, powers[x->i], Bcoeff[x->j]);
                add_sssaaaaaa(acc2, acc1, acc0,
                              acc2, acc1, acc0, WORD(0), pp1, pp0);
            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

        NMOD_RED3(Acoeff[Alen], acc2, acc1, acc0, ctx->ffinfo->mod);
        Alen += (Acoeff[Alen] != 0);

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

    TMP_END;
}


/* exponents of B are multiprecision */
void _nmod_mpoly_evaluate_one_ui_mp(nmod_mpoly_t A, const nmod_mpoly_t B,
                             slong var, ulong val, const nmod_mpoly_ctx_t ctx)
{
    int new;
    ulong acc0, acc1, acc2, pp0, pp1;
    slong i, j, N, bits;
    slong main_off;
    fmpz_t main_exp;
    ulong * cmpmask, * main_one;
    slong Aalloc, Alen, Blen;
    mp_limb_t * Acoeff, * Bcoeff;
    ulong * Aexp, * Bexp;
    ulong * main_exps;
    mp_limb_t * powers;
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
    FLINT_ASSERT(bits > FLINT_BITS);
    FLINT_ASSERT(Blen > 0);

    TMP_START;

    N = mpoly_words_per_exp(bits, ctx->minfo);
    main_one = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_gen_oneexp_offset_mp(main_one, &main_off, var, N, bits, ctx->minfo);
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
        node = mpoly_rbtree_get_fmpz(&new, tree, main_exp);
        if (new)
        {
            node->data = (void *) i;
        } else
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
    powers = (mp_limb_t *) TMP_ALLOC(tree->size*sizeof(mp_limb_t));
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
    powers[i] = nmod_pow_fmpz(val, (fmpz*)(&root->key), ctx->ffinfo->mod);

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
    nmod_mpoly_fit_bits(A, bits, ctx);
    A->bits = bits;
    Acoeff = A->coeffs;
    Aexp = A->exps;
    Aalloc = A->alloc;
    Alen = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        _nmod_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, N);

        mpoly_monomial_set(Aexp + Alen*N, exp, N);
        
        acc0 = acc1 = acc2 = 0;
        do {
            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
            do {
                *store++ = x->i;
                *store++ = x->j;
                umul_ppmm(pp1, pp0, powers[x->i], Bcoeff[x->j]);
                add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);
            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

        NMOD_RED3(Acoeff[Alen], acc2, acc1, acc0, ctx->ffinfo->mod);
        Alen += (Acoeff[Alen] != 0);

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

    TMP_END;
}


void nmod_mpoly_evaluate_one_ui(nmod_mpoly_t A, const nmod_mpoly_t B,
                              slong var, ulong val, const nmod_mpoly_ctx_t ctx)
{
    if (B->length == 0)
    {
        nmod_mpoly_zero(A, ctx);
        return;
    }

    if (A == B)
    {
        nmod_mpoly_t T;
        nmod_mpoly_init(T, ctx);
        nmod_mpoly_evaluate_one_ui(T, B, var, val, ctx);
        nmod_mpoly_swap(A, T, ctx);
        nmod_mpoly_clear(T, ctx);
        return;
    }

    NMOD_RED(val, val, ctx->ffinfo->mod);

    if (B->bits <= FLINT_BITS)
    {
        _nmod_mpoly_evaluate_one_ui_sp(A, B, var, val, ctx);
    } else
    {
        _nmod_mpoly_evaluate_one_ui_mp(A, B, var, val, ctx);
    }
}

