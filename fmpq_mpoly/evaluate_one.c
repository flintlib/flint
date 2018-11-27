/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

/* exponents of B are not multiprecision */
void _fmpq_mpoly_evaluate_one_fmpq_sp(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                       slong var, const fmpq_t val, const fmpq_mpoly_ctx_t ctx)
{
    int new;
    slong i, j, N;
    mp_bitcnt_t bits;
    slong main_shift, main_off;
    ulong main_exp, emin, emax;
    ulong * cmpmask, * one;
    slong Aalloc, Alen, Blen;
    fmpz * Acoeff, * Bcoeff;
    ulong * Aexp, * Bexp;
    ulong * main_exps;
    slong poweralloc;
    fmpz * powers;
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
    fmpq_t u;
    TMP_INIT;

    Blen = B->zpoly->length;
    Bcoeff = B->zpoly->coeffs;
    Bexp = B->zpoly->exps;
    bits = B->zpoly->bits;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(Blen > 0);

    TMP_START;

    fmpq_init(u);

    N = mpoly_words_per_exp(bits, ctx->zctx->minfo);
    one = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_gen_oneexp_offset_shift(one, &main_off, &main_shift,
                                               var, N, bits, ctx->zctx->minfo);
    mpoly_get_cmpmask(cmpmask, N, bits, ctx->zctx->minfo);

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
    emin = -UWORD(1);
    emax = UWORD(0);
    
    for (i = 0; i < Blen; i++)
    {
        ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
        main_exp = (Bexp[N*i + main_off] >> main_shift) & mask;
        emin = FLINT_MIN(emin, main_exp);
        emax = FLINT_MAX(emax, main_exp);
        node = mpoly_rbtree_get(&new, tree, main_exp);
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
    FLINT_ASSERT(emin <= emax);

    fmpz_pow_ui(fmpq_numref(u), fmpq_numref(val), emin);
    fmpz_pow_ui(fmpq_denref(u), fmpq_denref(val), emax);
    fmpq_mul(A->content, B->content, u);

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
    poweralloc = tree->size;
    powers = _fmpz_vec_init(poweralloc);
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

    FLINT_ASSERT(root->key >= emin);
    FLINT_ASSERT(root->key <= emax);
    mpoly_monomial_mul_ui(main_exps + N*i, one, N, root->key);
    
    fmpz_pow_ui(fmpq_numref(u), fmpq_numref(val), root->key - emin);
    fmpz_pow_ui(fmpq_denref(u), fmpq_denref(val), emax - root->key);
    fmpz_mul(powers + i, fmpq_numref(u), fmpq_denref(u));

    x = chain + i;
    x->i = i;
    x->j = (slong) root->data;

    x->next = NULL;
    mpoly_monomial_sub(exp_array + N*i, Bexp + N*x->j, main_exps + N*i, N);
    _mpoly_heap_insert(heap, exp_array + i*N, x,
                                      &next_loc, &heap_len, N, cmpmask);

    i++;    
    node = root->right;
    flint_free(root);
    root = node;

    goto looper;
done:
    FLINT_ASSERT(i == tree->size);

    /* take from heap and put into A */
    fmpz_mpoly_fit_bits(A->zpoly, bits, ctx->zctx);
    A->zpoly->bits = bits;
    Acoeff = A->zpoly->coeffs;
    Aexp = A->zpoly->exps;
    Aalloc = A->zpoly->alloc;
    Alen = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, N);

        mpoly_monomial_set(Aexp + Alen*N, exp, N);
        
        fmpz_zero(Acoeff + Alen);
        do {
            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
            do {
                *store++ = x->i;
                *store++ = x->j;
                fmpz_addmul(Acoeff + Alen, powers + x->i, Bcoeff + x->j);
            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

        Alen += !fmpz_is_zero(Acoeff + Alen);

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
                mpoly_monomial_sub(exp_array + N*i, Bexp + x->j*N,
                                                       main_exps + N*i, N);
                _mpoly_heap_insert(heap, exp_array + i*N, x,
                                      &next_loc, &heap_len, N, cmpmask);
            }
        }
    }

    A->zpoly->coeffs = Acoeff;
    A->zpoly->exps = Aexp;
    A->zpoly->alloc = Aalloc;
    _fmpz_mpoly_set_length(A->zpoly, Alen, ctx->zctx);

    fmpq_mpoly_reduce(A, ctx);

    _fmpz_vec_clear(powers, poweralloc);

    fmpq_clear(u);

    TMP_END;
}

/* exponents of B are multiprecision */
void _fmpq_mpoly_evaluate_one_fmpq_mp(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                       slong var, const fmpq_t val, const fmpq_mpoly_ctx_t ctx)
{
    int new;
    slong i, j, N, bits;
    slong main_off;
    fmpz_t main_exp, emin, emax;
    ulong * cmpmask, * main_one;
    slong Aalloc, Alen, Blen;
    fmpz * Acoeff, * Bcoeff;
    ulong * Aexp, * Bexp;
    ulong * main_exps;
    slong poweralloc;
    fmpz * powers;
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
    fmpz_t t;
    fmpq_t u;
    TMP_INIT;

    Blen = B->zpoly->length;
    Bcoeff = B->zpoly->coeffs;
    Bexp = B->zpoly->exps;
    bits = B->zpoly->bits;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(bits > FLINT_BITS);
    FLINT_ASSERT(Blen > 0);

    TMP_START;

    fmpz_init(t);
    fmpq_init(u);
    fmpz_init(emin);
    fmpz_init(emax);

    N = mpoly_words_per_exp(bits, ctx->zctx->minfo);
    main_one = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_gen_oneexp_offset_mp(main_one, &main_off, var, N, bits, ctx->zctx->minfo);
    mpoly_get_cmpmask(cmpmask, N, bits, ctx->zctx->minfo);

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
        if (i == 0)
        {
            fmpz_set(emin, main_exp);
            fmpz_set(emax, main_exp);
        }
        else
        {
            if (fmpz_cmp(emin, main_exp) > 0)
                fmpz_set(emin, main_exp);
            if (fmpz_cmp(emax, main_exp) < 0)
                fmpz_set(emax, main_exp);
        }
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
    FLINT_ASSERT(fmpz_cmp(emax, emin) >= 0);

    fmpz_pow_fmpz(fmpq_numref(u), fmpq_numref(val), emin);
    fmpz_pow_fmpz(fmpq_denref(u), fmpq_denref(val), emax);
    fmpq_mul(A->content, B->content, u);

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
    poweralloc = tree->size;
    powers = _fmpz_vec_init(poweralloc);
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

    fmpz_sub(t, &root->key, emin);
    fmpz_pow_fmpz(fmpq_numref(u), fmpq_numref(val), t);
    fmpz_sub(t, emax, &root->key);
    fmpz_pow_fmpz(fmpq_denref(u), fmpq_denref(val), t);
    fmpz_mul(powers + i, fmpq_numref(u), fmpq_denref(u));


    x = chain + i;
    x->i = i;
    x->j = (slong) root->data;

    x->next = NULL;
    mpoly_monomial_sub_mp(exp_array + N*i, Bexp + N*x->j, main_exps + N*i, N);
    _mpoly_heap_insert(heap, exp_array + i*N, x,
                                      &next_loc, &heap_len, N, cmpmask);

    i++;    
    node = root->right;
    fmpz_clear(&root->key);
    flint_free(root);
    root = node;

    goto looper;
done:
    FLINT_ASSERT(i == tree->size);

    /* take from heap and put into A */
    fmpz_mpoly_fit_bits(A->zpoly, bits, ctx->zctx);
    A->zpoly->bits = bits;
    Acoeff = A->zpoly->coeffs;
    Aexp = A->zpoly->exps;
    Aalloc = A->zpoly->alloc;
    Alen = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, N);

        mpoly_monomial_set(Aexp + Alen*N, exp, N);
        
        fmpz_zero(Acoeff + Alen);
        do {
            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
            do {
                *store++ = x->i;
                *store++ = x->j;
                fmpz_addmul(Acoeff + Alen, powers + x->i, Bcoeff + x->j);
            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

        Alen += !fmpz_is_zero(Acoeff + Alen);

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

    A->zpoly->coeffs = Acoeff;
    A->zpoly->exps = Aexp;
    A->zpoly->alloc = Aalloc;
    _fmpz_mpoly_set_length(A->zpoly, Alen, ctx->zctx);

    _fmpz_vec_clear(powers, poweralloc);
    fmpq_clear(u);
    fmpz_clear(t);
    fmpz_clear(emin);
    fmpz_clear(emax);

    fmpq_mpoly_reduce(A, ctx);

    TMP_END;
}

void fmpq_mpoly_evaluate_one_fmpq(fmpq_mpoly_t A,
                           const fmpq_mpoly_t B, slong var, const fmpq_t val,
                                                   const fmpq_mpoly_ctx_t ctx)
{
    if (B->zpoly->length == 0)
    {
        fmpq_mpoly_zero(A, ctx);
        return;
    }

    if (A == B)
    {
        fmpq_mpoly_t T;
        fmpq_mpoly_init(T, ctx);
        fmpq_mpoly_evaluate_one_fmpq(T, B, var, val, ctx);
        fmpq_mpoly_swap(A, T, ctx);
        fmpq_mpoly_clear(T, ctx);
        return;
    }

    if (B->zpoly->bits <= FLINT_BITS)
    {
        _fmpq_mpoly_evaluate_one_fmpq_sp(A, B, var, val, ctx);
    } else
    {
        _fmpq_mpoly_evaluate_one_fmpq_mp(A, B, var, val, ctx);
    }
}

