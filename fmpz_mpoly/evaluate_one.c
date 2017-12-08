/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"
#include "assert.h"


/*
    evaluate a f(xbar) at x_var = val,
*/
void fmpz_mpoly_evaluate_one_fmpz(fmpz_mpoly_t poly1, fmpz_mpoly_t poly2,
                                   slong var, fmpz_t val, fmpz_mpoly_ctx_t ctx)
{
    int new;
    slong i, j, N, bits;
    slong main_exp, main_var = var, main_shift, main_off;
    ulong * cmpmask, * one;
    slong p1_alloc, p1_len, p2_len;
    fmpz * p1_coeff, * p2_coeff;
    ulong * p1_exp, * p2_exp;
    slong * main_exps;
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
    TMP_INIT;

    p2_len = poly2->length;
    p2_coeff = poly2->coeffs;
    p2_exp = poly2->exps;
    bits = poly2->bits;

    if (p2_len == 0)
    {
        fmpz_mpoly_zero(poly1, ctx);
        return;
    }

    if (poly1 == poly2)
    {
        fmpz_mpoly_t temp;
        fmpz_mpoly_init(temp, ctx);
        fmpz_mpoly_evaluate_one_fmpz(temp, poly2, var, val, ctx);
        fmpz_mpoly_swap(poly1, temp, ctx);
        fmpz_mpoly_clear(temp, ctx);
        return;
    }

    TMP_START;

    N = mpoly_words_per_exp(bits, ctx->minfo);
    one = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_gen_oneexp_offset_shift(one, &main_off, &main_shift,
                                                main_var, N, bits, ctx->minfo);
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
    inds = (slong *) TMP_ALLOC(p2_len*sizeof(slong));
    mpoly_rbtree_init(tree);
    for (i = 0; i < p2_len; i++)
    {
        ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
        main_exp = (p2_exp[N*i + main_off] >> main_shift) & mask;
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

    /* manually traverse tree and add node data to heap */
    heap_len = 1;
    next_loc = tree->size + 4;
    heap = (mpoly_heap_s *) TMP_ALLOC((tree->size + 1)*sizeof(mpoly_heap_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(tree->size*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*tree->size*sizeof(slong));
    exp_array = (ulong *) TMP_ALLOC(tree->size*N*sizeof(ulong));

    i = 0;
    stack_size = 0;
    main_exps = (slong *) TMP_ALLOC(tree->size*sizeof(fmpz));
    powers = (fmpz *) TMP_ALLOC(tree->size*sizeof(fmpz));
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

    main_exps[i] = root->key;
    fmpz_init(powers + i);
    fmpz_pow_ui(powers + i, val, root->key);

    x = chain + i;
    x->i = i;
    x->j = (slong) root->data;

    x->next = NULL;
    mpoly_monomial_msub(exp_array + i*N, p2_exp + x->j*N, root->key, one, N);
    _mpoly_heap_insert(heap, exp_array + i*N, x,
                                      &next_loc, &heap_len, N, cmpmask);

    i++;    
    node = root->right;
    flint_free(root);
    root = node;

    goto looper;
done:
    assert(i == tree->size);

    /* take from heap and put into poly1 */
    fmpz_mpoly_fit_bits(poly1, bits, ctx);
    poly1->bits = bits;
    p1_coeff = poly1->coeffs;
    p1_exp = poly1->exps;
    p1_alloc = poly1->alloc;
    p1_len = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        _fmpz_mpoly_fit_length(&p1_coeff, &p1_exp, &p1_alloc, p1_len + 1, N);

        mpoly_monomial_set(p1_exp + p1_len*N, exp, N);
        fmpz_zero(p1_coeff + p1_len);
        do
        {
            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
            do
            {
                *store++ = x->i;
                *store++ = x->j;
                fmpz_addmul(p1_coeff + p1_len, powers + x->i, p2_coeff + x->j);
            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

        p1_len += !fmpz_is_zero(p1_coeff + p1_len);

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
                mpoly_monomial_msub(exp_array + i*N, p2_exp + x->j*N,
                                                         main_exps[i], one, N);
                _mpoly_heap_insert(heap, exp_array + i*N, x,
                                      &next_loc, &heap_len, N, cmpmask);
            }
        }
    }

    poly1->coeffs = p1_coeff;
    poly1->exps = p1_exp;
    poly1->alloc = p1_alloc;
    _fmpz_mpoly_set_length(poly1, p1_len, ctx);

    for (i = 0; i < tree->size; i++)
        fmpz_clear(powers + i);

    TMP_END;
}
