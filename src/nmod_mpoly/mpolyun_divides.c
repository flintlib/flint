/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

int _nmod_mpolyn_divides(
    nmod_mpolyn_t Q,
    const nmod_mpolyn_t A,
    const nmod_mpolyn_t B,
    slong N,
    const ulong * cmpmask,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    int lt_divides;
    flint_bitcnt_t bits = Q->bits;
    slong i, j, Qlen, s;
    slong next_loc;
    slong heap_len = 2; /* heap zero index unused */
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    ulong mask;
    slong * hind;
    n_poly_t r, acc;
    TMP_INIT;

    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(bits == Q->bits);
    FLINT_ASSERT(bits == A->bits);
    FLINT_ASSERT(bits == B->bits);

    TMP_START;

    n_poly_init(r);
    n_poly_init(acc);

    /* alloc array of heap nodes which can be chained together */
    next_loc = B->length + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap_s *) TMP_ALLOC((B->length + 1)*sizeof(mpoly_heap_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(B->length*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*B->length*sizeof(slong));

    /* array of exponent vectors, each of "N" words */
    exps = (ulong *) TMP_ALLOC(B->length*N*sizeof(ulong));
    /* list of pointers to available exponent vectors */
    exp_list = (ulong **) TMP_ALLOC(B->length*sizeof(ulong *));
    /* set up list of available exponent vectors */
    exp_next = 0;
    for (i = 0; i < B->length; i++)
        exp_list[i] = exps + i*N;

    /* space for flagged heap indices */
    hind = (slong *) TMP_ALLOC(B->length*sizeof(slong));
    for (i = 0; i < B->length; i++)
        hind[i] = 1;

    mask = mpoly_overflow_mask_sp(bits);

    /* s is the number of terms * (latest quotient) we should put into heap */
    s = B->length;

    /* insert (-1, 0, exp2[0]) into heap */
    x = chain + 0;
    x->i = -WORD(1);
    x->j = 0;
    x->next = NULL;
    heap[1].next = x;
    heap[1].exp = exp_list[exp_next++];
    mpoly_monomial_set(heap[1].exp, A->exps + N*0, N);

    Qlen = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        if (mpoly_monomial_overflows(exp, N, mask))
            goto not_exact_division;

        nmod_mpolyn_fit_length(Q, Qlen + 1, ctx);

        lt_divides = mpoly_monomial_divides(Q->exps + Qlen*N,
                                                  exp, B->exps + N*0, N, mask);

        n_poly_zero(acc);
        do
        {
            exp_list[--exp_next] = heap[1].exp;
            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
            do
            {
                *store++ = x->i;
                *store++ = x->j;
                if (x->i != -WORD(1))
                    hind[x->i] |= WORD(1);

                if (x->i == -WORD(1))
                {
                    n_poly_mod_add(acc, acc, A->coeffs + x->j, ctx->mod);
                }
                else
                {
                    n_poly_mod_mul(r, B->coeffs + x->i, Q->coeffs + x->j, ctx->mod);
                    n_poly_mod_sub(acc, acc, r, ctx->mod);
                }
            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

        /* process nodes taken from the heap */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            if (i == -WORD(1))
            {
                /* take next dividend term */
                if (j + 1 < A->length)
                {
                    x = chain + 0;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    mpoly_monomial_set(exp_list[exp_next], A->exps + x->j*N, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
            }
            else
            {
                /* should we go right? */
                if (  (i + 1 < B->length)
                   && (hind[i + 1] == 2*j + 1)
                   )
                {
                    x = chain + i + 1;
                    x->i = i + 1;
                    x->j = j;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;

                    mpoly_monomial_add(exp_list[exp_next], B->exps + x->i*N,
                                                           Q->exps + x->j*N, N);

                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
                /* should we go up? */
                if (j + 1 == Qlen)
                {
                    s++;
                } else if (  ((hind[i] & 1) == 1)
                          && ((i == 1) || (hind[i - 1] >= 2*(j + 2) + 1))
                          )
                {
                    x = chain + i;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;

                    mpoly_monomial_add(exp_list[exp_next], B->exps + x->i*N,
                                                           Q->exps + x->j*N, N);

                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
            }
        }

        if (n_poly_is_zero(acc))
            continue;

        n_poly_mod_divrem(Q->coeffs + Qlen, r, acc, B->coeffs + 0, ctx->mod);
        if (!n_poly_is_zero(r))
            goto not_exact_division;

        if (!lt_divides)
            goto not_exact_division;

        if (s > 1)
        {
            i = 1;
            x = chain + i;
            x->i = i;
            x->j = Qlen;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;

            mpoly_monomial_add(exp_list[exp_next], B->exps + x->i*N,
                                                   Q->exps + x->j*N, N);

            exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
        }
        s = 1;
        Qlen++;
    }

    success = 1;

cleanup:

    n_poly_clear(r);
    n_poly_clear(acc);

    Q->length = Qlen;

    TMP_END;

    return success;

not_exact_division:

    success = 0;
    goto cleanup;
}

int nmod_mpolyn_divides(
    nmod_mpolyn_t Q,
    const nmod_mpolyn_t A,
    const nmod_mpolyn_t B,
    const nmod_mpoly_ctx_t ctx)
{
    slong N;
    ulong * cmpmask;
    flint_bitcnt_t bits = Q->bits;
    int success;

    TMP_INIT;
    TMP_START;

    N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, bits, ctx->minfo);

    success = _nmod_mpolyn_divides(Q, A, B, N, cmpmask, ctx);

    TMP_END;

    return success;
}


/* The following functions are currently untested and unused. */

void _nmod_mpolyn_add(
    nmod_mpolyn_t A,
    const nmod_mpolyn_t B,
    const nmod_mpolyn_t C,
    slong N,
    const nmod_mpoly_ctx_t ctx)
{
    slong i = 0, j = 0, Alen = 0;

    FLINT_ASSERT(N == mpoly_words_per_exp(A->bits, ctx->minfo));
    FLINT_ASSERT(N == mpoly_words_per_exp(B->bits, ctx->minfo));
    FLINT_ASSERT(N == mpoly_words_per_exp(C->bits, ctx->minfo));

    nmod_mpolyn_fit_length(A, B->length + C->length, ctx);

    while (i < B->length && j < C->length)
    {
        int cmp = mpoly_monomial_cmp_nomask(B->exps + i*N, C->exps + j*N, N);

        if (cmp > 0)
        {
            n_poly_set(A->coeffs + Alen, B->coeffs + i);
            mpoly_monomial_set(A->exps + Alen*N, B->exps + i*N, N);
            i++;
            Alen++;
        }
        else if (cmp == 0)
        {
            n_poly_mod_add(A->coeffs + Alen, B->coeffs + i, C->coeffs + j, ctx->mod);
            mpoly_monomial_set(A->exps + Alen*N, B->exps + i*N, N);
            i++;
            j++;
            Alen += !n_poly_is_zero(A->coeffs + Alen);
        }
        else
        {
            n_poly_set(A->coeffs + Alen, C->coeffs + j);
            mpoly_monomial_set(A->exps + Alen*N, C->exps + j*N, N);
            j++;
            Alen++;
        }
    }

    while (i < B->length)
    {
        n_poly_set(A->coeffs + Alen, B->coeffs + i);
        mpoly_monomial_set(A->exps + Alen*N, B->exps + i*N, N);
        i++;
        Alen++;
    }

    while (j < C->length)
    {
        n_poly_set(A->coeffs + Alen, C->coeffs + j);
        mpoly_monomial_set(A->exps + Alen*N, C->exps + j*N, N);
        j++;
        Alen++;
    }

    A->length = Alen;
}


/* A = D - B*C, D may be modified if saveD == 0 */
void _nmod_mpolyn_mulsub(
    nmod_mpolyn_t A,
    const nmod_mpolyn_t D, int saveD,
    const nmod_mpolyn_t B,
    const nmod_mpolyn_t C,
    slong N,
    const ulong * cmpmask,
    const nmod_mpoly_ctx_t ctx)
{
    flint_bitcnt_t bits = A->bits;
    slong i, j;
    slong next_loc;
    slong heap_len = 2; /* heap zero index unused */
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    slong Di;
    slong Alen;
    slong Dlen = D->length;
    slong Blen = B->length;
    slong Clen = C->length;
    ulong * Dexp = D->exps;
    ulong * Bexp = B->exps;
    ulong * Cexp = C->exps;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    slong * hind;
    n_poly_t t;
    TMP_INIT;

    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(bits == A->bits);
    FLINT_ASSERT(bits == D->bits);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == C->bits);

    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(C->length > 0);

    TMP_START;

    n_poly_init(t);

    next_loc = B->length + 4; /* something bigger than heap can ever be */
    heap = (mpoly_heap_s *) TMP_ALLOC((B->length + 1)*sizeof(mpoly_heap_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(B->length*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*B->length*sizeof(slong));
    exps = (ulong *) TMP_ALLOC(B->length*N*sizeof(ulong));
    exp_list = (ulong **) TMP_ALLOC(B->length*sizeof(ulong *));
    for (i = 0; i < B->length; i++)
        exp_list[i] = exps + i*N;

    /* space for heap indices */
    hind = (slong *) TMP_ALLOC(B->length*sizeof(slong));
    for (i = 0; i < B->length; i++)
        hind[i] = 1;

    /* start with no heap nodes and no exponent vectors in use */
    exp_next = 0;

    /* put (0, 0, Bexp[0] + Cexp[0]) on heap */
    x = chain + 0;
    x->i = 0;
    x->j = 0;
    x->next = NULL;

    heap[1].next = x;
    heap[1].exp = exp_list[exp_next++];

    mpoly_monomial_add(heap[1].exp, B->exps + N*0, C->exps + N*0, N);

    hind[0] = 2*1 + 0;

    Alen = 0;
    Di = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        while (Di < Dlen && mpoly_monomial_gt(Dexp + N*Di, exp, N, cmpmask))
        {
            nmod_mpolyn_fit_length(A, Alen + 1, ctx);
            mpoly_monomial_set(A->exps + N*Alen, D->exps + N*Di, N);
            if (saveD)
                n_poly_set(A->coeffs + Alen, D->coeffs + Di);
            else
                n_poly_swap(A->coeffs + Alen, D->coeffs + Di);
            Alen++;
            Di++;
        }

        nmod_mpolyn_fit_length(A, Alen + 1, ctx);

        mpoly_monomial_set(A->exps + N*Alen, exp, N);

        if (Di < Dlen && mpoly_monomial_equal(Dexp + N*Di, exp, N))
        {
            if (saveD)
                n_poly_set(A->coeffs + Alen, D->coeffs + Di);
            else
                n_poly_swap(A->coeffs + Alen, D->coeffs + Di);
            Di++;
        }
        else
        {
            n_poly_zero(A->coeffs + Alen);
        }

        do
        {
            exp_list[--exp_next] = heap[1].exp;
            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
            do
            {
                hind[x->i] |= WORD(1);
                *store++ = x->i;
                *store++ = x->j;
                n_poly_mod_mul(t, B->coeffs + x->i, C->coeffs + x->j, ctx->mod);
                n_poly_mod_sub(A->coeffs + Alen, A->coeffs + Alen, t, ctx->mod);
            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

        Alen += !n_poly_is_zero(A->coeffs + Alen);

        /* process nodes taken from the heap */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            /* should we go right? */
            if (  (i + 1 < Blen)
                && (hind[i + 1] == 2*j + 1)
               )
            {
                x = chain + i + 1;
                x->i = i + 1;
                x->j = j;
                x->next = NULL;

                hind[x->i] = 2*(x->j+1) + 0;

                if (bits <= FLINT_BITS)
                    mpoly_monomial_add(exp_list[exp_next], Bexp + N*x->i,
                                                           Cexp + N*x->j, N);
                else
                    mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i,
                                                              Cexp + N*x->j, N);

                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
            }

            /* should we go up? */
            if (  (j + 1 < Clen)
               && ((hind[i] & 1) == 1)
               && ((i == 0) || (hind[i - 1] >= 2*(j + 2) + 1))
               )
            {
                x = chain + i;
                x->i = i;
                x->j = j + 1;
                x->next = NULL;

                hind[x->i] = 2*(x->j+1) + 0;

                if (bits <= FLINT_BITS)
                    mpoly_monomial_add(exp_list[exp_next], Bexp + N*x->i,
                                                           Cexp + N*x->j, N);
                else
                    mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i,
                                                              Cexp + N*x->j, N);

                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
            }
        }
    }

    FLINT_ASSERT(Di <= Dlen);

    nmod_mpolyn_fit_length(A, Alen + Dlen - Di, ctx);

    for (i = 0; i < Dlen - Di; i++)
        if (saveD)
            n_poly_set(A->coeffs + Alen + i, D->coeffs + Di + i);
        else
            n_poly_swap(A->coeffs + Alen + i, D->coeffs + Di + i);

    mpoly_copy_monomials(A->exps + N*Alen, Dexp + N*Di, Dlen - Di, N);
    Alen += Dlen - Di;

    A->length = Alen;

    n_poly_clear(t);

    TMP_END;
}

/*
    Q = A/B return 1 if division was exact.
    nmainvars is the number of main vars.
*/
int nmod_mpolyun_divides(
    nmod_mpolyun_t Q,
    const nmod_mpolyun_t A,
    const nmod_mpolyun_t B,
    const nmod_mpoly_ctx_t ctx)
{
    slong nmainvars = 1;
    flint_bitcnt_t bits = A->bits;
    nmod_mpolyn_struct * Bcoeff = B->coeffs;
    slong Blen = B->length;
    ulong * Bexp = B->exps;
    nmod_mpolyn_struct * Acoeff = A->coeffs;
    slong Alen = A->length;
    ulong * Aexp = A->exps;
    slong N;
    ulong * cmpmask;    /* cmp mask for lesser variables */
    nmod_mpolyn_t T, S;
    int success;
    ulong maskhi = 0;   /* main variables are in lex */
    int lt_divides;
    slong i, j, s;
    slong next_loc, heap_len;
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    slong * hind;
    ulong mask, exp, maxexp = Aexp[Alen - 1];

    FLINT_ASSERT(Alen > 0);
    FLINT_ASSERT(Blen > 0);
    FLINT_ASSERT(bits == A->bits);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == Q->bits);

    N = mpoly_words_per_exp(bits, ctx->minfo);
    cmpmask = (ulong *) flint_malloc(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, bits, ctx->minfo);

    /* alloc array of heap nodes which can be chained together */
    next_loc = Blen + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap1_s *) flint_malloc((Blen + 1)*sizeof(mpoly_heap1_s));
    chain = (mpoly_heap_t *) flint_malloc(Blen*sizeof(mpoly_heap_t));
    store = store_base = (slong *) flint_malloc(2*Blen*sizeof(mpoly_heap_t *));

    /* space for flagged heap indices */
    hind = (slong *) flint_malloc(Blen*sizeof(slong));
    for (i = 0; i < B->length; i++)
        hind[i] = 1;

    /* mask with high bit set in each field of main exponent vector */
    mask = 0;
    for (i = 0; i < nmainvars; i++)
        mask = (mask << (FLINT_BITS/nmainvars)) + (UWORD(1) << (FLINT_BITS/nmainvars - 1));

    /* s is the number of terms * (latest quotient) we should put into heap */
    s = Blen;

    /* insert (-1, 0, exp2[0]) into heap */
    heap_len = 2;
    x = chain + 0;
    x->i = -WORD(1);
    x->j = 0;
    x->next = NULL;
    HEAP_ASSIGN(heap[1], Aexp[0], x);

    Q->length = 0;

    nmod_mpolyn_init(T, bits, ctx);
    nmod_mpolyn_init(S, bits, ctx);

    while (heap_len > 1)
    {
        exp = heap[1].exp;

        if (mpoly_monomial_overflows1(exp, mask))
        {
            goto not_exact_division;
        }

        nmod_mpolyun_fit_length(Q, Q->length + 1, ctx);
        lt_divides = mpoly_monomial_divides1(Q->exps + Q->length,
                                                           exp, Bexp[0], mask);

        T->length = 0;
        do
        {
            x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
            do
            {
                *store++ = x->i;
                *store++ = x->j;
                if (x->i != -WORD(1))
                    hind[x->i] |= WORD(1);

                if (x->i == -WORD(1))
                {
                    _nmod_mpolyn_add(S, T, Acoeff + x->j, N, ctx);
                }
                else
                {
                    _nmod_mpolyn_mulsub(S, T, 0,
                             Bcoeff + x->i, Q->coeffs + x->j, N, cmpmask, ctx);
                }
                nmod_mpolyn_swap(S, T);
            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && heap[1].exp == exp);

        /* process nodes taken from the heap */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            if (i == -WORD(1))
            {
                /* take next dividend term */
                if (j + 1 < Alen)
                {
                    x = chain + 0;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    _mpoly_heap_insert1(heap, Aexp[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
                }
            }
            else
            {
                /* should we go right? */
                if (  (i + 1 < Blen)
                   && (hind[i + 1] == 2*j + 1)
                   )
                {
                    x = chain + i + 1;
                    x->i = i + 1;
                    x->j = j;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;
                    _mpoly_heap_insert1(heap, Bexp[x->i] + Q->exps[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
                }
                /* should we go up? */
                if (j + 1 == Q->length)
                {
                    s++;
                } else if (  ((hind[i] & 1) == 1)
                          && ((i == 1) || (hind[i - 1] >= 2*(j + 2) + 1))
                          )
                {
                    x = chain + i;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;
                    _mpoly_heap_insert1(heap, Bexp[x->i] + Q->exps[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
                }
            }
        }

        if (T->length == 0)
        {
            continue;
        }

        if (mpoly_monomials_overflow_test(T->exps, T->length, bits, ctx->minfo))
        {
            goto not_exact_division;
        }

        if (!_nmod_mpolyn_divides(Q->coeffs + Q->length,
                                               T, Bcoeff + 0, N, cmpmask, ctx))
        {
            goto not_exact_division;
        }

        if (!lt_divides || (exp^maskhi) < (maxexp^maskhi))
        {
            goto not_exact_division;
        }

        /* put newly generated quotient term back into the heap if necessary */
        if (s > 1)
        {
            i = 1;
            x = chain + i;
            x->i = i;
            x->j = Q->length;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;
            _mpoly_heap_insert1(heap, Bexp[x->i] + Q->exps[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
        }
        s = 1;
        Q->length++;
    }

    success = 1;

cleanup:

    nmod_mpolyn_clear(T, ctx);
    nmod_mpolyn_clear(S, ctx);

    flint_free(cmpmask);
    flint_free(heap);
    flint_free(chain);
    flint_free(store);
    flint_free(hind);

    return success;

not_exact_division:

    success = 0;
    Q->length = 0;
    goto cleanup;
}
