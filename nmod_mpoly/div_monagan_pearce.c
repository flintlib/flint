/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

static int _nmod_mpoly_div_monagan_pearce1(
    nmod_mpoly_t Q,
    const mp_limb_t * Acoeffs, const ulong * Aexps, slong Alen,
    const mp_limb_t * Bcoeffs, const ulong * Bexps, slong Blen,
    flint_bitcnt_t bits,
    ulong maskhi,
    nmod_t fctx)
{
    slong i, j, q_len, s;
    slong next_loc, heap_len = 2;
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    mp_limb_t * Qcoeffs = Q->coeffs;
    ulong * Qexps = Q->exps;
    slong * hind;
    ulong mask, exp;
    int lt_divides;
    mp_limb_t lc_minus_inv, acc0, acc1, acc2, pp1, pp0;
    TMP_INIT;

    TMP_START;

    /* alloc array of heap nodes which can be chained together */
    next_loc = Blen + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap1_s *) TMP_ALLOC((Blen + 1)*sizeof(mpoly_heap1_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(Blen*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*Blen*sizeof(slong));

    /* space for flagged heap indicies */
    hind = (slong *) TMP_ALLOC(Blen*sizeof(slong));
    for (i = 0; i < Blen; i++)
        hind[i] = 1;

    mask = mpoly_overflow_mask_sp(bits);

    q_len = WORD(0);

    /* s is the number of terms * (latest quotient) we should put into heap */
    s = Blen;

    /* insert (-1, 0, Aexps[0]) into heap */
    x = chain + 0;
    x->i = -WORD(1);
    x->j = 0;
    x->next = NULL;
    HEAP_ASSIGN(heap[1], Aexps[0], x);

    /* precompute leading cofficient info */
    lc_minus_inv = fctx.n - nmod_inv(Bcoeffs[0], fctx);

    while (heap_len > 1)
    {
        exp = heap[1].exp;

        if (mpoly_monomial_overflows1(exp, mask))
            goto exp_overflow;

        _nmod_mpoly_fit_length(&Qcoeffs, &Q->coeffs_alloc,
                               &Qexps, &Q->exps_alloc, 1, q_len + 1);

        lt_divides = mpoly_monomial_divides1(Qexps + q_len, exp, Bexps[0], mask);

        acc0 = acc1 = acc2 = 0;
        if (!lt_divides)
        {
            /* optimation: coeff arithmetic not needed */

            if (mpoly_monomial_gt1(Bexps[0], exp, maskhi))
            {
                /* optimization: no more quotient terms possible */
                goto cleanup;
            }

            do {
                x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
                do {
                    *store++ = x->i;
                    *store++ = x->j;
                    if (x->i != -WORD(1))
                        hind[x->i] |= WORD(1);

                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);
        }
        else
        {
            do {
                x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
                do {
                    *store++ = x->i;
                    *store++ = x->j;

                    if (x->i == -WORD(1))
                    {
                        add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), WORD(0), fctx.n - Acoeffs[x->j]);
                    }
                    else
                    {
                        hind[x->i] |= WORD(1);
                        umul_ppmm(pp1, pp0, Bcoeffs[x->i], Qcoeffs[x->j]);
                        add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);
                    }
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);

            NMOD_RED3(acc0, acc2, acc1, acc0, fctx);
        }

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
                    _mpoly_heap_insert1(heap, Aexps[x->j], x,
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
                    _mpoly_heap_insert1(heap, Bexps[x->i] + Qexps[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
                }
                /* should we go up? */
                if (j + 1 == q_len)
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
                    _mpoly_heap_insert1(heap, Bexps[x->i] + Qexps[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
                }
            }
        }

        /* try to divide accumulated term by leading term */
        if (!lt_divides)
            continue;

        if (acc0 == 0)
            continue;

        Qcoeffs[q_len] = nmod_mul(acc0, lc_minus_inv, fctx);

        /* put newly generated quotient term back into the heap if neccesary */
        if (s > 1)
        {
            i = 1;
            x = chain + i;
            x->i = i;
            x->j = q_len;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;
            _mpoly_heap_insert1(heap, Bexps[x->i] + Qexps[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
        }
        s = 1;
        q_len++;
    }

cleanup:

    Q->coeffs = Qcoeffs;
    Q->exps = Qexps;
    Q->length = q_len;

    TMP_END;

    return 1;

exp_overflow:

    Q->coeffs = Qcoeffs;
    Q->exps = Qexps;
    Q->length = 0;

    TMP_END;

    return 0;
}




static int _nmod_mpoly_div_monagan_pearce(
    nmod_mpoly_t Q,
    const mp_limb_t * Acoeffs, const ulong * Aexps, slong Alen,
    const mp_limb_t * Bcoeffs, const ulong * Bexps, slong Blen,
    flint_bitcnt_t bits,
    slong N,
    const ulong * cmpmask,
    nmod_t fctx)
{
    slong i, j, q_len, s;
    slong next_loc;
    slong heap_len = 2; /* heap zero index unused */
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    mp_limb_t * Qcoeffs = Q->coeffs;
    ulong * Qexps = Q->exps;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    ulong mask;
    slong * hind;
    int lt_divides;
    mp_limb_t lc_minus_inv, acc0, acc1, acc2, pp1, pp0;
    TMP_INIT;

    if (N == 1)
    {
        return _nmod_mpoly_div_monagan_pearce1(Q, Acoeffs, Aexps, Alen,
                                   Bcoeffs, Bexps, Blen, bits, cmpmask[0], fctx);
    }

    TMP_START;

    /* alloc array of heap nodes which can be chained together */
    next_loc = Blen + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap_s *) TMP_ALLOC((Blen + 1)*sizeof(mpoly_heap_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(Blen*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*Blen*sizeof(slong));

    /* array of exponent vectors, each of "N" words */
    exps = (ulong *) TMP_ALLOC(Blen*N*sizeof(ulong));
    /* list of pointers to available exponent vectors */
    exp_list = (ulong **) TMP_ALLOC(Blen*sizeof(ulong *));
    /* space to save copy of current exponent vector */
    exp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    /* set up list of available exponent vectors */
    exp_next = 0;
    for (i = 0; i < Blen; i++)
        exp_list[i] = exps + i*N;

    /* space for flagged heap indicies */
    hind = (slong *) TMP_ALLOC(Blen*sizeof(slong));
    for (i = 0; i < Blen; i++)
        hind[i] = 1;

    mask = bits <= FLINT_BITS ? mpoly_overflow_mask_sp(bits) : 0;

    q_len = WORD(0);
   
    /* s is the number of terms * (latest quotient) we should put into heap */
    s = Blen;
   
    /* insert (-1, 0, Aexps[0]) into heap */
    x = chain + 0;
    x->i = -WORD(1);
    x->j = 0;
    x->next = NULL;
    heap[1].next = x;
    heap[1].exp = exp_list[exp_next++];
    mpoly_monomial_set(heap[1].exp, Aexps, N);

    /* precompute leading cofficient info */
    lc_minus_inv = fctx.n - nmod_inv(Bcoeffs[0], fctx);
   
    while (heap_len > 1)
    {
        _nmod_mpoly_fit_length(&Qcoeffs, &Q->coeffs_alloc,
                               &Qexps, &Q->exps_alloc, N, q_len + 1);

        mpoly_monomial_set(exp, heap[1].exp, N);

        if (bits <= FLINT_BITS)
        {
            if (mpoly_monomial_overflows(exp, N, mask))
                goto exp_overflow;

            lt_divides = mpoly_monomial_divides(Qexps + q_len*N, exp, Bexps, N, mask);
        }
        else
        {
            if (mpoly_monomial_overflows_mp(exp, N, bits))
                goto exp_overflow;

            lt_divides = mpoly_monomial_divides_mp(Qexps + q_len*N, exp, Bexps, N, bits);
        }

        acc0 = acc1 = acc2 = 0;

        if (!lt_divides)
        {
            /* optimation: coeff arithmetic not needed */

            if (mpoly_monomial_gt(Bexps + 0, exp, N, cmpmask))
            {
                /* optimization: no more quotient terms possible */
                goto cleanup;
            }

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

                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));
        }
        else
        {
            do {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
                do {
                    *store++ = x->i;
                    *store++ = x->j;

                    if (x->i == -WORD(1))
                    {
                        add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0,
                                     WORD(0), WORD(0), fctx.n - Acoeffs[x->j]);
                    }
                    else
                    {
                        hind[x->i] |= WORD(1);
                        umul_ppmm(pp1, pp0, Bcoeffs[x->i], Qcoeffs[x->j]);
                        add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);
                    }
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));
        }

        NMOD_RED3(acc0, acc2, acc1, acc0, fctx);

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
                    mpoly_monomial_set(exp_list[exp_next], Aexps + x->j*N, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
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

                    if (bits <= FLINT_BITS)
                        mpoly_monomial_add(exp_list[exp_next], Bexps + x->i*N,
                                                           Qexps + x->j*N, N);
                    else
                        mpoly_monomial_add_mp(exp_list[exp_next], Bexps + x->i*N,
                                                           Qexps + x->j*N, N);

                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
                /* should we go up? */
                if (j + 1 == q_len)
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

                    if (bits <= FLINT_BITS)
                        mpoly_monomial_add(exp_list[exp_next], Bexps + x->i*N,
                                                           Qexps + x->j*N, N);
                    else
                        mpoly_monomial_add_mp(exp_list[exp_next], Bexps + x->i*N,
                                                           Qexps + x->j*N, N);

                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
            }
        }

        /* try to divide accumulated term by leading term */

        if (acc0 == 0)
            continue;

        if (!lt_divides)
            continue;

        Qcoeffs[q_len] = nmod_mul(acc0, lc_minus_inv, fctx);

        /* put newly generated quotient term back into the heap if neccesary */
        if (s > 1)
        {
            i = 1;
            x = chain + i;
            x->i = i;
            x->j = q_len;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;
            mpoly_monomial_add_mp(exp_list[exp_next], Bexps + x->i*N,
                                                      Qexps + x->j*N, N);
            exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
        }
        s = 1;
        q_len++;
    }

cleanup:

    Q->coeffs = Qcoeffs;
    Q->exps = Qexps;
    Q->length = q_len;

    TMP_END;

    return 1;

exp_overflow:

    Q->coeffs = Qcoeffs;
    Q->exps = Qexps;
    Q->length = 0;

    TMP_END;

    return 0;
}


void nmod_mpoly_div_monagan_pearce(
    nmod_mpoly_t Q,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx)
{
    slong N;
    flint_bitcnt_t Qbits;
    ulong * Aexps = A->exps, * Bexps = B->exps;
    ulong * cmpmask;
    int freeAexps = 0, freeBexps = 0;
    nmod_mpoly_t TQ;
    nmod_mpoly_struct * q;

    if (nmod_mpoly_is_zero(B, ctx))
    {
        flint_throw(FLINT_DIVZERO, "nmod_mpoly_div_monagan_pearce: divide by zero");
    }

    if (nmod_mpoly_is_zero(A, ctx))
    {
        nmod_mpoly_zero(Q, ctx);
        return;
    }

    nmod_mpoly_init(TQ, ctx);

    Qbits = FLINT_MAX(A->bits, B->bits);
    Qbits = mpoly_fix_bits(Qbits, ctx->minfo);

    N = mpoly_words_per_exp(Qbits, ctx->minfo);
    cmpmask = (ulong *) flint_malloc(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, Qbits, ctx->minfo);

    /* ensure input exponents packed to same size as output exponents */
    if (Qbits > A->bits)
    {
        freeAexps = 1;
        Aexps = (ulong *) flint_malloc(N*A->length*sizeof(ulong));
        mpoly_repack_monomials(Aexps, Qbits, A->exps, A->bits, A->length, ctx->minfo);
    }

    if (Qbits > B->bits)
    {
        freeBexps = 1;
        Bexps = (ulong *) flint_malloc(N*B->length*sizeof(ulong));
        mpoly_repack_monomials(Bexps, Qbits, B->exps, B->bits, B->length, ctx->minfo);
    }

    /* check divisor leading monomial is at most that of the dividend */
    if (mpoly_monomial_lt(Aexps, Bexps, N, cmpmask))
    {
        nmod_mpoly_zero(Q, ctx);
        goto cleanup;
    }

    /* take care of aliasing */
    if (Q == A || Q == B)
        q = TQ;
    else
        q = Q;

    /* do division with remainder */
    while (1)
    {
        nmod_mpoly_fit_length_reset_bits(q, A->length/B->length + 1, Qbits, ctx);

        if (_nmod_mpoly_div_monagan_pearce(q, A->coeffs, Aexps, A->length,
                     B->coeffs, Bexps, B->length, Qbits, N, cmpmask, ctx->mod))
        {
            break;
        }

        Qbits = mpoly_fix_bits(Qbits + 1, ctx->minfo);

        N = mpoly_words_per_exp(Qbits, ctx->minfo);
        cmpmask = (ulong *) flint_realloc(cmpmask, N*sizeof(ulong));
        mpoly_get_cmpmask(cmpmask, N, Qbits, ctx->minfo);

        if (freeAexps)
            flint_free(Aexps);
        Aexps = (ulong *) flint_malloc(N*A->length*sizeof(ulong));
        mpoly_repack_monomials(Aexps, Qbits, A->exps, A->bits, A->length, ctx->minfo);
        freeAexps = 1; 

        if (freeBexps)
            flint_free(Bexps);
        Bexps = (ulong *) flint_malloc(N*B->length*sizeof(ulong));
        mpoly_repack_monomials(Bexps, Qbits, B->exps, B->bits, B->length, ctx->minfo);
        freeBexps = 1; 
    }

    /* deal with aliasing */
    if (Q == A || Q == B)
        nmod_mpoly_swap(Q, TQ, ctx);

cleanup:

    nmod_mpoly_clear(TQ, ctx);

    if (freeAexps)
        flint_free(Aexps);

    if (freeBexps)
        flint_free(Bexps);

    flint_free(cmpmask);
}
