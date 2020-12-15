/*
    Copyright (C) 2017-2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

static int _nmod_mpoly_divrem_monagan_pearce1_binomial(
    nmod_mpoly_t Q,
    nmod_mpoly_t R,
    const mp_limb_t * Acoeffs, const ulong * Aexps, slong Alen,
    const mp_limb_t * Bcoeffs, const ulong * Bexps,
    flint_bitcnt_t bits,
    ulong maskhi,
    nmod_t mod)
{
    mp_limb_t * Qcoeffs = Q->coeffs;
    mp_limb_t * Rcoeffs = R->coeffs;
    ulong * Qexps = Q->exps;
    ulong * Rexps = R->exps;
    ulong lexp, mask = mpoly_overflow_mask_sp(bits);
    mp_limb_t lcoeff;
    mp_limb_t lc_inv = nmod_inv(Bcoeffs[0], mod);
    mp_limb_t mBcoeff1 = mod.n - Bcoeffs[1];
    slong Qlen = 0;
    slong Rlen = 0;
    slong Aidx = 0;
    slong Qidx = 0;

    while (1)
    {
        FLINT_ASSERT(0 <= Aidx && Aidx <= Alen);
        FLINT_ASSERT(0 <= Qidx && Qidx <= Qlen);

        if (Aidx < Alen)
        {
            lexp = Aexps[Aidx];

            if (Qidx < Qlen)
            {
                ulong thisexp = Bexps[1] + Qexps[Qidx];
                int cmp = mpoly_monomial_cmp1(lexp, thisexp, maskhi);
                if (cmp < 0)
                {
                    lexp = thisexp;
                    lcoeff = nmod_mul(mBcoeff1, Qcoeffs[Qidx], mod);
                    Qidx++;
                }
                else if (cmp == 0)
                {
                    lcoeff = Acoeffs[Aidx];
                    NMOD_ADDMUL(lcoeff, mBcoeff1, Qcoeffs[Qidx], mod);
                    Aidx++;
                    Qidx++;
                }
                else
                {
                    lcoeff = Acoeffs[Aidx];
                    Aidx++;
                }
            }
            else
            {
                lcoeff = Acoeffs[Aidx];
                Aidx++;
            }
        }
        else if (Qidx < Qlen)
        {
            lexp = Bexps[1] + Qexps[Qidx];
            lcoeff = nmod_mul(mBcoeff1, Qcoeffs[Qidx], mod);
            Qidx++;
        }
        else
        {
            break;
        }

        if (mpoly_monomial_overflows1(lexp, mask))
            goto exp_overflow;

        if (lcoeff == 0)
            continue;

        _nmod_mpoly_fit_length(&Qcoeffs, &Q->coeffs_alloc,
                               &Qexps, &Q->exps_alloc, 1, Qlen + 1);

        if (!mpoly_monomial_divides1(Qexps + Qlen, lexp, Bexps[0], mask))
        {

            _nmod_mpoly_fit_length(&Rcoeffs, &R->coeffs_alloc,
                                   &Rexps, &R->exps_alloc, 1, Rlen + 1);
            Rcoeffs[Rlen] = lcoeff;
            Rexps[Rlen] = lexp;
            Rlen++;
            continue;
        }

        if (lc_inv == 1)
            Qcoeffs[Qlen] = lcoeff;
        else
            Qcoeffs[Qlen] = nmod_mul(lcoeff, lc_inv, mod);
        Qlen++;
    }

    Q->coeffs = Qcoeffs;
    Q->exps = Qexps;
    Q->length = Qlen;

    R->coeffs = Rcoeffs;
    R->exps = Rexps;
    R->length = Rlen;

    return 1;

exp_overflow:

    Q->coeffs = Qcoeffs;
    Q->exps = Qexps;
    Q->length = 0;

    R->coeffs = Rcoeffs;
    R->exps = Rexps;
    R->length = 0;

    return 0;
}



static int _nmod_mpoly_divrem_monagan_pearce1(
    nmod_mpoly_t Q,
    nmod_mpoly_t R,
    const mp_limb_t * Acoeffs, const ulong * Aexps, slong Alen,
    const mp_limb_t * Bcoeffs, const ulong * Bexps, slong Blen,
    flint_bitcnt_t bits,
    ulong maskhi,
    nmod_t fctx)
{
    slong i, j, Qlen, Rlen, s;
    slong next_loc, heap_len = 2;
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    mp_limb_t * Qcoeffs = Q->coeffs;
    mp_limb_t * Rcoeffs = R->coeffs;
    ulong * Qexps = Q->exps;
    ulong * Rexps = R->exps;
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

    /* quotient and remainder poly indices start at -1 */
    Qlen = WORD(0);
    Rlen = WORD(0);

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
                               &Qexps, &Q->exps_alloc, 1, Qlen + 1);

        lt_divides = mpoly_monomial_divides1(Qexps + Qlen, exp, Bexps[0], mask);

        acc0 = acc1 = acc2 = 0;
        do {
            x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
            do {
                *store++ = x->i;
                *store++ = x->j;

                if (x->i == -WORD(1))
                {
                    add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0,
                                 UWORD(0), UWORD(0), fctx.n - Acoeffs[x->j]);
                }
                else
                {
                    hind[x->i] |= WORD(1);
                    umul_ppmm(pp1, pp0, Bcoeffs[x->i], Qcoeffs[x->j]);
                    add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, UWORD(0), pp1, pp0);
                }

            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && heap[1].exp == exp);

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
                    _mpoly_heap_insert1(heap, Aexps[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
                }
            } else
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
                    _mpoly_heap_insert1(heap, Bexps[x->i] + Qexps[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
                }
            }
        }

        /* try to divide accumulated term by leading term */
        if (acc0 == 0)
        {
            continue;
        }
        if (!lt_divides)
        {
            _nmod_mpoly_fit_length(&Rcoeffs, &R->coeffs_alloc,
                                   &Rexps, &R->exps_alloc, 1, Rlen + 1);
            Rcoeffs[Rlen] = fctx.n - acc0;
            Rexps[Rlen] = exp;
            Rlen++;
            continue;
        }

        Qcoeffs[Qlen] = nmod_mul(acc0, lc_minus_inv, fctx);

        /* put newly generated quotient term back into the heap if neccesary */
        if (s > 1)
        {
            i = 1;
            x = chain + i;
            x->i = i;
            x->j = Qlen;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;
            _mpoly_heap_insert1(heap, Bexps[x->i] + Qexps[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
        }
        s = 1;
        Qlen++;
    }

    Q->coeffs = Qcoeffs;
    Q->exps = Qexps;
    Q->length = Qlen;

    R->coeffs = Rcoeffs;
    R->exps = Rexps;
    R->length = Rlen;

    TMP_END;

    return 1;

exp_overflow:

    Q->coeffs = Qcoeffs;
    Q->exps = Qexps;
    Q->length = 0;

    R->coeffs = Rcoeffs;
    R->exps = Rexps;
    R->length = 0;

    TMP_END;

    return 0;
}



static int _nmod_mpoly_divrem_monagan_pearce(
    nmod_mpoly_t Q,
    nmod_mpoly_t R,
    const mp_limb_t * Acoeffs, const ulong * Aexps, slong Alen,
    const mp_limb_t * Bcoeffs, const ulong * Bexps, slong Blen,
    slong bits,
    slong N,
    const ulong * cmpmask,
    nmod_t fctx)
{
    slong i, j, Qlen, Rlen, s;
    slong next_loc;
    slong heap_len = 2; /* heap zero index unused */
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    mp_limb_t * Qcoeffs = Q->coeffs;
    mp_limb_t * Rcoeffs = R->coeffs;
    ulong * Qexps = Q->exps;
    ulong * Rexps = R->exps;
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
        if (Blen == 2)
            return _nmod_mpoly_divrem_monagan_pearce1_binomial(Q, R,
                 Acoeffs, Aexps, Alen, Bcoeffs, Bexps, bits, cmpmask[0], fctx);
        else
            return _nmod_mpoly_divrem_monagan_pearce1(Q, R, Acoeffs, Aexps, Alen,
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

    Qlen = WORD(0);
    Rlen = WORD(0);
   
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
                               &Qexps, &Q->exps_alloc, N, Qlen + 1);

        mpoly_monomial_set(exp, heap[1].exp, N);

        if (bits <= FLINT_BITS)
        {
            if (mpoly_monomial_overflows(exp, N, mask))
                goto exp_overflow2;

            lt_divides = mpoly_monomial_divides(Qexps + Qlen*N, exp, Bexps, N, mask);
        }
        else
        {
            if (mpoly_monomial_overflows_mp(exp, N, bits))
                goto exp_overflow2;

            lt_divides = mpoly_monomial_divides_mp(Qexps + Qlen*N, exp, Bexps, N, bits);
        }

        acc0 = acc1 = acc2 = 0;
        do {
            exp_list[--exp_next] = heap[1].exp;
            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
            do {
                *store++ = x->i;
                *store++ = x->j;

                if (x->i == -WORD(1))
                {
                    add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0,
                                   UWORD(0), UWORD(0), fctx.n - Acoeffs[x->j]);
                }
                else
                {
                    hind[x->i] |= WORD(1);
                    umul_ppmm(pp1, pp0, Bcoeffs[x->i], Qcoeffs[x->j]);
                    add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, UWORD(0), pp1, pp0);
                }
            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

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
                    mpoly_monomial_add_mp(exp_list[exp_next], Bexps + x->i*N,
                                                            Qexps + x->j*N, N);
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
        {
            _nmod_mpoly_fit_length(&Rcoeffs, &R->coeffs_alloc,
                                   &Rexps, &R->exps_alloc, N, Rlen + 1);
            Rcoeffs[Rlen] = fctx.n - acc0;
            mpoly_monomial_set(Rexps + Rlen*N, exp, N);
            Rlen++;
            continue;
        }

        Qcoeffs[Qlen] = nmod_mul(acc0, lc_minus_inv, fctx);

        /* put newly generated quotient term back into the heap if neccesary */
        if (s > 1)
        {
            i = 1;
            x = chain + i;
            x->i = i;
            x->j = Qlen;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;
            mpoly_monomial_add_mp(exp_list[exp_next], Bexps + x->i*N,
                                                     Qexps + x->j*N, N);
            exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
        }
        s = 1;
        Qlen++;
    }

    Q->coeffs = Qcoeffs;
    Q->exps = Qexps;
    Q->length = Qlen;

    R->coeffs = Rcoeffs;
    R->exps = Rexps;
    R->length = Rlen;

    TMP_END;

    return 1;

exp_overflow2:

    Q->coeffs = Qcoeffs;
    Q->exps = Qexps;
    Q->length = 0;

    R->coeffs = Rcoeffs;
    R->exps = Rexps;
    R->length = 0;

    TMP_END;

    return 0;
}

void nmod_mpoly_divrem_monagan_pearce(
    nmod_mpoly_t Q,
    nmod_mpoly_t R,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx)
{
    slong N;
    flint_bitcnt_t QRbits;
    ulong * Aexps = A->exps, * Bexps = B->exps;
    ulong * cmpmask;
    int freeAexps = 0, freeBexps = 0;
    nmod_mpoly_t TQ, TR;
    nmod_mpoly_struct * q, * r;

    if (nmod_mpoly_is_zero(B, ctx))
    {
        if (nmod_mpoly_ctx_modulus(ctx) == 1)
        {
            nmod_mpoly_zero(Q, ctx);
            nmod_mpoly_zero(R, ctx);
            return;
        }
        else
        {
            flint_throw(FLINT_DIVZERO, "nmod_mpoly_divrem_monagan_pearce: divide by zero");
        }
    }

    if (nmod_mpoly_is_zero(A, ctx))
    {
        nmod_mpoly_zero(Q, ctx);
        nmod_mpoly_zero(R, ctx);
        return;
    }

    nmod_mpoly_init(TQ, ctx);
    nmod_mpoly_init(TR, ctx);

    QRbits = FLINT_MAX(A->bits, B->bits);
    QRbits = mpoly_fix_bits(QRbits, ctx->minfo);

    N = mpoly_words_per_exp(QRbits, ctx->minfo);
    cmpmask = (ulong *) flint_malloc(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, QRbits, ctx->minfo);

    /* ensure input exponents packed to same size as output exponents */
    if (QRbits != A->bits)
    {
        freeAexps = 1;
        Aexps = (ulong *) flint_malloc(N*A->length*sizeof(ulong));
        mpoly_repack_monomials(Aexps, QRbits, A->exps, A->bits, A->length, ctx->minfo);
    }

    if (QRbits != B->bits)
    {
        freeBexps = 1;
        Bexps = (ulong *) flint_malloc(N*B->length*sizeof(ulong));
        mpoly_repack_monomials(Bexps, QRbits, B->exps, B->bits, B->length, ctx->minfo);
    }

    /* check divisor leading monomial is at most that of the dividend */
    if (mpoly_monomial_lt(Aexps, Bexps, N, cmpmask))
    {
        nmod_mpoly_set(R, A, ctx);
        nmod_mpoly_zero(Q, ctx);
        goto cleanup;
    }

    /* take care of aliasing */
    if (Q == A || Q == B)
        q = TQ;
    else
        q = Q;

    if (R == A || R == B)
        r = TR;
    else
        r = R;

    /* do division with remainder */
    while (1)
    {
        nmod_mpoly_fit_length_reset_bits(q, A->length/B->length + 1, QRbits, ctx);
        nmod_mpoly_fit_length_reset_bits(r, B->length, QRbits, ctx);

        if (_nmod_mpoly_divrem_monagan_pearce(q, r,
                     A->coeffs, Aexps, A->length, B->coeffs, Bexps, B->length,
                                                 QRbits, N, cmpmask, ctx->mod))
        {
            break;
        }

        QRbits = mpoly_fix_bits(QRbits + 1, ctx->minfo);

        N = mpoly_words_per_exp(QRbits, ctx->minfo);
        cmpmask = (ulong *) flint_realloc(cmpmask, N*sizeof(ulong));
        mpoly_get_cmpmask(cmpmask, N, QRbits, ctx->minfo);

        if (freeAexps)
            flint_free(Aexps);
        Aexps = (ulong *) flint_malloc(N*A->length*sizeof(ulong));
        mpoly_repack_monomials(Aexps, QRbits, A->exps, A->bits, A->length, ctx->minfo);
        freeAexps = 1; 

        if (freeBexps)
            flint_free(Bexps);
        Bexps = (ulong *) flint_malloc(N*B->length*sizeof(ulong));
        mpoly_repack_monomials(Bexps, QRbits, B->exps, B->bits, B->length, ctx->minfo);
        freeBexps = 1; 
    }

    /* deal with aliasing */
    if (Q == A || Q == B)
        nmod_mpoly_swap(Q, TQ, ctx);

    if (R == A || R == B)
        nmod_mpoly_swap(R, TR, ctx);

cleanup:

    nmod_mpoly_clear(TQ, ctx);
    nmod_mpoly_clear(TR, ctx);

    if (freeAexps)
        flint_free(Aexps);

    if (freeBexps)
        flint_free(Bexps);

    flint_free(cmpmask);
}
