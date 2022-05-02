/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

static int _fq_nmod_mpoly_divrem_monagan_pearce1_binomial(
    fq_nmod_mpoly_t Q,
    fq_nmod_mpoly_t R,
    const mp_limb_t * Acoeffs, const ulong * Aexps, slong Alen,
    const mp_limb_t * Bcoeffs, const ulong * Bexps,
    flint_bitcnt_t bits,
    ulong maskhi,
    const fq_nmod_ctx_t fqctx)
{
    slong d = fq_nmod_ctx_degree(fqctx);
    mp_limb_t * Qcoeffs = Q->coeffs;
    mp_limb_t * Rcoeffs = R->coeffs;
    ulong * Qexps = Q->exps;
    ulong * Rexps = R->exps;
    ulong lexp, mask = mpoly_overflow_mask_sp(bits);
    mp_limb_t * tmp, * lc_inv, * mBcoeff1;
    int lc_inv_is_one;
    slong Qlen = 0;
    slong Rlen = 0;
    slong Aidx = 0;
    slong Qidx = 0;
    TMP_INIT;

    TMP_START;

    tmp = (mp_limb_t *) TMP_ALLOC(d*(2 + FLINT_MAX(N_FQ_MUL_ITCH,
                                            N_FQ_INV_ITCH))*sizeof(mp_limb_t));
    lc_inv = tmp + d*FLINT_MAX(N_FQ_MUL_ITCH, N_FQ_INV_ITCH);
    mBcoeff1 = lc_inv + d;

    _n_fq_inv(lc_inv, Bcoeffs + d*0, fqctx, tmp);
    _n_fq_neg(mBcoeff1, Bcoeffs + d*1, d, fqctx->mod);

    lc_inv_is_one = _n_fq_is_one(lc_inv, d);

    while (1)
    {
        FLINT_ASSERT(0 <= Aidx && Aidx <= Alen);
        FLINT_ASSERT(0 <= Qidx && Qidx <= Qlen);

        _fq_nmod_mpoly_fit_length(&Qcoeffs, &Q->coeffs_alloc, d,
                                  &Qexps, &Q->exps_alloc, 1, Qlen + 1);

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
                    _n_fq_mul(Qcoeffs + d*Qlen, mBcoeff1,
                                                Qcoeffs + d*Qidx, fqctx, tmp);
                    Qidx++;
                }
                else if (cmp == 0)
                {
                    _n_fq_mul(Qcoeffs + d*Qlen, mBcoeff1,
                                                Qcoeffs + d*Qidx, fqctx, tmp);
                    _n_fq_add(Qcoeffs + d*Qlen, Qcoeffs + d*Qlen,
                                                Acoeffs + d*Aidx, d, fqctx->mod);
                    Aidx++;
                    Qidx++;
                }
                else
                {
                    _n_fq_set(Qcoeffs + d*Qlen, Acoeffs + d*Aidx, d);
                    Aidx++;
                }
            }
            else
            {
                _n_fq_set(Qcoeffs + d*Qlen, Acoeffs + d*Aidx, d);
                Aidx++;
            }
        }
        else if (Qidx < Qlen)
        {
            lexp = Bexps[1] + Qexps[Qidx];
            _n_fq_mul(Qcoeffs + d*Qlen, mBcoeff1, Qcoeffs + d*Qidx, fqctx, tmp);
            Qidx++;
        }
        else
        {
            break;
        }

        if (mpoly_monomial_overflows1(lexp, mask))
            goto exp_overflow;

        if (_n_fq_is_zero(Qcoeffs + d*Qlen, d))
            continue;

        if (!mpoly_monomial_divides1(Qexps + Qlen, lexp, Bexps[0], mask))
        {
            _fq_nmod_mpoly_fit_length(&Rcoeffs, &R->coeffs_alloc, d,
                                      &Rexps, &R->exps_alloc, 1, Rlen + 1);
            _n_fq_set(Rcoeffs + d*Rlen, Qcoeffs + d*Qlen, d);
            Rexps[Rlen] = lexp;
            Rlen++;
            continue;
        }

        if (!lc_inv_is_one)
            _n_fq_mul(Qcoeffs + d*Qlen, Qcoeffs + d*Qlen, lc_inv, fqctx, tmp);
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

static int _fq_nmod_mpoly_divrem_monagan_pearce1(
    fq_nmod_mpoly_t Q,
    fq_nmod_mpoly_t R,
    const mp_limb_t * Acoeffs, const ulong * Aexps, slong Alen,
    const mp_limb_t * Bcoeffs, const ulong * Bexps, slong Blen,
    flint_bitcnt_t bits,
    ulong maskhi,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
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
    mp_limb_t * lc_minus_inv, * t;
    int lazy_size = _n_fq_dot_lazy_size(Blen, ctx);
    TMP_INIT;

    TMP_START;

    t = (mp_limb_t *) TMP_ALLOC(6*d*sizeof(mp_limb_t));
    lc_minus_inv = (mp_limb_t *) TMP_ALLOC(d*sizeof(mp_limb_t));

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
    _n_fq_inv(lc_minus_inv, Bcoeffs + d*0, ctx, t);
    _n_fq_neg(lc_minus_inv, lc_minus_inv, d, ctx->mod);

    while (heap_len > 1)
    {
        exp = heap[1].exp;

        if (mpoly_monomial_overflows1(exp, mask))
            goto exp_overflow;

        _fq_nmod_mpoly_fit_length(&Qcoeffs, &Q->coeffs_alloc, d,
                                  &Qexps, &Q->exps_alloc, 1, Qlen + 1);

        lt_divides = mpoly_monomial_divides1(Qexps + Qlen, exp, Bexps[0], mask);

        _n_fq_zero(Qcoeffs + d*Qlen, d);
        _nmod_vec_zero(t, 6*d);

        switch (lazy_size)
        {
        case 1:
            do {
                x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
                do {
                    *store++ = x->i;
                    *store++ = x->j;
                    if (x->i == -WORD(1))
                    {
                        n_fq_sub(Qcoeffs + d*Qlen, Qcoeffs + d*Qlen,
                                                   Acoeffs + d*x->j, ctx);
                    }
                    else
                    {
                        hind[x->i] |= WORD(1);
                        _n_fq_madd2_lazy1(t, Bcoeffs + d*x->i,
                                             Qcoeffs + d*x->j, d);
                    }

                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);
            _n_fq_reduce2_lazy1(t, d, ctx->mod);
            break;

        case 2:
            do {
                x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
                do {
                    *store++ = x->i;
                    *store++ = x->j;
                    if (x->i == -WORD(1))
                    {
                        n_fq_sub(Qcoeffs + d*Qlen, Qcoeffs + d*Qlen,
                                                   Acoeffs + d*x->j, ctx);
                    }
                    else
                    {
                        hind[x->i] |= WORD(1);
                        _n_fq_madd2_lazy2(t, Bcoeffs + d*x->i,
                                             Qcoeffs + d*x->j, d);
                    }

                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);
            _n_fq_reduce2_lazy2(t, d, ctx->mod);
            break;

        case 3:
            do {
                x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
                do {
                    *store++ = x->i;
                    *store++ = x->j;
                    if (x->i == -WORD(1))
                    {
                        n_fq_sub(Qcoeffs + d*Qlen, Qcoeffs + d*Qlen,
                                                   Acoeffs + d*x->j, ctx);
                    }
                    else
                    {
                        hind[x->i] |= WORD(1);
                        _n_fq_madd2_lazy3(t, Bcoeffs + d*x->i,
                                             Qcoeffs + d*x->j, d);
                    }

                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);
            _n_fq_reduce2_lazy3(t, d, ctx->mod);
            break;

        default:
            do {
                x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
                do {
                    *store++ = x->i;
                    *store++ = x->j;
                    if (x->i == -WORD(1))
                    {
                        n_fq_sub(Qcoeffs + d*Qlen, Qcoeffs + d*Qlen,
                                                   Acoeffs + d*x->j, ctx);
                    }
                    else
                    {
                        hind[x->i] |= WORD(1);
                        _n_fq_madd2(t, Bcoeffs + d*x->i,
                                       Qcoeffs + d*x->j, ctx, t + 2*d);
                    }

                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);
            break;
        }

        _nmod_vec_add(t, t, Qcoeffs + d*Qlen, d, ctx->mod);
        _n_fq_reduce2(Qcoeffs + d*Qlen, t, ctx, t + 2*d);

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
        if (_n_fq_is_zero(Qcoeffs + d*Qlen, d))
            continue;

        if (!lt_divides)
        {
            _fq_nmod_mpoly_fit_length(&Rcoeffs, &R->coeffs_alloc, d,
                                      &Rexps, &R->exps_alloc, 1, Rlen + 1);
            _n_fq_neg(Rcoeffs + d*Rlen, Qcoeffs + d*Qlen, d, ctx->mod);
            Rexps[Rlen] = exp;
            Rlen++;
            continue;
        }

        _n_fq_mul(Qcoeffs + d*Qlen, Qcoeffs + d*Qlen, lc_minus_inv, ctx, t);

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


static int _fq_nmod_mpoly_divrem_monagan_pearce(
    fq_nmod_mpoly_t Q,
    fq_nmod_mpoly_t R,
    const mp_limb_t * Acoeffs, const ulong * Aexps, slong Alen,
    const mp_limb_t * Bcoeffs, const ulong * Bexps, slong Blen,
    flint_bitcnt_t bits,
    slong N,
    const ulong * cmpmask,
    const fq_nmod_ctx_t fqctx)
{
    slong d = fq_nmod_ctx_degree(fqctx);
    slong i, j, s;
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
    slong Qlen;
    slong Rlen;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    ulong mask;
    slong * hind;
    int lt_divides;
    mp_limb_t * lc_minus_inv, * pp;
    TMP_INIT;

    if (N == 1)
    {
        if (Blen == 2)
            return _fq_nmod_mpoly_divrem_monagan_pearce1_binomial(Q, R,
                Acoeffs, Aexps, Alen, Bcoeffs, Bexps, bits, cmpmask[0], fqctx);
        else        
            return _fq_nmod_mpoly_divrem_monagan_pearce1(Q, R, Acoeffs, Aexps,
                          Alen, Bcoeffs, Bexps, Blen, bits, cmpmask[0], fqctx);
    }

    TMP_START;

    pp = (mp_limb_t *) TMP_ALLOC(d*sizeof(mp_limb_t));
    lc_minus_inv = (mp_limb_t *) TMP_ALLOC(d*sizeof(mp_limb_t));

    /* alloc array of heap nodes which can be chained together */
    next_loc = Blen + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap_s *) TMP_ALLOC((Blen + 1)*sizeof(mpoly_heap_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(Blen*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*Blen*sizeof(slong));

    exps = (ulong *) TMP_ALLOC(Blen*N*sizeof(ulong));
    exp_list = (ulong **) TMP_ALLOC(Blen*sizeof(ulong *));
    exp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    exp_next = 0;
    for (i = 0; i < Blen; i++)
        exp_list[i] = exps + i*N;

    hind = (slong *) TMP_ALLOC(Blen*sizeof(slong));
    for (i = 0; i < Blen; i++)
        hind[i] = 1;

    mask = bits <= FLINT_BITS ? mpoly_overflow_mask_sp(bits) : 0;

    Qlen = 0;
    Rlen = 0;
   
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
    n_fq_inv(lc_minus_inv, Bcoeffs + d*0, fqctx);
    _n_fq_neg(lc_minus_inv, lc_minus_inv, d, fqctx->mod);
   
    while (heap_len > 1)
    {
        _fq_nmod_mpoly_fit_length(&Qcoeffs, &Q->coeffs_alloc, d,
                                  &Qexps, &Q->exps_alloc, N, Qlen + 1);

        mpoly_monomial_set(exp, heap[1].exp, N);

        if (bits <= FLINT_BITS)
        {
            if (mpoly_monomial_overflows(exp, N, mask))
                goto exp_overflow;

            lt_divides = mpoly_monomial_divides(Qexps + N*Qlen, exp, Bexps, N, mask);
        }
        else
        {
            if (mpoly_monomial_overflows_mp(exp, N, bits))
                goto exp_overflow;

            lt_divides = mpoly_monomial_divides_mp(Qexps + N*Qlen, exp, Bexps, N, bits);
        }

        _n_fq_zero(Qcoeffs + d*Qlen, d);
        do {
            exp_list[--exp_next] = heap[1].exp;
            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
            do {
                *store++ = x->i;
                *store++ = x->j;

                if (x->i == -WORD(1))
                {
                    n_fq_sub(Qcoeffs + d*Qlen, Qcoeffs + d*Qlen, Acoeffs + d*x->j, fqctx);
                }
                else
                {
                    hind[x->i] |= WORD(1);
                    n_fq_mul(pp, Bcoeffs + d*x->i, Qcoeffs + d*x->j, fqctx);
                    n_fq_add(Qcoeffs + d*Qlen, Qcoeffs + d*Qlen, pp, fqctx);
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
                if (j + 1 < Alen)
                {
                    x = chain + 0;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    mpoly_monomial_set(exp_list[exp_next], Aexps + N*x->j, N);
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
                    mpoly_monomial_add_mp(exp_list[exp_next], Bexps + N*x->i,
                                                             Qexps + N*x->j, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
                /* should we go up? */
                if (j + 1 == Qlen)
                {
                    s++;
                }
                else if (  ((hind[i] & 1) == 1)
                          && ((i == 1) || (hind[i - 1] >= 2*(j + 2) + 1))  )
                {
                    x = chain + i;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;
                    mpoly_monomial_add_mp(exp_list[exp_next], Bexps + N*x->i,
                                                             Qexps + N*x->j, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
            }
        }

        if (_n_fq_is_zero(Qcoeffs + d*Qlen, d))
            continue;

        /* try to divide accumulated term by leading term */
        if (!lt_divides)
        {
            _fq_nmod_mpoly_fit_length(&Rcoeffs, &R->coeffs_alloc, d,
                                      &Rexps, &R->exps_alloc, N, Rlen + 1);
            _n_fq_neg(Rcoeffs + d*Rlen, Qcoeffs + d*Qlen, d, fqctx->mod);
            mpoly_monomial_set(Rexps + N*Rlen, exp, N);
            Rlen++;
            continue;
        }

        n_fq_mul(Qcoeffs + d*Qlen, Qcoeffs + d*Qlen, lc_minus_inv, fqctx);

        /* put newly generated quotient term back into the heap if neccesary */
        if (s > 1)
        {
            i = 1;
            x = chain + i;
            x->i = i;
            x->j = Qlen;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;
            mpoly_monomial_add_mp(exp_list[exp_next], Bexps + N*x->i,
                                                     Qexps + N*x->j, N);
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

void fq_nmod_mpoly_divrem_monagan_pearce(
    fq_nmod_mpoly_t Q,
    fq_nmod_mpoly_t R,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong N;
    flint_bitcnt_t QRbits;
    ulong * Aexps = A->exps, * Bexps = B->exps;
    ulong * cmpmask;
    int freeAexps = 0, freeBexps = 0;
    fq_nmod_mpoly_t TQ, TR;
    fq_nmod_mpoly_struct * q, * r;

    if (fq_nmod_mpoly_is_zero(B, ctx))
    {
        flint_throw(FLINT_DIVZERO, "fq_nmod_mpoly_divrem_monagan_pearce: divide by zero");
    }

    if (fq_nmod_mpoly_is_zero(A, ctx))
    {
        fq_nmod_mpoly_zero(Q, ctx);
        fq_nmod_mpoly_zero(R, ctx);
        return;
    }

    fq_nmod_mpoly_init(TQ, ctx);
    fq_nmod_mpoly_init(TR, ctx);

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
        fq_nmod_mpoly_set(R, A, ctx);
        fq_nmod_mpoly_zero(Q, ctx);
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
        fq_nmod_mpoly_fit_length_reset_bits(q, A->length/B->length + 1, QRbits, ctx);
        fq_nmod_mpoly_fit_length_reset_bits(r, B->length, QRbits, ctx);

        if (_fq_nmod_mpoly_divrem_monagan_pearce(q, r,
                     A->coeffs, Aexps, A->length, B->coeffs, Bexps, B->length,
                                               QRbits, N, cmpmask, ctx->fqctx))
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
        fq_nmod_mpoly_swap(Q, TQ, ctx);

    if (R == A || R == B)
        fq_nmod_mpoly_swap(R, TR, ctx);

cleanup:

    fq_nmod_mpoly_clear(TQ, ctx);
    fq_nmod_mpoly_clear(TR, ctx);

    if (freeAexps)
        flint_free(Aexps);

    if (freeBexps)
        flint_free(Bexps);

    flint_free(cmpmask);
}
