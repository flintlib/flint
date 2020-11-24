/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

/* A = D - B*C */
slong _nmod_mpoly_mulsub1(nmod_mpoly_t A,
                 const mp_limb_t * Dcoeff, const ulong * Dexp, slong Dlen,
                 const mp_limb_t * Bcoeff, const ulong * Bexp, slong Blen,
                 const mp_limb_t * Ccoeff, const ulong * Cexp, slong Clen,
                                         ulong maskhi, nmod_t fctx)
{
    slong i, j;
    slong next_loc;
    slong heap_len = 2; /* heap zero index unused */
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    slong Di;
    slong Alen;
    mp_limb_t * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    ulong exp;
    slong * hind;
    mp_limb_t acc0, acc1, acc2, pp1, pp0;
    TMP_INIT;

    FLINT_ASSERT(Blen > 0);
    FLINT_ASSERT(Clen > 0);

    TMP_START;

    next_loc = Blen + 4; /* something bigger than heap can ever be */
    heap = (mpoly_heap1_s *) TMP_ALLOC((Blen + 1)*sizeof(mpoly_heap1_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(Blen*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*Blen*sizeof(slong));

    /* space for heap indices */
    hind = (slong *) TMP_ALLOC(Blen*sizeof(slong));
    for (i = 0; i < Blen; i++)
        hind[i] = 1;

    /* put (0, 0, Bexp[0] + Cexp[0]) on heap */
    x = chain + 0;
    x->i = 0;
    x->j = 0;
    x->next = NULL;
    HEAP_ASSIGN(heap[1], Bexp[0] + Cexp[0], x);
    hind[0] = 2*1 + 0;

    Alen = 0;
    Di = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        while (Di < Dlen && mpoly_monomial_gt1(Dexp[Di], exp, maskhi))
        {
            _nmod_mpoly_fit_length(&Acoeff, &A->coeffs_alloc,
                                   &Aexp, &A->exps_alloc, 1, Alen + 1);
            Aexp[Alen] = Dexp[Di];
            Acoeff[Alen] = Dcoeff[Di];
            Alen++;
            Di++;
        }

        _nmod_mpoly_fit_length(&Acoeff, &A->coeffs_alloc,
                               &Aexp, &A->exps_alloc, 1, Alen + 1);

        Aexp[Alen] = exp;

        acc0 = acc1 = acc2 = 0;
        do
        {
            x = _mpoly_heap_pop1(heap, &heap_len, maskhi);

            hind[x->i] |= WORD(1);
            *store++ = x->i;
            *store++ = x->j;
            umul_ppmm(pp1, pp0, Bcoeff[x->i], Ccoeff[x->j]);
            add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, 0, pp1, pp0);

            while ((x = x->next) != NULL)
            {
                hind[x->i] |= WORD(1);
                *store++ = x->i;
                *store++ = x->j;
                umul_ppmm(pp1, pp0, Bcoeff[x->i], Ccoeff[x->j]);
                add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, 0, pp1, pp0);
            }
        } while (heap_len > 1 && heap[1].exp == exp);

        NMOD_RED3(pp0, acc2, acc1, acc0, fctx);

        pp1 = 0;
        if (Di < Dlen && Dexp[Di] == exp)
        {
            pp1 = Dcoeff[Di];
            Di++;
        }

        Acoeff[Alen] = nmod_sub(pp1, pp0, fctx);

        Alen += Acoeff[Alen] != 0;

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

                hind[x->i] = 2*(x->j + 1) + 0;
                _mpoly_heap_insert1(heap, Bexp[x->i] + Cexp[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
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

                hind[x->i] = 2*(x->j + 1) + 0;
                _mpoly_heap_insert1(heap, Bexp[x->i] + Cexp[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
            }
        }
    }

    FLINT_ASSERT(Di <= Dlen);
    if (Di < Dlen)
    {
        _nmod_mpoly_fit_length(&Acoeff, &A->coeffs_alloc,
                               &Aexp, &A->exps_alloc, 1, Alen + Dlen - Di);
        _nmod_vec_set(Acoeff + Alen, Dcoeff + Di, Dlen - Di);
        mpoly_copy_monomials(Aexp + 1*Alen, Dexp + 1*Di, Dlen - Di, 1);
        Alen += Dlen - Di;
    }

    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->length = Alen;

    TMP_END;

    return Alen;
}

/* A = D - B*C */
void _nmod_mpoly_mulsub(nmod_mpoly_t A,
                 const mp_limb_t * Dcoeff, const ulong * Dexp, slong Dlen,
                 const mp_limb_t * Bcoeff, const ulong * Bexp, slong Blen,
                 const mp_limb_t * Ccoeff, const ulong * Cexp, slong Clen,
              flint_bitcnt_t bits, slong N, const ulong * cmpmask, nmod_t fctx)
{
    slong i, j;
    slong next_loc;
    slong heap_len = 2; /* heap zero index unused */
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    slong Di;
    slong Alen;
    mp_limb_t * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    slong * hind;
    mp_limb_t acc0, acc1, acc2, pp1, pp0;
    TMP_INIT;

    FLINT_ASSERT(Blen > 0);
    FLINT_ASSERT(Clen > 0);

    if (N == 1)
    {
        _nmod_mpoly_mulsub1(A, Dcoeff, Dexp, Dlen,
                     Bcoeff, Bexp, Blen, Ccoeff, Cexp, Clen, cmpmask[0], fctx);
        return;
    }

    TMP_START;

    next_loc = Blen + 4; /* something bigger than heap can ever be */
    heap = (mpoly_heap_s *) TMP_ALLOC((Blen + 1)*sizeof(mpoly_heap_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(Blen*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*Blen*sizeof(slong));
    exps = (ulong *) TMP_ALLOC(Blen*N*sizeof(ulong));
    exp_list = (ulong **) TMP_ALLOC(Blen*sizeof(ulong *));
    for (i = 0; i < Blen; i++)
        exp_list[i] = exps + i*N;

    /* space for heap indices */
    hind = (slong *) TMP_ALLOC(Blen*sizeof(slong));
    for (i = 0; i < Blen; i++)
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

    if (bits <= FLINT_BITS)
        mpoly_monomial_add(heap[1].exp, Bexp + N*0, Cexp + N*0, N);
    else
        mpoly_monomial_add_mp(heap[1].exp, Bexp + N*0, Cexp + N*0, N);

    hind[0] = 2*1 + 0;

    Alen = 0;
    Di = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        while (Di < Dlen && mpoly_monomial_gt(Dexp + N*Di, exp, N, cmpmask))
        {
            _nmod_mpoly_fit_length(&Acoeff, &A->coeffs_alloc,
                                   &Aexp, &A->exps_alloc, N, Alen + 1);
            mpoly_monomial_set(Aexp + N*Alen, Dexp + N*Di, N);
            Acoeff[Alen] = Dcoeff[Di];
            Alen++;
            Di++;
        }

        _nmod_mpoly_fit_length(&Acoeff, &A->coeffs_alloc,
                               &Aexp, &A->exps_alloc, N, Alen + 1);

        mpoly_monomial_set(Aexp + N*Alen, exp, N);

        acc0 = acc1 = acc2 = 0;
        do
        {
            exp_list[--exp_next] = heap[1].exp;
            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);

            hind[x->i] |= WORD(1);
            *store++ = x->i;
            *store++ = x->j;
            umul_ppmm(pp1, pp0, Bcoeff[x->i], Ccoeff[x->j]);
            add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, 0, pp1, pp0);

            while ((x = x->next) != NULL)
            {
                hind[x->i] |= WORD(1);
                *store++ = x->i;
                *store++ = x->j;
                umul_ppmm(pp1, pp0, Bcoeff[x->i], Ccoeff[x->j]);
                add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, 0, pp1, pp0);
            }
        } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

        NMOD_RED3(pp0, acc2, acc1, acc0, fctx);

        pp1 = 0;
        if (Di < Dlen && mpoly_monomial_equal(Dexp + N*Di, exp, N))
        {
            pp1 = Dcoeff[Di];
            Di++;
        }

        Acoeff[Alen] = nmod_sub(pp1, pp0, fctx);

        Alen += Acoeff[Alen] != 0;

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
    if (Di < Dlen)
    {
        _nmod_mpoly_fit_length(&Acoeff, &A->coeffs_alloc,
                               &Aexp, &A->exps_alloc, N, Alen + Dlen - Di);
        _nmod_vec_set(Acoeff + Alen, Dcoeff + Di, Dlen - Di);
        mpoly_copy_monomials(Aexp + N*Alen, Dexp + N*Di, Dlen - Di, N);
        Alen += Dlen - Di;
    }

    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->length = Alen;

    TMP_END;
}


/*
    Q = A/B return 1 if division was exact.
    nmainvars is the number of main vars.
*/
int nmod_mpolyuu_divides(
    nmod_mpolyu_t Q,
    const nmod_mpolyu_t A,
    const nmod_mpolyu_t B,
    slong nmainvars,
    const nmod_mpoly_ctx_t ctx)
{
    flint_bitcnt_t bits = A->bits;
    nmod_mpoly_struct * Bcoeff = B->coeffs;
    slong Blen = B->length;
    ulong * Bexp = B->exps;
    nmod_mpoly_struct * Acoeff = A->coeffs;
    slong Alen = A->length;
    ulong * Aexp = A->exps;
    nmod_mpoly_struct * a, * b, * q;
    slong N;
    ulong * cmpmask;    /* cmp mask for lesser variables */
    nmod_mpoly_t T, S;
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
    TMP_INIT;

    FLINT_ASSERT(Alen > 0);
    FLINT_ASSERT(Blen > 0);
    FLINT_ASSERT(bits == A->bits);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == Q->bits);

    TMP_START;

    N = mpoly_words_per_exp(bits, ctx->minfo);
    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, bits, ctx->minfo);

    /* alloc array of heap nodes which can be chained together */
    next_loc = Blen + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap1_s *) TMP_ALLOC((Blen + 1)*sizeof(mpoly_heap1_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(Blen*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*Blen*sizeof(slong));

    /* space for flagged heap indicies */
    hind = (slong *) TMP_ALLOC(Blen*sizeof(slong));
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

    nmod_mpoly_init3(T, 16, bits, ctx);
    nmod_mpoly_init3(S, 16, bits, ctx);

    while (heap_len > 1)
    {
        exp = heap[1].exp;

        if (mpoly_monomial_overflows1(exp, mask))
        {
            goto not_exact_division;
        }

        nmod_mpolyu_fit_length(Q, Q->length + 1, ctx);
        lt_divides = mpoly_monomial_divides1(Q->exps + Q->length, exp, Bexp[0], mask);

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
                    a = Acoeff + x->j;
                    nmod_mpoly_fit_length(S, T->length + a->length, ctx);
                    S->length = _nmod_mpoly_add(
                                    S->coeffs, S->exps,
                                    T->coeffs, T->exps, T->length,
                                    a->coeffs, a->exps, a->length,
                                                      N, cmpmask, ctx->mod);
                }
                else
                {
                    b = Bcoeff + x->i;
                    q = Q->coeffs + x->j;
                    _nmod_mpoly_mulsub(S, T->coeffs, T->exps, T->length,
                                          b->coeffs, b->exps, b->length,
                                          q->coeffs, q->exps, q->length,
                                                bits, N, cmpmask, ctx->mod);
                }
                nmod_mpoly_swap(S, T, ctx);

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

        q = Q->coeffs + Q->length;
        FLINT_ASSERT(q->bits == bits);
        b = Bcoeff + 0;
        if (!_nmod_mpoly_divides_monagan_pearce(q,
                                            T->coeffs, T->exps, T->length,
                                            b->coeffs, b->exps, b->length,
                                               bits, N, cmpmask, ctx->mod))
        {
            goto not_exact_division;
        }

        if (!lt_divides || (exp^maskhi) < (maxexp^maskhi))
        {
            goto not_exact_division;
        }

        /* put newly generated quotient term back into the heap if neccesary */
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

    nmod_mpoly_clear(T, ctx);
    nmod_mpoly_clear(S, ctx);

    TMP_END;

    return success;

not_exact_division:

    success = 0;
    Q->length = 0;
    goto cleanup;
}

