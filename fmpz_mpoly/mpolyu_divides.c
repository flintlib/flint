/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

/* A = D - B*C, D may be modified if saveD == 0 */
slong _fmpz_mpoly_mulsub1(
                fmpz ** A_coeff, ulong ** A_exp, slong * A_alloc,
                 const fmpz * Dcoeff, const ulong * Dexp, slong Dlen, int saveD,
                 const fmpz * Bcoeff, const ulong * Bexp, slong Blen,
                 const fmpz * Ccoeff, const ulong * Cexp, slong Clen,
                              ulong maskhi)
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
    slong Aalloc = *A_alloc;
    fmpz * Acoeff = *A_coeff;
    ulong * Aexp = *A_exp;
    ulong exp;
    slong * hind;
    int small;
    ulong acc_sm[3], pp0, pp1;
    TMP_INIT;

    FLINT_ASSERT(Blen > 0);
    FLINT_ASSERT(Clen > 0);

    TMP_START;

    small =   _fmpz_mpoly_fits_small(Bcoeff, Blen)
           && _fmpz_mpoly_fits_small(Ccoeff, Clen);

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
            _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, 1);
            Aexp[Alen] = Dexp[Di];
            if (saveD)
                fmpz_set(Acoeff + Alen, Dcoeff + Di);
            else
                fmpz_swap(Acoeff + Alen, (fmpz *)(Dcoeff + Di));
            Alen++;
            Di++;
        }

        _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, 1);

        Aexp[Alen] = exp;

        if (small)
        {
            acc_sm[0] = acc_sm[1] = acc_sm[2] = 0;

            do
            {
                x = _mpoly_heap_pop1(heap, &heap_len, maskhi);

                hind[x->i] |= WORD(1);
                *store++ = x->i;
                *store++ = x->j;
                FLINT_ASSERT(!COEFF_IS_MPZ(Bcoeff[x->i]));
                FLINT_ASSERT(!COEFF_IS_MPZ(Ccoeff[x->j]));
                smul_ppmm(pp1, pp0, Bcoeff[x->i], Ccoeff[x->j]);
                sub_dddmmmsss(acc_sm[2], acc_sm[1], acc_sm[0],
                              acc_sm[2], acc_sm[1], acc_sm[0],
                              FLINT_SIGN_EXT(pp1), pp1, pp0);

                while ((x = x->next) != NULL)
                {
                    hind[x->i] |= WORD(1);
                    *store++ = x->i;
                    *store++ = x->j;
                    FLINT_ASSERT(!COEFF_IS_MPZ(Bcoeff[x->i]));
                    FLINT_ASSERT(!COEFF_IS_MPZ(Ccoeff[x->j]));
                    smul_ppmm(pp1, pp0, Bcoeff[x->i], Ccoeff[x->j]);
                    sub_dddmmmsss(acc_sm[2], acc_sm[1], acc_sm[0],
                                  acc_sm[2], acc_sm[1], acc_sm[0],
                                  FLINT_SIGN_EXT(pp1), pp1, pp0);
                }
            } while (heap_len > 1 && heap[1].exp == exp);

            fmpz_set_signed_uiuiui(Acoeff + Alen,
                                             acc_sm[2], acc_sm[1], acc_sm[0]);

            if (Di < Dlen && Dexp[Di] == exp)
            {
                fmpz_add(Acoeff + Alen, Acoeff + Alen, Dcoeff + Di);
                Di++;
            }
        }
        else
        {
            if (Di < Dlen && Dexp[Di] == exp)
            {
                fmpz_set(Acoeff + Alen, Dcoeff + Di);
                Di++;
            }
            else
            {
                fmpz_zero(Acoeff + Alen);
            }

            do
            {
                x = _mpoly_heap_pop1(heap, &heap_len, maskhi);

                hind[x->i] |= WORD(1);
                *store++ = x->i;
                *store++ = x->j;
                fmpz_submul(Acoeff + Alen, Bcoeff + x->i, Ccoeff + x->j);

                while ((x = x->next) != NULL)
                {
                    hind[x->i] |= WORD(1);
                    *store++ = x->i;
                    *store++ = x->j;
                    fmpz_submul(Acoeff + Alen, Bcoeff + x->i, Ccoeff + x->j);
                }
            } while (heap_len > 1 && heap[1].exp == exp);
        }

        Alen += !fmpz_is_zero(Acoeff + Alen);

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
    _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + Dlen - Di, 1);
    if (saveD)
        _fmpz_vec_set(Acoeff + Alen, Dcoeff + Di, Dlen - Di);
    else
        _fmpz_vec_swap(Acoeff + Alen, (fmpz *)(Dcoeff + Di), Dlen - Di);
    mpoly_copy_monomials(Aexp + 1*Alen, Dexp + 1*Di, Dlen - Di, 1);
    Alen += Dlen - Di;

    *A_coeff = Acoeff;
    *A_exp = Aexp;
    *A_alloc = Aalloc;

    TMP_END;

    return Alen;
}

/* A = D - B*C, D may be modified if saveD == 0 */
slong _fmpz_mpoly_mulsub(
                fmpz ** A_coeff, ulong ** A_exp, slong * A_alloc,
                 const fmpz * Dcoeff, const ulong * Dexp, slong Dlen, int saveD,
                 const fmpz * Bcoeff, const ulong * Bexp, slong Blen,
                 const fmpz * Ccoeff, const ulong * Cexp, slong Clen,
                              flint_bitcnt_t bits, slong N, const ulong * cmpmask)
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
    slong Aalloc = *A_alloc;
    fmpz * Acoeff = *A_coeff;
    ulong * Aexp = *A_exp;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    slong * hind;
    int small;
    ulong acc_sm[3], pp0, pp1;
    TMP_INIT;

    FLINT_ASSERT(Blen > 0);
    FLINT_ASSERT(Clen > 0);

    if (N == 1)
    {
        return _fmpz_mpoly_mulsub1(A_coeff, A_exp, A_alloc,
                                    Dcoeff, Dexp, Dlen, saveD,
                                    Bcoeff, Bexp, Blen,
                                    Ccoeff, Cexp, Clen, cmpmask[0]);
    }

    TMP_START;

    /* whether input coeffs are small, thus accumulation fit in three words */
    small =   _fmpz_mpoly_fits_small(Bcoeff, Blen)
           && _fmpz_mpoly_fits_small(Ccoeff, Clen);

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
            _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, N);
            mpoly_monomial_set(Aexp + N*Alen, Dexp + N*Di, N);
            if (saveD)
                fmpz_set(Acoeff + Alen, Dcoeff + Di);
            else
                fmpz_swap(Acoeff + Alen, (fmpz *)(Dcoeff + Di));
            Alen++;
            Di++;
        }

        _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, N);

        mpoly_monomial_set(Aexp + N*Alen, exp, N);

        if (small)
        {
            acc_sm[0] = acc_sm[1] = acc_sm[2] = 0;

            do
            {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);

                hind[x->i] |= WORD(1);
                *store++ = x->i;
                *store++ = x->j;
                FLINT_ASSERT(!COEFF_IS_MPZ(Bcoeff[x->i]));
                FLINT_ASSERT(!COEFF_IS_MPZ(Ccoeff[x->j]));
                smul_ppmm(pp1, pp0, Bcoeff[x->i], Ccoeff[x->j]);
                sub_dddmmmsss(acc_sm[2], acc_sm[1], acc_sm[0],
                              acc_sm[2], acc_sm[1], acc_sm[0],
                              FLINT_SIGN_EXT(pp1), pp1, pp0);

                while ((x = x->next) != NULL)
                {
                    hind[x->i] |= WORD(1);
                    *store++ = x->i;
                    *store++ = x->j;
                    FLINT_ASSERT(!COEFF_IS_MPZ(Bcoeff[x->i]));
                    FLINT_ASSERT(!COEFF_IS_MPZ(Ccoeff[x->j]));
                    smul_ppmm(pp1, pp0, Bcoeff[x->i], Ccoeff[x->j]);
                    sub_dddmmmsss(acc_sm[2], acc_sm[1], acc_sm[0],
                                  acc_sm[2], acc_sm[1], acc_sm[0],
                                  FLINT_SIGN_EXT(pp1), pp1, pp0);
                }
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

            fmpz_set_signed_uiuiui(Acoeff + Alen,
                                             acc_sm[2], acc_sm[1], acc_sm[0]);

            if (Di < Dlen && mpoly_monomial_equal(Dexp + N*Di, exp, N))
            {
                fmpz_add(Acoeff + Alen, Acoeff + Alen, Dcoeff + Di);
                Di++;
            }
        }
        else
        {
            if (Di < Dlen && mpoly_monomial_equal(Dexp + N*Di, exp, N))
            {
                fmpz_set(Acoeff + Alen, Dcoeff + Di);
                Di++;
            }
            else
            {
                fmpz_zero(Acoeff + Alen);
            }

            do
            {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);

                hind[x->i] |= WORD(1);
                *store++ = x->i;
                *store++ = x->j;
                fmpz_submul(Acoeff + Alen, Bcoeff + x->i, Ccoeff + x->j);

                while ((x = x->next) != NULL)
                {
                    hind[x->i] |= WORD(1);
                    *store++ = x->i;
                    *store++ = x->j;
                    fmpz_submul(Acoeff + Alen, Bcoeff + x->i, Ccoeff + x->j);
                }
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));
        }

        Alen += !fmpz_is_zero(Acoeff + Alen);

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
    _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + Dlen - Di, N);
    if (saveD)
        _fmpz_vec_set(Acoeff + Alen, Dcoeff + Di, Dlen - Di);
    else
        _fmpz_vec_swap(Acoeff + Alen, (fmpz *)(Dcoeff + Di), Dlen - Di);
    mpoly_copy_monomials(Aexp + N*Alen, Dexp + N*Di, Dlen - Di, N);
    Alen += Dlen - Di;

    *A_coeff = Acoeff;
    *A_exp = Aexp;
    *A_alloc = Aalloc;

    TMP_END;

    return Alen;
}


/*
    Q = A/B return 1 if division was exact.
    nmainvars is the number of main vars.
*/
int fmpz_mpolyuu_divides(
    fmpz_mpolyu_t Q,
    const fmpz_mpolyu_t A,
    const fmpz_mpolyu_t B,
    slong nmainvars,
    const fmpz_mpoly_ctx_t ctx)
{
    flint_bitcnt_t bits = A->bits;
    fmpz_mpoly_struct * Bcoeff = B->coeffs;
    slong Blen = B->length;
    ulong * Bexp = B->exps;
    fmpz_mpoly_struct * Acoeff = A->coeffs;
    slong Alen = A->length;
    ulong * Aexp = A->exps;
    fmpz_mpoly_struct * a, * b, * q;
    slong N;
    ulong * cmpmask;    /* cmp mask for lesser variables */
    fmpz_mpoly_t T, S;
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

    fmpz_mpoly_init3(T, 16, bits, ctx);
    fmpz_mpoly_init3(S, 16, bits, ctx);

    while (heap_len > 1)
    {
        exp = heap[1].exp;

        if (mpoly_monomial_overflows1(exp, mask))
        {
            goto not_exact_division;
        }

        fmpz_mpolyu_fit_length(Q, Q->length + 1, ctx);
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
                    fmpz_mpoly_fit_length(S, T->length + a->length, ctx);
                    S->length = _fmpz_mpoly_add(
                                    S->coeffs, S->exps,
                                    T->coeffs, T->exps, T->length,
                                    a->coeffs, a->exps, a->length,
                                                              N, cmpmask);
                }
                else
                {
                    b = Bcoeff + x->i;
                    q = Q->coeffs + x->j;
                    S->length = _fmpz_mpoly_mulsub(
                                    &S->coeffs, &S->exps, &S->alloc,
                                    T->coeffs, T->exps, T->length, 0,
                                    b->coeffs, b->exps, b->length,
                                    q->coeffs, q->exps, q->length,
                                                         bits, N, cmpmask);
                }
                fmpz_mpoly_swap(S, T, ctx);

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
        q->length = _fmpz_mpoly_divides_monagan_pearce(
                            &q->coeffs, &q->exps, &q->alloc,
                            T->coeffs, T->exps, T->length,
                            b->coeffs, b->exps, b->length,
                                              bits, N, cmpmask);
        if (q->length == 0)
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

    fmpz_mpoly_clear(T, ctx);
    fmpz_mpoly_clear(S, ctx);

    TMP_END;

    return success;

not_exact_division:

    success = 0;
    Q->length = 0;
    goto cleanup;
}

