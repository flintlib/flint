/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

static int _nmod_mpoly_divides_monagan_pearce1(
    nmod_mpoly_t Q,
    const mp_limb_t * coeff2, const ulong * exp2, slong len2,
    const mp_limb_t * coeff3, const ulong * exp3, slong len3,
    slong bits,
    ulong maskhi,
    nmod_t fctx)
{
    int lt_divides;
    slong i, j, q_len, s;
    slong next_loc, heap_len;
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    mp_limb_t * q_coeff = Q->coeffs;
    ulong * q_exp = Q->exps;
    slong * hind;
    ulong mask, exp, maxexp = exp2[len2 - 1];
    mp_limb_t lc_minus_inv, acc0, acc1, acc2, pp1, pp0;
    TMP_INIT;

    TMP_START;

    /* alloc array of heap nodes which can be chained together */
    next_loc = len3 + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap1_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap1_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(len3*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*len3*sizeof(slong));

    /* space for flagged heap indicies */
    hind = (slong *) TMP_ALLOC(len3*sizeof(slong));
    for (i = 0; i < len3; i++)
        hind[i] = 1;

    mask = mpoly_overflow_mask_sp(bits);

    q_len = WORD(0);

    /* s is the number of terms * (latest quotient) we should put into heap */
    s = len3;

    /* insert (-1, 0, exp2[0]) into heap */
    heap_len = 2;
    x = chain + 0;
    x->i = -WORD(1);
    x->j = 0;
    x->next = NULL;
    HEAP_ASSIGN(heap[1], exp2[0], x);

    /* precompute leading cofficient info */
    lc_minus_inv = fctx.n - nmod_inv(coeff3[0], fctx);

    while (heap_len > 1)
    {
        exp = heap[1].exp;

        if (mpoly_monomial_overflows1(exp, mask))
            goto not_exact_division;

        _nmod_mpoly_fit_length(&q_coeff, &Q->coeffs_alloc,
                               &q_exp, &Q->exps_alloc, 1, q_len + 1);

        lt_divides = mpoly_monomial_divides1(q_exp + q_len, exp, exp3[0], mask);

        acc0 = acc1 = acc2 = 0;
        do {
            x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
            do {
                *store++ = x->i;
                *store++ = x->j;

                if (x->i == -WORD(1))
                {
                    add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0,
                                 WORD(0), WORD(0), fctx.n - coeff2[x->j]);
                }
                else
                {
                    hind[x->i] |= WORD(1);
                    umul_ppmm(pp1, pp0, coeff3[x->i], q_coeff[x->j]);
                    add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);
                }
            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && heap[1].exp == exp);

        NMOD_RED3(q_coeff[q_len], acc2, acc1, acc0, fctx);

        /* process nodes taken from the heap */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            if (i == -WORD(1))
            {
                /* take next dividend term */
                if (j + 1 < len2)
                {
                    x = chain + 0;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    _mpoly_heap_insert1(heap, exp2[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
                }

            } else
            {
                /* should we go right? */
                if (  (i + 1 < len3)
                   && (hind[i + 1] == 2*j + 1)
                   )
                {
                    x = chain + i + 1;
                    x->i = i + 1;
                    x->j = j;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;
                    _mpoly_heap_insert1(heap, exp3[x->i] + q_exp[x->j], x,
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
                    _mpoly_heap_insert1(heap, exp3[x->i] + q_exp[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
                }
            }
        }

        q_coeff[q_len] = nmod_mul(q_coeff[q_len], lc_minus_inv, fctx);

        if (q_coeff[q_len] == 0)
        {
            continue;
        }

        if (!lt_divides || (exp^maskhi) < (maxexp^maskhi))
            goto not_exact_division;

        /* put newly generated quotient term back into the heap if neccesary */
        if (s > 1)
        {
            i = 1;
            x = chain + i;
            x->i = i;
            x->j = q_len;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;
            _mpoly_heap_insert1(heap, exp3[x->i] + q_exp[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
        }
        s = 1;
        q_len++;
    }

    Q->coeffs = q_coeff;
    Q->exps = q_exp;
    Q->length = q_len;

    TMP_END;

    return 1;

not_exact_division:

    Q->coeffs = q_coeff;
    Q->exps = q_exp;
    Q->length = 0;

    TMP_END;

    return 0;
}


int _nmod_mpoly_divides_monagan_pearce(
    nmod_mpoly_t Q,
    const mp_limb_t * coeff2, const ulong * exp2, slong len2,
    const mp_limb_t * coeff3, const ulong * exp3, slong len3,
    flint_bitcnt_t bits,
    slong N,
    const ulong * cmpmask,
    nmod_t fctx)
{
    int lt_divides;
    slong i, j, q_len, s;
    slong next_loc, heap_len;
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    mp_limb_t * q_coeff = Q->coeffs;
    ulong * q_exp = Q->exps;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    mp_limb_t lc_minus_inv, acc0, acc1, acc2, pp1, pp0;
    ulong mask;
    slong * hind;
    TMP_INIT;

    if (N == 1)
        return _nmod_mpoly_divides_monagan_pearce1(Q, coeff2, exp2, len2,
                                   coeff3, exp3, len3, bits, cmpmask[0], fctx);

    TMP_START;

    /* alloc array of heap nodes which can be chained together */
    next_loc = len3 + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(len3*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*len3*sizeof(slong));

    /* array of exponent vectors, each of "N" words */
    exps = (ulong *) TMP_ALLOC(len3*N*sizeof(ulong));
    /* list of pointers to available exponent vectors */
    exp_list = (ulong **) TMP_ALLOC(len3*sizeof(ulong *));
    /* space to save copy of current exponent vector */
    exp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    /* set up list of available exponent vectors */
    exp_next = 0;
    for (i = 0; i < len3; i++)
        exp_list[i] = exps + i*N;

    /* space for flagged heap indicies */
    hind = (slong *) TMP_ALLOC(len3*sizeof(slong));
    for (i = 0; i < len3; i++)
        hind[i] = 1;

    mask = bits <= FLINT_BITS ? mpoly_overflow_mask_sp(bits) : 0;

    q_len = WORD(0);

    /* s is the number of terms * (latest quotient) we should put into heap */
    s = len3;

    /* insert (-1, 0, exp2[0]) into heap */
    heap_len = 2;
    x = chain + 0;
    x->i = -WORD(1);
    x->j = 0;
    x->next = NULL;
    heap[1].next = x;
    heap[1].exp = exp_list[exp_next++];
    mpoly_monomial_set(heap[1].exp, exp2, N);

    /* precompute leading cofficient info */
    lc_minus_inv = fctx.n - nmod_inv(coeff3[0], fctx);

    while (heap_len > 1)
    {
        _nmod_mpoly_fit_length(&q_coeff, &Q->coeffs_alloc,
                               &q_exp, &Q->exps_alloc, N, q_len + 1);

        mpoly_monomial_set(exp, heap[1].exp, N);

        if (bits <= FLINT_BITS)
        {
            if (mpoly_monomial_overflows(exp, N, mask))
                goto not_exact_division;
            lt_divides = mpoly_monomial_divides(q_exp + q_len*N, exp, exp3, N, mask);
        }
        else
        {
            if (mpoly_monomial_overflows_mp(exp, N, bits))
                goto not_exact_division;
            lt_divides = mpoly_monomial_divides_mp(q_exp + q_len*N, exp, exp3, N, bits);
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
                                 WORD(0), WORD(0), fctx.n - coeff2[x->j]);
                }
                else
                {
                    hind[x->i] |= WORD(1);
                    umul_ppmm(pp1, pp0, coeff3[x->i], q_coeff[x->j]);
                    add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);                    
                }
            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

        NMOD_RED3(q_coeff[q_len], acc2, acc1, acc0, fctx);

        /* process nodes taken from the heap */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            if (i == -WORD(1))
            {
                /* take next dividend term */
                if (j + 1 < len2)
                {
                    x = chain + 0;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    mpoly_monomial_set(exp_list[exp_next], exp2 + x->j*N, N);
                    if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                      &next_loc, &heap_len, N, cmpmask))
                        exp_next--;
                }
            } else
            {
                /* should we go up */
                if (  (i + 1 < len3)
                   && (hind[i + 1] == 2*j + 1)
                   )
                {
                    x = chain + i + 1;
                    x->i = i + 1;
                    x->j = j;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;

                    if (bits <= FLINT_BITS)
                        mpoly_monomial_add(exp_list[exp_next], exp3 + x->i*N,
                                                              q_exp + x->j*N, N);
                    else
                        mpoly_monomial_add_mp(exp_list[exp_next], exp3 + x->i*N,
                                                                 q_exp + x->j*N, N);

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
                        mpoly_monomial_add(exp_list[exp_next], exp3 + x->i*N,
                                                              q_exp + x->j*N, N);
                    else
                        mpoly_monomial_add_mp(exp_list[exp_next], exp3 + x->i*N,
                                                                 q_exp + x->j*N, N);

                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
            }
        }

        q_coeff[q_len] = nmod_mul(q_coeff[q_len], lc_minus_inv, fctx);

        if (q_coeff[q_len] == 0)
        {
            continue;
        }

        if (!lt_divides ||
            mpoly_monomial_gt(exp2 + N*(len2 - 1), exp, N, cmpmask))
        {
            goto not_exact_division;
        }

        if (s > 1)
        {
            i = 1;
            x = chain + i;
            x->i = i;
            x->j = q_len;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;

            if (bits <= FLINT_BITS)
                mpoly_monomial_add(exp_list[exp_next], exp3 + x->i*N, q_exp + x->j*N, N);
            else
                mpoly_monomial_add_mp(exp_list[exp_next], exp3 + x->i*N, q_exp + x->j*N, N);

            exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
        }
        s = 1;      
        q_len++;
    }

    Q->coeffs = q_coeff;
    Q->exps = q_exp;
    Q->length = q_len;

    TMP_END;

    return 1;

not_exact_division:

    Q->coeffs = q_coeff;
    Q->exps = q_exp;
    Q->length = 0;

    TMP_END;

    return 0;
}

/* return 1 if quotient is exact */
int nmod_mpoly_divides_monagan_pearce(
    nmod_mpoly_t Q,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, N;
    flint_bitcnt_t Qbits;
    fmpz * Amaxfields, * Bmaxfields;
    ulong * cmpmask;
    ulong * exp2 = A->exps, * exp3 = B->exps, * expq;
    int divides, easy_exit, free2 = 0, free3 = 0;
    ulong mask = 0;
    TMP_INIT;

    if (B->length == 0)
    {
        if (A->length == 0 || nmod_mpoly_ctx_modulus(ctx) == 1)
        {
            nmod_mpoly_set(Q, A, ctx);
            return 1;
        }
        else
        {
            flint_throw(FLINT_DIVZERO, "nmod_mpoly_divides_monagan_pearce: divide by zero");
        }
    }

    if (A->length == 0)
    {
        nmod_mpoly_zero(Q, ctx);
        return 1;
    }

    TMP_START;

    Amaxfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    Bmaxfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_init(Amaxfields + i);
        fmpz_init(Bmaxfields + i);
    }

    mpoly_max_fields_fmpz(Amaxfields, A->exps, A->length,
                                                      A->bits, ctx->minfo);
    mpoly_max_fields_fmpz(Bmaxfields, B->exps, B->length,
                                                      B->bits, ctx->minfo);
    easy_exit = 0;
    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        /*
            cannot be exact division if any max field from A
            is less than corresponding max field from B
        */
        if (fmpz_cmp(Amaxfields + i, Bmaxfields + i) < 0)
            easy_exit = 1;
    }

    Qbits = 1 + _fmpz_vec_max_bits(Amaxfields, ctx->minfo->nfields);
    Qbits = FLINT_MAX(Qbits, A->bits);
    Qbits = FLINT_MAX(Qbits, B->bits);
    Qbits = mpoly_fix_bits(Qbits, ctx->minfo);

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_clear(Amaxfields + i);
        fmpz_clear(Bmaxfields + i);
    }

    if (easy_exit)
    {
        nmod_mpoly_zero(Q, ctx);
        divides = 0;
        goto cleanup;
    }

    N = mpoly_words_per_exp(Qbits, ctx->minfo);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, Qbits, ctx->minfo);

    /* temporary space to check leading monomials divide */
    expq = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    /* quick check for easy case of inexact division of leading monomials */
    if (Qbits == A->bits && Qbits == B->bits && A->exps[N - 1] < B->exps[N - 1])
    {
        nmod_mpoly_zero(Q, ctx);
        divides = 0;
        goto cleanup;
    }

    /* ensure input exponents packed to same size as output exponents */
    if (Qbits > A->bits)
    {
        free2 = 1;
        exp2 = (ulong *) flint_malloc(N*A->length*sizeof(ulong));
        mpoly_repack_monomials(exp2, Qbits, A->exps, A->bits,
                                                    A->length, ctx->minfo);
    }

    if (Qbits > B->bits)
    {
        free3 = 1;
        exp3 = (ulong *) flint_malloc(N*B->length*sizeof(ulong));
        mpoly_repack_monomials(exp3, Qbits, B->exps, B->bits,
                                                    B->length, ctx->minfo);
    }

    /* check leading monomial divides exactly */
    if (Qbits <= FLINT_BITS)
    {
        /* mask with high bit of each exponent vector field set */
        for (i = 0; i < FLINT_BITS/Qbits; i++)
            mask = (mask << Qbits) + (UWORD(1) << (Qbits - 1));

        if (!mpoly_monomial_divides(expq, exp2, exp3, N, mask))
        {
            nmod_mpoly_zero(Q, ctx);
            divides = 0;
            goto cleanup;
        }
    }
    else
    {
        if (!mpoly_monomial_divides_mp(expq, exp2, exp3, N, Qbits))
        {
            nmod_mpoly_zero(Q, ctx);
            divides = 0;
            goto cleanup;
        }
    }

    /* deal with aliasing and divide polynomials */
    if (Q == A || Q == B)
    {
        nmod_mpoly_t temp;
        nmod_mpoly_init3(temp, A->length/B->length + 1, Qbits, ctx);
        divides = _nmod_mpoly_divides_monagan_pearce(temp,
                                      A->coeffs, exp2, A->length,
                                      B->coeffs, exp3, B->length,
                                            Qbits, N, cmpmask, ctx->mod);
        nmod_mpoly_swap(temp, Q, ctx);
        nmod_mpoly_clear(temp, ctx);
    }
    else
    {
        nmod_mpoly_fit_length_reset_bits(Q, A->length/B->length + 1, Qbits, ctx);

        divides = _nmod_mpoly_divides_monagan_pearce(Q,
                                    A->coeffs, exp2, A->length,
                                    B->coeffs, exp3, B->length,
                                            Qbits, N, cmpmask, ctx->mod);
    }

cleanup:

    if (free2)
        flint_free(exp2);

    if (free3)
        flint_free(exp3);

    TMP_END;

    return divides;
}
