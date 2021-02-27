/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

static int _fmpz_mod_mpoly_divrem_monagan_pearce1_binomial(
    fmpz_mod_mpoly_t Q,
    fmpz_mod_mpoly_t R,
    const fmpz * Acoeffs, const ulong * Aexps, slong Alen,
    const fmpz * Bcoeffs, const ulong * Bexps,
    flint_bitcnt_t bits,
    ulong cmpmask,
    const fmpz_mod_ctx_t fctx)
{
    fmpz * Qcoeffs = Q->coeffs;
    fmpz * Rcoeffs = R->coeffs;
    ulong * Qexps = Q->exps;
    ulong * Rexps = R->exps;
    ulong lexp, mask = mpoly_overflow_mask_sp(bits);
    fmpz_t lc_inv, mBcoeff1;
    slong Qlen = 0;
    slong Rlen = 0;
    slong Aidx = 0;
    slong Qidx = 0;

    fmpz_init(lc_inv);
    fmpz_init(mBcoeff1);
    
    fmpz_mod_inv(lc_inv, Bcoeffs + 0, fctx);
    fmpz_mod_neg(mBcoeff1, Bcoeffs + 1, fctx);

    while (1)
    {
        FLINT_ASSERT(0 <= Aidx && Aidx <= Alen);
        FLINT_ASSERT(0 <= Qidx && Qidx <= Qlen);

        _fmpz_mod_mpoly_fit_length(&Qcoeffs, &Q->coeffs_alloc,
                                   &Qexps, &Q->exps_alloc, 1, Qlen + 1);

        if (Aidx < Alen)
        {
            lexp = Aexps[Aidx];

            if (Qidx < Qlen)
            {
                ulong thisexp = Bexps[1] + Qexps[Qidx];
                int cmp = mpoly_monomial_cmp1(lexp, thisexp, cmpmask);
                if (cmp < 0)
                {
                    lexp = thisexp;
                    fmpz_mod_mul(Qcoeffs + Qlen, mBcoeff1, Qcoeffs + Qidx, fctx);
                    Qidx++;
                }
                else if (cmp == 0)
                {
                    fmpz_mod_mul(Qcoeffs + Qlen, mBcoeff1, Qcoeffs + Qidx, fctx);
                    fmpz_mod_add(Qcoeffs + Qlen, Qcoeffs + Qlen, Acoeffs + Aidx, fctx);
                    Aidx++;
                    Qidx++;
                }
                else
                {
                    fmpz_set(Qcoeffs + Qlen, Acoeffs + Aidx);
                    Aidx++;
                }
            }
            else
            {
                fmpz_set(Qcoeffs + Qlen, Acoeffs + Aidx);
                Aidx++;
            }
        }
        else if (Qidx < Qlen)
        {
            lexp = Bexps[1] + Qexps[Qidx];
            fmpz_mod_mul(Qcoeffs + Qlen, mBcoeff1, Qcoeffs + Qidx, fctx);
            Qidx++;
        }
        else
        {
            break;
        }

        if (mpoly_monomial_overflows1(lexp, mask))
            goto exp_overflow;

        if (fmpz_is_zero(Qcoeffs + Qlen))
            continue;

        if (!mpoly_monomial_divides1(Qexps + Qlen, lexp, Bexps[0], mask))
        {
            _fmpz_mod_mpoly_fit_length(&Rcoeffs, &R->coeffs_alloc,
                                       &Rexps, &R->exps_alloc, 1, Rlen + 1);
            fmpz_swap(Rcoeffs + Rlen, Qcoeffs + Qlen);
            Rexps[Rlen] = lexp;
            Rlen++;
            continue;
        }

        if (!fmpz_is_one(lc_inv))
            fmpz_mod_mul(Qcoeffs + Qlen, Qcoeffs + Qlen, lc_inv, fctx);
        Qlen++;
    }

    Q->coeffs = Qcoeffs;
    Q->exps = Qexps;
    Q->length = Qlen;

    R->coeffs = Rcoeffs;
    R->exps = Rexps;
    R->length = Rlen;

    fmpz_clear(lc_inv);
    fmpz_clear(mBcoeff1);

    return 1;

exp_overflow:

    Q->coeffs = Qcoeffs;
    Q->exps = Qexps;
    Q->length = 0;

    R->coeffs = Rcoeffs;
    R->exps = Rexps;
    R->length = 0;

    fmpz_clear(lc_inv);
    fmpz_clear(mBcoeff1);

    return 0;
}


static int _fmpz_mod_mpoly_divrem_monagan_pearce1(
    fmpz_mod_mpoly_t Q,
    fmpz_mod_mpoly_t R,
    const fmpz * Acoeffs, const ulong * Aexps, slong Alen,
    const fmpz * Bcoeffs, const ulong * Bexps, slong Blen,
    flint_bitcnt_t bits,
    ulong cmpmask,
    const fmpz_mod_ctx_t fctx)
{
    slong i, j, Qlen, Rlen, s;
    slong next_loc, heap_len = 2;
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    fmpz * Qcoeffs = Q->coeffs;
    fmpz * Rcoeffs = R->coeffs;
    ulong * Qexps = Q->exps;
    ulong * Rexps = R->exps;
    slong * hind;
    ulong mask, exp;
    int lt_divides;
    mpz_t t, acc, modulus;
    ulong acc_sm[3];
    fmpz_t lc_minus_inv;
    TMP_INIT;         

    mpz_init(t);
    mpz_init(acc);
    mpz_init(modulus);
    fmpz_get_mpz(modulus, fmpz_mod_ctx_modulus(fctx));

    fmpz_init(lc_minus_inv);
    fmpz_mod_inv(lc_minus_inv, Bcoeffs + 0, fctx);
    fmpz_mod_neg(lc_minus_inv, lc_minus_inv, fctx);

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

    while (heap_len > 1)
    {
        exp = heap[1].exp;

        if (mpoly_monomial_overflows1(exp, mask))
            goto exp_overflow;

        _fmpz_mod_mpoly_fit_length(&Qcoeffs, &Q->coeffs_alloc,
                                   &Qexps, &Q->exps_alloc, 1, Qlen + 1);

        lt_divides = mpoly_monomial_divides1(Qexps + Qlen, exp, Bexps[0], mask);

        mpz_set_ui(acc, 0);
        acc_sm[2] = acc_sm[1] = acc_sm[0] = 0;

        do {
            x = _mpoly_heap_pop1(heap, &heap_len, cmpmask);
            do {
                *store++ = x->i;
                *store++ = x->j;

                if (x->i == -UWORD(1))
                {
                    fmpz Aj = Acoeffs[x->j];

                    if (COEFF_IS_MPZ(Aj))
                        mpz_sub(acc, acc, COEFF_TO_PTR(Aj));
                    else
                        flint_mpz_sub_ui(acc, acc, Aj);
                }
                else
                {
                    fmpz Bi = Bcoeffs[x->i];
                    fmpz Qj = Qcoeffs[x->j];

                    hind[x->i] |= WORD(1);

                    if (COEFF_IS_MPZ(Bi) && COEFF_IS_MPZ(Qj))
                    {
                        mpz_addmul(acc, COEFF_TO_PTR(Bi), COEFF_TO_PTR(Qj));
                    }
                    else if (COEFF_IS_MPZ(Bi) && !COEFF_IS_MPZ(Qj))
                    {
                        flint_mpz_addmul_ui(acc, COEFF_TO_PTR(Bi), Qj);
                    }
                    else if (!COEFF_IS_MPZ(Bi) && COEFF_IS_MPZ(Qj))
                    {
                        flint_mpz_addmul_ui(acc, COEFF_TO_PTR(Qj), Bi);
                    }
                    else
                    {
                        ulong pp1, pp0;
                        umul_ppmm(pp1, pp0, Bi, Qj);
                        add_sssaaaaaa(acc_sm[2], acc_sm[1], acc_sm[0],
                                      acc_sm[2], acc_sm[1], acc_sm[0],
                                      0, pp1, pp0);
                    }
                }
            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && heap[1].exp == exp);

        flint_mpz_add_uiuiui(acc, acc, acc_sm[2], acc_sm[1], acc_sm[0]);
        if (mpz_sgn(acc) < 0)
            mpz_add(acc, acc, modulus);

        mpz_tdiv_qr(t, _fmpz_promote(Qcoeffs + Qlen), acc, modulus);
        _fmpz_demote_val(Qcoeffs + Qlen);

        /* process nodes taken from the heap */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            if (i == -UWORD(1))
            {
                /* take next dividend term */
                if (j + 1 < Alen)
                {
                    x = chain + 0;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    _mpoly_heap_insert1(heap, Aexps[x->j], x,
                                                 &next_loc, &heap_len, cmpmask);
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
                                                 &next_loc, &heap_len, cmpmask);
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
                                                 &next_loc, &heap_len, cmpmask);
                }
            }
        }

        if (fmpz_is_zero(Qcoeffs + Qlen))
            continue;

        if (!lt_divides)
        {
            _fmpz_mod_mpoly_fit_length(&Rcoeffs, &R->coeffs_alloc,
                                       &Rexps, &R->exps_alloc, 1, Rlen + 1);
            fmpz_sub(Rcoeffs + Rlen, fmpz_mod_ctx_modulus(fctx), Qcoeffs + Qlen);
            Rexps[Rlen] = exp;
            Rlen++;
            continue;
        }

        fmpz_mod_mul(Qcoeffs + Qlen, Qcoeffs + Qlen, lc_minus_inv, fctx);

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
                                                 &next_loc, &heap_len, cmpmask);
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

    fmpz_clear(lc_minus_inv);
    mpz_clear(t);
    mpz_clear(acc);
    mpz_clear(modulus);

    return 1;

exp_overflow:

    Q->coeffs = Qcoeffs;
    Q->exps = Qexps;
    Q->length = 0;

    R->coeffs = Rcoeffs;
    R->exps = Rexps;
    R->length = 0;

    TMP_END;

    fmpz_clear(lc_minus_inv);
    mpz_clear(t);
    mpz_clear(acc);
    mpz_clear(modulus);

    return 0;
}



static int _fmpz_mod_mpoly_divrem_monagan_pearce(
    fmpz_mod_mpoly_t Q,
    fmpz_mod_mpoly_t R,
    const fmpz * Acoeffs, const ulong * Aexps, slong Alen,
    const fmpz * Bcoeffs, const ulong * Bexps, slong Blen,
    slong bits,
    slong N,
    const ulong * cmpmask,
    const fmpz_mod_ctx_t fctx)
{
    slong i, j, Qlen, Rlen, s;
    slong next_loc;
    slong heap_len = 2; /* heap zero index unused */
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    fmpz * Qcoeffs = Q->coeffs;
    fmpz * Rcoeffs = R->coeffs;
    ulong * Qexps = Q->exps;
    ulong * Rexps = R->exps;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    ulong mask;
    slong * hind;
    int lt_divides;
    mpz_t t, acc, modulus;
    ulong acc_sm[3];
    fmpz_t lc_minus_inv;
    TMP_INIT;

    if (N == 1)
    {
        if (Blen == 2)
            return _fmpz_mod_mpoly_divrem_monagan_pearce1_binomial(Q, R,
                 Acoeffs, Aexps, Alen, Bcoeffs, Bexps, bits, cmpmask[0], fctx);
        else
            return _fmpz_mod_mpoly_divrem_monagan_pearce1(Q, R, Acoeffs, Aexps, Alen,
                                   Bcoeffs, Bexps, Blen, bits, cmpmask[0], fctx);
    }

    mpz_init(t);
    mpz_init(acc);
    mpz_init(modulus);
    fmpz_get_mpz(modulus, fmpz_mod_ctx_modulus(fctx));

    fmpz_init(lc_minus_inv);
    fmpz_mod_inv(lc_minus_inv, Bcoeffs + 0, fctx);
    fmpz_mod_neg(lc_minus_inv, lc_minus_inv, fctx);

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

    mask = (bits <= FLINT_BITS) ? mpoly_overflow_mask_sp(bits) : 0;

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
   
    while (heap_len > 1)
    {
        _fmpz_mod_mpoly_fit_length(&Qcoeffs, &Q->coeffs_alloc,
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

        mpz_set_ui(acc, 0);
        acc_sm[2] = acc_sm[1] = acc_sm[0] = 0;

        do {
            exp_list[--exp_next] = heap[1].exp;
            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
            do {
                *store++ = x->i;
                *store++ = x->j;

                if (x->i == -WORD(1))
                {
                    fmpz Aj = Acoeffs[x->j];

                    if (COEFF_IS_MPZ(Aj))
                        mpz_sub(acc, acc, COEFF_TO_PTR(Aj));
                    else
                        flint_mpz_sub_ui(acc, acc, Aj);
                }
                else
                {
                    fmpz Bi = Bcoeffs[x->i];
                    fmpz Qj = Qcoeffs[x->j];

                    hind[x->i] |= WORD(1);

                    if (COEFF_IS_MPZ(Bi) && COEFF_IS_MPZ(Qj))
                    {
                        mpz_addmul(acc, COEFF_TO_PTR(Bi), COEFF_TO_PTR(Qj));
                    }
                    else if (COEFF_IS_MPZ(Bi) && !COEFF_IS_MPZ(Qj))
                    {
                        flint_mpz_addmul_ui(acc, COEFF_TO_PTR(Bi), Qj);
                    }
                    else if (!COEFF_IS_MPZ(Bi) && COEFF_IS_MPZ(Qj))
                    {
                        flint_mpz_addmul_ui(acc, COEFF_TO_PTR(Qj), Bi);
                    }
                    else
                    {
                        ulong pp1, pp0;
                        umul_ppmm(pp1, pp0, Bi, Qj);
                        add_sssaaaaaa(acc_sm[2], acc_sm[1], acc_sm[0],
                                      acc_sm[2], acc_sm[1], acc_sm[0],
                                      0, pp1, pp0);
                    }
                }
            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

        flint_mpz_add_uiuiui(acc, acc, acc_sm[2], acc_sm[1], acc_sm[0]);
        if (mpz_sgn(acc) < 0)
            mpz_add(acc, acc, modulus);

        mpz_tdiv_qr(t, _fmpz_promote(Qcoeffs + Qlen), acc, modulus);
        _fmpz_demote_val(Qcoeffs + Qlen);

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

        if (fmpz_is_zero(Qcoeffs + Qlen))
            continue;

        if (!lt_divides)
        {
            _fmpz_mod_mpoly_fit_length(&Rcoeffs, &R->coeffs_alloc,
                                       &Rexps, &R->exps_alloc, N, Rlen + 1);
            fmpz_sub(Rcoeffs + Rlen, fmpz_mod_ctx_modulus(fctx), Qcoeffs + Qlen);
            mpoly_monomial_set(Rexps + Rlen*N, exp, N);
            Rlen++;
            continue;
        }

        fmpz_mod_mul(Qcoeffs + Qlen, Qcoeffs + Qlen, lc_minus_inv, fctx);

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

    fmpz_clear(lc_minus_inv);
    mpz_clear(t);
    mpz_clear(acc);
    mpz_clear(modulus);

    return 1;

exp_overflow2:

    Q->coeffs = Qcoeffs;
    Q->exps = Qexps;
    Q->length = 0;

    R->coeffs = Rcoeffs;
    R->exps = Rexps;
    R->length = 0;

    TMP_END;

    fmpz_clear(lc_minus_inv);
    mpz_clear(t);
    mpz_clear(acc);
    mpz_clear(modulus);

    return 0;
}

void fmpz_mod_mpoly_divrem_monagan_pearce(
    fmpz_mod_mpoly_t Q,
    fmpz_mod_mpoly_t R,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong N;
    flint_bitcnt_t QRbits;
    ulong * Aexps = A->exps, * Bexps = B->exps;
    ulong * cmpmask;
    int freeAexps = 0, freeBexps = 0;
    fmpz_mod_mpoly_t TQ, TR;
    fmpz_mod_mpoly_struct * q, * r;

    if (fmpz_mod_mpoly_is_zero(B, ctx))
    {
        if (fmpz_is_one(fmpz_mod_ctx_modulus(ctx->ffinfo)))
        {
            fmpz_mod_mpoly_zero(Q, ctx);
            fmpz_mod_mpoly_zero(R, ctx);
            return;
        }
        else
        {
            flint_throw(FLINT_DIVZERO, "fmpz_mod_mpoly_divrem_monagan_pearce: divide by zero");
        }
    }

    if (fmpz_mod_mpoly_is_zero(A, ctx))
    {
        fmpz_mod_mpoly_zero(Q, ctx);
        fmpz_mod_mpoly_zero(R, ctx);
        return;
    }

    fmpz_mod_mpoly_init(TQ, ctx);
    fmpz_mod_mpoly_init(TR, ctx);

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
        fmpz_mod_mpoly_set(R, A, ctx);
        fmpz_mod_mpoly_zero(Q, ctx);
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
        fmpz_mod_mpoly_fit_length_reset_bits(q, A->length/B->length + 1, QRbits, ctx);
        fmpz_mod_mpoly_fit_length_reset_bits(r, B->length, QRbits, ctx);

        if (_fmpz_mod_mpoly_divrem_monagan_pearce(q, r,
                     A->coeffs, Aexps, A->length, B->coeffs, Bexps, B->length,
                                              QRbits, N, cmpmask, ctx->ffinfo))
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
        fmpz_mod_mpoly_swap(Q, TQ, ctx);

    if (R == A || R == B)
        fmpz_mod_mpoly_swap(R, TR, ctx);

cleanup:

    fmpz_mod_mpoly_clear(TQ, ctx);
    fmpz_mod_mpoly_clear(TR, ctx);

    if (freeAexps)
        flint_free(Aexps);

    if (freeBexps)
        flint_free(Bexps);

    flint_free(cmpmask);
}
