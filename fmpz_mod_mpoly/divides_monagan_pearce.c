/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

static int _fmpz_mod_mpoly_divides_monagan_pearce1(
    fmpz_mod_mpoly_t Q,
    const fmpz * Acoeffs, const ulong * Aexps, slong Alen,
    const fmpz * Bcoeffs, const ulong * Bexps, slong Blen,
    slong bits,
    ulong cmpmask,
    const fmpz_mod_ctx_t fctx)
{
    int lt_divides;
    slong i, j, Qlen, s;
    slong next_loc, heap_len;
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    fmpz * Qcoeffs = Q->coeffs;
    ulong * Qexps = Q->exps;
    slong * hind;
    ulong mask, exp, maxexp = Aexps[Alen - 1];
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

    Qlen = 0;

    /* s is the number of terms * (latest quotient) we should put into heap */
    s = Blen;

    /* insert (-1, 0, Aexps[0]) into heap */
    heap_len = 2;
    x = chain + 0;
    x->i = -WORD(1);
    x->j = 0;
    x->next = NULL;
    HEAP_ASSIGN(heap[1], Aexps[0], x);

    while (heap_len > 1)
    {
        exp = heap[1].exp;

        if (mpoly_monomial_overflows1(exp, mask))
            goto not_exact_division;

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
                                                 &next_loc, &heap_len, cmpmask);
                }
            }
            else
            {
                /* should we go right? */
                if ((i + 1 < Blen) && (hind[i + 1] == 2*j + 1))
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
                }
                else if (((hind[i] & 1) == 1) &&
                         ((i == 1) || (hind[i - 1] >= 2*(j + 2) + 1)))
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

        fmpz_mod_mul(Qcoeffs + Qlen, Qcoeffs + Qlen, lc_minus_inv, fctx);
        if (fmpz_is_zero(Qcoeffs + Qlen))
            continue;

        if (!lt_divides || (exp^cmpmask) < (maxexp^cmpmask))
            goto not_exact_division;

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

cleanup:

    Q->coeffs = Qcoeffs;
    Q->exps = Qexps;
    Q->length = Qlen;

    TMP_END;

    fmpz_clear(lc_minus_inv);
    mpz_clear(t);
    mpz_clear(acc);
    mpz_clear(modulus);

    return Qlen > 0;

not_exact_division:

    Qlen = 0;
    goto cleanup;
}


int _fmpz_mod_mpoly_divides_monagan_pearce(
    fmpz_mod_mpoly_t Q,
    const fmpz * Acoeffs, const ulong * Aexps, slong Alen,
    const fmpz * Bcoeffs, const ulong * Bexps, slong Blen,
    flint_bitcnt_t bits,
    slong N,
    const ulong * cmpmask,
    const fmpz_mod_ctx_t fctx)
{
    int lt_divides;
    slong i, j, Qlen, s;
    slong next_loc, heap_len;
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    fmpz * Qcoeffs = Q->coeffs;
    ulong * Qexps = Q->exps;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    slong * hind;
    ulong mask;
    mpz_t t, acc, modulus;
    ulong acc_sm[3];
    fmpz_t lc_minus_inv;
    TMP_INIT;

    if (N == 1)
        return _fmpz_mod_mpoly_divides_monagan_pearce1(Q, Acoeffs, Aexps, Alen,
                                 Bcoeffs, Bexps, Blen, bits, cmpmask[0], fctx);

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

    Qlen = WORD(0);

    /* s is the number of terms * (latest quotient) we should put into heap */
    s = Blen;

    /* insert (-1, 0, Aexps[0]) into heap */
    heap_len = 2;
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
                goto not_exact_division;

            lt_divides = mpoly_monomial_divides(Qexps + Qlen*N, exp, Bexps,
                                                                      N, mask);
        }
        else
        {
            if (mpoly_monomial_overflows_mp(exp, N, bits))
                goto not_exact_division;

            lt_divides = mpoly_monomial_divides_mp(Qexps + Qlen*N, exp, Bexps,
                                                                      N, bits);
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
                /* should we go up */
                if ((i + 1 < Blen) && (hind[i + 1] == 2*j + 1))
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
                else if (((hind[i] & 1) == 1) &&
                         ((i == 1) || (hind[i - 1] >= 2*(j + 2) + 1)))
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

        fmpz_mod_mul(Qcoeffs + Qlen, Qcoeffs + Qlen, lc_minus_inv, fctx);
        if (fmpz_is_zero(Qcoeffs + Qlen))
            continue;

        if (!lt_divides ||
            mpoly_monomial_gt(Aexps + N*(Alen - 1), exp, N, cmpmask))
        {
            goto not_exact_division;
        }

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

cleanup:

    Q->coeffs = Qcoeffs;
    Q->exps = Qexps;
    Q->length = Qlen;

    TMP_END;

    fmpz_clear(lc_minus_inv);
    mpz_clear(t);
    mpz_clear(acc);
    mpz_clear(modulus);

    return Qlen > 0;

not_exact_division:

    Qlen = 0;
    goto cleanup;
}


int _fmpz_mod_mpoly_divides_monagan_pearce_maxfields(
    fmpz_mod_mpoly_t Q,
    const fmpz_mod_mpoly_t A, fmpz * maxAfields,
    const fmpz_mod_mpoly_t B, fmpz * maxBfields,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, N;
    flint_bitcnt_t Qbits;
    ulong * cmpmask;
    ulong * Aexps = A->exps, * Bexps = B->exps, * expq;
    int divides, freeAexps = 0, freeBexps = 0;
    TMP_INIT;

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        /*
            cannot be exact division if any max field from A
            is less than corresponding max field from B
        */
        if (fmpz_cmp(maxAfields + i, maxBfields + i) < 0)
        {
            fmpz_mod_mpoly_zero(Q, ctx);
            return 0;
        }
    }

    Qbits = 1 + _fmpz_vec_max_bits(maxAfields, ctx->minfo->nfields);
    Qbits = FLINT_MAX(Qbits, A->bits);
    Qbits = FLINT_MAX(Qbits, B->bits);
    Qbits = mpoly_fix_bits(Qbits, ctx->minfo);

    TMP_START;

    N = mpoly_words_per_exp(Qbits, ctx->minfo);
    cmpmask = TMP_ARRAY_ALLOC(2*N, ulong);
    expq = cmpmask + N; /* temp space to check leading monomials divide */
    mpoly_get_cmpmask(cmpmask, N, Qbits, ctx->minfo);

    /* quick check for easy case of inexact division of leading monomials */
    if (Qbits == A->bits && Qbits == B->bits && A->exps[N - 1] < B->exps[N - 1])
    {
        fmpz_mod_mpoly_zero(Q, ctx);
        divides = 0;
        goto cleanup;
    }

    /* ensure input exponents packed to same size as output exponents */
    if (Qbits != A->bits)
    {
        freeAexps = 1;
        Aexps = (ulong *) flint_malloc(N*A->length*sizeof(ulong));
        mpoly_repack_monomials(Aexps, Qbits, A->exps, A->bits,
                                                    A->length, ctx->minfo);
    }

    if (Qbits != B->bits)
    {
        freeBexps = 1;
        Bexps = (ulong *) flint_malloc(N*B->length*sizeof(ulong));
        mpoly_repack_monomials(Bexps, Qbits, B->exps, B->bits,
                                                    B->length, ctx->minfo);
    }

    /* check leading monomial divides exactly */
    if (Qbits > FLINT_BITS)
        divides = mpoly_monomial_divides_mp(expq, Aexps, Bexps, N, Qbits);
    else
        divides = mpoly_monomial_divides(expq, Aexps, Bexps, N,
                                                mpoly_overflow_mask_sp(Qbits));
    if (!divides)
    {
        fmpz_mod_mpoly_zero(Q, ctx);
        goto cleanup;
    }

    /* deal with aliasing and divide polynomials */
    if (Q == A || Q == B)
    {
        fmpz_mod_mpoly_t T;
        fmpz_mod_mpoly_init3(T, A->length/B->length + 1, Qbits, ctx);
        divides = _fmpz_mod_mpoly_divides_monagan_pearce(T,
                                      A->coeffs, Aexps, A->length,
                                      B->coeffs, Bexps, B->length,
                                            Qbits, N, cmpmask, ctx->ffinfo);
        fmpz_mod_mpoly_swap(T, Q, ctx);
        fmpz_mod_mpoly_clear(T, ctx);
    }
    else
    {
        fmpz_mod_mpoly_fit_length_reset_bits(Q, A->length/B->length + 1, Qbits, ctx);

        divides = _fmpz_mod_mpoly_divides_monagan_pearce(Q,
                                    A->coeffs, Aexps, A->length,
                                    B->coeffs, Bexps, B->length,
                                            Qbits, N, cmpmask, ctx->ffinfo);
    }

cleanup:

    if (freeAexps)
        flint_free(Aexps);

    if (freeBexps)
        flint_free(Bexps);

    TMP_END;

    return divides;
}

/* return 1 if quotient is exact */
int fmpz_mod_mpoly_divides_monagan_pearce(
    fmpz_mod_mpoly_t Q,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int divides;
    slong i;
    fmpz * maxAfields, * maxBfields;
    TMP_INIT;

    if (fmpz_mod_mpoly_is_zero(B, ctx))
    {
        if (!fmpz_mod_mpoly_is_zero(A, ctx) &&
            !fmpz_is_one(fmpz_mod_ctx_modulus(ctx->ffinfo)))
        {
            flint_throw(FLINT_DIVZERO, "fmpz_mod_mpoly_divides_monagan_pearce: divide by zero");
        }

        fmpz_mod_mpoly_zero(Q, ctx);
        return 1;
    }

    if (fmpz_mod_mpoly_is_zero(A, ctx))
    {
        fmpz_mod_mpoly_zero(Q, ctx);
        return 1;
    }

    TMP_START;

    maxAfields = TMP_ARRAY_ALLOC(2*ctx->minfo->nfields, fmpz);
    maxBfields = maxAfields + ctx->minfo->nfields;
    for (i = 0; i < 2*ctx->minfo->nfields; i++)
        fmpz_init(maxAfields + i);

    mpoly_max_fields_fmpz(maxAfields, A->exps, A->length, A->bits, ctx->minfo);
    mpoly_max_fields_fmpz(maxBfields, B->exps, B->length, B->bits, ctx->minfo);

    divides = _fmpz_mod_mpoly_divides_monagan_pearce_maxfields(Q,
                                            A, maxAfields, B, maxBfields, ctx);

    for (i = 0; i < 2*ctx->minfo->nfields; i++)
        fmpz_clear(maxAfields + i);

    TMP_END;
    return divides;
}
