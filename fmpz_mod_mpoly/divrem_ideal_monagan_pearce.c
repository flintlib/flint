/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

int _fmpz_mod_mpoly_divrem_ideal_monagan_pearce(
    fmpz_mod_mpoly_struct ** Q,
    fmpz_mod_mpoly_t R,
    const fmpz * Acoeffs, const ulong * Aexps, slong Alen,
    fmpz_mod_mpoly_struct * const * Bs, ulong * const * Bexps, slong Blen,
    slong N,
    flint_bitcnt_t bits,
    const fmpz_mod_mpoly_ctx_t ctx,
    const ulong * cmpmask)
{
    int overflows, divides, success;
    slong i, j, p, w;
    slong next_loc;
    slong * store, * store_base;
    slong len3;
    slong heap_len = 2; /* heap zero index unused */
    mpoly_heap_s * heap;
    mpoly_nheap_t ** chains;
    slong ** hinds;
    mpoly_nheap_t * x;
    fmpz * r_coeff = R->coeffs;
    ulong * r_exp = R->exps;
    slong r_len;
    ulong * exp, * exps, * texp;
    ulong ** exp_list;
    slong exp_next;
    ulong mask;
    slong * q_len, * s;
    fmpz * lc_inv;
    fmpz_t acc;
    TMP_INIT;

    TMP_START;

    fmpz_init(acc);

    /* TODO figure out how to avoid TMP_ALLOC in a loop */
    chains = TMP_ARRAY_ALLOC(Blen, mpoly_nheap_t *);
    hinds = TMP_ARRAY_ALLOC(Blen, slong *);
    len3 = 0;
    for (w = 0; w < Blen; w++)
    {
        chains[w] = TMP_ARRAY_ALLOC(Bs[w]->length, mpoly_nheap_t);
        hinds[w] = TMP_ARRAY_ALLOC(Bs[w]->length, slong);
        len3 += Bs[w]->length;
        for (i = 0; i < Bs[w]->length; i++)
            hinds[w][i] = 1;
    }

    next_loc = len3 + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap_s));
    store = store_base = (slong *) TMP_ALLOC(3*len3*sizeof(slong *));

    exps = (ulong *) TMP_ALLOC(len3*N*sizeof(ulong));
    exp_list = (ulong **) TMP_ALLOC(len3*sizeof(ulong *));
    texp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    exp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    q_len = (slong *) TMP_ALLOC(Blen*sizeof(slong));
    s = (slong *) TMP_ALLOC(Blen*sizeof(slong));

    exp_next = 0;
    for (i = 0; i < len3; i++)
        exp_list[i] = exps + i*N;

    mask = bits <= FLINT_BITS ? mpoly_overflow_mask_sp(bits) : 0;

    for (w = 0; w < Blen; w++)
    {
        q_len[w] = WORD(0);
        s[w] = Bs[w]->length;
    }
    r_len = WORD(0);
   
    x = chains[0] + 0;
    x->i = -WORD(1);
    x->j = 0;
    x->p = -WORD(1);
    x->next = NULL;
    heap[1].next = x;
    heap[1].exp = exp_list[exp_next++];
    mpoly_monomial_set(heap[1].exp, Aexps, N);

    /* precompute leading coeff info */
    lc_inv = TMP_ARRAY_ALLOC(Blen, fmpz);
    for (w = 0; w < Blen; w++)
    {
        fmpz_init(lc_inv + w);
        fmpz_mod_inv(lc_inv + w, Bs[w]->coeffs + 0, ctx->ffinfo);
    }

    while (heap_len > 1)
    {
        mpoly_monomial_set(exp, heap[1].exp, N);

        overflows =  (bits <= FLINT_BITS) ?
                                    mpoly_monomial_overflows(exp, N, mask) :
                                    mpoly_monomial_overflows_mp(exp, N, bits);
        if (overflows)
            goto exp_overflow;

        fmpz_zero(acc);
        do {
            exp_list[--exp_next] = heap[1].exp;
            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
            do {
                *store++ = x->i;
                *store++ = x->j;
                *store++ = x->p;

                if (x->i == -WORD(1))
                {
                    fmpz_add(acc, acc, Acoeffs + x->j);
                }
                else
                {
                    hinds[x->p][x->i] |= WORD(1);
                    fmpz_submul(acc, Bs[x->p]->coeffs + x->i, Q[x->p]->coeffs + x->j);
                }

            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

        fmpz_mod_set_fmpz(acc, acc, ctx->ffinfo);

        while (store > store_base)
        {
            p = *--store;
            j = *--store;
            i = *--store;
            if (i == -WORD(1))
            {
                if (j + 1 < Alen)
                {
                    x = chains[0] + 0;
                    x->i = -WORD(1);
                    x->j = j + 1;
                    x->p = p;
                    x->next = NULL;
                    mpoly_monomial_set(exp_list[exp_next], Aexps + x->j*N, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
            }
            else
            {
                /* should we go right? */
                if ((i + 1 < Bs[p]->length) &&
                    (hinds[p][i + 1] == 2*j + 1))
                {
                    x = chains[p] + i + 1;
                    x->i = i + 1;
                    x->j = j;
                    x->p = p;
                    x->next = NULL;
                    hinds[p][x->i] = 2*(x->j + 1) + 0;
                    mpoly_monomial_add_mp(exp_list[exp_next],
                              Bexps[x->p] + N*x->i, Q[x->p]->exps + N*x->j, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }

                /* should we go up? */
                if (j + 1 == q_len[p])
                {
                    s[p]++;
                }
                else if (((hinds[p][i] & 1) == 1) &&
                         ((i == 1) || (hinds[p][i - 1] >= 2*(j + 2) + 1)))
                {
                    x = chains[p] + i;
                    x->i = i;
                    x->j = j + 1;
                    x->p = p;
                    x->next = NULL;
                    hinds[p][x->i] = 2*(x->j + 1) + 0;
                    mpoly_monomial_add_mp(exp_list[exp_next],
                              Bexps[x->p] + N*x->i, Q[x->p]->exps + N*x->j, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
            }
        }

        if (fmpz_is_zero(acc))
            continue;

        for (w = 0; w < Blen; w++)
        {
            divides = (bits <= FLINT_BITS) ?
                 mpoly_monomial_divides(texp, exp, Bexps[w] + N*0, N, mask) :
                 mpoly_monomial_divides_mp(texp, exp, Bexps[w] + N*0, N, bits);

            if (!divides)
                continue;

            fmpz_mod_mpoly_fit_length(Q[w], q_len[w] + 1, ctx);
            mpoly_monomial_set(Q[w]->exps + N*q_len[w], texp, N);
            fmpz_mod_mul(Q[w]->coeffs + q_len[w], acc, lc_inv + w, ctx->ffinfo);

            if (s[w] > 1)
            {
                i = 1;
                x = chains[w] + i;
                x->i = i;
                x->j = q_len[w];
                x->p = w;
                x->next = NULL;
                hinds[w][x->i] = 2*(x->j + 1) + 0;
                mpoly_monomial_add_mp(exp_list[exp_next], Bexps[w] + N*x->i, 
                                                       Q[w]->exps + N*x->j, N);
                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
            }
            s[w] = 1;
            q_len[w]++;
            goto break_continue; /* continue in heap loop */
        }

        /* if get here, no leading terms divided */
        _fmpz_mod_mpoly_fit_length(&r_coeff, &R->coeffs_alloc,
                                   &r_exp, &R->exps_alloc, N, r_len + 1);
        fmpz_set(r_coeff + r_len, acc);
        mpoly_monomial_set(r_exp + N*r_len, exp, N);
        r_len++;

break_continue:;
    }

    R->coeffs = r_coeff;
    R->exps = r_exp;
    R->length = r_len;

    for (i = 0; i < Blen; i++)
        Q[i]->length = q_len[i];

    success = 1;

cleanup:

    for (w = 0; w < Blen; w++)
        fmpz_clear(lc_inv + w);

    fmpz_clear(acc);

    TMP_END;

    return success;

exp_overflow:

    R->coeffs = r_coeff;
    R->exps = r_exp;
    R->length = 0;

    for (i = 0; i < Blen; i++)
        Q[i]->length = 0;

    success = 1;
    goto cleanup;
}

/* Assumes divisor polys don't alias any output polys */
void fmpz_mod_mpoly_divrem_ideal_monagan_pearce(
    fmpz_mod_mpoly_struct ** Q,
    fmpz_mod_mpoly_t R,
    const fmpz_mod_mpoly_t A,
    fmpz_mod_mpoly_struct * const * B,
    slong Blen,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, N;
    flint_bitcnt_t QRbits;
    slong len3 = 0;
    ulong * cmpmask;
    ulong * Aexps;
    ulong ** Bexps;
    int freeAexps, * freeBexps;
    fmpz_mod_mpoly_t TR;
    fmpz_mod_mpoly_struct * r;
    TMP_INIT;

    for (i = 0; i < Blen; i++)
    {  
        len3 = FLINT_MAX(len3, B[i]->length);
        if (fmpz_mod_mpoly_is_zero(B[i], ctx))
            flint_throw(FLINT_DIVZERO,
                 "fmpz_mod_mpoly_divrem_ideal_monagan_pearce: divide by zero");
    }

    /* dividend is zero, write out quotients and remainder */
    if (fmpz_mod_mpoly_is_zero(A, ctx))
    {
        fmpz_mod_mpoly_zero(R, ctx);
        for (i = 0; i < Blen; i++)
            fmpz_mod_mpoly_zero(Q[i], ctx);
        return;
    }

    TMP_START;

    fmpz_mod_mpoly_init(TR, ctx);

    freeBexps = TMP_ARRAY_ALLOC(Blen, int);
    Bexps = TMP_ARRAY_ALLOC(Blen, ulong *);

    /* compute maximum degrees that can occur in any input or output polys */
    QRbits = A->bits;
    for (i = 0; i < Blen; i++)
        QRbits = FLINT_MAX(QRbits, B[i]->bits);
    QRbits = mpoly_fix_bits(QRbits, ctx->minfo);

    N = mpoly_words_per_exp(QRbits, ctx->minfo);
    cmpmask = FLINT_ARRAY_ALLOC(N, ulong);
    mpoly_get_cmpmask(cmpmask, N, QRbits, ctx->minfo);

    /* ensure input exponents packed to same size as output exponents */
    Aexps = A->exps;
    freeAexps = 0;
    if (QRbits > A->bits)
    {
        freeAexps = 1;
        Aexps = FLINT_ARRAY_ALLOC(N*A->length, ulong);
        mpoly_repack_monomials(Aexps, QRbits, A->exps, A->bits,
                                                        A->length, ctx->minfo);
    }

    for (i = 0; i < Blen; i++)
    {
        Bexps[i] = B[i]->exps;
        freeBexps[i] = 0;
        if (QRbits > B[i]->bits)
        {
            freeBexps[i] = 1;
            Bexps[i] = FLINT_ARRAY_ALLOC(N*B[i]->length, ulong);
            mpoly_repack_monomials(Bexps[i], QRbits, B[i]->exps, B[i]->bits,
                                                     B[i]->length, ctx->minfo);
        }
    }

    /* check leading mon. of at least one divisor is at most that of dividend */
    for (i = 0; i < Blen; i++)
    {
        if (!mpoly_monomial_lt(Aexps + N*0, Bexps[i] + N*0, N, cmpmask))
            break;
    }

    if (i == Blen)
    {
        fmpz_mod_mpoly_set(R, A, ctx);
        for (i = 0; i < Blen; i++)
            fmpz_mod_mpoly_zero(Q[i], ctx);
        goto cleanup;
    }

    /* take care of aliasing */
    if (R == A)
        r = TR;
    else
        r = R;

    /* do division with remainder */
    while (1)
    {
        fmpz_mod_mpoly_fit_length_reset_bits(r, len3, QRbits, ctx);        
        for (i = 0; i < Blen; i++)
            fmpz_mod_mpoly_fit_length_reset_bits(Q[i], 1, QRbits, ctx);

        if (_fmpz_mod_mpoly_divrem_ideal_monagan_pearce(Q, r, A->coeffs, Aexps,
                           A->length, B, Bexps, Blen, N, QRbits, ctx, cmpmask))
        {
            break;
        }

        QRbits = mpoly_fix_bits(QRbits + 1, ctx->minfo);
        N = mpoly_words_per_exp(QRbits, ctx->minfo);
        cmpmask = FLINT_ARRAY_REALLOC(cmpmask, N, ulong);
        mpoly_get_cmpmask(cmpmask, N, QRbits, ctx->minfo);

        if (freeAexps)
            flint_free(Aexps);
        Aexps = FLINT_ARRAY_ALLOC(N*A->length, ulong);
        mpoly_repack_monomials(Aexps, QRbits, A->exps, A->bits,
                                                        A->length, ctx->minfo);
        freeAexps = 1; 

        for (i = 0; i < Blen; i++)
        {
            if (freeBexps[i])
                flint_free(Bexps[i]);

            Bexps[i] = FLINT_ARRAY_ALLOC(N*B[i]->length, ulong);
            mpoly_repack_monomials(Bexps[i], QRbits, B[i]->exps, B[i]->bits,
                                                     B[i]->length, ctx->minfo);
            freeBexps[i] = 1;
        }
    }

    /* take care of aliasing */
    if (R == A)
        fmpz_mod_mpoly_swap(R, TR, ctx);

cleanup:

    fmpz_mod_mpoly_clear(TR, ctx);

    if (freeAexps)
        flint_free(Aexps);

    for (i = 0; i < Blen; i++)
    {
        if (freeBexps[i])
            flint_free(Bexps[i]);
    }

    flint_free(cmpmask);

    TMP_END;
}

