/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "nmod_mpoly.h"

/*
   As for divrem_monagan_pearce1 except that an array of divisor polynomials is
   passed and an array of quotient polynomials is returned. These are not in
   low level format.
*/
int _nmod_mpoly_divrem_ideal_monagan_pearce1(
    nmod_mpoly_struct ** Q,
    nmod_mpoly_t R,
    const mp_limb_t * Acoeffs, const ulong * Aexps, slong Alen,
    nmod_mpoly_struct * const * Bs, ulong * const * Bexps, slong Blen,
    flint_bitcnt_t bits,
    const nmod_mpoly_ctx_t ctx,
    ulong maskhi)
{
    slong i, j, p, w;
    slong next_loc;
    slong * store, * store_base;
    slong len3;
    slong heap_len = 2; /* heap zero index unused */
    mpoly_heap1_s * heap;
    mpoly_nheap_t ** chains;
    slong ** hinds;
    mpoly_nheap_t * x;
    mp_limb_t * r_coeff = R->coeffs;
    ulong * r_exp = R->exps;
    slong r_len;
    ulong exp, texp;
    ulong mask;
    slong * q_len, * s;
    mp_limb_t * lc_minus_inv, acc0, acc1, acc2, pp1, pp0;
    TMP_INIT;

    TMP_START;

    chains = (mpoly_nheap_t **) TMP_ALLOC(Blen*sizeof(mpoly_nheap_t *));
    hinds = (slong **) TMP_ALLOC(Blen*sizeof(slong *));

    len3 = 0;
    for (w = 0; w < Blen; w++)
    {
        chains[w] = (mpoly_nheap_t *) TMP_ALLOC((Bs[w]->length)*sizeof(mpoly_nheap_t));
        hinds[w] = (slong *) TMP_ALLOC((Bs[w]->length)*sizeof(slong));
        len3 += Bs[w]->length;
        for (i = 0; i < Bs[w]->length; i++)
            hinds[w][i] = 1;
    }

    next_loc = len3 + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap1_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap1_s));
    store = store_base = (slong *) TMP_ALLOC(3*len3*sizeof(slong *));

    q_len = (slong *) TMP_ALLOC(Blen*sizeof(slong));
    s = (slong *) TMP_ALLOC(Blen*sizeof(slong));

    mask = mpoly_overflow_mask_sp(bits);

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
    HEAP_ASSIGN(heap[1], Aexps[0], x);

    /* precompute leading coeff info */
    lc_minus_inv = (mp_limb_t *) TMP_ALLOC(Blen*sizeof(mp_limb_t));
    for (w = 0; w < Blen; w++)
        lc_minus_inv[w] = ctx->mod.n - nmod_inv(Bs[w]->coeffs[0], ctx->mod);

    while (heap_len > 1)
    {
        exp = heap[1].exp;
        if (mpoly_monomial_overflows1(exp, mask))
            goto exp_overflow;

        acc0 = acc1 = acc2 = 0;
        do
        {
            x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
            do
            {
                *store++ = x->i;
                *store++ = x->j;
                *store++ = x->p;

                if (x->i != -WORD(1))
                    hinds[x->p][x->i] |= WORD(1);

                if (x->i == -WORD(1))
                {
                    add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0,
                           WORD(0), WORD(0), ctx->mod.n - Acoeffs[x->j]);
                } else
                {
                    umul_ppmm(pp1, pp0, Bs[x->p]->coeffs[x->i], Q[x->p]->coeffs[x->j]);
                    add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);
                }

            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && heap[1].exp == exp);

        NMOD_RED3(acc0, acc2, acc1, acc0, ctx->mod);

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
                    _mpoly_heap_insert1(heap, Aexps[x->j], x, &next_loc, &heap_len, maskhi);
                }
            } else
            {
                if ( (i + 1 < Bs[p]->length)
                   && (hinds[p][i + 1] == 2*j + 1)
                   )
                {
                    x = chains[p] + i + 1;
                    x->i = i + 1;
                    x->j = j;
                    x->p = p;
                    x->next = NULL;
                    hinds[p][x->i] = 2*(x->j + 1) + 0;
                    _mpoly_heap_insert1(heap, Bexps[x->p][x->i] +
                                Q[x->p]->exps[x->j], x, &next_loc, &heap_len, maskhi);
                }
                if (j + 1 == q_len[p])
                {
                    s[p]++;
                } else if (  ((hinds[p][i] & 1) == 1)
                          && ((i == 1) || (hinds[p][i - 1] >= 2*(j + 2) + 1))
                          )
                {
                    x = chains[p] + i;
                    x->i = i;
                    x->j = j + 1;
                    x->p = p;
                    x->next = NULL;
                    hinds[p][x->i] = 2*(x->j + 1) + 0;
                    _mpoly_heap_insert1(heap, Bexps[x->p][x->i] +
                                Q[x->p]->exps[x->j], x, &next_loc, &heap_len, maskhi);
                }
            }
        }

        if (acc0 == 0)
            continue;

        for (w = 0; w < Blen; w++)
        {
            if (mpoly_monomial_divides1(&texp, exp, Bexps[w][0], mask))
            {
                nmod_mpoly_fit_length(Q[w], q_len[w] + 1, ctx);
                Q[w]->coeffs[q_len[w]] = nmod_mul(acc0, lc_minus_inv[w], ctx->mod);
                Q[w]->exps[q_len[w]] = texp;
                if (s[w] > 1)
                {
                    i = 1;
                    x = chains[w] + i;
                    x->i = i;
                    x->j = q_len[w];
                    x->p = w;
                    x->next = NULL;
                    hinds[w][x->i] = 2*(x->j + 1) + 0;
                    _mpoly_heap_insert1(heap, Bexps[w][x->i] + Q[w]->exps[x->j],
                                               x, &next_loc, &heap_len, maskhi);
                }
                s[w] = 1;
                q_len[w]++;
                goto break_continue; /* break out of w for loop and continue in heap loop */
            }
        }

        /* if get here, no leading terms divided */
        _nmod_mpoly_fit_length(&r_coeff, &R->coeffs_alloc,
                               &r_exp, &R->exps_alloc, 1, r_len + 1);
        r_coeff[r_len] = ctx->mod.n - acc0;
        r_exp[r_len] = exp;
        r_len++;

break_continue:;
    }

    R->coeffs = r_coeff;
    R->exps = r_exp;
    R->length = r_len;

    for (i = 0; i < Blen; i++)
        Q[i]->length = q_len[i];

    TMP_END;

    return 1;

exp_overflow:

    R->coeffs = r_coeff;
    R->exps = r_exp;
    R->length = 0;

    for (i = 0; i < Blen; i++)
        Q[i]->length = 0;

    TMP_END;

    return 0;
}

/*
   As for divrem_monagan_pearce except that an array of divisor polynomials is
   passed and an array of quotient polynomials is returned. These are not in
   low level format.
*/
int _nmod_mpoly_divrem_ideal_monagan_pearce(
    nmod_mpoly_struct ** Q,
    nmod_mpoly_t R,
    const mp_limb_t * Acoeffs, const ulong * Aexps, slong Alen,
    nmod_mpoly_struct * const * Bs, ulong * const * Bexps, slong Blen,
    slong N,
    flint_bitcnt_t bits,
    const nmod_mpoly_ctx_t ctx,
    const ulong * cmpmask)
{
    slong i, j, p, w;
    slong next_loc;
    slong * store, * store_base;
    slong len3;
    slong heap_len = 2; /* heap zero index unused */
    mpoly_heap_s * heap;
    mpoly_nheap_t ** chains;
    slong ** hinds;
    mpoly_nheap_t * x;
    mp_limb_t * r_coeff = R->coeffs;
    ulong * r_exp = R->exps;
    slong r_len;
    ulong * exp, * exps, * texp;
    ulong ** exp_list;
    slong exp_next;
    ulong mask;
    slong * q_len, * s;
    mp_limb_t * lc_minus_inv, acc0, acc1, acc2, pp1, pp0;
    TMP_INIT;

    if (N == 1)
    {
        return _nmod_mpoly_divrem_ideal_monagan_pearce1(Q, R,
                   Acoeffs, Aexps, Alen, Bs, Bexps, Blen, bits, ctx, cmpmask[0]);
    }

    TMP_START;
   
    chains = (mpoly_nheap_t **) TMP_ALLOC(Blen*sizeof(mpoly_nheap_t *));
    hinds = (slong **) TMP_ALLOC(Blen*sizeof(mpoly_heap_t *));
    len3 = 0;
    for (w = 0; w < Blen; w++)
    {
        chains[w] = (mpoly_nheap_t *) TMP_ALLOC((Bs[w]->length)*sizeof(mpoly_nheap_t));
        hinds[w] = (slong *) TMP_ALLOC((Bs[w]->length)*sizeof(slong));
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
    lc_minus_inv = (mp_limb_t *) TMP_ALLOC(Blen*sizeof(mp_limb_t));
    for (w = 0; w < Blen; w++)
        lc_minus_inv[w] = ctx->mod.n - nmod_inv(Bs[w]->coeffs[0], ctx->mod);

    while (heap_len > 1)
    {
        mpoly_monomial_set(exp, heap[1].exp, N);
        if (bits <= FLINT_BITS)
        {
            if (mpoly_monomial_overflows(exp, N, mask))
                goto exp_overflow;
        }
        else
        {
            if (mpoly_monomial_overflows_mp(exp, N, bits))
                goto exp_overflow;
        }

        acc0 = acc1 = acc2 = 0;
        do
        {
            exp_list[--exp_next] = heap[1].exp;
            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
            do
            {
                *store++ = x->i;
                *store++ = x->j;
                *store++ = x->p;

                if (x->i == -WORD(1))
                {
                    add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0,
                         UWORD(0), UWORD(0), ctx->mod.n - Acoeffs[x->j]);
                }
                else
                {
                    hinds[x->p][x->i] |= WORD(1);
                    umul_ppmm(pp1, pp0, Bs[x->p]->coeffs[x->i], Q[x->p]->coeffs[x->j]);
                    add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, UWORD(0), pp1, pp0);
                }

            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

        NMOD_RED3(acc0, acc2, acc1, acc0, ctx->mod);

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
                if ( (i + 1 < Bs[p]->length)
                   && (hinds[p][i + 1] == 2*j + 1)
                   )
                {
                    x = chains[p] + i + 1;
                    x->i = i + 1;
                    x->j = j;
                    x->p = p;
                    x->next = NULL;
                    hinds[p][x->i] = 2*(x->j + 1) + 0;
                    mpoly_monomial_add_mp(exp_list[exp_next], Bexps[x->p] + x->i*N,
                                                Q[x->p]->exps + x->j*N, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }

                /* should we go up? */
                if (j + 1 == q_len[p])
                {
                    s[p]++;
                } else if (  ((hinds[p][i] & 1) == 1)
                          && ((i == 1) || (hinds[p][i - 1] >= 2*(j + 2) + 1))
                          )
                {
                    x = chains[p] + i;
                    x->i = i;
                    x->j = j + 1;
                    x->p = p;
                    x->next = NULL;
                    hinds[p][x->i] = 2*(x->j + 1) + 0;
                    mpoly_monomial_add_mp(exp_list[exp_next], Bexps[x->p] + x->i*N,
                                                Q[x->p]->exps + x->j*N, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
            }
        }

        if (acc0 == 0)
            continue;

        for (w = 0; w < Blen; w++)
        {
            int divides;

            if (bits <= FLINT_BITS)
            {
                divides = mpoly_monomial_divides(texp, exp, Bexps[w] + N*0, N, mask);
            }
            else
            {
                divides = mpoly_monomial_divides_mp(texp, exp, Bexps[w] + N*0, N, bits);
            }

            if (divides)
            {
                nmod_mpoly_fit_length(Q[w], q_len[w] + 1, ctx);
                Q[w]->coeffs[q_len[w]] = nmod_mul(acc0, lc_minus_inv[w], ctx->mod);
                mpoly_monomial_set(Q[w]->exps + N*q_len[w], texp, N);
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
                goto break_continue; /* break out of w for loop and continue in heap loop */
            }
        }

        /* if get here, no leading terms divided */
        _nmod_mpoly_fit_length(&r_coeff, &R->coeffs_alloc,
                               &r_exp, &R->exps_alloc, N, r_len + 1);
        r_coeff[r_len] = ctx->mod.n - acc0;
        mpoly_monomial_set(r_exp + r_len*N, exp, N);
        r_len++;

break_continue:;
    }

    R->coeffs = r_coeff;
    R->exps = r_exp;
    R->length = r_len;

    for (i = 0; i < Blen; i++)
        Q[i]->length = q_len[i];

    TMP_END;

    return 1;

exp_overflow:

    R->coeffs = r_coeff;
    R->exps = r_exp;
    R->length = 0;

    for (i = 0; i < Blen; i++)
        Q[i]->length = 0;

    TMP_END;

    return 0;
}

/* Assumes divisor polys don't alias any output polys */
void nmod_mpoly_divrem_ideal_monagan_pearce(
    nmod_mpoly_struct ** Q,
    nmod_mpoly_t R,
    const nmod_mpoly_t A,
    nmod_mpoly_struct * const * B,
    slong Blen,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, N;
    flint_bitcnt_t QRbits;
    slong len3 = 0;
    ulong * cmpmask;
    ulong * Aexps;
    ulong ** Bexps;
    int freeAexps, * freeBexps;
    nmod_mpoly_t TR;
    nmod_mpoly_struct * r;
    TMP_INIT;

    for (i = 0; i < Blen; i++)
    {  
        len3 = FLINT_MAX(len3, B[i]->length);
        if (nmod_mpoly_is_zero(B[i], ctx))
        {
            flint_throw(FLINT_DIVZERO, "nmod_mpoly_divrem_ideal_monagan_pearce: divide by zero");
        }
    }

    /* dividend is zero, write out quotients and remainder */
    if (nmod_mpoly_is_zero(A, ctx))
    {
        nmod_mpoly_zero(R, ctx);
        for (i = 0; i < Blen; i++)
            nmod_mpoly_zero(Q[i], ctx);
        return;
    }

    TMP_START;

    nmod_mpoly_init(TR, ctx);

    freeBexps = (int *) TMP_ALLOC(Blen*sizeof(int));
    Bexps = (ulong **) TMP_ALLOC(Blen*sizeof(ulong *));

    /* compute maximum degrees that can occur in any input or output polys */
    QRbits = A->bits;
    for (i = 0; i < Blen; i++)
        QRbits = FLINT_MAX(QRbits, B[i]->bits);
    QRbits = mpoly_fix_bits(QRbits, ctx->minfo);

    N = mpoly_words_per_exp(QRbits, ctx->minfo);
    cmpmask = (ulong *) flint_malloc(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, QRbits, ctx->minfo);

    /* ensure input exponents packed to same size as output exponents */
    Aexps = A->exps;
    freeAexps = 0;
    if (QRbits > A->bits)
    {
        freeAexps = 1;
        Aexps = (ulong *) flint_malloc(N*A->length*sizeof(ulong));
        mpoly_repack_monomials(Aexps, QRbits, A->exps, A->bits, A->length, ctx->minfo);
    }

    for (i = 0; i < Blen; i++)
    {
        Bexps[i] = B[i]->exps;
        freeBexps[i] = 0;
        if (QRbits > B[i]->bits)
        {
            freeBexps[i] = 1;
            Bexps[i] = (ulong *) flint_malloc(N*B[i]->length*sizeof(ulong));
            mpoly_repack_monomials(Bexps[i], QRbits, B[i]->exps, B[i]->bits, B[i]->length, ctx->minfo);
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
        nmod_mpoly_set(R, A, ctx);
        for (i = 0; i < Blen; i++)
            nmod_mpoly_zero(Q[i], ctx);
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
        nmod_mpoly_fit_length_reset_bits(r, len3, QRbits, ctx);        
        for (i = 0; i < Blen; i++)
            nmod_mpoly_fit_length_reset_bits(Q[i], 1, QRbits, ctx);

        if (_nmod_mpoly_divrem_ideal_monagan_pearce(Q, r, A->coeffs, Aexps,
                            A->length, B, Bexps, Blen, N, QRbits, ctx, cmpmask))
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

        for (i = 0; i < Blen; i++)
        {
            if (freeBexps[i])
                flint_free(Bexps[i]);

            Bexps[i] = (ulong *) flint_malloc(N*B[i]->length*sizeof(ulong));
            mpoly_repack_monomials(Bexps[i], QRbits, B[i]->exps, B[i]->bits, B[i]->length, ctx->minfo);
            freeBexps[i] = 1;
        }
    }

    /* take care of aliasing */
    if (R == A)
        nmod_mpoly_swap(R, TR, ctx);

cleanup:

    nmod_mpoly_clear(TR, ctx);

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

