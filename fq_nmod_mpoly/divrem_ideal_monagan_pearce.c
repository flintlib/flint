/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "fq_nmod_mpoly.h"

/*
   As for divrem_monagan_pearce except that an array of divisor polynomials is
   passed and an array of quotient polynomials is returned. These are not in
   low level format.
*/
int _fq_nmod_mpoly_divrem_ideal_monagan_pearce(
    fq_nmod_mpoly_struct ** Q,
    fq_nmod_mpoly_t R,
    const mp_limb_t * poly2, const ulong * exp2, slong len2,
    fq_nmod_mpoly_struct * const * poly3, ulong * const * exp3, slong len,
    slong N,
    flint_bitcnt_t bits,
    const fq_nmod_mpoly_ctx_t ctx,
    const ulong * cmpmask)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
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
    mp_limb_t * acc, * pp, * lc_minus_inv;
    TMP_INIT;

    TMP_START;
   
    chains = (mpoly_nheap_t **) TMP_ALLOC(len*sizeof(mpoly_nheap_t *));
    hinds = (slong **) TMP_ALLOC(len*sizeof(slong *));
    len3 = 0;
    for (w = 0; w < len; w++)
    {
        chains[w] = (mpoly_nheap_t *) TMP_ALLOC((poly3[w]->length)*sizeof(mpoly_nheap_t));
        hinds[w] = (slong *) TMP_ALLOC((poly3[w]->length)*sizeof(slong));
        len3 += poly3[w]->length;
        for (i = 0; i < poly3[w]->length; i++)
            hinds[w][i] = 1;
    }

    next_loc = len3 + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap_s));
    store = store_base = (slong *) TMP_ALLOC(3*len3*sizeof(slong));

    exps = (ulong *) TMP_ALLOC(len3*N*sizeof(ulong));
    exp_list = (ulong **) TMP_ALLOC(len3*sizeof(ulong *));
    texp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    exp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    q_len = (slong *) TMP_ALLOC(len*sizeof(slong));
    s = (slong *) TMP_ALLOC(len*sizeof(slong));

    exp_next = 0;
    for (i = 0; i < len3; i++)
        exp_list[i] = exps + i*N;

    mask = bits <= FLINT_BITS ? mpoly_overflow_mask_sp(bits) : 0;

    r_len = 0;
    for (w = 0; w < len; w++)
    {
        q_len[w] = 0;
        s[w] = poly3[w]->length;
    }
   
    x = chains[0] + 0;
    x->i = -WORD(1);
    x->j = 0;
    x->p = -WORD(1);
    x->next = NULL;
    heap[1].next = x;
    heap[1].exp = exp_list[exp_next++];
    mpoly_monomial_set(heap[1].exp, exp2, N);

    /* precompute leading coeff info */

    pp = (mp_limb_t *) TMP_ALLOC(d*sizeof(mp_limb_t));
    acc = (mp_limb_t *) TMP_ALLOC(d*sizeof(mp_limb_t));
    lc_minus_inv = (mp_limb_t *) TMP_ALLOC(d*len*sizeof(mp_limb_t));
    for (w = 0; w < len; w++)
    {
        n_fq_inv(lc_minus_inv + d*w, poly3[w]->coeffs + d*0, ctx->fqctx);
        _n_fq_neg(lc_minus_inv + d*w, lc_minus_inv + d*w, d, ctx->fqctx->mod);
    }

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

        _n_fq_zero(acc, d);
        do {
            exp_list[--exp_next] = heap[1].exp;
            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
            do {
                *store++ = x->i;
                *store++ = x->j;
                *store++ = x->p;

                if (x->i == -WORD(1))
                {
                    n_fq_sub(acc, acc, poly2 + d*x->j, ctx->fqctx);
                }
                else
                {
                    hinds[x->p][x->i] |= WORD(1);
                    n_fq_mul(pp, poly3[x->p]->coeffs + d*x->i,
                                 Q[x->p]->coeffs + d*x->j, ctx->fqctx);
                    n_fq_add(acc, acc, pp, ctx->fqctx);
                }
            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

        while (store > store_base)
        {
            p = *--store;
            j = *--store;
            i = *--store;
            if (i == -WORD(1))
            {
                if (j + 1 < len2)
                {
                    x = chains[0] + 0;
                    x->i = -WORD(1);
                    x->j = j + 1;
                    x->p = p;
                    x->next = NULL;
                    mpoly_monomial_set(exp_list[exp_next], exp2 + N*x->j, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
            }
            else
            {
                /* should we go right? */
                if ( (i + 1 < poly3[p]->length)
                   && (hinds[p][i + 1] == 2*j + 1)
                   )
                {
                    x = chains[p] + i + 1;
                    x->i = i + 1;
                    x->j = j;
                    x->p = p;
                    x->next = NULL;
                    hinds[p][x->i] = 2*(x->j + 1) + 0;
                    mpoly_monomial_add_mp(exp_list[exp_next], exp3[x->p] + N*x->i,
                                                Q[x->p]->exps + N*x->j, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }

                /* should we go up? */
                if (j + 1 == q_len[p])
                {
                    s[p]++;
                }
                else if (  ((hinds[p][i] & 1) == 1)
                          && ((i == 1) || (hinds[p][i - 1] >= 2*(j + 2) + 1)) )
                {
                    x = chains[p] + i;
                    x->i = i;
                    x->j = j + 1;
                    x->p = p;
                    x->next = NULL;
                    hinds[p][x->i] = 2*(x->j + 1) + 0;
                    mpoly_monomial_add_mp(exp_list[exp_next], exp3[x->p] + N*x->i,
                                                Q[x->p]->exps + N*x->j, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
            }
        }

        if (_n_fq_is_zero(acc, d))
            continue;

        for (w = 0; w < len; w++)
        {
            int divides;

            if (bits <= FLINT_BITS)
                divides = mpoly_monomial_divides(texp, exp, exp3[w] + N*0, N, mask);
            else
                divides = mpoly_monomial_divides_mp(texp, exp, exp3[w] + N*0, N, bits);

            if (divides)
            {
                fq_nmod_mpoly_fit_length(Q[w], q_len[w] + 1, ctx);
                n_fq_mul(Q[w]->coeffs + d*q_len[w], acc, lc_minus_inv + d*w, ctx->fqctx);
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
                    mpoly_monomial_add_mp(exp_list[exp_next], exp3[w] + N*x->i, 
                                                   Q[w]->exps + N*x->j, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
                s[w] = 1;
                q_len[w]++;
                /* break out of w for loop and continue in heap loop */
                goto break_continue;
            }
        }

        /* if get here, no leading terms divided */
        _fq_nmod_mpoly_fit_length(&r_coeff, &R->coeffs_alloc, d,
                                  &r_exp, &R->exps_alloc, N, r_len + 1);
        _n_fq_neg(r_coeff + d*r_len, acc, d, ctx->fqctx->mod);
        mpoly_monomial_set(r_exp + N*r_len, exp, N);
        r_len++;

break_continue:;
    }

    R->coeffs = r_coeff;
    R->exps = r_exp;
    R->length = r_len;

    for (i = 0; i < len; i++)
        Q[i]->length = q_len[i];

    TMP_END;

    return 1;

exp_overflow:

    R->coeffs = r_coeff;
    R->exps = r_exp;
    R->length = 0;

    for (i = 0; i < len; i++)
        Q[i]->length = 0;

    TMP_END;

    return 0;
}

/* Assumes divisor polys don't alias any output polys */
void fq_nmod_mpoly_divrem_ideal_monagan_pearce(
    fq_nmod_mpoly_struct ** Q,
    fq_nmod_mpoly_t R,
    const fq_nmod_mpoly_t A,
    fq_nmod_mpoly_struct * const * B,
    slong len,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, N;
    flint_bitcnt_t QRbits;
    slong len3 = 0;
    ulong * cmpmask;
    ulong * Aexps;
    ulong ** Bexps;
    int freeAexps, * freeBexps;
    fq_nmod_mpoly_t TR;
    fq_nmod_mpoly_struct * r;
    TMP_INIT;

    for (i = 0; i < len; i++)
    {  
        len3 = FLINT_MAX(len3, B[i]->length);
        if (fq_nmod_mpoly_is_zero(B[i], ctx))
        {
            flint_throw(FLINT_DIVZERO, "fq_nmod_mpoly_divrem_ideal_monagan_pearce: divide by zero");
        }
    }

    /* dividend is zero, write out quotients and remainder */
    if (fq_nmod_mpoly_is_zero(A, ctx))
    {
        fq_nmod_mpoly_zero(R, ctx);
        for (i = 0; i < len; i++)
            fq_nmod_mpoly_zero(Q[i], ctx);
        return;
    }

    TMP_START;

    fq_nmod_mpoly_init(TR, ctx);

    freeBexps = (int *) TMP_ALLOC(len*sizeof(int));
    Bexps = (ulong **) TMP_ALLOC(len*sizeof(ulong *));

    /* compute maximum degrees that can occur in any input or output polys */
    QRbits = A->bits;
    for (i = 0; i < len; i++)
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

    for (i = 0; i < len; i++)
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
    for (i = 0; i < len; i++)
    {
        if (!mpoly_monomial_lt(Aexps + N*0, Bexps[i] + N*0, N, cmpmask))
            break;
    }

    if (i == len)
    {
        fq_nmod_mpoly_set(R, A, ctx);
        for (i = 0; i < len; i++)
            fq_nmod_mpoly_zero(Q[i], ctx);
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
        fq_nmod_mpoly_fit_length_reset_bits(r, len3, QRbits, ctx);        
        for (i = 0; i < len; i++)
            fq_nmod_mpoly_fit_length_reset_bits(Q[i], 1, QRbits, ctx);

        if (_fq_nmod_mpoly_divrem_ideal_monagan_pearce(Q, r, A->coeffs, Aexps,
                            A->length, B, Bexps, len, N, QRbits, ctx, cmpmask))
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

        for (i = 0; i < len; i++)
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
        fq_nmod_mpoly_swap(R, TR, ctx);

cleanup:

    fq_nmod_mpoly_clear(TR, ctx);

    if (freeAexps)
        flint_free(Aexps);

    for (i = 0; i < len; i++)
    {
        if (freeBexps[i])
            flint_free(Bexps[i]);
    }

    flint_free(cmpmask);

    TMP_END;
}

