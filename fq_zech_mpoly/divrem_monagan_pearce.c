/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly.h"

static slong _fq_zech_mpoly_divrem_monagan_pearce(slong * lenr,
                  fq_zech_struct ** polyq,      ulong ** expq, slong * allocq,
                  fq_zech_struct ** polyr,      ulong ** expr, slong * allocr,
            const fq_zech_struct * coeff2, const ulong * exp2, slong len2,
            const fq_zech_struct * coeff3, const ulong * exp3, slong len3,
         slong bits, slong N, const ulong * cmpmask, const fq_zech_ctx_t fqctx)
{
    slong i, j, q_len, r_len, s;
    slong next_loc;
    slong heap_len = 2; /* heap zero index unused */
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    fq_zech_struct * q_coeff = *polyq;
    fq_zech_struct * r_coeff = *polyr;
    ulong * q_exp = *expq;
    ulong * r_exp = *expr;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    ulong mask;
    slong * hind;
    int lt_divides;
    fq_zech_t lc_minus_inv, pp;
    TMP_INIT;

    TMP_START;
    fq_zech_init(pp, fqctx);
    fq_zech_init(lc_minus_inv, fqctx);

    /* alloc array of heap nodes which can be chained together */
    next_loc = len3 + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(len3*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*len3*sizeof(mpoly_heap_t *));

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
    r_len = WORD(0);
   
    /* s is the number of terms * (latest quotient) we should put into heap */
    s = len3;
   
    /* insert (-1, 0, exp2[0]) into heap */
    x = chain + 0;
    x->i = -WORD(1);
    x->j = 0;
    x->next = NULL;
    heap[1].next = x;
    heap[1].exp = exp_list[exp_next++];
    mpoly_monomial_set(heap[1].exp, exp2, N);

    /* precompute leading cofficient info */
    fq_zech_inv(lc_minus_inv, coeff3 + 0, fqctx);
    fq_zech_neg(lc_minus_inv, lc_minus_inv, fqctx);
   
    while (heap_len > 1)
    {
        _fq_zech_mpoly_fit_length(&q_coeff, &q_exp, allocq, q_len + 1, N, fqctx);

        mpoly_monomial_set(exp, heap[1].exp, N);

        if (bits <= FLINT_BITS)
        {
            if (mpoly_monomial_overflows(exp, N, mask))
                goto exp_overflow2;
            lt_divides = mpoly_monomial_divides(q_exp + q_len*N, exp, exp3, N, mask);
        }
        else
        {
            if (mpoly_monomial_overflows_mp(exp, N, bits))
                goto exp_overflow2;
            lt_divides = mpoly_monomial_divides_mp(q_exp + q_len*N, exp, exp3, N, bits);
        }

        fq_zech_zero(q_coeff + q_len, fqctx);
        do {
            exp_list[--exp_next] = heap[1].exp;
            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
            do {
                *store++ = x->i;
                *store++ = x->j;
                if (x->i != -WORD(1))
                    hind[x->i] |= WORD(1);

                if (x->i == -WORD(1))
                {
                    fq_zech_sub(q_coeff + q_len, q_coeff + q_len, coeff2 + x->j, fqctx);
                }
                else
                {
                    fq_zech_mul(pp, coeff3 + x->i, q_coeff + x->j, fqctx);
                    fq_zech_add(q_coeff + q_len, q_coeff + q_len, pp, fqctx);
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
                if (j + 1 < len2)
                {
                    x = chain + 0;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    mpoly_monomial_set(exp_list[exp_next], exp2 + x->j*N, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
            }
            else
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
                    mpoly_monomial_add_mp(exp_list[exp_next], exp3 + x->i*N,
                                                             q_exp + x->j*N, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
                /* should we go up? */
                if (j + 1 == q_len)
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
                    mpoly_monomial_add_mp(exp_list[exp_next], exp3 + x->i*N,
                                                             q_exp + x->j*N, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
            }
        }

        /* try to divide accumulated term by leading term */
        if (fq_zech_is_zero(q_coeff + q_len, fqctx))
        {
            continue;
        }
        if (!lt_divides)
        {
            _fq_zech_mpoly_fit_length(&r_coeff, &r_exp, allocr, r_len + 1, N, fqctx);
            fq_zech_neg(r_coeff + r_len, q_coeff + q_len, fqctx);
            mpoly_monomial_set(r_exp + r_len*N, exp, N);
            r_len++;
            continue;
        }

        fq_zech_mul(q_coeff + q_len, q_coeff + q_len, lc_minus_inv, fqctx);

        /* put newly generated quotient term back into the heap if neccesary */
        if (s > 1)
        {
            i = 1;
            x = chain + i;
            x->i = i;
            x->j = q_len;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;
            mpoly_monomial_add_mp(exp_list[exp_next], exp3 + x->i*N,
                                                     q_exp + x->j*N, N);
            exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
        }
        s = 1;
        q_len++;
    }


cleanup2:

    (*polyq) = q_coeff;
    (*expq) = q_exp;
    (*polyr) = r_coeff;
    (*expr) = r_exp;

    /* set remainder poly length */
    (*lenr) = r_len;

    TMP_END;
    fq_zech_clear(pp, fqctx);
    fq_zech_clear(lc_minus_inv, fqctx);

    /* return quotient poly length */
    return q_len;

exp_overflow2:
    q_len = 0;
    r_len = 0;
    goto cleanup2;
}

void fq_zech_mpoly_divrem_monagan_pearce(fq_zech_mpoly_t q, fq_zech_mpoly_t r,
                      const fq_zech_mpoly_t poly2, const fq_zech_mpoly_t poly3,
                                                 const fq_zech_mpoly_ctx_t ctx)
{
    slong exp_bits, N, lenq = 0, lenr = 0;
    ulong * exp2 = poly2->exps, * exp3 = poly3->exps;
    ulong * cmpmask;
    int free2 = 0, free3 = 0;
    fq_zech_mpoly_t temp1, temp2;
    fq_zech_mpoly_struct * tq, * tr;

    if (poly3->length == 0)
    {
        flint_throw(FLINT_DIVZERO, "Divide by zero in nmod_mpoly_divrem_monagan_pearce");
    }

    if (poly2->length == 0)
    {
        fq_zech_mpoly_zero(q, ctx);
        fq_zech_mpoly_zero(r, ctx);
        return;
    }

    exp_bits = FLINT_MAX(poly2->bits, poly3->bits);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);

    N = mpoly_words_per_exp(exp_bits, ctx->minfo);
    cmpmask = (ulong *) flint_malloc(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, exp_bits, ctx->minfo);

    /* ensure input exponents packed to same size as output exponents */
    if (exp_bits > poly2->bits)
    {
        free2 = 1;
        exp2 = (ulong *) flint_malloc(N*poly2->length*sizeof(ulong));
        mpoly_repack_monomials(exp2, exp_bits, poly2->exps, poly2->bits,
                                                    poly2->length, ctx->minfo);
    }

    if (exp_bits > poly3->bits)
    {
        free3 = 1;
        exp3 = (ulong *) flint_malloc(N*poly3->length*sizeof(ulong));
        mpoly_repack_monomials(exp3, exp_bits, poly3->exps, poly3->bits,
                                                    poly3->length, ctx->minfo);
    }

    /* check divisor leading monomial is at most that of the dividend */
    if (mpoly_monomial_lt(exp2, exp3, N, cmpmask))
    {
        fq_zech_mpoly_set(r, poly2, ctx);
        fq_zech_mpoly_zero(q, ctx);
        goto cleanup3;
    }

    /* take care of aliasing */
    if (q == poly2 || q == poly3)
    {
        fq_zech_mpoly_init2(temp1, poly2->length/poly3->length + 1,                                                                          ctx);
        fq_zech_mpoly_fit_bits(temp1, exp_bits, ctx);
        temp1->bits = exp_bits;
        tq = temp1;
    }
    else
    {
        fq_zech_mpoly_fit_length(q, poly2->length/poly3->length + 1, ctx);
        fq_zech_mpoly_fit_bits(q, exp_bits, ctx);
        q->bits = exp_bits;
        tq = q;
    }

    if (r == poly2 || r == poly3)
    {
        fq_zech_mpoly_init2(temp2, poly3->length, ctx);
        fq_zech_mpoly_fit_bits(temp2, exp_bits, ctx);
        temp2->bits = exp_bits;
        tr = temp2;
    }
    else
    {
        fq_zech_mpoly_fit_length(r, poly3->length, ctx);
        fq_zech_mpoly_fit_bits(r, exp_bits, ctx);
        r->bits = exp_bits;
        tr = r;
    }

    /* do division with remainder */
    while ((lenq = _fq_zech_mpoly_divrem_monagan_pearce(&lenr, &tq->coeffs, &tq->exps,
                      &tq->alloc, &tr->coeffs, &tr->exps, &tr->alloc, poly2->coeffs, exp2, 
                       poly2->length, poly3->coeffs, exp3, poly3->length, exp_bits,
                                         N, cmpmask, ctx->fqctx)) == 0
          && lenr == 0)
   {
        ulong * old_exp2 = exp2, * old_exp3 = exp3;
        slong old_exp_bits = exp_bits;

        exp_bits = mpoly_fix_bits(exp_bits + 1, ctx->minfo);

        N = mpoly_words_per_exp(exp_bits, ctx->minfo);
        cmpmask = (ulong *) flint_realloc(cmpmask, N*sizeof(ulong));
        mpoly_get_cmpmask(cmpmask, N, exp_bits, ctx->minfo);

        exp2 = (ulong *) flint_malloc(N*poly2->length*sizeof(ulong));
        mpoly_repack_monomials(exp2, exp_bits, old_exp2, old_exp_bits,
                                                    poly2->length, ctx->minfo);

        exp3 = (ulong *) flint_malloc(N*poly3->length*sizeof(ulong));
        mpoly_repack_monomials(exp3, exp_bits, old_exp3, old_exp_bits,
                                                    poly3->length, ctx->minfo);

        if (free2)
            flint_free(old_exp2);

        if (free3)
            flint_free(old_exp3);

        free2 = free3 = 1; 

        fq_zech_mpoly_fit_bits(tq, exp_bits, ctx);
        tq->bits = exp_bits;

        fq_zech_mpoly_fit_bits(tr, exp_bits, ctx);
        tr->bits = exp_bits;
    }

    /* deal with aliasing */
    if (q == poly2 || q == poly3)
    {
        fq_zech_mpoly_swap(temp1, q, ctx);
        fq_zech_mpoly_clear(temp1, ctx);
    }

    if (r == poly2 || r == poly3)
    {
        fq_zech_mpoly_swap(temp2, r, ctx);
        fq_zech_mpoly_clear(temp2, ctx);
    }

    _fq_zech_mpoly_set_length(q, lenq, ctx);
    _fq_zech_mpoly_set_length(r, lenr, ctx);

cleanup3:

    if (free2)
        flint_free(exp2);

    if (free3)
        flint_free(exp3);

    flint_free(cmpmask);
}
