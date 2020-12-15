/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly.h"

static slong _fq_zech_mpoly_divides_monagan_pearce(
                     fq_zech_struct ** coeff1,      ulong ** exp1, slong * alloc,
                const fq_zech_struct * coeff2, const ulong * exp2, slong len2,
                const fq_zech_struct * coeff3, const ulong * exp3, slong len3,
     flint_bitcnt_t bits, slong N, const ulong * cmpmask, const fq_zech_ctx_t fqctx)
{
    int lt_divides;
    slong i, j, q_len, s;
    slong next_loc, heap_len;
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    fq_zech_struct * q_coeff = * coeff1;
    ulong * q_exp = * exp1;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    fq_zech_t lc_minus_inv, pp;
    ulong mask;
    slong * hind;
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
    fq_zech_inv(lc_minus_inv, coeff3 + 0, fqctx);
    fq_zech_neg(lc_minus_inv, lc_minus_inv, fqctx);

    while (heap_len > 1)
    {
        _fq_zech_mpoly_fit_length(&q_coeff, &q_exp, alloc, q_len + 1, N, fqctx);

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

        fq_zech_zero(q_coeff + q_len, fqctx);
        do {
            exp_list[--exp_next] = heap[1].exp;
            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
            do {
                *store++ = x->i;
                *store++ = x->j;

                if (x->i == -WORD(1))
                {
                    fq_zech_sub(q_coeff + q_len, q_coeff + q_len, coeff2 + x->j, fqctx);
                }
                else
                {
                    hind[x->i] |= WORD(1);
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

        fq_zech_mul(q_coeff + q_len, q_coeff + q_len, lc_minus_inv, fqctx);

        if (fq_zech_is_zero(q_coeff + q_len, fqctx))
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


cleanup:

    *coeff1 = q_coeff;
    *exp1 = q_exp;

    TMP_END;
    fq_zech_clear(pp, fqctx);
    fq_zech_clear(lc_minus_inv, fqctx);

    return q_len;

not_exact_division:
    q_len = 0;
    goto cleanup;
}

/* return 1 if quotient is exact */
int fq_zech_mpoly_divides_monagan_pearce(fq_zech_mpoly_t poly1,
                  const fq_zech_mpoly_t poly2, const fq_zech_mpoly_t poly3,
                                                    const fq_zech_mpoly_ctx_t ctx)
{
    slong i, N, len = 0;
    flint_bitcnt_t exp_bits;
    fmpz * max_fields2, * max_fields3;
    ulong * cmpmask;
    ulong * exp2 = poly2->exps, * exp3 = poly3->exps, * expq;
    int easy_exit, free2 = 0, free3 = 0;
    ulong mask = 0;
    TMP_INIT;

    if (poly3->length == 0)
        flint_throw(FLINT_DIVZERO, "Divide by zero in fq_zech_mpoly_divides_monagan_pearce");

    if (poly2->length == 0)
    {
        fq_zech_mpoly_zero(poly1, ctx);
        return 1;
    }

    TMP_START;

    max_fields2 = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    max_fields3 = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_init(max_fields2 + i);
        fmpz_init(max_fields3 + i);
    }

    mpoly_max_fields_fmpz(max_fields2, poly2->exps, poly2->length,
                                                      poly2->bits, ctx->minfo);
    mpoly_max_fields_fmpz(max_fields3, poly3->exps, poly3->length,
                                                      poly3->bits, ctx->minfo);
    easy_exit = 0;
    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        /*
            cannot be exact division if any max field from poly2
            is less than corresponding max field from poly3
        */
        if (fmpz_cmp(max_fields2 + i, max_fields3 + i) < 0)
            easy_exit = 1;
    }

    exp_bits = _fmpz_vec_max_bits(max_fields2, ctx->minfo->nfields);
    exp_bits = FLINT_MAX(MPOLY_MIN_BITS, exp_bits + 1);
    exp_bits = FLINT_MAX(exp_bits, poly2->bits);
    exp_bits = FLINT_MAX(exp_bits, poly3->bits);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_clear(max_fields2 + i);
        fmpz_clear(max_fields3 + i);
    }

    if (easy_exit)
    {
        len = 0;
        goto cleanup;
    }

    N = mpoly_words_per_exp(exp_bits, ctx->minfo);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, exp_bits, ctx->minfo);

    /* temporary space to check leading monomials divide */
    expq = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    /* quick check for easy case of inexact division of leading monomials */
    if (poly2->bits == poly3->bits && N == 1 && 
       poly2->exps[0] < poly3->exps[0])
    {
        goto cleanup;
    }

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

    /* check leading monomial divides exactly */
    if (exp_bits <= FLINT_BITS)
    {
        /* mask with high bit of each exponent vector field set */
        for (i = 0; i < FLINT_BITS/exp_bits; i++)
            mask = (mask << exp_bits) + (UWORD(1) << (exp_bits - 1));

        if (!mpoly_monomial_divides(expq, exp2, exp3, N, mask))
        {
            len = 0;
            goto cleanup;
        }
    } else
    {
        if (!mpoly_monomial_divides_mp(expq, exp2, exp3, N, exp_bits))
        {
            len = 0;
            goto cleanup;
        }
    }

    /* deal with aliasing and divide polynomials */
    if (poly1 == poly2 || poly1 == poly3)
    {
      fq_zech_mpoly_t temp;

      fq_zech_mpoly_init2(temp, poly2->length/poly3->length + 1, ctx);
      fq_zech_mpoly_fit_bits(temp, exp_bits, ctx);
      temp->bits = exp_bits;

      len = _fq_zech_mpoly_divides_monagan_pearce(&temp->coeffs, &temp->exps,
                            &temp->alloc, poly2->coeffs, exp2, poly2->length,
                              poly3->coeffs, exp3, poly3->length, exp_bits, N,
                                                  cmpmask, ctx->fqctx);

      fq_zech_mpoly_swap(temp, poly1, ctx);

      fq_zech_mpoly_clear(temp, ctx);
   } else
   {
      fq_zech_mpoly_fit_length(poly1, poly2->length/poly3->length + 1, ctx);
      fq_zech_mpoly_fit_bits(poly1, exp_bits, ctx);
      poly1->bits = exp_bits;

      len = _fq_zech_mpoly_divides_monagan_pearce(&poly1->coeffs, &poly1->exps,
                            &poly1->alloc, poly2->coeffs, exp2, poly2->length,
                              poly3->coeffs, exp3, poly3->length, exp_bits, N,
                                                  cmpmask, ctx->fqctx);
   }

cleanup:

   _fq_zech_mpoly_set_length(poly1, len, ctx);

   if (free2)
      flint_free(exp2);

   if (free3)
      flint_free(exp3);

   TMP_END;

   /* division is exact if len is nonzero */
   return (len != 0);
}
