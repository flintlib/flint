/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"


slong _nmod_mpoly_mul_johnson1(mp_limb_t ** coeff1, ulong ** exp1, slong * alloc,
              const mp_limb_t * coeff2, const ulong * exp2, slong len2,
              const mp_limb_t * coeff3, const ulong * exp3, slong len3,
                                          ulong maskhi, const nmodf_ctx_t fctx)
{
    slong i, j;
    slong next_loc;
    slong Q_len = 0, heap_len = 2; /* heap zero index unused */
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * Q;
    mpoly_heap_t * x;
    slong len1;
    mp_limb_t * p1 = * coeff1;
    ulong * e1 = * exp1;
    slong * hind;
    ulong exp;
    ulong acc0, acc1, acc2, pp0, pp1;
    TMP_INIT;

    TMP_START;

    next_loc = len2 + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap1_s *) TMP_ALLOC((len2 + 1)*sizeof(mpoly_heap1_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(len2*sizeof(mpoly_heap_t));
    Q = (slong *) TMP_ALLOC(2*len2*sizeof(slong));

    /* space for heap indices */
    hind = (slong *) TMP_ALLOC(len2*sizeof(slong));
    for (i = 0; i < len2; i++)
        hind[i] = 1;

    /* put (0, 0, exp2[0] + exp3[0]) on heap */
    x = chain + 0;
    x->i = 0;
    x->j = 0;
    x->next = NULL;

    HEAP_ASSIGN(heap[1], exp2[0] + exp3[0], x);
    hind[0] = 2*1 + 0;

    len1 = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        _nmod_mpoly_fit_length(&p1, &e1, alloc, len1 + 1, 1);

        e1[len1] = exp;

        acc0 = acc1 = acc2 = 0;
        do
        {
            x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
         
            hind[x->i] |= WORD(1);
            Q[Q_len++] = x->i;
            Q[Q_len++] = x->j;
            umul_ppmm(pp1, pp0, coeff2[x->i], coeff3[x->j]);
            add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);

            while ((x = x->next) != NULL)
            {
                hind[x->i] |= WORD(1);
                Q[Q_len++] = x->i;
                Q[Q_len++] = x->j;
                umul_ppmm(pp1, pp0, coeff2[x->i], coeff3[x->j]);
                add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);
            }
        } while (heap_len > 1 && heap[1].exp == exp);

        NMOD_RED3(p1[len1], acc2, acc1, acc0, fctx->mod);
        len1 += (p1[len1] != 0);
      
        while (Q_len > 0)
        {
            j = Q[--Q_len];
            i = Q[--Q_len];

            /* should we go right? */
            if (  (i + 1 < len2)
               && (hind[i + 1] == 2*j + 1)
               )
            {
                x = chain + i + 1;
                x->i = i + 1;
                x->j = j;
                x->next = NULL;

                hind[x->i] = 2*(x->j + 1) + 0;
                _mpoly_heap_insert1(heap, exp2[x->i] + exp3[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
            }

            /* should we go up? */
            if (  (j + 1 < len3)
               && ((hind[i] & 1) == 1)
               && (  (i == 0)
                  || (hind[i - 1] >= 2*(j + 2) + 1)
                  )
               )
            {
                x = chain + i;
                x->i = i;
                x->j = j + 1;
                x->next = NULL;

                hind[x->i] = 2*(x->j + 1) + 0;
                _mpoly_heap_insert1(heap, exp2[x->i] + exp3[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
            }
        }
    }

    (* coeff1) = p1;
    (* exp1) = e1;
   
    TMP_END;

    return len1;
}


slong _nmod_mpoly_mul_johnson(mp_limb_t ** coeff1, ulong ** exp1, slong * alloc,
                 const mp_limb_t * coeff2, const ulong * exp2, slong len2,
                 const mp_limb_t * coeff3, const ulong * exp3, slong len3,
      mp_bitcnt_t bits, slong N, const ulong * cmpmask, const nmodf_ctx_t fctx)
{
    slong i, j;
    slong next_loc;
    slong Q_len = 0, heap_len = 2; /* heap zero index unused */
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * Q;
    mpoly_heap_t * x;
    slong len1;
    mp_limb_t * p1 = * coeff1;
    ulong * e1 = *exp1;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    slong * hind;
    ulong acc0, acc1, acc2, pp0, pp1;
    TMP_INIT;

    if (N == 1)
        return _nmod_mpoly_mul_johnson1(coeff1, exp1, alloc,
                     coeff2, exp2, len2, coeff3, exp3, len3, cmpmask[0], fctx);

    TMP_START;

    next_loc = len2 + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap_s *) TMP_ALLOC((len2 + 1)*sizeof(mpoly_heap_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(len2*sizeof(mpoly_heap_t));
    Q = (slong *) TMP_ALLOC(2*len2*sizeof(slong));
    exps = (ulong *) TMP_ALLOC(len2*N*sizeof(ulong));
    exp_list = (ulong **) TMP_ALLOC(len2*sizeof(ulong *));
    for (i = 0; i < len2; i++)
        exp_list[i] = exps + i*N;

    hind = (slong *) TMP_ALLOC(len2*sizeof(slong));
    for (i = 0; i < len2; i++)
        hind[i] = 1;

    /* start with no heap nodes and no exponent vectors in use */
    exp_next = 0;

    /* put (0, 0, exp2[0] + exp3[0]) on heap */
    x = chain + 0;
    x->i = 0;
    x->j = 0;
    x->next = NULL;

    heap[1].next = x;
    heap[1].exp = exp_list[exp_next++];

    if (bits <= FLINT_BITS)
        mpoly_monomial_add(heap[1].exp, exp2, exp3, N);
    else
        mpoly_monomial_add_mp(heap[1].exp, exp2, exp3, N);

    hind[0] = 2*1 + 0;

    len1 = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        _nmod_mpoly_fit_length(&p1, &e1, alloc, len1 + 1, N);

        mpoly_monomial_set(e1 + len1*N, exp, N);

        acc0 = acc1 = acc2 = 0;
        do
        {
            exp_list[--exp_next] = heap[1].exp;

            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);

            hind[x->i] |= WORD(1);
            Q[Q_len++] = x->i;
            Q[Q_len++] = x->j;
            umul_ppmm(pp1, pp0, coeff2[x->i], coeff3[x->j]);
            add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);

            while ((x = x->next) != NULL)
            {
                hind[x->i] |= WORD(1);
                Q[Q_len++] = x->i;
                Q[Q_len++] = x->j;
                umul_ppmm(pp1, pp0, coeff2[x->i], coeff3[x->j]);
                add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);
            }
        } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

        NMOD_RED3(p1[len1], acc2, acc1, acc0, fctx->mod);
        len1 += (p1[len1] != 0);

        while (Q_len > 0)
        {
            /* take node from store */
            j = Q[--Q_len];
            i = Q[--Q_len];

            /* should we go right? */
            if (  (i + 1 < len2)
               && (hind[i + 1] == 2*j + 1)
               )
            {
                x = chain + i + 1;
                x->i = i + 1;
                x->j = j;
                x->next = NULL;

                hind[x->i] = 2*(x->j+1) + 0;

                if (bits <= FLINT_BITS)
                    mpoly_monomial_add(exp_list[exp_next], exp2 + x->i*N, exp3 + x->j*N, N);
                else
                    mpoly_monomial_add_mp(exp_list[exp_next], exp2 + x->i*N, exp3 + x->j*N, N);

                if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                      &next_loc, &heap_len, N, cmpmask))
                    exp_next--;
            }

            /* should we go up? */
            if (  (j + 1 < len3)
               && ((hind[i] & 1) == 1)
               && (  (i == 0)
                  || (hind[i - 1] >= 2*(j + 2) + 1)
                  )
               )
            {
                x = chain + i;
                x->i = i;
                x->j = j + 1;
                x->next = NULL;

                hind[x->i] = 2*(x->j+1) + 0;

                if (bits <= FLINT_BITS)
                    mpoly_monomial_add(exp_list[exp_next], exp2 + x->i*N, exp3 + x->j*N, N);
                else
                    mpoly_monomial_add_mp(exp_list[exp_next], exp2 + x->i*N, exp3 + x->j*N, N);

                if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                      &next_loc, &heap_len, N, cmpmask))
                    exp_next--;
            }
        }
    }

    (* coeff1) = p1;
    (* exp1) = e1;

    TMP_END;

    return len1;
}

void nmod_mpoly_mul_johnson(nmod_mpoly_t poly1, const nmod_mpoly_t poly2,
                          const nmod_mpoly_t poly3, const nmod_mpoly_ctx_t ctx)
{
    slong i, N, len1 = 0;
    mp_bitcnt_t exp_bits;
    fmpz * max_fields2, * max_fields3;
    ulong * cmpmask;
    ulong * exp2 = poly2->exps, * exp3 = poly3->exps;
    int free2 = 0, free3 = 0;
    TMP_INIT;

    if (poly2->length == 0 || poly3->length == 0)
    {
        nmod_mpoly_zero(poly1, ctx);
        return;
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
    _fmpz_vec_add(max_fields2, max_fields2, max_fields3, ctx->minfo->nfields);

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

    N = mpoly_words_per_exp(exp_bits, ctx->minfo);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, exp_bits, ctx->minfo);

    /* ensure input exponents are packed into same sized fields as output */
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

    /* deal with aliasing and do multiplication */
    if (poly1 == poly2 || poly1 == poly3)
    {
        nmod_mpoly_t temp;

        nmod_mpoly_init2(temp, poly2->length + poly3->length - 1, ctx);
        nmod_mpoly_fit_bits(temp, exp_bits, ctx);
        temp->bits = exp_bits;

        if (poly2->length >= poly3->length)
        {
            len1 = _nmod_mpoly_mul_johnson(&temp->coeffs, &temp->exps, &temp->alloc,
                                      poly3->coeffs, exp3, poly3->length,
                                      poly2->coeffs, exp2, poly2->length,
                                          exp_bits, N, cmpmask, ctx->ffinfo);
        } else
        {
            len1 = _nmod_mpoly_mul_johnson(&temp->coeffs, &temp->exps, &temp->alloc,
                                      poly2->coeffs, exp2, poly2->length,
                                      poly3->coeffs, exp3, poly3->length,
                                          exp_bits, N, cmpmask, ctx->ffinfo);
        }

        nmod_mpoly_swap(temp, poly1, ctx);
        nmod_mpoly_clear(temp, ctx);
   } else
   {
        nmod_mpoly_fit_length(poly1, poly2->length + poly3->length - 1, ctx);
        nmod_mpoly_fit_bits(poly1, exp_bits, ctx);
        poly1->bits = exp_bits;

        if (poly2->length > poly3->length)
        {
            len1 = _nmod_mpoly_mul_johnson(&poly1->coeffs, &poly1->exps, &poly1->alloc,
                                      poly3->coeffs, exp3, poly3->length,
                                      poly2->coeffs, exp2, poly2->length,
                                          exp_bits, N, cmpmask, ctx->ffinfo);
        } else
        {
            len1 = _nmod_mpoly_mul_johnson(&poly1->coeffs, &poly1->exps, &poly1->alloc,
                                      poly2->coeffs, exp2, poly2->length,
                                      poly3->coeffs, exp3, poly3->length,
                                          exp_bits, N, cmpmask, ctx->ffinfo);
        }
   }

   if (free2)
      flint_free(exp2);

   if (free3)
      flint_free(exp3);

   _nmod_mpoly_set_length(poly1, len1, ctx);

   TMP_END;
}
