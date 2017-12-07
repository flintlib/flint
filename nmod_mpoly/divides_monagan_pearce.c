/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

slong _nmod_mpoly_divides_monagan_pearce1(
                     mp_limb_t ** coeff1,      ulong ** exp1, slong * alloc,
                const mp_limb_t * coeff2, const ulong * exp2, slong len2,
                const mp_limb_t * coeff3, const ulong * exp3, slong len3,
                              slong bits, ulong maskhi, const nmodf_ctx_t fctx)
{
    int lt_divides;
    slong i, j, k, s;
    slong next_loc, heap_len;
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    mp_limb_t * p1 = * coeff1;
    ulong * e1 = * exp1;
    slong * hind;
    ulong mask, exp, maxexp = exp2[len2 - 1];
    mp_limb_t lc_minus_inv, acc0, acc1, acc2, pp1, pp0;
    TMP_INIT;

    TMP_START;

    /* alloc array of heap nodes which can be chained together */
    next_loc = len3 + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap1_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap1_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(len3*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*len3*sizeof(mpoly_heap_t *));

    /* space for flagged heap indicies */
    hind = (slong *) TMP_ALLOC(len3*sizeof(slong));
    for (i = 0; i < len3; i++)
        hind[i] = 1;

    /* mask with high bit set in each field of exponent vector */
    mask = 0;
    for (i = 0; i < FLINT_BITS/bits; i++)
        mask = (mask << bits) + (UWORD(1) << (bits - 1));

    /* output poly index starts at -1, will be immediately updated to 0 */
    k = -WORD(1);

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
    lc_minus_inv = fctx->mod.n - nmod_inv(coeff3[0], fctx->mod);

    while (heap_len > 1)
    {
        exp = heap[1].exp;

        if (mpoly_monomial_overflows1(exp, mask))
            goto not_exact_division;

        k++;
        _nmod_mpoly_fit_length(&p1, &e1, alloc, k + 1, 1);

        lt_divides = mpoly_monomial_divides1(e1 + k, exp, exp3[0], mask);

        acc0 = acc1 = acc2 = 0;
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
                    add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), WORD(0), fctx->mod.n - coeff2[x->j]);
                } else
                {
                    umul_ppmm(pp1, pp0, coeff3[x->i], p1[x->j]);
                    add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);
                }
            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && heap[1].exp == exp);

        NMOD_RED3(p1[k], acc2, acc1, acc0, fctx->mod);
        p1[k] = nmod_mul(p1[k], lc_minus_inv, fctx->mod);

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
                    _mpoly_heap_insert1(heap, exp3[x->i] + e1[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
                }
                /* should we go up? */
                if (j + 1 == k)
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
                    _mpoly_heap_insert1(heap, exp3[x->i] + e1[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
                }
            }
        }

        if (p1[k] == 0)
        {
            k--;
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
            x->j = k;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;
            _mpoly_heap_insert1(heap, exp3[x->i] + e1[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
        }
        s = 1;
    }

    k++;

cleanup:

    * coeff1 = p1;
    * exp1 = e1;

    TMP_END;

    return k;

not_exact_division:
    k = 0;
    goto cleanup;
}


slong _nmod_mpoly_divides_monagan_pearce(
                     mp_limb_t ** coeff1,      ulong ** exp1, slong * alloc,
                const mp_limb_t * coeff2, const ulong * exp2, slong len2,
                const mp_limb_t * coeff3, const ulong * exp3, slong len3,
       slong bits, slong N, ulong maskhi, ulong masklo, const nmodf_ctx_t fctx)
{
    int lt_divides;
    slong i, j, k, s;
    slong next_loc, heap_len;
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    mp_limb_t * p1 = * coeff1;
    ulong * e1 = * exp1;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    mp_limb_t lc_minus_inv, acc0, acc1, acc2, pp1, pp0;
    ulong mask;
    slong * hind;
    TMP_INIT;

    /* if exponent vectors are all one word, call specialised version */
    if (N == 1)
        return _nmod_mpoly_divides_monagan_pearce1(coeff1, exp1, alloc,
                   coeff2, exp2, len2, coeff3, exp3, len3, bits, maskhi, fctx);

    TMP_START;

    /* alloc array of heap nodes which can be chained together */
    next_loc = len3 + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(len3*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*len3*sizeof(mpoly_heap_t *));

    /* array of exponent vectors, each of "N" words */
    exps = (ulong *) TMP_ALLOC(len3*N*sizeof(ulong));
    /* list of pointers to available exponent vectors */
    exp_list = (ulong **) TMP_ALLOC(len3*sizeof(ulong *));
    /* set up list of available exponent vectors */
    exp_next = 0;
    for (i = 0; i < len3; i++)
        exp_list[i] = exps + i*N;

    /* space for flagged heap indicies */
    hind = (slong *) TMP_ALLOC(len3*sizeof(slong));
    for (i = 0; i < len3; i++)
        hind[i] = 1;

    /* mask with high bit set in each word of each field of exponent vector */
    mask = 0;
    for (i = 0; i < FLINT_BITS/bits; i++)
        mask = (mask << bits) + (UWORD(1) << (bits - 1));

    /* output poly index starts at -1, will be immediately updated to 0 */
    k = -WORD(1);

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
    lc_minus_inv = fctx->mod.n - nmod_inv(coeff3[0], fctx->mod);

    while (heap_len > 1)
    {
        exp = heap[1].exp;

        if (mpoly_monomial_overflows(exp, N, mask))
            goto not_exact_division;
      
        k++;
        _nmod_mpoly_fit_length(&p1, &e1, alloc, k + 1, N);

        lt_divides = mpoly_monomial_divides(e1 + k*N, exp, exp3, N, mask);

        acc0 = acc1 = acc2 = 0;
        do
        {
            exp_list[--exp_next] = heap[1].exp;
            x = _mpoly_heap_pop(heap, &heap_len, N, maskhi, masklo);
            do
            {
                *store++ = x->i;
                *store++ = x->j;
                if (x->i != -WORD(1))
                    hind[x->i] |= WORD(1);

                if (x->i == -WORD(1))
                {
                    add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), WORD(0), fctx->mod.n - coeff2[x->j]);
                } else
                {
                    umul_ppmm(pp1, pp0, coeff3[x->i], p1[x->j]);
                    add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);                    
                }
            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

        NMOD_RED3(p1[k], acc2, acc1, acc0, fctx->mod);
        p1[k] = nmod_mul(p1[k], lc_minus_inv, fctx->mod);

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
                                      &next_loc, &heap_len, N, maskhi, masklo))
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
                    mpoly_monomial_add(exp_list[exp_next], exp3 + x->i*N,
                                                           e1   + x->j*N, N);
                    if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                      &next_loc, &heap_len, N, maskhi, masklo))
                        exp_next--;
                }
                /* should we go up? */
                if (j + 1 == k)
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
                    mpoly_monomial_add(exp_list[exp_next], exp3 + x->i*N,
                                                           e1   + x->j*N, N);
                    if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                      &next_loc, &heap_len, N, maskhi, masklo))
                        exp_next--;
                }
            }
        }

        if (p1[k] == 0)
        {
            k--;
            continue;
        }

        if (!lt_divides ||
                mpoly_monomial_gt(exp, exp2 + N*(len2 - 1), N, maskhi, masklo))
            goto not_exact_division;

        if (s > 1)
        {
            i = 1;
            x = chain + i;
            x->i = i;
            x->j = k;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;
            mpoly_monomial_add(exp_list[exp_next], exp3 + x->i*N,
                                                   e1   + x->j*N, N);
            if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                  &next_loc, &heap_len, N, maskhi, masklo))
                exp_next--;
        }
        s = 1;      
    }

    k++;

cleanup:

    * coeff1 = p1;
    * exp1 = e1;

    TMP_END;

    return k;

not_exact_division:
    k = 0;
    goto cleanup;
}

/* return 1 if quotient is exact */
int nmod_mpoly_divides_monagan_pearce(nmod_mpoly_t poly1,
                  const nmod_mpoly_t poly2, const nmod_mpoly_t poly3,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i, bits, exp_bits, N, len = 0;
    ulong max, * max_fields2, * max_fields3;
    ulong maskhi, masklo;
    ulong * exp2 = poly2->exps, * exp3 = poly3->exps, * expq;
    int free2 = 0, free3 = 0;
    ulong mask = 0;
    TMP_INIT;

    if (poly3->length == 0)
        flint_throw(FLINT_DIVZERO, "Divide by zero in nmod_mpoly_divides_monagan_pearce");

    if (poly2->length == 0)
    {
        nmod_mpoly_zero(poly1, ctx);
        return 1;
    }

    TMP_START;

    max_fields2 = (ulong *) TMP_ALLOC(ctx->n*sizeof(ulong));
    max_fields3 = (ulong *) TMP_ALLOC(ctx->n*sizeof(ulong));
    mpoly_max_fields_ui(max_fields2, poly2->exps, poly2->length,
                                                      poly2->bits, ctx->minfo);
    mpoly_max_fields_ui(max_fields3, poly3->exps, poly3->length,
                                                      poly3->bits, ctx->minfo);
    max = 0;
    for (i = 0; i < ctx->n; i++)
    {
        if (max_fields2[i] > max)
            max = max_fields2[i];

        /* cannot be exact division if poly2 degrees less than those of poly3 */
        if (max_fields2[i] < max_fields3[i])
        {
            len = 0;
            goto cleanup;
        }
    }

    /* compute number of bits required for exponent fields */
    bits = FLINT_BIT_COUNT(max);
    exp_bits = FLINT_MAX(WORD(8), bits + 1); /* extra bit required for signs */
    exp_bits = FLINT_MAX(exp_bits, poly2->bits);
    exp_bits = FLINT_MAX(exp_bits, poly3->bits);
    exp_bits = mpoly_optimize_bits(exp_bits, ctx->n);

    masks_from_bits_ord(maskhi, masklo, exp_bits, ctx->ord);
    N = words_per_exp(ctx->n, exp_bits);

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
        mpoly_unpack_monomials(exp2, exp_bits, poly2->exps, poly2->bits,
                                                        poly2->length, ctx->n);
    }

    if (exp_bits > poly3->bits)
    {
        free3 = 1;
        exp3 = (ulong *) flint_malloc(N*poly3->length*sizeof(ulong));
        mpoly_unpack_monomials(exp3, exp_bits, poly3->exps, poly3->bits,
                                                        poly3->length, ctx->n);
    }

    /* mask with high bit of each exponent vector field set */
    for (i = 0; i < FLINT_BITS/exp_bits; i++)
        mask = (mask << exp_bits) + (UWORD(1) << (exp_bits - 1));

    /* check leading monomial divides exactly */
    if (!mpoly_monomial_divides(expq, exp2, exp3, N, mask))
    {
        len = 0;
        goto cleanup;
    }

    /* deal with aliasing and divide polynomials */
    if (poly1 == poly2 || poly1 == poly3)
    {
      nmod_mpoly_t temp;

      nmod_mpoly_init2(temp, poly2->length/poly3->length + 1, ctx);
      nmod_mpoly_fit_bits(temp, exp_bits, ctx);
      temp->bits = exp_bits;

      len = _nmod_mpoly_divides_monagan_pearce(&temp->coeffs, &temp->exps,
                            &temp->alloc, poly2->coeffs, exp2, poly2->length,
                              poly3->coeffs, exp3, poly3->length, exp_bits, N,
                                                  maskhi, masklo, ctx->ffinfo);

      nmod_mpoly_swap(temp, poly1, ctx);

      nmod_mpoly_clear(temp, ctx);
   } else
   {
      nmod_mpoly_fit_length(poly1, poly2->length/poly3->length + 1, ctx);
      nmod_mpoly_fit_bits(poly1, exp_bits, ctx);
      poly1->bits = exp_bits;

      len = _nmod_mpoly_divides_monagan_pearce(&poly1->coeffs, &poly1->exps,
                            &poly1->alloc, poly2->coeffs, exp2, poly2->length,
                              poly3->coeffs, exp3, poly3->length, exp_bits, N,
                                                  maskhi, masklo, ctx->ffinfo);
   }

cleanup:

   _nmod_mpoly_set_length(poly1, len, ctx);

   if (free2)
      flint_free(exp2);

   if (free3)
      flint_free(exp3);

   TMP_END;

   /* division is exact if len is nonzero */
   return (len != 0);
}
