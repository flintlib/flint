/*
    Copyright (C) 2016 William Hart
    Copyright (C) 2017-2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void fmpz_mpoly_div(fmpz_mpoly_t Q, const fmpz_mpoly_t A,
                              const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)
{
    /* TODO !!! */
    fmpz_mpoly_div_monagan_pearce(Q, A, B, ctx);
}

/*
   Set polyq to the quotient poly2 by poly3 discarding remainder (with notional
   remainder coeffs reduced modulo the leading coeff of poly3), and return
   the length of the quotient. This version of the function assumes the
   exponent vectors all fit in a single word. The exponent vectors are
   assumed to have fields with the given number of bits. Assumes input polys
   are nonzero. Implements "Polynomial division using dynamic arrays, heaps
   and packed exponents" by Michael Monagan and Roman Pearce [1]. The word
   "maxhi" is set to a mask for the degree field of the exponent vector (in
   the case of a degree ordering) to facilitate the reverse ordering in
   degrevlex. Division is from left to right with a heap with largest
   exponent at the head. Quotient poly is written in order.
   [1] http://www.cecm.sfu.ca/~rpearcea/sdmp/sdmp_paper.pdf 
*/
slong _fmpz_mpoly_div_monagan_pearce1(fmpz ** polyq, ulong ** expq,
                slong * allocq, const fmpz * poly2, const ulong * exp2,
            slong len2, const fmpz * poly3, const ulong * exp3, slong len3,
                                                      slong bits, ulong maskhi)
{
    slong i, j, q_len, s;
    slong next_loc, heap_len = 2;
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    fmpz * q_coeff = *polyq;
    ulong * q_exp = *expq;
    slong * hind;
    ulong mask, exp;
    fmpz_t r, acc_lg;
    ulong acc_sm[3];
    int lt_divides, small;
    slong bits2, bits3;
    ulong lc_norm = 0, lc_abs = 0, lc_sign = 0, lc_n = 0, lc_i = 0;
    TMP_INIT;

    TMP_START;

    fmpz_init(acc_lg);
    fmpz_init(r);

   /* whether intermediate computations q - a*b will fit in three words */
   bits2 = _fmpz_vec_max_bits(poly2, len2);
   bits3 = _fmpz_vec_max_bits(poly3, len3);
   /* allow one bit for sign, one bit for subtraction */
   small = FLINT_ABS(bits2) <= (FLINT_ABS(bits3) + FLINT_BIT_COUNT(len3) +
           SMALL_FMPZ_BITCOUNT_MAX) && FLINT_ABS(bits3) <= SMALL_FMPZ_BITCOUNT_MAX;

    /* alloc array of heap nodes which can be chained together */
    next_loc = len3 + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap1_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap1_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(len3*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*len3*sizeof(slong));

    /* space for flagged heap indices */
    hind = (slong *) TMP_ALLOC(len3*sizeof(slong));
    for (i = 0; i < len3; i++)
        hind[i] = 1;

    mask = mpoly_overflow_mask_sp(bits);

    q_len = WORD(0);

    /* see description of divisor heap division in paper */
    s = len3;
   
    /* insert (-1, 0, exp2[0]) into heap */
    x = chain + 0;
    x->i = -WORD(1);
    x->j = 0;
    x->next = NULL;
    HEAP_ASSIGN(heap[1], exp2[0], x);

    /* precompute leading coefficient info assuming "small" case */
    if (small)
    {
        lc_abs = FLINT_ABS(poly3[0]);
        lc_sign = FLINT_SIGN_EXT(poly3[0]);
        count_leading_zeros(lc_norm, lc_abs);
        lc_n = lc_abs << lc_norm;
        invert_limb(lc_i, lc_n);
    }

    while (heap_len > 1)
    {
        exp = heap[1].exp;

        if (mpoly_monomial_overflows1(exp, mask))
            goto exp_overflow;

        _fmpz_mpoly_fit_length(&q_coeff, &q_exp, allocq, q_len + 1, 1);
        lt_divides = mpoly_monomial_divides1(q_exp + q_len, exp, exp3[0], mask);

        /* take nodes from heap with exponent matching exp */

        if (!lt_divides)
        {
            /* optimization: coeff arithmetic not needed */

            if (mpoly_monomial_gt1(exp3[0], exp, maskhi))
            {
                /* optimization: no more quotient terms possible */
                goto cleanup;
            }

            do
            {
                x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
                do
                {
                    *store++ = x->i;
                    *store++ = x->j;
                    if (x->i != -WORD(1))
                        hind[x->i] |= WORD(1);

                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);
        }
        else if (small)
        {
            /* optimization: small coeff arithmetic, acc_sm used below */

            acc_sm[0] = acc_sm[1] = acc_sm[2] = 0;
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
                        _fmpz_mpoly_add_uiuiui_fmpz(acc_sm, poly2 + x->j);
                    else
                        _fmpz_mpoly_submul_uiuiui_fmpz(acc_sm, poly3[x->i], q_coeff[x->j]);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);
        }
        else
        {
            /* general coeff arithmetic */

            fmpz_zero(acc_lg);  
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
                        fmpz_add(acc_lg, acc_lg, poly2 + x->j);
                    else
                        fmpz_submul(acc_lg, poly3 + x->i, q_coeff + x->j);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);
        }

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

        /* try to divide accumulated term by leading term */
        if (!lt_divides)
            continue;

        if (small)
        {
            ulong d0, d1, ds = acc_sm[2];

            /* d1:d0 = abs(acc_sm[1:0]) assuming ds is sign extension of acc_sm[1] */
            sub_ddmmss(d1, d0, acc_sm[1]^ds, acc_sm[0]^ds, ds, ds);
            
            if ((acc_sm[0] | acc_sm[1] | acc_sm[2]) == 0)
                continue;

            if (ds == FLINT_SIGN_EXT(acc_sm[1]) && d1 < lc_abs)
            {
                ulong qq, rr, nhi, nlo;
                FLINT_ASSERT(0 < lc_norm && lc_norm < FLINT_BITS);
                nhi = (d1 << lc_norm) | (d0 >> (FLINT_BITS - lc_norm));
                nlo = d0 << lc_norm;
                udiv_qrnnd_preinv(qq, rr, nhi, nlo, lc_n, lc_i);
                (void) rr; /* silence compiler warning */

                if (qq == 0)
                    continue;

                if (qq <= COEFF_MAX)
                {
                    _fmpz_demote(q_coeff + q_len);
                    q_coeff[q_len] = qq;
                    if (ds != lc_sign)
                        q_coeff[q_len] = -q_coeff[q_len];
                }
                else
                {
                    small = 0;
                    fmpz_set_ui(q_coeff + q_len, qq);
                    if (ds != lc_sign)
                        fmpz_neg(q_coeff + q_len, q_coeff + q_len);
                }
            }
            else
            {
                small = 0;
                fmpz_set_signed_uiuiui(acc_lg, acc_sm[2], acc_sm[1], acc_sm[0]);
                goto large_lt_divides;
            }
        }
        else
        {
            if (fmpz_is_zero(acc_lg))
                continue;

large_lt_divides:

            fmpz_fdiv_qr(q_coeff + q_len, r, acc_lg, poly3 + 0);
            if (fmpz_is_zero(q_coeff + q_len))
                continue;
        }

        /* put newly generated quotient term back into the heap if necessary */
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


cleanup:

    fmpz_clear(acc_lg);
    fmpz_clear(r);

    (*polyq) = q_coeff;
    (*expq) = q_exp;

    TMP_END;

    /* return quotient poly length */
    return q_len;

exp_overflow:
    for (i = 0; i < q_len; i++)
        _fmpz_demote(q_coeff + i);
    q_len = -WORD(1);
    goto cleanup;
}


slong _fmpz_mpoly_div_monagan_pearce(fmpz ** polyq,
           ulong ** expq, slong * allocq, const fmpz * poly2,
   const ulong * exp2, slong len2, const fmpz * poly3, const ulong * exp3, 
                   slong len3, slong bits, slong N, const ulong * cmpmask)
{
    slong i, j, q_len, s;
    slong next_loc;
    slong heap_len = 2; /* heap zero index unused */
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    fmpz * q_coeff = *polyq;
    ulong * q_exp = *expq;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    ulong mask;
    fmpz_t r, acc_lg;
    ulong acc_sm[3];
    slong * hind;
    int lt_divides, small;
    slong bits2, bits3;
    ulong lc_norm = 0, lc_abs = 0, lc_sign = 0, lc_n = 0, lc_i = 0;
    TMP_INIT;

    if (N == 1)
        return _fmpz_mpoly_div_monagan_pearce1(polyq, expq, allocq,
                       poly2, exp2, len2, poly3, exp3, len3, bits, cmpmask[0]);

    TMP_START;

    fmpz_init(acc_lg);
    fmpz_init(r);

    /* whether intermediate computations q - a*b will fit in three words */
    bits2 = _fmpz_vec_max_bits(poly2, len2);
    bits3 = _fmpz_vec_max_bits(poly3, len3);
    /* allow one bit for sign, one bit for subtraction */
    small = FLINT_ABS(bits2) <= (FLINT_ABS(bits3) + FLINT_BIT_COUNT(len3) + SMALL_FMPZ_BITCOUNT_MAX)
         && FLINT_ABS(bits3) <= SMALL_FMPZ_BITCOUNT_MAX;


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

    /* space for flagged heap indices */
    hind = (slong *) TMP_ALLOC(len3*sizeof(slong));
    for (i = 0; i < len3; i++)
        hind[i] = 1;

    mask = bits <= FLINT_BITS ? mpoly_overflow_mask_sp(bits) : 0;

    q_len = WORD(0);
   
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

    /* precompute leading coefficient info in "small" case */
    if (small)
    {
        lc_abs = FLINT_ABS(poly3[0]);
        lc_sign = FLINT_SIGN_EXT(poly3[0]);
        count_leading_zeros(lc_norm, lc_abs);
        lc_n = lc_abs << lc_norm;
        invert_limb(lc_i, lc_n);
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

        _fmpz_mpoly_fit_length(&q_coeff, &q_exp, allocq, q_len + 1, N);

        if (bits <= FLINT_BITS)
            lt_divides = mpoly_monomial_divides(q_exp + q_len*N, exp, exp3, N, mask);
        else
            lt_divides = mpoly_monomial_divides_mp(q_exp + q_len*N, exp, exp3, N, bits);

        /* take nodes from heap with exponent matching exp */

        if (!lt_divides)
        {
            /* optimization: coeff arithmetic not needed */

            if (mpoly_monomial_gt(exp3 + 0, exp, N, cmpmask))
            {
                /* optimization: no more quotient terms possible */
                goto cleanup;
            }

            do
            {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
                do
                {
                    *store++ = x->i;
                    *store++ = x->j;
                    if (x->i != -WORD(1))
                        hind[x->i] |= WORD(1);

                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));
        }
        else if (small)
        {
            /* optimization: small coeff arithmetic, acc_sm used below */

            acc_sm[0] = acc_sm[1] = acc_sm[2] = 0;
            do
            {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
                do
                {
                    *store++ = x->i;
                    *store++ = x->j;
                    if (x->i != -WORD(1))
                        hind[x->i] |= WORD(1);

                    if (x->i == -WORD(1))
                        _fmpz_mpoly_add_uiuiui_fmpz(acc_sm, poly2 + x->j);
                    else
                        _fmpz_mpoly_submul_uiuiui_fmpz(acc_sm, poly3[x->i], q_coeff[x->j]);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));
        }
        else
        {
            /* general coeff arithmetic */

            fmpz_zero(acc_lg);
            do
            {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
                do
                {
                    *store++ = x->i;
                    *store++ = x->j;
                    if (x->i != -WORD(1))
                        hind[x->i] |= WORD(1);

                    if (x->i == -WORD(1))
                        fmpz_add(acc_lg, acc_lg, poly2 + x->j);
                    else
                        fmpz_submul(acc_lg, poly3 + x->i, q_coeff + x->j);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));
        }

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

        /* try to divide accumulated term by leading term */
        if (!lt_divides)
        {
            continue;
        }
        if (small)
        {
            ulong d0, d1, ds = acc_sm[2];

            /* d1:d0 = abs(acc_sm[1:0]) assuming ds is sign extension of acc_sm[1] */
            sub_ddmmss(d1, d0, acc_sm[1]^ds, acc_sm[0]^ds, ds, ds);
            
            if ((acc_sm[0] | acc_sm[1] | acc_sm[2]) == 0)
            {
                continue;
            }

            if (ds == FLINT_SIGN_EXT(acc_sm[1]) && d1 < lc_abs)
            {
                ulong qq, rr, nhi, nlo;
                FLINT_ASSERT(0 < lc_norm && lc_norm < FLINT_BITS);
                nhi = (d1 << lc_norm) | (d0 >> (FLINT_BITS - lc_norm));
                nlo = d0 << lc_norm;
                udiv_qrnnd_preinv(qq, rr, nhi, nlo, lc_n, lc_i);
                (void) rr;

                if (qq == 0)
                    continue;

                if (qq <= COEFF_MAX)
                {
                    _fmpz_demote(q_coeff + q_len);
                    q_coeff[q_len] = qq;
                    if (ds != lc_sign)
                        q_coeff[q_len] = -q_coeff[q_len];
                }
                else
                {
                    small = 0;
                    fmpz_set_ui(q_coeff + q_len, qq);
                    if (ds != lc_sign)
                        fmpz_neg(q_coeff + q_len, q_coeff + q_len);
                }
            }
            else
            {
                small = 0;
                fmpz_set_signed_uiuiui(acc_lg, acc_sm[2], acc_sm[1], acc_sm[0]);
                goto large_lt_divides;
            }
        }
        else
        {
            if (fmpz_is_zero(acc_lg))
            {
                continue;
            }
large_lt_divides:
            fmpz_fdiv_qr(q_coeff + q_len, r, acc_lg, poly3 + 0);
            if (fmpz_is_zero(q_coeff + q_len))
            {
                continue;
            }
        }

        /* put newly generated quotient term back into the heap if necessary */
        if (s > 1)
        {
            i = 1;
            x = chain + i;
            x->i = i;
            x->j = q_len;
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
        s = 1;
        q_len++;
    }


cleanup:

    fmpz_clear(acc_lg);
    fmpz_clear(r);

    (*polyq) = q_coeff;
    (*expq) = q_exp;

    TMP_END;

    /* return quotient poly length */
    return q_len;

exp_overflow:
    for (i = 0; i < q_len; i++)
        _fmpz_demote(q_coeff + i);
    q_len = -WORD(1);
    goto cleanup;

}

void fmpz_mpoly_div_monagan_pearce(fmpz_mpoly_t q, const fmpz_mpoly_t poly2, 
                          const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)
{
    slong exp_bits, N, lenq = 0;
    ulong * exp2 = poly2->exps, * exp3 = poly3->exps;
    ulong * cmpmask;
    int free2 = 0, free3 = 0;
    fmpz_mpoly_t temp1;
    fmpz_mpoly_struct * tq;

   /* check divisor is nonzero */
   if (poly3->length == 0)
      flint_throw(FLINT_DIVZERO, "Divide by zero in fmpz_mpoly_div_monagan_pearce");

   /* dividend zero, write out quotient */
   if (poly2->length == 0)
   {
      fmpz_mpoly_zero(q, ctx);
      return;
   }

   /* maximum bits in quotient exps and inputs is max for poly2 and poly3 */
   exp_bits = FLINT_MAX(poly2->bits, poly3->bits);

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
      fmpz_mpoly_zero(q, ctx);
      goto cleanup3;
   }

   /* take care of aliasing */
   if (q == poly2 || q == poly3)
   {
      fmpz_mpoly_init2(temp1, poly2->length/poly3->length + 1, ctx);
      fmpz_mpoly_fit_bits(temp1, exp_bits, ctx);
      temp1->bits = exp_bits;
      tq = temp1;
   } else
   {
      fmpz_mpoly_fit_length(q, poly2->length/poly3->length + 1, ctx);
      fmpz_mpoly_fit_bits(q, exp_bits, ctx);
      q->bits = exp_bits;
      tq = q;
   }

   /* do division without remainder */
   while ((lenq = _fmpz_mpoly_div_monagan_pearce(&tq->coeffs, &tq->exps,
                         &tq->alloc, poly2->coeffs, exp2, poly2->length, 
                                     poly3->coeffs, exp3, poly3->length,
                                      exp_bits, N, cmpmask)) == -WORD(1)  )
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

      fmpz_mpoly_fit_bits(tq, exp_bits, ctx);
      tq->bits = exp_bits;
   }

   /* take care of aliasing */
   if (q == poly2 || q == poly3)
   {
      fmpz_mpoly_swap(temp1, q, ctx);
      fmpz_mpoly_clear(temp1, ctx);
   } 

   _fmpz_mpoly_set_length(q, lenq, ctx);

cleanup3:

   if (free2)
      flint_free(exp2);

   if (free3)
      flint_free(exp3);

   flint_free(cmpmask);
}

void fmpz_mpoly_divrem(fmpz_mpoly_t Q, fmpz_mpoly_t R, const fmpz_mpoly_t A,
                              const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)
{
    /* TODO !!! */
    fmpz_mpoly_divrem_monagan_pearce(Q, R, A, B, ctx);
}

/* improve locality */
#define BLOCK 128
#define MAX_ARRAY_SIZE (WORD(300000))

/*
   Set polyq to the quotient and polyr to the remainder of poly2 divided
   by poly3, and return the length of the quotient. The polynomials have
   their exponents tightly packed, with mixed bases equal to the largest
   exponent for each variable, e.g. the input polys have exponents of the
   form a_0 + a_1*b1 + a_2*b_2*b_2 + .... where b_0, b_1, b_2, etc, are
   the bases, which are equal to the largest possible exponent for
   each of the respective variables in the exponent. The dividend poly3
   is assumed to be nonzero. There are assumed to be "num" variables and
   the bases b_i are passed in the array "mults". The function reallocates
   its output. The quotient and remainder are written out in reverse order.
   The quotient and remainder poly are not assumed to be zero on input.
   The quotient and remainder terms are appended to the existing terms in
   those polys. 
*/
slong _fmpz_mpoly_divrem_array_tight(slong * lenr,
 fmpz ** polyq, ulong ** expq, slong * allocq, slong len0,
       fmpz ** polyr, ulong ** expr, slong * allocr, slong len1,
                  const fmpz * poly2, const ulong * exp2, slong len2,
                        const fmpz * poly3, const ulong * exp3, slong len3,
                                                      slong * mults, slong num)
{
   slong i, j, q, r, prod, bits1, bits2, bits3, k = len0, l = len1;
   slong max3 = (slong) exp3[0]; /* largest exponent in poly3 */
   slong * prods;
   fmpz c3 = poly3[0];
   /* abs val of leading coeff of poly3 */
   ulong u3 = ((ulong) FLINT_ABS(c3)) >> 1;
   fmpz * p1 = *polyq, * p2 = *polyr;
   ulong * e1 = *expq, * e2 = *expr;
   int small;
   TMP_INIT;

   TMP_START;

   prods = (slong *) TMP_ALLOC((num + 1)*sizeof(slong));

   /*
      compute products 1, b0, b0*b1, b0*b1*b2 ...
      from list of bases b0, b1, b2, ...
   */
   prods[0] = 1;
   for (i = 1; i <= num; i++)
      prods[i] = mults[i - 1]*prods[i - 1];

   prod = prods[num];

   /* compute bound on poly2 - q*poly3 assuming quotient remains small */
   bits2 = _fmpz_vec_max_bits(poly2, len2);
   bits3 = _fmpz_vec_max_bits(poly3, len3);
   /* we assume a bound of SMALL_FMPZ_BITCOUNT_MAX for coefficients of the quotient */
   bits1 = FLINT_ABS(bits3) + FLINT_BITS + FLINT_BIT_COUNT(len3) - 2;

   small = FLINT_ABS(bits2) <= bits1 && FLINT_ABS(bits3) <= SMALL_FMPZ_BITCOUNT_MAX;
   bits1 += 2; /* incr. so poly2 - q*poly3 doesn't overflow and for sign */
 
   /* input coeffs small and intermediate computations fit two words */
   if (small && bits1 <= 2*FLINT_BITS)
   {
      ulong * t2 = (ulong *) TMP_ALLOC(2*prod*sizeof(ulong));

      for (i = 0; i < 2*prod; i++)
         t2[i] = 0;

      /* poly2 to array format */
      _fmpz_mpoly_to_ulong_array2(t2, poly2, exp2, len2);

      /* for each term of poly2 array relevant to quotient */
      for (i = prod - 1; i >= max3; i--)
      {
         ulong * ptr = t2 + 2*i;
         ulong p[2];

         /* if coeff is nonzero */
         if (ptr[0] != 0 || ptr[1] != 0)
         {
            if (0 > (slong) ptr[1])
               mpn_neg(p, ptr, 2);
            else
               flint_mpn_copyi(p, ptr, 2);

            /* not exact quotient monomial, thus remainder monomial */
            if (!mpoly_monomial_divides_tight(i, max3, prods, num))
            {
               /* realloc remainder poly */
               _fmpz_mpoly_fit_length(&p2, &e2, allocr, l + 1, 1);

               /* set remainder coeff... */
               fmpz_set_signed_uiui(p2 + l, ptr[1], ptr[0]);

               /* ...and exponent */
               e2[l++] = i;
            } else /* monomials can be divided exactly */
            {
               /* check quotient won't overflow a word */
               if (u3 <= p[1] || (u3 == 0 && 0 > (slong) p[0])) /* quotient too large */
               {
                  for (j = len0; j < k; j++)
                     _fmpz_demote(p1 + j);
                  for (j = len1; j < l; j++)
                     _fmpz_demote(p2 + j);
                  k = len0;
                  l = len1;

                  goto big;
               }

               /* quotient and remainder of coeffs */
               sdiv_qrnnd(q, r, ptr[1], ptr[0], c3);

               /* check coefficient is small, else restart with multiprec code */
               if (COEFF_IS_MPZ(FLINT_ABS(q))) /* quotient too large */
               {
                  for (j = len0; j < k; j++)
                     _fmpz_demote(p1 + j);
                  for (j = len1; j < l; j++)
                     _fmpz_demote(p2 + j);
                  k = len0;
                  l = len1;

                  goto big;
               }

               /* check coeff quotient was exact */
               if (r != 0) /* not an exact division */
               {
                  /* realloc remainder poly */
                  _fmpz_mpoly_fit_length(&p2, &e2, allocr, l + 1, 1);

                  /* set remainder coeff... */
                  fmpz_set_si(p2 + l, (slong) r);  

                  /* ... and exponent */
                  e2[l++] = i;
               }

               if (q != 0)
               {
                  /* submul a - q*b */
                  _fmpz_mpoly_submul_array1_slong2_1(t2, q, i - max3,
                                                            poly3, exp3, len3);

                  /* realloc quotient poly */
                  _fmpz_mpoly_fit_length(&p1, &e1, allocq, k + 1, 1);
                  /* set quotient coeff and exponent */
                  fmpz_set_si(p1 + k, q);
                  e1[k++] = i - max3;
               }
            }         
         }
      }

      /* all remaining terms are remainder terms */
      for ( ; i >= 0; i--)
      {
         ulong * ptr = t2 + 2*i;

         /* if coeff nonzero */
         if (ptr[0] != 0 || ptr[1] != 0)  /* not an exact division */
         {
            /* realloc remainder poly */
            _fmpz_mpoly_fit_length(&p2, &e2, allocr, l + 1, 1);

            /* set remainder coeff... */
            fmpz_set_signed_uiui(p2 + l, ptr[1], ptr[0]);

            /* and exponent */
            e2[l++] = i;
         }
      }
   }

   /* not done, coeffs small and intermediate computations fit three words */
   if (k == len0 && l == len1 && small)
   {
      ulong * t2 = (ulong *) TMP_ALLOC(3*prod*sizeof(ulong));

      for (i = 0; i < 3*prod; i++)
         t2[i] = 0;

      /* poly2 to array format */
      _fmpz_mpoly_to_ulong_array(t2, poly2, exp2, len2);

      /* for each term of poly2 array relevant to exact quotient */
      for (i = prod - 1; i >= max3; i--)
      {
         ulong * ptr = t2 + 3*i;
         ulong p[3];

         /* if coeff is nonzero */
         if (ptr[0] != 0 || ptr[1] != 0 || ptr[2] != 0)
         {
            if (0 > (slong) ptr[2])
               mpn_neg(p, ptr, 3);
            else
               flint_mpn_copyi(p, ptr, 3);

            /* not exact quotient monomial, thus remainder monomial */
            if (!mpoly_monomial_divides_tight(i, max3, prods, num))
            {
               /* realloc remainder poly */
               _fmpz_mpoly_fit_length(&p2, &e2, allocr, l + 1, 1);

               /* set remainder coeff... */
               fmpz_set_signed_uiuiui(p2 + l, ptr[2], ptr[1], ptr[0]);

               /* ... and exponent */
               e2[l++] = i;
            } else /* monomials can be divided exact */
            {
               /* check quotient won't overflow a word */
               if (p[2] > 0 || u3 <= p[1] || (u3 == 0 && 0 > (slong) p[0])) /* quotient too large */
               {
                  for (j = len0; j < k; j++)
                     _fmpz_demote(p1 + j);
                  for (j = len1; j < l; j++)
                     _fmpz_demote(p2 + j);
                  k = len0;
                  l = len1;

                  goto big;
               }

               /* quotient and remainder of coeffs */
               sdiv_qrnnd(q, r, ptr[1], ptr[0], c3);

               /* check coefficient is small, else restart with multiprec code */
               if (COEFF_IS_MPZ(FLINT_ABS(q))) /* quotient too large */
               {
                  for (j = len0; j < k; j++)
                     _fmpz_demote(p1 + j);
                  for (j = len1; j < l; j++)
                     _fmpz_demote(p2 + j);
                  k = len0;
                  l = len1;

                  goto big;
               }

               /* check if coeff quotient was exact */
               if (r != 0) /* else remainder term */ 
               {
                  /* reallocate remainder poly */
                  _fmpz_mpoly_fit_length(&p2, &e2, allocr, l + 1, 1);

                  /* set remainder coeff... */
                  fmpz_set_si(p2 + l, (slong) r);  

                  /* and exponent */
                  e2[l++] = i;
               }

               /* if nonzero quotient */
               if (q != 0)
               {
                  /* submul a - q*b */
                  _fmpz_mpoly_submul_array1_slong_1(t2, q, i - max3,
                                                            poly3, exp3, len3);

                  /* realloc quotient poly */
                  _fmpz_mpoly_fit_length(&p1, &e1, allocq, k + 1, 1);
                  
                  /* set quotient coeff and exponent */
                  fmpz_set_si(p1 + k, q);

                  e1[k++] = i - max3;
               }
            }
         }
      }

      /* all remaining terms are remainder terms */
      for ( ; i >= 0; i--)
      {
         ulong * ptr = t2 + 3*i;

         /* if coeff nonzero */
         if (ptr[0] != 0 || ptr[1] != 0 || ptr[2] != 0) 
         {
            /* not an exact division */

            /* realloc remainder poly */
            _fmpz_mpoly_fit_length(&p2, &e2, allocr, l + 1, 1);

            /* set remainder coeff... */
            fmpz_set_signed_uiuiui(p2 + l, ptr[2], ptr[1], ptr[0]);

            /* ...and exponent */
            e2[l++] = i;
         }
      }
   }

big:

   /* if still not done, use multiprecision coeffs instead */
   if (k == len0 && l == len1)
   {
      fmpz * t2 = (fmpz *) TMP_ALLOC(prod*sizeof(fmpz));
      fmpz_t fq, fr;

      fmpz_init(fq);
      fmpz_init(fr);

      for (i = 0; i < prod; i++)
         fmpz_init(t2 + i);

      /* poly2 to array format */
      _fmpz_mpoly_to_fmpz_array(t2, poly2, exp2, len2);
      
      /* for each term of poly2 array relevant to exact quotient */
      for (i = prod - 1; i >= max3; i--)
      {
         /* if coeff is nonzero */
         if (!fmpz_is_zero(t2 + i))
         {
            /* not exact quotient monomial, thus remainder monomial */
            if (!mpoly_monomial_divides_tight(i, max3, prods, num))
            {
               /* realloc remainder poly */
               _fmpz_mpoly_fit_length(&p2, &e2, allocr, l + 1, 1);

               /* set remainder coeff... */
               fmpz_set(p2 + l, t2 + i);
               
               /* ... and exponent */
               e2[l++] = i;
            } else /* monomials can be divided exactly */
            {
               /* quotient and remainder of coeffs */
               fmpz_fdiv_qr(fq, fr, t2 + i, poly3);

               /* check if coeff quotient was exact */
               if (!fmpz_is_zero(fr)) /* else remainder term */
               {
                  /* realloc remainder poly */
                  _fmpz_mpoly_fit_length(&p2, &e2, allocr, l + 1, 1);

                  /* set remainder coeff... */
                  fmpz_set(p2 + l, fr);  

                  /* and exponent */
                  e2[l++] = i;
               }

               /* if nonzero quotient */
               if (!fmpz_is_zero(fq))
               {
                  /* submul a - q*b */
                  _fmpz_mpoly_submul_array1_fmpz_1(t2, fq, i - max3,
                                                            poly3, exp3, len3);

                  /* realloc quotient poly */
                  _fmpz_mpoly_fit_length(&p1, &e1, allocq, k + 1, 1);
            
                  /* set quotient coeff and exponent */
                  fmpz_set(p1 + k, fq);
                  e1[k++] = i - max3;
               }
            }
         }
      }

      /* all remaining terms are remainder terms */
      for ( ; i >= 0; i--)
      {
         /* if coeff nonzero */
         if (!fmpz_is_zero(t2 + i))
         {
            /* remainder */

            /* realloc remainder poly */
            _fmpz_mpoly_fit_length(&p2, &e2, allocr, l + 1, 1);

            /* set remainder coeff... */
            fmpz_set(p2 + l, t2 + i);

            /* ... and exponent */
            e2[l++] = i;
         }
      }

      fmpz_clear(fq);
      fmpz_clear(fr);

      for (i = 0; i < prod; i++)
         fmpz_clear(t2 + i);
   }

   (*polyq) = p1;
   (*expq) = e1;
   (*polyr) = p2;
   (*expr) = e2;

   /* set remainder poly length */
   (*lenr) = l - len1;

   TMP_END;

   /* return quotient poly length */
   return k - len0;
}

/*
   Use dense array division to set polyq, polyr to poly2/poly3 in num + 1
   variables, given a list of multipliers to tightly pack exponents and a
   number of bits for the fields of the exponents of the result, assuming
   no aliasing. classical exact division in main variable, array
   multiplication (submul) for multivariate coefficients in remaining num
   variables. The function reallocates its output and returns the length
   of the quotient poly. It is assumed that poly2 is not zero. The
   quotient and remainder are written in reverse order.
*/
slong _fmpz_mpoly_divrem_array_chunked(slong * lenr,
            fmpz ** polyq, ulong ** expq, slong * allocq,
                 fmpz ** polyr, ulong ** expr, slong * allocr, 
                const fmpz * poly2, const ulong * exp2, slong len2, 
        const fmpz * poly3, const ulong * exp3, slong len3, slong * mults, 
                                                         slong num, slong bits)
{
   slong i, j, k, l = 0, m, prod, len = 0, l1, l2, l3;
   slong bits1, bits2, bits3 = 0, tlen, talloc;
   slong shift = bits*(num);
   slong * i1, * i2, * i3, * n1, * n2, * n3, * prods;
   slong * b1, * b3, * maxb1, * maxb3, * max_exp1, * max_exp3;
   ulong * e2, * e3, * texp, * p2;
   fmpz * temp;
   int small;
   TMP_INIT;

   TMP_START;
   
   /*
      compute products 1, b0, b0*b1, b0*b1*b2 ...
      from list of bases b0, b1, b2, ...
   */
   prods = (slong *) TMP_ALLOC((num + 1)*sizeof(slong));

   prods[0] = 1;
   for (i = 0; i < num; i++)
      prods[i + 1] = prods[i]*mults[i];
   prod = prods[num];

   /* lengths of poly2, poly3 and polyq in chunks */
   l2 = 1 + (slong) (exp2[0] >> shift);
   l3 = 1 + (slong) (exp3[0] >> shift);

   l1 = FLINT_MAX(l2 - l3 + 1, 0);

   /* compute indices and lengths of coefficients of polys in main variable */

   i1 = (slong *) TMP_ALLOC(l1*sizeof(slong));
   n1 = (slong *) TMP_ALLOC(l1*sizeof(slong));
   i2 = (slong *) TMP_ALLOC(l2*sizeof(slong));
   n2 = (slong *) TMP_ALLOC(l2*sizeof(slong));
   i3 = (slong *) TMP_ALLOC(l3*sizeof(slong));
   n3 = (slong *) TMP_ALLOC(l3*sizeof(slong));
   b1 = (slong *) TMP_ALLOC(l1*sizeof(slong));
   maxb1 = (slong *) TMP_ALLOC(l1*sizeof(slong));
   b3 = (slong *) TMP_ALLOC(l3*sizeof(slong));
   maxb3 = (slong *) TMP_ALLOC(l3*sizeof(slong));
   max_exp1 = (slong *) TMP_ALLOC(l1*num*sizeof(slong));
   max_exp3 = (slong *) TMP_ALLOC(l3*num*sizeof(slong));

   mpoly_main_variable_terms1(i2, n2, exp2, l2, len2, num + 1, num + 1, bits);
   mpoly_main_variable_terms1(i3, n3, exp3, l3, len3, num + 1, num + 1, bits);

   /* work out max bits for each coeff and optimal bits */

   for (i = 0; i < l3; i++)
   {
      _fmpz_vec_sum_max_bits(&b3[i], &maxb3[i], poly3+i3[i], n3[i]);

      if (bits3 < maxb3[i])
         bits3 = maxb3[i];
   }

   /* pack input coefficients tightly */

   e2 = (ulong *) TMP_ALLOC(len2*sizeof(ulong));
   e3 = (ulong *) TMP_ALLOC(len3*sizeof(ulong));

   mpoly_pack_monomials_tight(e2, exp2, len2, mults, num, bits);
   mpoly_pack_monomials_tight(e3, exp3, len3, mults, num, bits);

   /* work out maximum exponents for each chunk */
   for (i = 0; i < l3; i++)
      mpoly_max_degrees_tight(max_exp3 + i*num, e3 + i3[i], n3[i], prods, num);
   
   /* bound poly2 coeffs and check input/output coeffs likely small */
   bits2 = _fmpz_vec_max_bits(poly2, len2);
   /* we assume a bound of SMALL_FMPZ_BITCOUNT_MAX for coefficients of the quotient */
   bits1 = FLINT_ABS(bits3) + FLINT_BITS + FLINT_BIT_COUNT(len3) - 2;

   small = FLINT_ABS(bits2) <= bits1 && FLINT_ABS(bits3) <= SMALL_FMPZ_BITCOUNT_MAX;

   /* alloc space for copy of coeff/chunk of poly2 */

   temp = (fmpz *) flint_calloc(n2[0] + 1, sizeof(fmpz));
   texp = (ulong *) flint_malloc((n2[0] + 1)*sizeof(ulong));
   talloc = n2[0] + 1; /* plus one so doubling always increases size */

   /* enough space for three words per coeff, even if only one or two needed */
   p2 = (ulong *) TMP_ALLOC(3*prod*sizeof(ulong));

   /* coefficients likely to be small */
   if (small)
   {
      /* for each chunk of poly2 */
      for (i = 0; i < l2; i++)
      {
         slong num1 = 0;
         bits1 = 0;

         /* if there are already quotient terms */
         if (i != 0)
         {
            /* compute bound on intermediate computations a - q*b */
            for (j = 0; j < i && j < l1; j++)
            {
               k = i - j;

               if (k < l3 && k >= 0)
               {
                  for (m = 0; m < num; m++)
                  {
                     if (max_exp1[j*num + m] + max_exp3[k*num + m] >= mults[m])
                     {
                        for (j = 0; j < len; j++)
                           _fmpz_demote((*polyq) + j);
                        for (j = 0; j < l; j++)
                           _fmpz_demote((*polyr) + j);
                        len = 0;
                        l = 0;
                        goto cleanup3;
                     }
                  }

                  bits1 = FLINT_MAX(bits1, FLINT_MIN(b1[j] + maxb3[k], maxb1[j] + b3[k]));
                  num1++;
               }
            }

            bits1 += FLINT_BIT_COUNT(num1);
            bits1 = FLINT_MAX(FLINT_ABS(bits2), bits1);

            bits1 += 2; /* bit for sign and so a - q*b doesn't overflow */
         } else
            bits1 = FLINT_ABS(bits2) + 1; /* extra bit for sign */

         /* intermediate computations fit in one word */
         if (bits1 <= FLINT_BITS)
         {
            for (j = 0; j < prod; j++)
               p2[j] = 0;

            /* convert relevant coeff/chunk of poly2 to array format */
            _fmpz_mpoly_to_ulong_array1(p2, poly2 + i2[i], e2 + i2[i], n2[i]);

            /* submuls */

            for (j = 0; j < i && j < l1; j++)
            {
               k = i - j;

               if (k < l3 && k >= 0)
               {
                  _fmpz_mpoly_submul_array1_slong1(p2, (*polyq) + i1[j],
                     (*expq) + i1[j], n1[j], poly3 + i3[k], e3 + i3[k], n3[k]);
               }
            }

            /* convert chunk from array format */
            tlen = _fmpz_mpoly_from_ulong_array1(&temp, &texp, &talloc, 
                                                      p2, mults, num, bits, 0);
         } else if (bits1 <= 2*FLINT_BITS) /* intermed comps fit two words */
         {
            for (j = 0; j < 2*prod; j++)
               p2[j] = 0;

            /* convert relevant coeff/chunk of poly2 to array format */
            _fmpz_mpoly_to_ulong_array2(p2, poly2 + i2[i], e2 + i2[i], n2[i]);

            /* submuls */

            for (j = 0; j < i && j < l1; j++)
            {
               k = i - j;

               if (k < l3 && k >= 0)
               {
                  _fmpz_mpoly_submul_array1_slong2(p2, (*polyq) + i1[j],
                     (*expq) + i1[j], n1[j], poly3 + i3[k], e3 + i3[k], n3[k]);
               }
            }

            /* convert chunk from array format */
            tlen = _fmpz_mpoly_from_ulong_array2(&temp, &texp, &talloc, 
                                                      p2, mults, num, bits, 0);
         } else /* intermed comps fit three words */
         {
            for (j = 0; j < 3*prod; j++)
               p2[j] = 0;

            /* convert relevant coeff/chunk of poly2 to array format */
            _fmpz_mpoly_to_ulong_array(p2, poly2 + i2[i], e2 + i2[i], n2[i]);

            /* submuls */

            for (j = 0; j < i && j < l1; j++)
            {
               k = i - j;

               if (k < l3 && k >= 0)
                  _fmpz_mpoly_submul_array1_slong(p2, (*polyq) + i1[j],
                     (*expq) + i1[j], n1[j], poly3 + i3[k], e3 + i3[k], n3[k]);
            }

            /* convert chunk from array format */
            tlen = _fmpz_mpoly_from_ulong_array(&temp, &texp, &talloc, 
                                                      p2, mults, num, bits, 0);
         }

         if (tlen != 0) /* nonzero coeff/chunk */
         {
            if (i < l1) /* potentially a quotient with remainder */
            {
               /* tightly pack chunk exponents */
               mpoly_pack_monomials_tight(texp, texp, tlen, mults, num, bits);

               /* set starting index for quotient chunk we are about to compute */
               i1[i] = len;

               /* compute quotient chunk and set length of quotient chunk */
	       n1[i] = _fmpz_mpoly_divrem_array_tight(lenr, polyq,
                 expq, allocq, len, polyr, expr, allocr, l, temp, texp,
                     tlen, poly3 + i3[0], e3 + i3[0], n3[0], mults, num);
																   
	       /* unpack remainder exponents */
	       mpoly_unpack_monomials_tight(*expr + l,
			                            *expr + l, *lenr, mults, num, bits);
														  
	       /* insert main variable */
	       for (j = 0; j < *lenr; j++)
	          (*expr)[l + j] += ((l2 - i - 1) << shift);
				  
	       /* work out maximum exponents for chunk */
               mpoly_max_degrees_tight(max_exp1 + i*num,
	                                   (*expq) + i1[i], n1[i], prods, num);
						 
	       /* check there were no exponent overflows */
	       for (m = 0; m < num; m++)
	       {
	          if (max_exp3[m] + max_exp1[i*num + m] >= mults[m])
                  {
                     for (j = 0; j < len; j++)
                        _fmpz_demote((*polyq) + j);
                     for (j = 0; j < l; j++)
                        _fmpz_demote((*polyr) + j);
                     len = 0;
		     l = 0;

	             goto cleanup3;
                  }				  
	       }

               /* check the quotient didn't have large coefficients */
               if (FLINT_ABS(_fmpz_vec_max_bits((*polyq) + len,
                                                      n1[i])) > SMALL_FMPZ_BITCOUNT_MAX)
               {
                  for (j = 0; j < len; j++)
                     _fmpz_demote((*polyq) + j);
                  for (j = 0; j < l; j++)
                     _fmpz_demote((*polyr) + j);
                  len = 0;
                  l = 0;

                  goto big;
               }

               /* abs bound and sum of abs vals of coeffs of quotient chunk */
               _fmpz_vec_sum_max_bits(&b1[i], &maxb1[i], (*polyq)+i1[i], n1[i]);

               /* update length of output quotient and remainder polys */
               len += n1[i];
               l += *lenr;
            } else /* remainder terms only */
            {
               /* realloc remainder poly */
               _fmpz_mpoly_fit_length(polyr, expr, allocr, l + tlen, 1);

               /* for each term in remainder chunk */
               for (j = 0; j < tlen; j++)
               {
                  /* set remainder coeff and exponent */
                  fmpz_set(*polyr + l + j, temp + j);
                  (*expr)[l + j] = (texp[j]) + ((l2 - i - 1) << shift);
               }

               l += tlen;
            }
         } else if (i < l1) /* zero chunk, no quotient or remainder */
         {
            /* set index and length of quotient chunk */
            i1[i] = len;
            n1[i] = 0;
            b1[i] = 0;
            maxb1[i] = 0;

            /* write out maximum exponents in chunk */
            mpoly_max_degrees_tight(max_exp1 + i*num,
	                                   (*expq) + i1[i], n1[i], prods, num);
         }
      }
   }

big:

   /* if not done, use multiprecision coeffs instead */
   if (len == 0 && l == 0)
   {
      fmpz * p2 = (fmpz *) TMP_ALLOC(prod*sizeof(fmpz));

      for (j = 0; j < prod; j++)
            fmpz_init(p2 + j);
      
      /* for each chunk of poly2 */
      for (i = 0; i < l2; i++)
      {
         for (j = 0; j < prod; j++)
            fmpz_zero(p2 + j);

         /* convert relevant coeff/chunk of poly2 to array format */
         _fmpz_mpoly_to_fmpz_array(p2, poly2 + i2[i], e2 + i2[i], n2[i]);

         /* submuls */

         for (j = 0; j < i && j < l1; j++)
         {
            k = i - j;

            if (k < l3 && k >= 0)
            {
	       for (m = 0; m < num; m++)
	       {
	          if (max_exp1[j*num + m] + max_exp3[k*num + m] >= mults[m])
	          {
                     for (j = 0; j < len; j++)
                        _fmpz_demote((*polyq) + j);
                     for (j = 0; j < l; j++)
                        _fmpz_demote((*polyr) + j);
		     len = 0;
	             l = 0;
						
                for (j = 0; j < prod; j++)
                   fmpz_clear(p2 + j);
		     goto cleanup3;
	          }
	       }

               _fmpz_mpoly_submul_array1_fmpz(p2, (*polyq) + i1[j],
                    (*expq) + i1[j], n1[j], poly3 + i3[k], e3 + i3[k], n3[k]);
	    }
         }

         /* convert chunk from array format */
         tlen = _fmpz_mpoly_from_fmpz_array(&temp, &texp, &talloc, 
                                                      p2, mults, num, bits, 0);

         if (tlen != 0) /* nonzero coeff/chunk */
         {
            if (i < l1) /* potentially a quotient with remainder */
            {
               /* tightly pack chunk exponents */
               mpoly_pack_monomials_tight(texp, texp, tlen, mults, num, bits);

               /* set start index of quotient chunk we are about to compute */
               i1[i] = len;
            
               /* compute quotient chunk and set length of quotient chunk */
               n1[i] = _fmpz_mpoly_divrem_array_tight(lenr, polyq,
                      expq, allocq, len, polyr, expr, allocr, l, temp, texp, 
                           tlen, poly3 + i3[0], e3 + i3[0], n3[0], mults, num);

	       /* unpack remainder exponents */
	       mpoly_unpack_monomials_tight(*expr + l,
			                *expr + l, *lenr, mults, num, bits);
														  
	       /* insert main variable */
	       for (j = 0; j < *lenr; j++)
			             (*expr)[l + j] += ((l2 - i - 1) << shift);

	       /* work out maximum exponents for chunk */
               mpoly_max_degrees_tight(max_exp1 + i*num,
	                                   (*expq) + i1[i], n1[i], prods, num);

	       /* check there were no exponent overflows */
	       for (m = 0; m < num; m++)
	       {
	          if (max_exp3[m] + max_exp1[i*num + m] >= mults[m])
                  {
                     for (j = 0; j < len; j++)
                        _fmpz_demote((*polyq) + j);
                     for (j = 0; j < l; j++)
                        _fmpz_demote((*polyr) + j);
                     len = 0;
	             l = 0;

                for (j = 0; j < prod; j++)
                   fmpz_clear(p2 + j);
	             goto cleanup3;
                  }				  
	       }

               /* abs bound and sum of abs vals of coeffs of quotient chunk */
               _fmpz_vec_sum_max_bits(&b1[i], &maxb1[i], (*polyq)+i1[i], n1[i]);

               /* update length of output quotient and remainder polys */
               len += n1[i];
               l += *lenr;
            } else /* remainder terms only */
            {
               /* realloc remainder poly */
               _fmpz_mpoly_fit_length(polyr, expr, allocr, l + tlen, 1);

               /* for each term in chunk */
               for (j = 0; j < tlen; j++)
               {
                  /* set remainder coeff and exponent */
                  fmpz_set(*polyr + l + j, temp + j);
                  (*expr)[l + j] = (texp[j]) + ((l2 - i - 1) << shift);
               }

               /* update length of output remainder poly */
               l += tlen;
            }
         } else if (i < l1) /* zero chunk, no quotient or remainder */
         {
            /* set index and length of quotient chunk */
            i1[i] = len;
            n1[i] = 0;
            b1[i] = 0;
            maxb1[i] = 0;

            /* write out maximum exponents in chunk */
            mpoly_max_degrees_tight(max_exp1 + i*num,
	                                   (*expq) + i1[i], n1[i], prods, num);
         }
      }

      for (j = 0; j < prod; j++)
            fmpz_clear(p2 + j);
   }

   /* if there were quotient terms */
   if (len != 0)
   {
      /* unpack monomials of quotient */
      mpoly_unpack_monomials_tight((*expq), (*expq), len, mults, num, bits);

      /* put main variable back in quotient */
      for (i = 0; i < l1; i++)
      {
         for (j = 0; j < n1[i]; j++)
         {
            (*expq)[i1[i] + j] += ((l1 - i - 1) << shift);
         }
      }
   }

cleanup3:

   for (j = 0; j < talloc; j++)
      fmpz_clear(temp + j);

   flint_free(temp);
   flint_free(texp);

   TMP_END;

   /* set remainder length */
   *lenr = l;

   /* return quotient length */
   return len;
}

/*
   Use dense array division to set polyq, polyr to poly2/poly3 in num variables,
   given a list of multipliers to tightly pack exponents and a number of bits
   for the fields of the exponents of the result, assuming no aliasing. The
   array "mults" is a list of bases to be used in encoding the array indices
   from the exponents. The function reallocates its output and returns the
   length of the quotient. It is assumed that poly2 is not zero. The quotient
   and remainder are written in reverse order.
*/
slong _fmpz_mpoly_divrem_array(slong * lenr,
       fmpz ** polyq, ulong ** expq, slong * allocq,
              fmpz ** polyr, ulong ** expr, slong * allocr, 
                const fmpz * poly2, const ulong * exp2, slong len2, 
        const fmpz * poly3, const ulong * exp3, slong len3, slong * mults, 
                                                         slong num, slong bits)
{
   slong i;
   ulong * e2, * e3;
   slong len, prod;
   slong * prods, * max_exp1, * max_exp3;
   TMP_INIT;

   TMP_START;

   /*
      compute products 1, b0, b0*b1, b0*b1*b2 ...
      from list of bases b0, b1, b2, ...
   */
   prods = (slong *) TMP_ALLOC((num + 1)*sizeof(slong));
  
   prods[0] = 1;
   for (i = 0; i < num; i++)
      prods[i + 1] = prods[i]*mults[i];
   prod = prods[num];

   /* if array size will be too large, chunk the polynomials */
   if (prod > MAX_ARRAY_SIZE)
   {
      TMP_END;
	  
	  return _fmpz_mpoly_divrem_array_chunked(lenr, polyq, expq, allocq,
                                     polyr, expr, allocr, poly2, exp2, len2,
                                      poly3, exp3, len3, mults, num - 1, bits);
   }
									  
   e2 = (ulong *) TMP_ALLOC(len2*sizeof(ulong));
   e3 = (ulong *) TMP_ALLOC(len3*sizeof(ulong));
   max_exp1 = (slong *) TMP_ALLOC(num*sizeof(slong));
   max_exp3 = (slong *) TMP_ALLOC(num*sizeof(slong));

   /* pack input exponents tightly with mixed bases specified by "mults" */

   mpoly_pack_monomials_tight(e2, exp2, len2, mults, num, bits);
   mpoly_pack_monomials_tight(e3, exp3, len3, mults, num, bits);

   /* do divrem on tightly packed polys */
   len = _fmpz_mpoly_divrem_array_tight(lenr, polyq, expq, allocq, 0,
                                  polyr, expr, allocr, 0, poly2, e2, len2,
                                                  poly3, e3, len3, mults, num);

   /* check for overflows */
   mpoly_max_degrees_tight(max_exp3, e3, len3, prods, num);
   mpoly_max_degrees_tight(max_exp1, *expq, len, prods, num);
   
   for (i = 0; i < num; i++)
   {
      if (max_exp3[i] + max_exp1[i] >= mults[i])
      {
	 len = 0;
         *lenr = 0;
		 
         break;
      }
   }
   /* unpack output quotient and remainder exponents */
   mpoly_unpack_monomials_tight((*expq), (*expq), len, mults, num, bits);
   mpoly_unpack_monomials_tight((*expr), (*expr), *lenr, mults, num, bits);

   TMP_END;

   return len;
}

/*
   Return 1 if q, r can be set to quotient and remainder of poly2 by poly3,
   else return 0 if array division not able to be performed.
*/
int fmpz_mpoly_divrem_array(fmpz_mpoly_t q, fmpz_mpoly_t r,
                    const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, 
                                                    const fmpz_mpoly_ctx_t ctx)
{
   slong i, exp_bits, N, lenq = 0, lenr = 0, array_size;
   ulong * max_fields, * max_fields2, * max_fields3;
   ulong * exp2 = poly2->exps, * exp3 = poly3->exps;
   int free2 = 0, free3 = 0;
   fmpz_mpoly_t temp1, temp2;
   fmpz_mpoly_struct * tq, * tr;
   int res = 0;

   TMP_INIT;

   /* check divisor is not zero */
   if (poly3->length == 0)
      flint_throw(FLINT_DIVZERO, "Divide by zero in fmpz_mpoly_divrem_array");

   /* dividend is zero */
   if (poly2->length == 0)
   {
      fmpz_mpoly_zero(q, ctx);
      fmpz_mpoly_zero(r, ctx);

      return 1;
   }

   TMP_START;


   /* compute maximum exponents for each variable */
   max_fields = (ulong *) TMP_ALLOC(ctx->minfo->nfields*sizeof(ulong));
   max_fields2 = (ulong *) TMP_ALLOC(ctx->minfo->nfields*sizeof(ulong));
   max_fields3 = (ulong *) TMP_ALLOC(ctx->minfo->nfields*sizeof(ulong));
   mpoly_max_fields_ui_sp(max_fields2, poly2->exps, poly2->length,
                                                      poly2->bits, ctx->minfo);
   mpoly_max_fields_ui_sp(max_fields3, poly3->exps, poly3->length,
                                                      poly3->bits, ctx->minfo);
   for (i = 0; i < ctx->minfo->nfields; i++)
      max_fields[i] = FLINT_MAX(max_fields2[i], max_fields3[i]);

   /* compute number of bits required for output exponents */
   exp_bits = FLINT_MAX(poly2->bits, poly3->bits);
   N = mpoly_words_per_exp(exp_bits, ctx->minfo);

   /* array division expects each exponent vector in one word */
   /* current code is wrong for reversed orderings */
   if (N != 1 || mpoly_ordering_isrev(ctx->minfo))
      goto cleanup;

   /* compute bounds on output exps, used as mixed bases for packing exps */
   array_size = 1;
   for (i = 0; i < ctx->minfo->nfields - 1; i++)
   {
      max_fields2[i] = max_fields[i] + 1;
      array_size *= max_fields2[i];
   }  
   max_fields2[ctx->minfo->nfields - 1] = max_fields[ctx->minfo->nfields - 1] + 1;
   
   /* if exponents too large for array multiplication, exit silently */
   if (array_size > MAX_ARRAY_SIZE)
      goto cleanup;

   /* expand input exponents to same number of bits as output */
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

   if (exp2[0] < exp3[0])
   {
      fmpz_mpoly_set(r, poly2, ctx);
      fmpz_mpoly_zero(q, ctx);
	  
      res = 1;
	  
      goto cleanup2;
   }

   /* handle aliasing and do array division */

   if (q == poly2 || q == poly3)
   {
      fmpz_mpoly_init2(temp1, poly2->length/poly3->length + 1, ctx);
      fmpz_mpoly_fit_bits(temp1, exp_bits, ctx);
      temp1->bits = exp_bits;

      tq = temp1;
   } else
   {
      fmpz_mpoly_fit_length(q, poly2->length/poly3->length + 1, ctx);
      fmpz_mpoly_fit_bits(q, exp_bits, ctx);
      q->bits = exp_bits;

      tq = q;
   }

   if (r == poly2 || r == poly3)
   {
      fmpz_mpoly_init2(temp2, poly3->length, ctx);
      fmpz_mpoly_fit_bits(temp2, exp_bits, ctx);
      temp2->bits = exp_bits;

      tr = temp2;
   } else
   {
      fmpz_mpoly_fit_length(r, poly3->length, ctx);
      fmpz_mpoly_fit_bits(r, exp_bits, ctx);
      r->bits = exp_bits;

      tr = r;
   }

   lenq = _fmpz_mpoly_divrem_array(&lenr, &tq->coeffs, &tq->exps,
        &tq->alloc, &tr->coeffs, &tr->exps, &tr->alloc, poly2->coeffs,
                  exp2, poly2->length, poly3->coeffs, exp3, poly3->length,
                         (slong *) max_fields2, ctx->minfo->nfields, exp_bits);

    res = (lenq != 0 || lenr != 0);

    if (res)
    {
        if (q == poly2 || q == poly3)
        {
            fmpz_mpoly_swap(temp1, q, ctx);
            fmpz_mpoly_clear(temp1, ctx);
        }

        if (r == poly2 || r == poly3)
        {
            fmpz_mpoly_swap(temp2, r, ctx);
            fmpz_mpoly_clear(temp2, ctx);
        }
    }
    else
    {
        if (q == poly2 || q == poly3)
        {
            fmpz_mpoly_clear(temp1, ctx);
        }

        if (r == poly2 || r == poly3)
        {
            fmpz_mpoly_clear(temp2, ctx);
        }

        for (i = q->length; i < q->alloc; i++)
        {
            _fmpz_demote(q->coeffs + i);
        }
        for (i = r->length; i < r->alloc; i++)
        {
            _fmpz_demote(r->coeffs + i);
        }
    }

    _fmpz_mpoly_set_length(q, lenq, ctx);
    _fmpz_mpoly_set_length(r, lenr, ctx);


cleanup2:

   if (free2)
      flint_free(exp2);

   if (free3)
      flint_free(exp3);

cleanup:

   TMP_END;

   return res;
}

void fmpz_mpoly_divrem_ideal(fmpz_mpoly_struct ** Q,
     fmpz_mpoly_t R, const fmpz_mpoly_t A, fmpz_mpoly_struct * const * B,
                                        slong len, const fmpz_mpoly_ctx_t ctx)
{
    /* TODO !!! */
    fmpz_mpoly_divrem_ideal_monagan_pearce(Q, R, A, B, len, ctx);
}

/* ensure rounding is towards -\infty like fmpz_fdiv_qr */
#define fdiv_qrnnd(qxx, rxx, nhixx, nloxx, dxx)        \
   do {                                                \
      slong __dxx = (dxx);                             \
      sdiv_qrnnd(qxx, rxx, nhixx, nloxx, dxx);         \
      if (((slong) (rxx) < 0 && (slong) (dxx) > 0) ||  \
          ((slong) (rxx) > 0 && (slong) (dxx) < 0))    \
      {                                                \
         (rxx) += __dxx;                               \
         (qxx)--;                                      \
      }                                                \
   } while (0)

/*
   As for divrem_monagan_pearce1 except that an array of divisor polynomials is
   passed and an array of quotient polynomials is returned. These are not in
   low level format.
*/
slong _fmpz_mpoly_divrem_ideal_monagan_pearce1(fmpz_mpoly_struct ** polyq,
  fmpz ** polyr, ulong ** expr, slong * allocr, const fmpz * poly2,
     const ulong * exp2, slong len2, fmpz_mpoly_struct * const * poly3,
       ulong * const * exp3, slong len, slong bits, const fmpz_mpoly_ctx_t ctx,
                                                                  ulong maskhi)
{
    slong i, j, p, l, w;
    slong next_loc;
    slong * store, * store_base;
    slong len3;
    slong heap_len = 2; /* heap zero index unused */
    mpoly_heap1_s * heap;
    mpoly_nheap_t ** chains;
    slong ** hinds;
    mpoly_nheap_t * x;
    fmpz * p2 = *polyr;
    ulong * e2 = *expr;
    ulong exp, texp;
    ulong c[3]; /* for accumulating coefficients */
    ulong mask;
    ulong * ub;
    slong * k, * s;
    fmpz_t qc, q;
    fmpz * mb;
    int small;
    slong bits2, bits3;
    int d1, d2, div_flag;
    TMP_INIT;

    TMP_START;

    fmpz_init(q);
    fmpz_init(qc);

    bits2 = _fmpz_vec_max_bits(poly2, len2);

    chains = (mpoly_nheap_t **) TMP_ALLOC(len*sizeof(mpoly_nheap_t *));
    hinds = (slong **) TMP_ALLOC(len*sizeof(slong *));

    bits3 = 0;
    len3 = 0;
    for (w = 0; w < len; w++)
    {
        chains[w] = (mpoly_nheap_t *) TMP_ALLOC((poly3[w]->length)*sizeof(mpoly_nheap_t));
        hinds[w] = (slong *) TMP_ALLOC((poly3[w]->length)*sizeof(slong));
        bits3 = FLINT_MAX(bits3, FLINT_ABS(fmpz_mpoly_max_bits(poly3[w])));
        len3 += poly3[w]->length;
        for (i = 0; i < poly3[w]->length; i++)
            hinds[w][i] = 1;
    }
      
    /* allow one bit for sign, one bit for subtraction */
    small = FLINT_ABS(bits2) <= (FLINT_ABS(bits3) + FLINT_BIT_COUNT(len3) + SMALL_FMPZ_BITCOUNT_MAX)
          && FLINT_ABS(bits3) <= SMALL_FMPZ_BITCOUNT_MAX;

    next_loc = len3 + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap1_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap1_s));
    store = store_base = (slong *) TMP_ALLOC(3*len3*sizeof(slong));

    k = (slong *) TMP_ALLOC(len*sizeof(slong));
    s = (slong *) TMP_ALLOC(len*sizeof(slong));
    ub = (ulong *) TMP_ALLOC(len*sizeof(ulong));
    mb = (fmpz * ) TMP_ALLOC(len*sizeof(fmpz));

    mask = mpoly_overflow_mask_sp(bits);

    for (w = 0; w < len; w++)
    {
        k[w] = -WORD(1);
        s[w] = poly3[w]->length;
    }
    l = -WORD(1);
   
    x = chains[0] + 0;
    x->i = -WORD(1);
    x->j = 0;
    x->p = -WORD(1);
    x->next = NULL;
    HEAP_ASSIGN(heap[1], exp2[0], x);

    for (i = 0; i < len; i++)
    {
        fmpz_init(mb + i);
        fmpz_neg(mb + i, poly3[i]->coeffs);
        ub[i] = ((ulong) FLINT_ABS(mb[i])) >> 1; /* abs(poly3[0])/2 */
    }

    while (heap_len > 1)
    {
        exp = heap[1].exp;
        if (mpoly_monomial_overflows1(exp, mask))
            goto exp_overflow;

        c[0] = c[1] = c[2] = 0;
        fmpz_zero(qc);
        while (heap_len > 1 && heap[1].exp == exp)
        {
            x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
            do
            {
                *store++ = x->i;
                *store++ = x->j;
                *store++ = x->p;
                if (x->i != -WORD(1))
                    hinds[x->p][x->i] |= WORD(1);

                if (small)
                {
                    if (x->i == -WORD(1))
                      _fmpz_mpoly_sub_uiuiui_fmpz(c, poly2 + x->j);
                    else
                      _fmpz_mpoly_addmul_uiuiui_fmpz(c,
                             poly3[x->p]->coeffs[x->i], polyq[x->p]->coeffs[x->j]);
                } else
                {
                    if (x->i == -WORD(1))
                        fmpz_sub(qc, qc, poly2 + x->j);
                    else
                        fmpz_addmul(qc, poly3[x->p]->coeffs + x->i,
                                                   polyq[x->p]->coeffs + x->j);
                }
            } while ((x = x->next) != NULL);
        }

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
                    _mpoly_heap_insert1(heap, exp2[x->j], x, &next_loc, &heap_len, maskhi);
                }
            } else
            {
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
                    _mpoly_heap_insert1(heap, exp3[x->p][x->i] +
                                polyq[x->p]->exps[x->j], x, &next_loc, &heap_len, maskhi);
                }
                if (j == k[p])
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
                    _mpoly_heap_insert1(heap, exp3[x->p][x->i] +
                                polyq[x->p]->exps[x->j], x, &next_loc, &heap_len, maskhi);
                }
            }
        }

        if ((small && (c[2] != 0 || c[1] != 0 || c[0] != 0)) ||
                                                 (!small && !fmpz_is_zero(qc)))
        {
            div_flag = 0;
            for (w = 0; w < len; w++)
            {
                d1 = mpoly_monomial_divides1(&texp, exp, exp3[w][0], mask);
                if (d1)
                {
                    d2 = 0;
                    if (small)
                    {
                        ulong d[3];
                        if (0 > (slong) c[2])
                            mpn_neg(d, c, 3);
                        else
                            flint_mpn_copyi(d, c, 3);

                        if (d[2] != 0 || ub[w] <= d[1] ||
                          (ub[w] == 0 && 0 > (slong) d[0])) /* quotient not a small */
                        {
                            fmpz_set_signed_uiuiui(qc, c[2], c[1], c[0]);
                            small = 0;
                        } else /* quotient fits a small */
                        {
                            slong r1;
                            slong tq;
                            fdiv_qrnnd(tq, r1, c[1], c[0], mb[w]);
                            if (tq > COEFF_MAX || tq < COEFF_MIN)
                            {
                                small = 0;
                                fmpz_set_signed_uiuiui(qc, c[2], c[1], c[0]);
                            } else
                            {
                                div_flag = (r1 == 0);
                                d2 = tq != 0;
                                if (d2)
                                {
                                    k[w]++;
                                    fmpz_mpoly_fit_length(polyq[w], k[w] + 1, ctx);
                                    fmpz_set_si(polyq[w]->coeffs + k[w], tq);
                                    polyq[w]->exps[k[w]] = texp;
                                }
                                c[0] = r1;
                                c[2] = c[1] = -(slong)(r1 < 0);
                            }
                        }
                    } 
                    /* quotient non-small case */
                    if (!small)
                    {
                        fmpz_fdiv_qr(q, qc, qc, mb + w);
                        div_flag = fmpz_is_zero(qc);
                        d2 = !fmpz_is_zero(q);
                        if (d2)
                        {
                            k[w]++;
                            fmpz_mpoly_fit_length(polyq[w], k[w] + 1, ctx);
                            fmpz_set(polyq[w]->coeffs + k[w], q);                     
                            polyq[w]->exps[k[w]] = texp;
                        }
                    }
                    if (d2)
                    {
                        if (s[w] > 1)
                        {
                            i = 1;
                            x = chains[w] + i;
                            x->i = i;
                            x->j = k[w];
                            x->p = w;
                            x->next = NULL;
                            hinds[w][x->i] = 2*(x->j + 1) + 0;
                            _mpoly_heap_insert1(heap, exp3[w][i] +
                                                  polyq[w]->exps[k[w]], x,
                                                 &next_loc, &heap_len, maskhi);
                        }
                        s[w] = 1;
                    }
                }
            }
            if (!div_flag)
            {
                l++;
                _fmpz_mpoly_fit_length(&p2, &e2, allocr, l + 1, 1);

                if (small)
                {
                    fmpz_set_signed_uiuiui(p2 + l, c[2], c[1], c[0]);
                    fmpz_neg(p2 + l, p2 + l);
                } else
                {
                    fmpz_neg(p2 + l, qc);
                }
                e2[l] = exp;
            }
        }
    }

cleanup:

   for (i = 0; i < len; i++)
      _fmpz_mpoly_set_length(polyq[i], k[i] + 1, ctx); 
   for (i = 0; i < len; i++)
      fmpz_clear(mb + i);
   fmpz_clear(qc);
   fmpz_clear(q);

   (*polyr) = p2;
   (*expr) = e2;
   
   TMP_END;
   return l + 1;

exp_overflow:
    for (i = 0; i < l; i++)
        _fmpz_demote(p2 + i);
    for (w = 0; w < len; w++)
    {
        for (i = 0; i < k[w]; i++)
            _fmpz_demote(polyq[w]->coeffs + i);
        k[w] = -WORD(1); /* we add 1 to this before exit */
    }
    l = -WORD(2); /* we add 1 to this upon return */
    goto cleanup;
}

/*
   As for divrem_monagan_pearce except that an array of divisor polynomials is
   passed and an array of quotient polynomials is returned. These are not in
   low level format.
*/
slong _fmpz_mpoly_divrem_ideal_monagan_pearce(fmpz_mpoly_struct ** polyq,
  fmpz ** polyr, ulong ** expr, slong * allocr, const fmpz * poly2,
     const ulong * exp2, slong len2, fmpz_mpoly_struct * const * poly3,
                     ulong * const * exp3, slong len, slong N, slong bits, 
                        const fmpz_mpoly_ctx_t ctx, const ulong * cmpmask)
{
    slong i, j, p, l, w;
    slong next_loc;
    slong * store, * store_base;
    slong len3;
    slong heap_len = 2; /* heap zero index unused */
    mpoly_heap_s * heap;
    mpoly_nheap_t ** chains;
    slong ** hinds;
    mpoly_nheap_t * x;
    fmpz * p2 = *polyr;
    ulong * e2 = *expr;
    ulong * exp, * exps, * texp;
    ulong ** exp_list;
    ulong c[3]; /* for accumulating coefficients */
    slong exp_next;
    ulong mask = 0;
    ulong * ub;
    slong * k, * s;
    fmpz_t qc, q;
    fmpz * mb;
    int small;
    slong bits2, bits3;
    int d1, d2, div_flag;
    TMP_INIT;

    if (N == 1)
        return _fmpz_mpoly_divrem_ideal_monagan_pearce1(polyq, polyr, expr,
           allocr, poly2, exp2, len2, poly3, exp3, len, bits, ctx, cmpmask[0]);

    TMP_START;

    fmpz_init(q);
    fmpz_init(qc);

    bits2 = _fmpz_vec_max_bits(poly2, len2);
   
    chains = (mpoly_nheap_t **) TMP_ALLOC(len*sizeof(mpoly_nheap_t *));
    hinds = (slong **) TMP_ALLOC(len*sizeof(slong *));
    bits3 = 0;
    len3 = 0;
    for (w = 0; w < len; w++)
    {
        chains[w] = (mpoly_nheap_t *) TMP_ALLOC((poly3[w]->length)*sizeof(mpoly_nheap_t));
        hinds[w] = (slong *) TMP_ALLOC((poly3[w]->length)*sizeof(slong));
        bits3 = FLINT_MAX(bits3, FLINT_ABS(fmpz_mpoly_max_bits(poly3[w])));
        len3 += poly3[w]->length;
        for (i = 0; i < poly3[w]->length; i++)
            hinds[w][i] = 1;
    }
      
    /* allow one bit for sign, one bit for subtraction */
    small = FLINT_ABS(bits2) <= (FLINT_ABS(bits3) +
           FLINT_BIT_COUNT(len3) + SMALL_FMPZ_BITCOUNT_MAX) &&
           FLINT_ABS(bits3) <= SMALL_FMPZ_BITCOUNT_MAX;

    next_loc = len3 + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap_s));
    store = store_base = (slong *) TMP_ALLOC(3*len3*sizeof(slong));

    exps = (ulong *) TMP_ALLOC(len3*N*sizeof(ulong));
    exp_list = (ulong **) TMP_ALLOC(len3*sizeof(ulong *));
    texp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    exp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    k = (slong *) TMP_ALLOC(len*sizeof(slong));
    s = (slong *) TMP_ALLOC(len*sizeof(slong));
    ub = (ulong *) TMP_ALLOC(len*sizeof(ulong));
    mb = (fmpz * ) TMP_ALLOC(len*sizeof(fmpz));

    exp_next = 0;
    for (i = 0; i < len3; i++)
        exp_list[i] = exps + i*N;

    mask = bits <= FLINT_BITS ? mpoly_overflow_mask_sp(bits) : 0;

    for (w = 0; w < len; w++)
    {
        k[w] = -WORD(1);
        s[w] = poly3[w]->length;
    }
    l = -WORD(1);
   
    x = chains[0] + 0;
    x->i = -WORD(1);
    x->j = 0;
    x->p = -WORD(1);
    x->next = NULL;
    heap[1].next = x;
    heap[1].exp = exp_list[exp_next++];
    mpoly_monomial_set(heap[1].exp, exp2, N);

    for (i = 0; i < len; i++)
    {
        fmpz_init(mb + i);
        fmpz_neg(mb + i, poly3[i]->coeffs);
        ub[i] = ((ulong) FLINT_ABS(mb[i])) >> 1; /* abs(poly3[0])/2 */
    }

    while (heap_len > 1)
    {
        mpoly_monomial_set(exp, heap[1].exp, N);

        c[0] = c[1] = c[2] = 0;
        fmpz_zero(qc);

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

        while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N))
        {
            exp_list[--exp_next] = heap[1].exp;
            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
            do
            {
                *store++ = x->i;
                *store++ = x->j;
                *store++ = x->p;
                if (x->i != -WORD(1))
                    hinds[x->p][x->i] |= WORD(1);

                if (small)
                {
                    if (x->i == -WORD(1))
                        _fmpz_mpoly_sub_uiuiui_fmpz(c, poly2 + x->j);
                    else
                        _fmpz_mpoly_addmul_uiuiui_fmpz(c,
                             poly3[x->p]->coeffs[x->i], polyq[x->p]->coeffs[x->j]);
                } else
                {
                    if (x->i == -WORD(1))
                        fmpz_sub(qc, qc, poly2 + x->j);
                    else
                        fmpz_addmul(qc, poly3[x->p]->coeffs + x->i,
                                                       polyq[x->p]->coeffs + x->j);
                }
            } while ((x = x->next) != NULL);
        }

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
                    mpoly_monomial_set(exp_list[exp_next], exp2 + x->j*N, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
            } else
            {
                /* should we go right */
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
                    mpoly_monomial_add_mp(exp_list[exp_next], exp3[x->p] + x->i*N,
                                                   polyq[x->p]->exps + x->j*N, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }

                /* should we go up? */
                if (j == k[p])
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
                    mpoly_monomial_add_mp(exp_list[exp_next], exp3[x->p] + x->i*N,
                                                   polyq[x->p]->exps + x->j*N, N);
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
            }
        }

        if ((small && (c[2] != 0 || c[1] != 0 || c[0] != 0)) ||
                                                 (!small && !fmpz_is_zero(qc)))
        {
            div_flag = 0;
            for (w = 0; w < len; w++)
            {
                if (bits <= FLINT_BITS)
                    d1 = mpoly_monomial_divides(texp, exp, exp3[w], N, mask);
                else
                    d1 = mpoly_monomial_divides_mp(texp, exp, exp3[w], N, bits);

                if (d1)
                {
                    d2 = 0;
                    if (small)
                    {
                        ulong d[3];
                        if (0 > (slong) c[2])
                            mpn_neg(d, c, 3);
                        else
                            flint_mpn_copyi(d, c, 3);

                        if (d[2] != 0 || ub[w] <= d[1] ||
                          (ub[w] == 0 && 0 > (slong) d[0])) /* quotient not a small */
                        {
                            small = 0;
                            fmpz_set_signed_uiuiui(qc, c[2], c[1], c[0]);
                        } else /* quotient fits a small */
                        {
                            slong r1;
                            slong tq;
                            fdiv_qrnnd(tq, r1, c[1], c[0], mb[w]);
                            if (tq > COEFF_MAX || tq < COEFF_MIN)
                            {
                                small = 0;
                                fmpz_set_signed_uiuiui(qc, c[2], c[1], c[0]);
                            } else
                            {
                                d2 = tq != 0;
                                div_flag = r1 == 0;
                                if (d2)
                                {
                                    k[w]++;
                                    fmpz_mpoly_fit_length(polyq[w], k[w] + 1, ctx);
                                    fmpz_set_si(polyq[w]->coeffs + k[w], tq);
                                    mpoly_monomial_set(polyq[w]->exps + k[w]*N, texp, N);
                                }
                                c[0] = r1;
                                c[2] = c[1] = -(slong)(r1 < 0);
                            }
                        }
                    } 
                    /* quotient non-small case */
                    if (!small)
                    {
                        fmpz_fdiv_qr(q, qc, qc, mb + w);
                        d2 = !fmpz_is_zero(q);
                        div_flag = fmpz_is_zero(qc);
                        if (d2)
                        {
                            k[w]++;
                            fmpz_mpoly_fit_length(polyq[w], k[w] + 1, ctx);
                            fmpz_set(polyq[w]->coeffs + k[w], q);
                            mpoly_monomial_set(polyq[w]->exps + k[w]*N, texp, N);
                        }
                    }
                    if (d2)
                    {
                        if (s[w] > 1)
                        {
                            i = 1;
                            x = chains[w] + i;
                            x->i = i;
                            x->j = k[w];
                            x->p = w;
                            x->next = NULL;
                            hinds[w][x->i] = 2*(x->j + 1) + 0;
                            mpoly_monomial_add_mp(exp_list[exp_next], exp3[w] + i*N, 
                                                   polyq[w]->exps + k[w]*N, N);
                            exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                        }
                        s[w] = 1;
                    }
                }
            }
            if (!div_flag)
            {
                l++;
                _fmpz_mpoly_fit_length(&p2, &e2, allocr, l + 1, N);
                if (small)
                {
                    fmpz_set_signed_uiuiui(p2 + l, c[2], c[1], c[0]);
                    fmpz_neg(p2 + l, p2 + l);
                } else
                {
                    fmpz_neg(p2 + l, qc);
                }
                mpoly_monomial_set(e2 + l*N, exp, N);
            }
        } 
    }

cleanup2:

   for (i = 0; i < len; i++)
      _fmpz_mpoly_set_length(polyq[i], k[i] + 1, ctx); 

   for (i = 0; i < len; i++)
      fmpz_clear(mb + i);
   fmpz_clear(qc);
   fmpz_clear(q);

   (*polyr) = p2;
   (*expr) = e2;
   
   TMP_END;

   return l + 1;

exp_overflow:
    for (i = 0; i < l; i++)
       _fmpz_demote(p2 + i);
    for (w = 0; w < len; w++)
    {
        for (i = 0; i < k[w]; i++)
           _fmpz_demote(polyq[w]->coeffs + i);
        k[w] = -WORD(1);
    }
    l = -WORD(2);
    goto cleanup2;

}

/* Assumes divisor polys don't alias any output polys */
void fmpz_mpoly_divrem_ideal_monagan_pearce(fmpz_mpoly_struct ** q, fmpz_mpoly_t r,
    const fmpz_mpoly_t poly2, fmpz_mpoly_struct * const * poly3, slong len,
                                                    const fmpz_mpoly_ctx_t ctx)
{
   slong i, exp_bits, N, lenr = 0;
   slong len3 = 0;
   ulong * cmpmask;
   ulong * exp2;
   ulong ** exp3;
   int free2 = 0;
   int * free3;
   fmpz_mpoly_t temp2;
   fmpz_mpoly_struct * tr;
   TMP_INIT;

   /* check none of the divisor polynomials is zero */
   for (i = 0; i < len; i++)
   {  
      if (poly3[i]->length == 0)
         flint_throw(FLINT_DIVZERO,
                   "Divide by zero in fmpz_mpoly_divrem_ideal_monagan_pearce");

      len3 = FLINT_MAX(len3, poly3[i]->length);
   }

   /* dividend is zero, write out quotients and remainder */
   if (poly2->length == 0)
   {
      for (i = 0; i < len; i++)
      {
         fmpz_mpoly_zero(q[i], ctx);
      }
      
      fmpz_mpoly_zero(r, ctx);

      return;
   }

   TMP_START;

   free3 = (int *) TMP_ALLOC(len*sizeof(int));
   exp3 = (ulong **) TMP_ALLOC(len*sizeof(ulong *));

   /* compute maximum degrees that can occur in any input or output polys */
   exp_bits = poly2->bits;
   for (i = 0; i < len; i++)
      exp_bits = FLINT_MAX(exp_bits, poly3[i]->bits);

    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);

    N = mpoly_words_per_exp(exp_bits, ctx->minfo);
    cmpmask = (ulong *) flint_malloc(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, exp_bits, ctx->minfo);

   /* ensure input exponents packed to same size as output exponents */
   exp2 = poly2->exps;
   free2 = 0;
   if (exp_bits > poly2->bits)
   {
      free2 = 1;
      exp2 = (ulong *) flint_malloc(N*poly2->length*sizeof(ulong));
      mpoly_repack_monomials(exp2, exp_bits, poly2->exps, poly2->bits,
                                                    poly2->length, ctx->minfo);
   }

   for (i = 0; i < len; i++)
   {
      exp3[i] = poly3[i]->exps;
      free3[i] = 0;
      if (exp_bits > poly3[i]->bits)
      {
         free3[i] = 1;
         exp3[i] = (ulong *) flint_malloc(N*poly3[i]->length*sizeof(ulong));
         mpoly_repack_monomials(exp3[i], exp_bits, poly3[i]->exps,
                                 poly3[i]->bits, poly3[i]->length, ctx->minfo);
      }
      fmpz_mpoly_fit_length(q[i], 1, ctx);
      fmpz_mpoly_fit_bits(q[i], exp_bits, ctx);
      q[i]->bits = exp_bits;
   }

   /* check leading mon. of at least one divisor is at most that of dividend */
   for (i = 0; i < len; i++)
   {
      if (!mpoly_monomial_lt(exp2, exp3[i], N, cmpmask))
         break;
   }

   if (i == len)
   {
      fmpz_mpoly_set(r, poly2, ctx);
      for (i = 0; i < len; i++)
         fmpz_mpoly_zero(q[i], ctx);

      goto cleanup3;
   }

   /* take care of aliasing */
   if (r == poly2)
   {
      fmpz_mpoly_init2(temp2, len3, ctx);
      fmpz_mpoly_fit_bits(temp2, exp_bits, ctx);
      temp2->bits = exp_bits;

      tr = temp2;
   } else
   {
      fmpz_mpoly_fit_length(r, len3, ctx);
      fmpz_mpoly_fit_bits(r, exp_bits, ctx);
      r->bits = exp_bits;

      tr = r;
   }

   /* do division with remainder */
   while (1)
   {
      slong old_exp_bits = exp_bits;
      ulong * old_exp2 = exp2, * old_exp3;

      lenr = _fmpz_mpoly_divrem_ideal_monagan_pearce(q, &tr->coeffs, &tr->exps,
                         &tr->alloc, poly2->coeffs, exp2, poly2->length,
                           poly3, exp3, len, N, exp_bits, ctx, cmpmask);

      if (lenr >= 0) /* check if division was successful */
         break;

      exp_bits = mpoly_fix_bits(exp_bits + 1, ctx->minfo);
      N = mpoly_words_per_exp(exp_bits, ctx->minfo);
      cmpmask = (ulong *) flint_realloc(cmpmask, N*sizeof(ulong));
      mpoly_get_cmpmask(cmpmask, N, exp_bits, ctx->minfo);

      exp2 = (ulong *) flint_malloc(N*poly2->length*sizeof(ulong));
      mpoly_repack_monomials(exp2, exp_bits, old_exp2, old_exp_bits,
                                                    poly2->length, ctx->minfo);

      if (free2)
         flint_free(old_exp2);

      free2 = 1;
 
      fmpz_mpoly_fit_bits(tr, exp_bits, ctx);
      tr->bits = exp_bits;

      for (i = 0; i < len; i++)
      {
         old_exp3 = exp3[i];

         exp3[i] = (ulong *) flint_malloc(N*poly3[i]->length*sizeof(ulong));
         mpoly_repack_monomials(exp3[i], exp_bits, old_exp3, old_exp_bits,
                                                 poly3[i]->length, ctx->minfo);
   
         if (free3[i])
            flint_free(old_exp3);

         free3[i] = 1; 

         fmpz_mpoly_fit_bits(q[i], exp_bits, ctx);
         q[i]->bits = exp_bits;
      }

   }

   /* take care of aliasing */
   if (r == poly2)
   {
      fmpz_mpoly_swap(temp2, r, ctx);
      fmpz_mpoly_clear(temp2, ctx);
   } 

   _fmpz_mpoly_set_length(r, lenr, ctx);

cleanup3:

   if (free2)
      flint_free(exp2);

   for (i = 0; i < len; i++)
   {
      if (free3[i])
         flint_free(exp3[i]);
   }

   flint_free(cmpmask);

   TMP_END;
}

/*
   Set polyq, polyr to the quotient and remainder of poly2 by poly3 (with
   remainder coeffs reduced modulo the leading coeff of poly3), and return
   the length of the quotient. This version of the function assumes the
   exponent vectors all fit in a single word. The exponent vectors are
   assumed to have fields with the given number of bits. Assumes input polys
   are nonzero. Implements "Polynomial division using dynamic arrays, heaps
   and packed exponents" by Michael Monagan and Roman Pearce [1], except that
   we use a heap with smallest exponent at head. Note that if a < b then
   (n - b) < (n - b) where n is the maximum value a and b can take. The word
   "maxn" is set to an exponent vector whose fields are all set to such a
   value n. This allows division from left to right with a heap with smallest
   exponent at the head. Quotient and remainder polys are written in reverse
   order.
   [1] http://www.cecm.sfu.ca/~rpearcea/sdmp/sdmp_paper.pdf 
*/

slong _fmpz_mpoly_divrem_monagan_pearce1(slong * lenr,
   fmpz ** polyq, ulong ** expq, slong * allocq, fmpz ** polyr,
  ulong ** expr, slong * allocr, const fmpz * poly2, const ulong * exp2,
            slong len2, const fmpz * poly3, const ulong * exp3, slong len3,
                                                      slong bits, ulong maskhi)
{
    slong i, j, q_len, r_len, s;
    slong next_loc, heap_len = 2;
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    fmpz * q_coeff = *polyq;
    fmpz * r_coeff = *polyr;
    ulong * q_exp = *expq;
    ulong * r_exp = *expr;
    slong * hind;
    ulong mask, exp;
    fmpz_t r, acc_lg;
    ulong acc_sm[3];
    int lt_divides, small;
    slong bits2, bits3;
    ulong lc_norm = 0, lc_abs = 0, lc_sign = 0, lc_n = 0, lc_i = 0;
    TMP_INIT;

    TMP_START;

    fmpz_init(acc_lg);
    fmpz_init(r);

    /* whether intermediate computations q - a*b will fit in three words */
    bits2 = _fmpz_vec_max_bits(poly2, len2);
    bits3 = _fmpz_vec_max_bits(poly3, len3);
    /* allow one bit for sign, one bit for subtraction */
    small = FLINT_ABS(bits2) <= (FLINT_ABS(bits3) + FLINT_BIT_COUNT(len3) + SMALL_FMPZ_BITCOUNT_MAX)
          && FLINT_ABS(bits3) <= SMALL_FMPZ_BITCOUNT_MAX;

    /* alloc array of heap nodes which can be chained together */
    next_loc = len3 + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap1_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap1_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(len3*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*len3*sizeof(slong));

    /* space for flagged heap indices */
    hind = (slong *) TMP_ALLOC(len3*sizeof(slong));
    for (i = 0; i < len3; i++)
        hind[i] = 1;

    mask = mpoly_overflow_mask_sp(bits);

    /* quotient and remainder poly indices start at -1 */
    q_len = WORD(0);
    r_len = WORD(0);

    /* s is the number of terms * (latest quotient) we should put into heap */
    s = len3;

    /* insert (-1, 0, exp2[0]) into heap */
    x = chain + 0;
    x->i = -WORD(1);
    x->j = 0;
    x->next = NULL;
    HEAP_ASSIGN(heap[1], exp2[0], x);

    /* precompute leading coefficient info in "small" case */
    if (small)
    {
        lc_abs = FLINT_ABS(poly3[0]);
        lc_sign = FLINT_SIGN_EXT(poly3[0]);
        count_leading_zeros(lc_norm, lc_abs);
        lc_n = lc_abs << lc_norm;
        invert_limb(lc_i, lc_n);
    }

    while (heap_len > 1)
    {
        exp = heap[1].exp;

        if (mpoly_monomial_overflows1(exp, mask))
            goto exp_overflow;

        _fmpz_mpoly_fit_length(&q_coeff, &q_exp, allocq, q_len + 1, 1);
        lt_divides = mpoly_monomial_divides1(q_exp + q_len, exp, exp3[0], mask);

        /* take nodes from heap with exponent matching exp */
        if (small)
        {
            acc_sm[0] = acc_sm[1] = acc_sm[2] = 0;
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
                        _fmpz_mpoly_add_uiuiui_fmpz(acc_sm, poly2 + x->j);
                    else
                        _fmpz_mpoly_submul_uiuiui_fmpz(acc_sm, poly3[x->i], q_coeff[x->j]);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);
        } else
        {
            fmpz_zero(acc_lg);  
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
                        fmpz_add(acc_lg, acc_lg, poly2 + x->j);
                    else
                        fmpz_submul(acc_lg, poly3 + x->i, q_coeff + x->j);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);
        }

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

        /* try to divide accumulated term by leading term */
        if (small)
        {
            ulong d0, d1, ds = acc_sm[2];

            /* d1:d0 = abs(acc_sm[1:0]) assuming ds is sign extension of acc_sm[1] */
            sub_ddmmss(d1, d0, acc_sm[1]^ds, acc_sm[0]^ds, ds, ds);
            
            if ((acc_sm[0] | acc_sm[1] | acc_sm[2]) == 0)
            {
                continue;
            }
            if (!lt_divides)
            {
                _fmpz_mpoly_fit_length(&r_coeff, &r_exp, allocr, r_len + 1, 1);
                fmpz_set_signed_uiuiui(r_coeff + r_len, acc_sm[2], acc_sm[1], acc_sm[0]);
                r_exp[r_len] = exp;
                r_len++;
                continue;
            }
            if (ds == FLINT_SIGN_EXT(acc_sm[1]) && d1 < lc_abs)
            {
                ulong qq, rr, nhi, nlo;
                FLINT_ASSERT(0 < lc_norm && lc_norm < FLINT_BITS);
                nhi = (d1 << lc_norm) | (d0 >> (FLINT_BITS - lc_norm));
                nlo = d0 << lc_norm;
                udiv_qrnnd_preinv(qq, rr, nhi, nlo, lc_n, lc_i);
                rr = rr >> lc_norm;
                if (rr != 0)
                {
                    _fmpz_mpoly_fit_length(&r_coeff, &r_exp, allocr, r_len + 1, 1);
                    if (ds == 0)
                        fmpz_set_si(r_coeff + r_len, rr);
                    else
                        fmpz_neg_ui(r_coeff + r_len, rr);
                    r_exp[r_len] = exp;
                    r_len++;
                }

                if (qq == 0)
                    continue;

                if (qq <= COEFF_MAX)
                {
                    _fmpz_demote(q_coeff + q_len);
                    q_coeff[q_len] = qq;
                    if (ds != lc_sign)
                        q_coeff[q_len] = -q_coeff[q_len];
                }
                else
                {
                    small = 0;
                    fmpz_set_ui(q_coeff + q_len, qq);
                    if (ds != lc_sign)
                        fmpz_neg(q_coeff + q_len, q_coeff + q_len);
                }
            }
            else
            {
                small = 0;
                fmpz_set_signed_uiuiui(acc_lg, acc_sm[2], acc_sm[1], acc_sm[0]);
                goto large_lt_divides;
            }
        }
        else
        {
            if (fmpz_is_zero(acc_lg))
            {
                continue;
            }
            if (!lt_divides)
            {
                _fmpz_mpoly_fit_length(&r_coeff, &r_exp, allocr, r_len + 1, 1);
                fmpz_set(r_coeff + r_len, acc_lg); 
                r_exp[r_len] = exp;
                r_len++;
                continue;
            }
large_lt_divides:
            fmpz_fdiv_qr(q_coeff + q_len, r, acc_lg, poly3 + 0);
            if (!fmpz_is_zero(r))
            {
                _fmpz_mpoly_fit_length(&r_coeff, &r_exp, allocr, r_len + 1, 1);
                fmpz_set(r_coeff + r_len, r);                     
                r_exp[r_len] = exp;
                r_len++;
            }
            if (fmpz_is_zero(q_coeff + q_len))
            {
                continue;
            }
        }

        /* put newly generated quotient term back into the heap if necessary */
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

cleanup:

    fmpz_clear(acc_lg);
    fmpz_clear(r);

   (*polyq) = q_coeff;
   (*expq) = q_exp;
   (*polyr) = r_coeff;
   (*expr) = r_exp;
   
   /* set remainder poly length */
   (*lenr) = r_len;

    TMP_END;

    return q_len;

exp_overflow:
    for (i = 0; i < q_len; i++)
        _fmpz_demote(q_coeff + i);
    for (i = 0; i < r_len; i++)
        _fmpz_demote(r_coeff + i);
    q_len = 0;
    r_len = 0;
    goto cleanup;
}




slong _fmpz_mpoly_divrem_monagan_pearce(slong * lenr,
  fmpz ** polyq, ulong ** expq, slong * allocq, fmpz ** polyr,
                  ulong ** expr, slong * allocr, const fmpz * poly2,
   const ulong * exp2, slong len2, const fmpz * poly3, const ulong * exp3, 
                   slong len3, slong bits, slong N, const ulong * cmpmask)
{
    slong i, j, q_len, r_len, s;
    slong next_loc;
    slong heap_len = 2; /* heap zero index unused */
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    fmpz * q_coeff = *polyq;
    fmpz * r_coeff = *polyr;
    ulong * q_exp = *expq;
    ulong * r_exp = *expr;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    ulong mask;
    fmpz_t r, acc_lg;
    ulong acc_sm[3];
    slong * hind;
    int lt_divides, small;
    slong bits2, bits3;
    ulong lc_norm = 0, lc_abs = 0, lc_sign = 0, lc_n = 0, lc_i = 0;
    TMP_INIT;

    if (N == 1)
        return _fmpz_mpoly_divrem_monagan_pearce1(lenr, polyq, expq, allocq,
                                     polyr, expr, allocr, poly2, exp2, len2,
                                          poly3, exp3, len3, bits, cmpmask[0]);

    TMP_START;

    fmpz_init(acc_lg);
    fmpz_init(r);

    /* whether intermediate computations q - a*b will fit in three words */
    bits2 = _fmpz_vec_max_bits(poly2, len2);
    bits3 = _fmpz_vec_max_bits(poly3, len3);
    /* allow one bit for sign, one bit for subtraction */
    small = FLINT_ABS(bits2) <= (FLINT_ABS(bits3) + FLINT_BIT_COUNT(len3) + SMALL_FMPZ_BITCOUNT_MAX)
         && FLINT_ABS(bits3) <= SMALL_FMPZ_BITCOUNT_MAX;


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

    /* space for flagged heap indices */
    hind = (slong *) TMP_ALLOC(len3*sizeof(slong));
    for (i = 0; i < len3; i++)
        hind[i] = 1;

    mask = bits <= FLINT_BITS ? mpoly_overflow_mask_sp(bits) : 0;

    /* quotient and remainder poly indices start at -1 */
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

    /* precompute leading coefficient info in "small" case */
    if (small)
    {
        lc_abs = FLINT_ABS(poly3[0]);
        lc_sign = FLINT_SIGN_EXT(poly3[0]);
        count_leading_zeros(lc_norm, lc_abs);
        lc_n = lc_abs << lc_norm;
        invert_limb(lc_i, lc_n);
    }
   
    while (heap_len > 1)
    {
        mpoly_monomial_set(exp, heap[1].exp, N);

        if (bits <= FLINT_BITS)
        {
            if (mpoly_monomial_overflows(exp, N, mask))
                goto exp_overflow2;
        }
        else
        {
            if (mpoly_monomial_overflows_mp(exp, N, bits))
                goto exp_overflow2;
        }
      
        _fmpz_mpoly_fit_length(&q_coeff, &q_exp, allocq, q_len + 1, N);

        if (bits <= FLINT_BITS)
            lt_divides = mpoly_monomial_divides(q_exp + q_len*N, exp, exp3, N, mask);
        else
            lt_divides = mpoly_monomial_divides_mp(q_exp + q_len*N, exp, exp3, N, bits);

        /* take nodes from heap with exponent matching exp */
        if (small) 
        {
            acc_sm[0] = acc_sm[1] = acc_sm[2] = 0;
            do
            {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
                do
                {
                    *store++ = x->i;
                    *store++ = x->j;
                    if (x->i != -WORD(1))
                        hind[x->i] |= WORD(1);

                    if (x->i == -WORD(1))
                        _fmpz_mpoly_add_uiuiui_fmpz(acc_sm, poly2 + x->j);
                    else
                        _fmpz_mpoly_submul_uiuiui_fmpz(acc_sm, poly3[x->i], q_coeff[x->j]);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));
        } else
        {
            fmpz_zero(acc_lg);
            do
            {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
                do
                {
                    *store++ = x->i;
                    *store++ = x->j;
                    if (x->i != -WORD(1))
                        hind[x->i] |= WORD(1);

                    if (x->i == -WORD(1))
                        fmpz_add(acc_lg, acc_lg, poly2 + x->j);
                    else
                        fmpz_submul(acc_lg, poly3 + x->i, q_coeff + x->j);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));
        }

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
                                                           q_exp   + x->j*N, N);
                    else
                        mpoly_monomial_add_mp(exp_list[exp_next], exp3 + x->i*N,
                                                           q_exp   + x->j*N, N);

                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
            }
        }

        /* try to divide accumulated term by leading term */
        if (small)
        {
            ulong d0, d1, ds = acc_sm[2];

            /* d1:d0 = abs(acc_sm[1:0]) assuming ds is sign extension of acc_sm[1] */
            sub_ddmmss(d1, d0, acc_sm[1]^ds, acc_sm[0]^ds, ds, ds);
            
            if ((acc_sm[0] | acc_sm[1] | acc_sm[2]) == 0)
            {
                continue;
            }
            if (!lt_divides)
            {
                _fmpz_mpoly_fit_length(&r_coeff, &r_exp, allocr, r_len + 1, N);
                fmpz_set_signed_uiuiui(r_coeff + r_len, acc_sm[2], acc_sm[1], acc_sm[0]);
                mpoly_monomial_set(r_exp + r_len*N, exp, N);
                r_len++;
                continue;
            }
            if (ds == FLINT_SIGN_EXT(acc_sm[1]) && d1 < lc_abs)
            {
                ulong qq, rr, nhi, nlo;
                FLINT_ASSERT(0 < lc_norm && lc_norm < FLINT_BITS);
                nhi = (d1 << lc_norm) | (d0 >> (FLINT_BITS - lc_norm));
                nlo = d0 << lc_norm;
                udiv_qrnnd_preinv(qq, rr, nhi, nlo, lc_n, lc_i);
                rr = rr >> lc_norm;
                if (rr != 0)
                {
                    _fmpz_mpoly_fit_length(&r_coeff, &r_exp, allocr, r_len + 1, N);
                    if (ds == 0)
                        fmpz_set_ui(r_coeff + r_len, rr);
                    else
                        fmpz_neg_ui(r_coeff + r_len, rr);
                    mpoly_monomial_set(r_exp + r_len*N, exp, N);
                    r_len++;
                }

                if (qq == 0)
                    continue;

                if (qq <= COEFF_MAX)
                {
                    _fmpz_demote(q_coeff + q_len);
                    q_coeff[q_len] = qq;
                    if (ds != lc_sign)
                        q_coeff[q_len] = -q_coeff[q_len];
                }
                else
                {
                    small = 0;
                    fmpz_set_ui(q_coeff + q_len, qq);
                    if (ds != lc_sign)
                        fmpz_neg(q_coeff + q_len, q_coeff + q_len);
                }
            }
            else
            {
                small = 0;
                fmpz_set_signed_uiuiui(acc_lg, acc_sm[2], acc_sm[1], acc_sm[0]);
                goto large_lt_divides;
            }
        }
        else
        {
            if (fmpz_is_zero(acc_lg))
            {
                continue;
            }
            if (!lt_divides)
            {
                _fmpz_mpoly_fit_length(&r_coeff, &r_exp, allocr, r_len + 1, N);
                fmpz_set(r_coeff + r_len, acc_lg); 
                mpoly_monomial_set(r_exp + r_len*N, exp, N);
                r_len++;
                continue;
            }
large_lt_divides:
            fmpz_fdiv_qr(q_coeff + q_len, r, acc_lg, poly3 + 0);
            if (!fmpz_is_zero(r))
            {
                _fmpz_mpoly_fit_length(&r_coeff, &r_exp, allocr, r_len + 1, N);
                fmpz_set(r_coeff + r_len, r);                     
                mpoly_monomial_set(r_exp + r_len*N, exp, N);
                r_len++;
            }
            if (fmpz_is_zero(q_coeff + q_len))
            {
                continue;
            }
        }

        /* put newly generated quotient term back into the heap if necessary */
        if (s > 1)
        {
            i = 1;
            x = chain + i;
            x->i = i;
            x->j = q_len;
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
        s = 1;
        q_len++;
    }

cleanup2:

    fmpz_clear(acc_lg);
    fmpz_clear(r);

    (*polyq) = q_coeff;
    (*expq) = q_exp;
    (*polyr) = r_coeff;
    (*expr) = r_exp;

    (*lenr) = r_len;

    TMP_END;

    return q_len;

exp_overflow2:
    for (i = 0; i < q_len; i++)
        _fmpz_demote(q_coeff + i);
    for (i = 0; i < r_len; i++)
        _fmpz_demote(r_coeff + i);
    q_len = 0;
    r_len = 0;
    goto cleanup2;
}

void fmpz_mpoly_divrem_monagan_pearce(fmpz_mpoly_t q, fmpz_mpoly_t r,
                  const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3,
                                                    const fmpz_mpoly_ctx_t ctx)
{
   slong exp_bits, N, lenq = 0, lenr = 0;
   ulong * exp2 = poly2->exps, * exp3 = poly3->exps;
   ulong * cmpmask;
   int free2 = 0, free3 = 0;
   fmpz_mpoly_t temp1, temp2;
   fmpz_mpoly_struct * tq, * tr;

   /* check divisor is nonzero */
   if (poly3->length == 0)
      flint_throw(FLINT_DIVZERO, "Divide by zero in fmpz_mpoly_divrem_monagan_pearce");

   /* dividend zero, write out quotient and remainder */
   if (poly2->length == 0)
   {
      fmpz_mpoly_zero(q, ctx);
      fmpz_mpoly_zero(r, ctx);

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
      fmpz_mpoly_set(r, poly2, ctx);
      fmpz_mpoly_zero(q, ctx);
      goto cleanup3;
   }

   /* take care of aliasing */
   if (q == poly2 || q == poly3)
   {
      fmpz_mpoly_init2(temp1, poly2->length/poly3->length + 1, ctx);
      fmpz_mpoly_fit_bits(temp1, exp_bits, ctx);
      temp1->bits = exp_bits;
      tq = temp1;
   } else
   {
      fmpz_mpoly_fit_length(q, poly2->length/poly3->length + 1, ctx);
      fmpz_mpoly_fit_bits(q, exp_bits, ctx);
      q->bits = exp_bits;
      tq = q;
   }

   if (r == poly2 || r == poly3)
   {
      fmpz_mpoly_init2(temp2, poly3->length, ctx);
      fmpz_mpoly_fit_bits(temp2, exp_bits, ctx);
      temp2->bits = exp_bits;
      tr = temp2;
   } else
   {
      fmpz_mpoly_fit_length(r, poly3->length, ctx);
      fmpz_mpoly_fit_bits(r, exp_bits, ctx);
      r->bits = exp_bits;
      tr = r;
   }

   /* do division with remainder */
   while ((lenq = _fmpz_mpoly_divrem_monagan_pearce(&lenr, &tq->coeffs, &tq->exps,
         &tq->alloc, &tr->coeffs, &tr->exps, &tr->alloc, poly2->coeffs, exp2, 
         poly2->length, poly3->coeffs, exp3, poly3->length, exp_bits,
                                                       N, cmpmask)) == 0
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

      fmpz_mpoly_fit_bits(tq, exp_bits, ctx);
      tq->bits = exp_bits;

      fmpz_mpoly_fit_bits(tr, exp_bits, ctx);
      tr->bits = exp_bits;
   }

   /* deal with aliasing */
   if (q == poly2 || q == poly3)
   {
      fmpz_mpoly_swap(temp1, q, ctx);
      fmpz_mpoly_clear(temp1, ctx);
   } 

   if (r == poly2 || r == poly3)
   {
      fmpz_mpoly_swap(temp2, r, ctx);
      fmpz_mpoly_clear(temp2, ctx);
   } 

   _fmpz_mpoly_set_length(q, lenq, ctx);
   _fmpz_mpoly_set_length(r, lenr, ctx);

cleanup3:

   if (free2)
      flint_free(exp2);

   if (free3)
      flint_free(exp3);

   flint_free(cmpmask);
}
