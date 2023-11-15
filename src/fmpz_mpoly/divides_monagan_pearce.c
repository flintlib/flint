/*
    Copyright (C) 2017 William Hart
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz_mpoly.h"

/*
   Set poly1 to poly2/poly3 if the division is exact, and return the length
   of the quotient. Otherwise return 0. This version of the function assumes
   the exponent vectors all fit in a single word. The exponent vectors are
   assumed to have fields with the given number of bits. Assumes input polys
   are nonzero. Implements "Polynomial division using dynamic arrays, heaps
   and packed exponents" by Michael Monagan and Roman Pearce [1], except that
   we divide from right to left and use a heap with smallest exponent at head.
   [1] http://www.cecm.sfu.ca/~rpearcea/sdmp/sdmp_paper.pdf
*/
slong _fmpz_mpoly_divides_monagan_pearce1(fmpz ** poly1, ulong ** exp1,
         slong * alloc, const fmpz * poly2, const ulong * exp2, slong len2,
                const fmpz * poly3, const ulong * exp3, slong len3, slong bits,
                                                                  ulong maskhi)
{
    slong i, j, k, s;
    slong next_loc, heap_len = 2;
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    fmpz * p1 = *poly1;
    ulong * e1 = *exp1;
    slong * hind;
    ulong mask, exp, maxexp = exp2[len2 - 1];
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

    /* output poly index starts at -1, will be immediately updated to 0 */
    k = -WORD(1);

    /* s is the number of terms * (latest quotient) we should put into heap */
    s = len3;

    /* insert (-1, 0, exp2[0]) into heap */
    x = chain + 0;
    x->i = -WORD(1);
    x->j = 0;
    x->next = NULL;
    HEAP_ASSIGN(heap[1], exp2[0], x);

    if (small)
    {
        lc_abs = FLINT_ABS(poly3[0]);
        lc_sign = FLINT_SIGN_EXT(poly3[0]);
        lc_norm = flint_clz(lc_abs);
        lc_n = lc_abs << lc_norm;
        lc_i = n_preinvert_limb_prenorm(lc_n);
    }

    while (heap_len > 1)
    {
        exp = heap[1].exp;

        if (mpoly_monomial_overflows1(exp, mask))
            goto not_exact_division;

        k++;
        _fmpz_mpoly_fit_length(&p1, &e1, alloc, k + 1, 1);

        lt_divides = mpoly_monomial_divides1(e1 + k, exp, exp3[0], mask);

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
                        _fmpz_mpoly_submul_uiuiui_fmpz(acc_sm, poly3[x->i], p1[x->j]);
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
                        fmpz_submul(acc_lg, poly3 + x->i, p1 + x->j);
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

        /* try to divide accumulated term by leading term */
        if (small)
        {
            ulong d0, d1, ds = acc_sm[2];

            /* d1:d0 = abs(acc_sm[1:0]) assuming ds is sign extension of acc_sm[1] */
            sub_ddmmss(d1, d0, acc_sm[1]^ds, acc_sm[0]^ds, ds, ds);

            if ((acc_sm[0] | acc_sm[1] | acc_sm[2]) == 0)
            {
                k--;
                continue;
            }

            if (ds == FLINT_SIGN_EXT(acc_sm[1]) && d1 < lc_abs)
            {
                ulong qq, rr, nhi, nlo;
                FLINT_ASSERT(0 < lc_norm && lc_norm < FLINT_BITS);
                nhi = (d1 << lc_norm) | (d0 >> (FLINT_BITS - lc_norm));
                nlo = d0 << lc_norm;
                udiv_qrnnd_preinv(qq, rr, nhi, nlo, lc_n, lc_i);
                if (rr != WORD(0))
                    goto not_exact_division;

                if ((qq & (WORD(3) << (SMALL_FMPZ_BITCOUNT_MAX))) == WORD(0))
                {
                    _fmpz_demote(p1 + k);
                    p1[k] = (qq^ds^lc_sign) - (ds^lc_sign);
                } else
                {
                    small = 0;
                    fmpz_set_ui(p1 + k, qq);
                    if (ds != lc_sign)
                        fmpz_neg(p1 + k, p1 + k);
                }
            } else
            {
                small = 0;
                fmpz_set_signed_uiuiui(acc_lg, acc_sm[2], acc_sm[1], acc_sm[0]);
                fmpz_fdiv_qr(p1 + k, r, acc_lg, poly3 + 0);
                if (!fmpz_is_zero(r))
                    goto not_exact_division;
            }

        } else
        {
            if (fmpz_is_zero(acc_lg))
            {
                k--;
                continue;
            }
            fmpz_fdiv_qr(p1 + k, r, acc_lg, poly3 + 0);
            if (!fmpz_is_zero(r))
                goto not_exact_division;
        }

        if (!lt_divides || (exp^maskhi) < (maxexp^maskhi))
            goto not_exact_division;

        /* put newly generated quotient term back into the heap if necessary */
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

    fmpz_clear(acc_lg);
    fmpz_clear(r);

    (*poly1) = p1;
    (*exp1) = e1;

    TMP_END;

    return k;

not_exact_division:
    for (i = 0; i <= k; i++)
        _fmpz_demote(p1 + i);
    k = 0;
    goto cleanup;
}


slong _fmpz_mpoly_divides_monagan_pearce(fmpz ** poly1, ulong ** exp1,
         slong * alloc, const fmpz * poly2, const ulong * exp2, slong len2,
       const fmpz * poly3, const ulong * exp3, slong len3, flint_bitcnt_t bits, slong N,
                                                         const ulong * cmpmask)
{
    slong i, j, k, s;
    slong next_loc;
    slong heap_len = 2; /* heap zero index unused */
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    fmpz * p1 = *poly1;
    ulong * e1 = *exp1;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    fmpz_t r, acc_lg;
    ulong acc_sm[3];
    ulong mask;
    slong * hind;
    int lt_divides, small;
    slong bits2, bits3;
    ulong lc_norm = 0, lc_abs = 0, lc_sign = 0, lc_n = 0, lc_i = 0;
    TMP_INIT;

    /* if exponent vectors are all one word, call specialised version */
    if (N == 1)
        return _fmpz_mpoly_divides_monagan_pearce1(poly1, exp1, alloc,
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

    /* output poly index starts at -1, will be immediately updated to 0 */
    k = -WORD(1);

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

    if (small)
    {
        lc_abs = FLINT_ABS(poly3[0]);
        lc_sign = FLINT_SIGN_EXT(poly3[0]);
        lc_norm = flint_clz(lc_abs);
        lc_n = lc_abs << lc_norm;
        lc_i = n_preinvert_limb_prenorm(lc_n);
    }

    while (heap_len > 1)
    {
        mpoly_monomial_set(exp, heap[1].exp, N);

        if (bits <= FLINT_BITS)
        {
            if (mpoly_monomial_overflows(exp, N, mask))
                goto not_exact_division;
        } else
        {
            if (mpoly_monomial_overflows_mp(exp, N, bits))
                goto not_exact_division;
        }

        k++;
        _fmpz_mpoly_fit_length(&p1, &e1, alloc, k + 1, N);

        if (bits <= FLINT_BITS)
            lt_divides = mpoly_monomial_divides(e1 + k*N, exp, exp3, N, mask);
        else
            lt_divides = mpoly_monomial_divides_mp(e1 + k*N, exp, exp3, N, bits);

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
                        _fmpz_mpoly_submul_uiuiui_fmpz(acc_sm, poly3[x->i], p1[x->j]);
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
                        fmpz_submul(acc_lg, poly3 + x->i, p1 + x->j);
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
                                                               e1 + x->j*N, N);
                    else
                        mpoly_monomial_add_mp(exp_list[exp_next], exp3 + x->i*N,
                                                               e1 + x->j*N, N);

                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
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

                    if (bits <= FLINT_BITS)
                        mpoly_monomial_add(exp_list[exp_next], exp3 + x->i*N,
                                                               e1 + x->j*N, N);
                    else
                        mpoly_monomial_add_mp(exp_list[exp_next], exp3 + x->i*N,
                                                               e1 + x->j*N, N);

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
                k--;
                continue;
            }

            if (ds == FLINT_SIGN_EXT(acc_sm[1]) && d1 < lc_abs)
            {
                ulong qq, rr, nhi, nlo;
                FLINT_ASSERT(0 < lc_norm && lc_norm < FLINT_BITS);
                nhi = (d1 << lc_norm) | (d0 >> (FLINT_BITS - lc_norm));
                nlo = d0 << lc_norm;
                udiv_qrnnd_preinv(qq, rr, nhi, nlo, lc_n, lc_i);
                if (rr != WORD(0))
                    goto not_exact_division;

                if ((qq & (WORD(3) << (SMALL_FMPZ_BITCOUNT_MAX))) == WORD(0))
                {
                    _fmpz_demote(p1 + k);
                    p1[k] = (qq^(ds^lc_sign)) - (ds^lc_sign);
                }
                else
                {
                    small = 0;
                    fmpz_set_ui(p1 + k, qq);
                    if (ds != lc_sign)
                        fmpz_neg(p1 + k, p1 + k);
                }
            } else
            {
                small = 0;
                fmpz_set_signed_uiuiui(acc_lg, acc_sm[2], acc_sm[1], acc_sm[0]);
                fmpz_fdiv_qr(p1 + k, r, acc_lg, poly3 + 0);
                if (!fmpz_is_zero(r))
                    goto not_exact_division;
            }

        } else
        {
            if (fmpz_is_zero(acc_lg))
            {
                k--;
                continue;
            }
            fmpz_fdiv_qr(p1 + k, r, acc_lg, poly3 + 0);
            if (!fmpz_is_zero(r))
                goto not_exact_division;
        }

        if (!lt_divides ||
                mpoly_monomial_gt(exp2 + (len2 - 1)*N, exp, N, cmpmask))
            goto not_exact_division;

        if (s > 1)
        {
            i = 1;
            x = chain + i;
            x->i = i;
            x->j = k;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;

            if (bits <= FLINT_BITS)
                mpoly_monomial_add(exp_list[exp_next], exp3 + x->i*N,
                                                               e1 + x->j*N, N);
            else
                mpoly_monomial_add_mp(exp_list[exp_next], exp3 + x->i*N,
                                                               e1 + x->j*N, N);

            exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
        }
        s = 1;
    }

    k++;

cleanup:

    fmpz_clear(acc_lg);
    fmpz_clear(r);

    (*poly1) = p1;
    (*exp1) = e1;

    TMP_END;

    return k;

not_exact_division:
    for (i = 0; i <= k; i++)
        _fmpz_demote(p1 + i);
    k = 0;
    goto cleanup;
}

/* return 1 if quotient is exact */
int fmpz_mpoly_divides_monagan_pearce(fmpz_mpoly_t poly1,
                  const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i, N, len = 0;
    flint_bitcnt_t exp_bits;
    fmpz * max_fields2, * max_fields3;
    ulong * cmpmask;
    ulong * exp2 = poly2->exps, * exp3 = poly3->exps, * expq;
    int easy_exit, free2 = 0, free3 = 0;
    ulong mask = 0;
    TMP_INIT;

   /* check divisor is nonzero */
   if (poly3->length == 0)
      flint_throw(FLINT_DIVZERO, "Divide by zero in fmpz_mpoly_divides_monagan_pearce");

   /* dividend zero, write out quotient */
   if (poly2->length == 0)
   {
      fmpz_mpoly_zero(poly1, ctx);

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
      fmpz_mpoly_t temp;

      fmpz_mpoly_init2(temp, poly2->length/poly3->length + 1, ctx);
      fmpz_mpoly_fit_bits(temp, exp_bits, ctx);
      temp->bits = exp_bits;

      len = _fmpz_mpoly_divides_monagan_pearce(&temp->coeffs, &temp->exps,
                            &temp->alloc, poly2->coeffs, exp2, poly2->length,
                              poly3->coeffs, exp3, poly3->length, exp_bits, N,
                                                              cmpmask);

      fmpz_mpoly_swap(temp, poly1, ctx);

      fmpz_mpoly_clear(temp, ctx);
   } else
   {
      fmpz_mpoly_fit_length(poly1, poly2->length/poly3->length + 1, ctx);
      fmpz_mpoly_fit_bits(poly1, exp_bits, ctx);
      poly1->bits = exp_bits;

      len = _fmpz_mpoly_divides_monagan_pearce(&poly1->coeffs, &poly1->exps,
                            &poly1->alloc, poly2->coeffs, exp2, poly2->length,
                              poly3->coeffs, exp3, poly3->length, exp_bits, N,
                                                               cmpmask);
   }

cleanup:

   _fmpz_mpoly_set_length(poly1, len, ctx);

   if (free2)
      flint_free(exp2);

   if (free3)
      flint_free(exp3);

   TMP_END;

   /* division is exact if len is nonzero */
   return (len != 0);
}
