/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz_mpoly.h"

slong _fmpz_mpoly_quasidiv_heap1(fmpz_t scale,
                              fmpz ** polyq, ulong ** expq, slong * allocq,
                         const fmpz * poly2, const ulong * exp2, slong len2,
                         const fmpz * poly3, const ulong * exp3, slong len3,
                                                      slong bits, ulong maskhi)
{
    slong i, j, s = len3;
    slong q_len = 0;
    slong next_loc;
    slong heap_len = 2; /* heap zero index unused */
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    slong * hind;
    fmpz * q_coeff = *polyq;
    ulong * q_exp = *expq;
    ulong acc_sm[3]; /* for accumulating coefficients */
    ulong mask, exp;
    ulong lc_norm = 0, lc_abs = 0, lc_sign = 0, lc_n = 0, lc_i = 0;
    fmpz_t lc_abs_lg, ns, gcd, acc_lg, r, tp;
    slong bits2, bits3;
    int lt_divides, scale_is_one, small;
    fmpz * qs;
    slong qs_alloc;
    TMP_INIT;

    TMP_START;

    fmpz_init(lc_abs_lg);
    fmpz_init(acc_lg);
    fmpz_init(r);
    fmpz_init(tp);
    fmpz_init(ns);
    fmpz_init(gcd);
    fmpz_set_ui(scale, 1);
    scale_is_one = 1;

    qs_alloc = 64;
    qs = (fmpz *) flint_calloc(qs_alloc, sizeof(fmpz));

    /* whether intermediate computations q - a*b will fit in three words */
    bits2 = _fmpz_vec_max_bits(poly2, len2);
    bits3 = _fmpz_vec_max_bits(poly3, len3);
    /* allow one bit for sign, one bit for subtraction */
    small = FLINT_ABS(bits2) <= (FLINT_ABS(bits3) + FLINT_BIT_COUNT(len3) +
           SMALL_FMPZ_BITCOUNT_MAX) && FLINT_ABS(bits3) <= SMALL_FMPZ_BITCOUNT_MAX;

    next_loc = len3 + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap1_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap1_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(len3*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*len3*sizeof(slong));

    /* space for flagged heap indices */
    hind = (slong *) TMP_ALLOC(len3*sizeof(slong));
    for (i = 0; i < len3; i++)
        hind[i] = 1;

    mask = mpoly_overflow_mask_sp(bits);

    /* insert (-1, 0, exp2[0]) into heap */
    x = chain + 0;
    x->i = -WORD(1);
    x->j = 0;
    x->next = NULL;
    HEAP_ASSIGN(heap[1], exp2[0], x);

    /* precompute leading coefficient info */
    fmpz_abs(lc_abs_lg, poly3 + 0);
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
        /* make sure quotient array has space for q_len + 1 entries */
        _fmpz_mpoly_fit_length(&q_coeff, &q_exp, allocq, q_len + 1, 1);
        if (q_len + 1 > qs_alloc)
        {
            slong len;
            len = FLINT_MAX(q_len + 1, 2*qs_alloc);
            qs = (fmpz *) flint_realloc(qs, len*sizeof(fmpz));
            if (len > qs_alloc)
                flint_mpn_zero((mp_ptr) (qs + qs_alloc), len - qs_alloc);
            qs_alloc = len;
        }

        exp = heap[1].exp;
        if (mpoly_monomial_overflows1(exp, mask))
            goto exp_overflow;

        lt_divides = mpoly_monomial_divides1(q_exp + q_len, exp, exp3[0], mask);

        /* accumulate terms from highest terms on heap */

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

            FLINT_ASSERT(scale_is_one);
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
        else if (scale_is_one)
        {
            /* general coeff arithmetic with scale 1 */

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
                        fmpz_addmul(acc_lg, scale, poly2 + x->j);
                    else
                    {
                        fmpz_divexact(tp, scale, qs + x->j);
                        fmpz_mul(tp, tp, q_coeff + x->j);
                        fmpz_submul(acc_lg, poly3 + x->i, tp);
                    }
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);
        }

        /* process popped nodes */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;
            if (i == -WORD(1))
            {
                if (j + 1 < len2)
                {
                    x = chain + 0;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    _mpoly_heap_insert1(heap, exp2[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
                }
            }
            else
            {
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

        if (!lt_divides)
            continue;

        /* try to divide accumulated term by leading term of divisor */
        if (small)
        {
            ulong d0, d1, ds = acc_sm[2];

            FLINT_ASSERT(scale_is_one);

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

                if (rr == 0 && (qq & (WORD(3) << (SMALL_FMPZ_BITCOUNT_MAX))) == 0)
                {
                    _fmpz_demote(q_coeff + q_len);
                    q_coeff[q_len] = (qq^ds^lc_sign) - (ds^lc_sign);
                    fmpz_set_ui(qs + q_len, 1);
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

            fmpz_gcd(gcd, acc_lg, poly3 + 0);
            fmpz_divexact(ns, lc_abs_lg, gcd);
            if (!fmpz_is_one(ns))
                scale_is_one = 0;
            fmpz_mul(q_coeff + q_len, ns, acc_lg);
            fmpz_divexact(q_coeff + q_len, q_coeff + q_len, poly3 + 0);
            fmpz_mul(qs + q_len, scale, ns);
            fmpz_set(scale, qs + q_len);
        }

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
        q_len++;
        s = 1;
    }

cleanup:

    (*polyq) = q_coeff;
    (*expq)  = q_exp;

    TMP_END;

    for (i = 0; i < q_len; i++)
    {
        fmpz_divexact(tp, scale, qs + i);
        fmpz_mul(q_coeff + i, q_coeff + i, tp);
    }
    for (i = 0; i < qs_alloc; i++)
        fmpz_clear(qs + i);

    flint_free(qs);

    fmpz_clear(lc_abs_lg);
    fmpz_clear(acc_lg);
    fmpz_clear(r);
    fmpz_clear(tp);
    fmpz_clear(ns);
    fmpz_clear(gcd);

    return q_len;

exp_overflow:
    for (i = 0; i <= q_len; i++)
        _fmpz_demote(q_coeff + i);
    q_len = -WORD(1);
    goto cleanup;
}


slong _fmpz_mpoly_quasidiv_heap(fmpz_t scale,
                               fmpz ** polyq, ulong ** expq, slong * allocq,
                           const fmpz * poly2, const ulong * exp2, slong len2,
                           const fmpz * poly3, const ulong * exp3, slong len3,
                                    slong bits, slong N, const ulong * cmpmask)
{
    slong i, j, s = len3;
    slong q_len = 0;
    slong next_loc;
    slong heap_len = 2; /* heap zero index unused */
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    slong * hind;

    fmpz * q_coeff = *polyq;
    ulong * q_exp = *expq;
    ulong * exp, * exps;
    ulong ** exp_list;
    ulong acc_sm[3]; /* for accumulating coefficients */
    slong exp_next;
    ulong mask;
    ulong lc_norm = 0, lc_abs = 0, lc_sign = 0, lc_n = 0, lc_i = 0;
    fmpz_t lc_abs_lg, ns, gcd, acc_lg, r, tp;
    slong bits2, bits3;
    int lt_divides, scale_is_one, small;
    fmpz * qs;
    slong qs_alloc;

    TMP_INIT;

    if (N == 1)
        return _fmpz_mpoly_quasidiv_heap1(scale, polyq, expq, allocq,
                                                 poly2, exp2, len2,
                                                 poly3, exp3, len3,
                                                             bits, cmpmask[0]);

    TMP_START;

    fmpz_init(lc_abs_lg);
    fmpz_init(acc_lg);
    fmpz_init(r);
    fmpz_init(tp);
    fmpz_init(ns);
    fmpz_init(gcd);
    fmpz_set_ui(scale, 1);
    scale_is_one = 1;

    qs_alloc = 64;
    qs = (fmpz *) flint_calloc(qs_alloc, sizeof(fmpz));

    /* whether intermediate computations q - a*b will fit in three words */
    bits2 = _fmpz_vec_max_bits(poly2, len2);
    bits3 = _fmpz_vec_max_bits(poly3, len3);
    /* allow one bit for sign, one bit for subtraction */
    small = FLINT_ABS(bits2) <= (FLINT_ABS(bits3) + FLINT_BIT_COUNT(len3) +
           SMALL_FMPZ_BITCOUNT_MAX) && FLINT_ABS(bits3) <= SMALL_FMPZ_BITCOUNT_MAX;

    next_loc = len3 + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(len3*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*len3*sizeof(slong));

    /* array of exponents of N words each */
    exps = (ulong *) TMP_ALLOC(len3*N*sizeof(ulong));
    /* array of pointers to unused exponent vectors */
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

    /* insert (-1, 0, exp2[0]) into heap */
    x = chain + 0;
    x->i = -WORD(1);
    x->j = 0;
    x->next = NULL;

    heap[1].next = x;
    heap[1].exp = exp_list[exp_next++];

    mpoly_monomial_set(heap[1].exp, exp2, N);

    /* precompute leading coefficient info */
    fmpz_abs(lc_abs_lg, poly3 + 0);
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
        /* make sure quotient array has space for q_len + 1 entries */
        _fmpz_mpoly_fit_length(&q_coeff, &q_exp, allocq, q_len + 1, N);
        if (q_len + 1 > qs_alloc)
        {
            slong len;
            len = FLINT_MAX(q_len + 1, 2*qs_alloc);
            qs = (fmpz *) flint_realloc(qs, len*sizeof(fmpz));
            if (len > qs_alloc)
                flint_mpn_zero((mp_ptr) (qs + qs_alloc), len - qs_alloc);
            qs_alloc = len;
        }

        mpoly_monomial_set(exp, heap[1].exp, N);

        if (bits <= FLINT_BITS)
        {
            if (mpoly_monomial_overflows(exp, N, mask))
                goto exp_overflow;

            lt_divides = mpoly_monomial_divides(q_exp + q_len*N, exp, exp3, N, mask);

        }
        else
        {
            if (mpoly_monomial_overflows_mp(exp, N, bits))
                goto exp_overflow;

            lt_divides = mpoly_monomial_divides_mp(q_exp + q_len*N, exp, exp3, N, bits);
        }

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

            FLINT_ASSERT(scale_is_one);
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
        else if (scale_is_one)
        {
            /* general coeff arithmetic with scale 1 */

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
                        fmpz_addmul(acc_lg, scale, poly2 + x->j);
                    else
                    {
                        fmpz_divexact(tp, scale, qs + x->j);
                        fmpz_mul(tp, tp, q_coeff + x->j);
                        fmpz_submul(acc_lg, poly3 + x->i, tp);
                    }
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));
        }

        /* process popped nodes */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;
            if (i == -WORD(1))
            {
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
                        mpoly_monomial_add(exp_list[exp_next], exp3  + x->i*N,
                                                           q_exp + x->j*N, N);
                    else
                        mpoly_monomial_add_mp(exp_list[exp_next], exp3  + x->i*N,
                                                           q_exp + x->j*N, N);

                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
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
                        mpoly_monomial_add(exp_list[exp_next], exp3  + x->i*N,
                                                           q_exp + x->j*N, N);
                    else
                        mpoly_monomial_add_mp(exp_list[exp_next], exp3  + x->i*N,
                                                           q_exp + x->j*N, N);

                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
            }
        }

        if (!lt_divides)
            continue;

        /* try to divide accumulated term by leading term of divisor */
        if (small)
        {
            ulong d0, d1, ds = acc_sm[2];

            FLINT_ASSERT(scale_is_one);

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

                if (rr ==0 && (qq & (WORD(3) << (SMALL_FMPZ_BITCOUNT_MAX))) == 0)
                {
                    _fmpz_demote(q_coeff + q_len);
                    q_coeff[q_len] = (qq^ds^lc_sign) - (ds^lc_sign);
                    fmpz_set_ui(qs + q_len, 1);
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

            fmpz_gcd(gcd, acc_lg, poly3 + 0);
            fmpz_divexact(ns, lc_abs_lg, gcd);
            if (!fmpz_is_one(ns))
                scale_is_one = 0;
            fmpz_mul(q_coeff + q_len, ns, acc_lg);
            fmpz_divexact(q_coeff + q_len, q_coeff + q_len, poly3 + 0);
            fmpz_mul(qs + q_len, scale, ns);
            fmpz_set(scale, qs + q_len);
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
                mpoly_monomial_add(exp_list[exp_next], exp3  + x->i*N,
                                                       q_exp + x->j*N, N);
            else
                mpoly_monomial_add_mp(exp_list[exp_next], exp3  + x->i*N,
                                                          q_exp + x->j*N, N);

            exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
        }
        q_len++;
        s = 1;
    }

cleanup:

    (*polyq) = q_coeff;
    (*expq)  = q_exp;

    TMP_END;

    for (i = 0; i < q_len; i++)
    {
        fmpz_divexact(tp, scale, qs + i);
        fmpz_mul(q_coeff + i, q_coeff + i, tp);
    }
    for (i = 0; i < qs_alloc; i++)
        fmpz_clear(qs + i);

    flint_free(qs);

    fmpz_clear(lc_abs_lg);
    fmpz_clear(acc_lg);
    fmpz_clear(r);
    fmpz_clear(tp);
    fmpz_clear(ns);
    fmpz_clear(gcd);

    return q_len;

exp_overflow:
    for (i = 0; i <= q_len; i++)
        _fmpz_demote(q_coeff + i);
    q_len = -WORD(1);
    goto cleanup;
}

void fmpz_mpoly_quasidiv_heap(fmpz_t scale, fmpz_mpoly_t q,
                  const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong exp_bits, N, lenq = 0;
    ulong * exp2 = poly2->exps, * exp3 = poly3->exps;
    ulong * cmpmask;
    int free2 = 0, free3 = 0;
    fmpz_mpoly_t temp1;
    fmpz_mpoly_struct * tq;

    /* check divisor is nonzero */
    if (poly3->length == 0)
        flint_throw(FLINT_DIVZERO, "Divide by zero in fmpz_mpoly_quasidiv_heap");

    fmpz_set_ui(scale, 1);

    /* dividend zero, write out quotient and remainder */
    if (poly2->length == 0)
    {
        fmpz_mpoly_zero(q, ctx);
        return;
    }

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

   /* do division with remainder */
   while ((lenq = _fmpz_mpoly_quasidiv_heap(scale, &tq->coeffs, &tq->exps,
                          &tq->alloc,  poly2->coeffs, exp2, poly2->length,
                                       poly3->coeffs, exp3, poly3->length,
                                            exp_bits, N, cmpmask)) == -WORD(1))
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
