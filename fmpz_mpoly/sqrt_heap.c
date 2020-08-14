/*
    Copyright (C) 2016, 2020 William Hart
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"
#include "longlong.h"

/*
   Set polyq to the square root of poly2 and return the length of the square
   root if it exists or zero otherwise. This version of the function assumes
   the exponent vectors all fit in a single word. The exponent vectors are
   assumed to have fields with the given number of bits. Assumes input poly
   is nonzero. Implements "Heap based multivariate square root" by William
   Hart. The word "maxhi" is set to a mask for the degree field of the
   exponent vector (in the case of a degree ordering) to facilitate the
   reverse ordering in degrevlex. Square root is from left to right with a
   heap with largest exponent at the head. Output poly is written in order.
*/
slong _fmpz_mpoly_sqrt_heap1(fmpz ** polyq, ulong ** expq,
                slong * allocq, const fmpz * poly2, const ulong * exp2,
                               slong len2, slong bits, ulong maskhi, int check)
{
    slong i, j, q_len;
    slong next_loc, heap_len = 1;
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    fmpz * q_coeff = *polyq;
    ulong * q_exp = *expq;
    ulong mask, exp, exp3 = 0;
    fmpz_t r, acc_lg, temp;
    ulong acc_sm[3]; /* three word accumulation for small coefficients */
    int lt_divides, small, process_coeffs;
    slong bits2;
    ulong lc_abs = 0; /* sqrt of lc (positive sign) */
    ulong lc_norm = 0; /* number of bits to shift sqrt(lc) to normalise */
    ulong lc_n = 0; /* sqrt of lc normalised */
    ulong lc_i = 0; /* precomputed inverse of sqrt(lc) */
    TMP_INIT;

    TMP_START;

    fmpz_init(acc_lg);
    fmpz_init(r);
    fmpz_init(temp);

   /* if intermediate computations a - c_i*c_j likely fit in three words */
   bits2 = _fmpz_vec_max_bits(poly2, len2);
   /* add 1 for sign, 1 for subtraction and 1 for multiplication by 2 */
   small = FLINT_ABS(bits2) + FLINT_BIT_COUNT(len2) + 3 <= 2*FLINT_BITS && 
           FLINT_ABS(bits2) <= 2*(FLINT_BITS - 3);

    /* alloc array of heap nodes which can be chained together */
    next_loc = len2 + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap1_s *) TMP_ALLOC((next_loc - 2)*sizeof(mpoly_heap1_s));
    chain = (mpoly_heap_t *) TMP_ALLOC((next_loc - 3)*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*(next_loc - 3)*sizeof(mpoly_heap_t *));

    /* mask with high bit set in each field of exponent vector */
    mask = 0;
    for (i = 0; i < FLINT_BITS/bits; i++)
        mask = (mask << bits) + (UWORD(1) << (bits - 1));

    q_len = 0;
   
    /* insert (-1, 1, exp2[1]) into heap */
    if (len2 > 1)
    {
       x = chain + 0;
       x->i = -WORD(1);
       x->j = 1;
       x->next = NULL;
       HEAP_ASSIGN(heap[1], exp2[1], x);
       heap_len++;
    }

    if (fmpz_sgn(poly2 + 0) < 0)/* not a perfect square coeff */
        goto cleanup; 

    _fmpz_mpoly_fit_length(&q_coeff, &q_exp, allocq, q_len + 1, 1);
   
    fmpz_sqrtrem(q_coeff + 0, r, poly2 + 0);

    q_len++;
       
    if (!fmpz_is_zero(r))
       goto not_sqrt;

    /* multiply by 2, we revert this at the end */
    fmpz_mul_2exp(q_coeff + 0, q_coeff + 0, 1);

    /* precompute leading coefficient info */
    if (small)
    {
        lc_abs = q_coeff[0];

        count_leading_zeros(lc_norm, lc_abs);
        lc_n = lc_abs << lc_norm;
        invert_limb(lc_i, lc_n);

        _fmpz_demote(q_coeff + 0);
        q_coeff[q_len] = lc_abs;            
    } 

    if (!mpoly_monomial_halves1(q_exp + 0, exp2[0], mask))
        goto not_sqrt; /* exponent is not square */

    /* optimisation, compute final exponent */
    if (!check)
    {
        if (!mpoly_monomial_halves1(&exp3, exp2[len2 - 1], mask))
            goto not_sqrt; /* exponent is not square */

        exp3 += q_exp[0]; /* overflow not possible */
    }

    while (heap_len > 1)
    {
        exp = heap[1].exp;

        if (mpoly_monomial_overflows1(exp, mask))
            goto exp_overflow;

        _fmpz_mpoly_fit_length(&q_coeff, &q_exp, allocq, q_len + 1, 1);

        lt_divides = mpoly_monomial_divides1(q_exp + q_len, exp, q_exp[0], mask);

        /* take nodes from heap with exponent matching exp */

        process_coeffs = check || !mpoly_monomial_gt1(exp, exp3, maskhi);

        if (!lt_divides && !check)
        {
            do
            {
                x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
                do
                {
                    *store++ = x->i;
                    *store++ = x->j;
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

                    if (process_coeffs)
                    {
                        if (x->i == -WORD(1))
                            _fmpz_mpoly_add_uiuiui_fmpz(acc_sm, poly2 + x->j);
                        else if (x->i == x->j)
                            _fmpz_mpoly_submul_uiuiui_fmpz(acc_sm, q_coeff[x->i], q_coeff[x->j]);
                        else
                        {
                            /* TODO: write function that performs this operation */
                            _fmpz_mpoly_submul_uiuiui_fmpz(acc_sm, q_coeff[x->i], q_coeff[x->j]);
                            _fmpz_mpoly_submul_uiuiui_fmpz(acc_sm, q_coeff[x->i], q_coeff[x->j]);
                        }
                    }
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

                    if (process_coeffs)
                    {
                        if (x->i == -WORD(1))
                            fmpz_add(acc_lg, acc_lg, poly2 + x->j);
                        else
                        {
                            fmpz_mul(temp, q_coeff + x->i, q_coeff + x->j);
                            if (x->i != x->j)
                                fmpz_mul_2exp(temp, temp, 1);
                            fmpz_sub(acc_lg, acc_lg, temp);
                        }
                    }
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
                /* take next input term */
                if (j + 1 < len2)
                {
                    x = chain + 0;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    if (check || !mpoly_monomial_gt1(exp3, exp2[x->j], maskhi))
                    {
                        _mpoly_heap_insert1(heap, exp2[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
                    }
                }
            } else
            {
                /* should we go right */
                if (j < i)
                {
                    x = chain + i;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;

                    if (check || !mpoly_monomial_gt1(exp3, q_exp[x->i] + q_exp[x->j], maskhi))
                    {
                        _mpoly_heap_insert1(heap, q_exp[x->i] + q_exp[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
                    }
                }
            }
        }

        /* try to divide accumulated term by leading term */
        if (!check && !lt_divides)
            continue;

        if (process_coeffs)
        {
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

                    if (rr != 0)
                        goto not_sqrt;
                
                    if (qq == 0)
                        continue;

                    if (qq <= COEFF_MAX)
                    {
                        _fmpz_demote(q_coeff + q_len);
                        q_coeff[q_len] = qq;
                        if (ds != 0)
                            q_coeff[q_len] = -q_coeff[q_len];
                    }
                    else
                    {
                        small = 0;
                        fmpz_set_ui(q_coeff + q_len, qq);
                        if (ds != 0)
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

                fmpz_fdiv_qr(q_coeff + q_len, r, acc_lg, q_coeff + 0);
            
                if (!fmpz_is_zero(r))
                    goto not_sqrt;
            }
        }

        /* put (q_len, 1) in heap */
        i = q_len;
        x = chain + i;
        x->i = i;
        x->j = 1;
        x->next = NULL;

        _mpoly_heap_insert1(heap, q_exp[i] + q_exp[1], x,
                                                 &next_loc, &heap_len, maskhi);

        q_len++;
    }

    /* divide extra factor of 2 back out of leading coefficient */
    fmpz_fdiv_q_2exp(q_coeff + 0, q_coeff + 0, 1);

cleanup:

    fmpz_clear(acc_lg);
    fmpz_clear(r);
    fmpz_clear(temp);

    (*polyq) = q_coeff;
    (*expq) = q_exp;

    TMP_END;

    /* return sqrt poly length, or zero if not a square root */
    return q_len;

not_sqrt:
    for (i = 0; i < q_len; i++)
        _fmpz_demote(q_coeff + i);
    q_len = 0;
    goto cleanup;

exp_overflow:
    for (i = 0; i < q_len; i++)
        _fmpz_demote(q_coeff + i);
    q_len = -WORD(1);
    goto cleanup;
}


slong _fmpz_mpoly_sqrt_heap(fmpz ** polyq,
           ulong ** expq, slong * allocq, const fmpz * poly2,
   const ulong * exp2, slong len2, 
                   slong bits, slong N, const ulong * cmpmask, int check)
{
    slong i, j, q_len;
    slong next_loc;
    slong heap_len = 1;
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    fmpz * q_coeff = *polyq;
    ulong * q_exp = *expq;
    ulong * exp, * exps, * exp3;
    ulong ** exp_list;
    slong exp_next;
    ulong mask;
    fmpz_t r, acc_lg, temp;
    ulong acc_sm[3];
    int lt_divides, small, halves, process_coeffs;
    slong bits2;
    ulong lc_abs = 0; /* sqrt of lc (positive sign) */
    ulong lc_norm = 0; /* number of bits to shift sqrt(lc) to normalise */
    ulong lc_n = 0; /* sqrt of lc normalised */
    ulong lc_i = 0; /* precomputed inverse of sqrt(lc) */
    TMP_INIT;

    if (N == 1)
        return _fmpz_mpoly_sqrt_heap1(polyq, expq, allocq,
                        poly2, exp2, len2, bits, cmpmask[0], check);

    TMP_START;

    fmpz_init(acc_lg);
    fmpz_init(r);
    fmpz_init(temp);

    /* if intermediate computations a - c_i*c_j likely fit in three words */
    bits2 = _fmpz_vec_max_bits(poly2, len2);
    /* add 1 for sign, 1 for subtraction and 1 for multiplication by 2 */
    small = FLINT_ABS(bits2) + FLINT_BIT_COUNT(len2) + 3 <= 2*FLINT_BITS && 
           FLINT_ABS(bits2) <= 2*(FLINT_BITS - 3);

    /* alloc array of heap nodes which can be chained together */
    next_loc = len2 + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap_s *) TMP_ALLOC((next_loc - 2)*sizeof(mpoly_heap_s));
    chain = (mpoly_heap_t *) TMP_ALLOC((next_loc - 3)*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*(next_loc - 3)*sizeof(mpoly_heap_t *));

    /* array of exponent vectors, each of "N" words */
    exps = (ulong *) TMP_ALLOC((next_loc - 3)*N*sizeof(ulong));
    /* list of pointers to available exponent vectors */
    exp_list = (ulong **) TMP_ALLOC((next_loc - 3)*sizeof(ulong *));
    /* space to save copy of current exponent vector */
    exp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    /* final exponent */
    exp3 = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    /* set up list of available exponent vectors */
    exp_next = 0;
    for (i = 0; i < next_loc - 3; i++)
        exp_list[i] = exps + i*N;

    /* mask with high bit set in each word of each field of exponent vector */
    mask = 0;
    for (i = 0; i < FLINT_BITS/bits; i++)
        mask = (mask << bits) + (UWORD(1) << (bits - 1));

    q_len = 0;
   
    /* insert (-1, 1, exp2[0]) into heap */
    if (len2 > 1)
    {
        x = chain + 0;
        x->i = -WORD(1);
        x->j = 1;
        x->next = NULL;
        heap[1].next = x;
        heap[1].exp = exp_list[exp_next++];
        mpoly_monomial_set(heap[1].exp, exp2 + N, N);
        heap_len++;
    }

    if (fmpz_sgn(poly2 + 0) < 0)/* not a perfect square coeff */
        goto not_sqrt; 

    _fmpz_mpoly_fit_length(&q_coeff, &q_exp, allocq, q_len + 1, 1);

    fmpz_sqrtrem(q_coeff + 0, r, poly2 + 0);

    q_len++;

    if (!fmpz_is_zero(r))
        goto not_sqrt;

    /* multiply by 2, we revert this at the end */
    fmpz_mul_2exp(q_coeff + 0, q_coeff + 0, 1);

    /* precompute leading coefficient info */
    if (small)
    {
        lc_abs = q_coeff[0];

        count_leading_zeros(lc_norm, lc_abs);
        lc_n = lc_abs << lc_norm;
        invert_limb(lc_i, lc_n);

        _fmpz_demote(q_coeff + 0);
        q_coeff[q_len] = lc_abs;            
    }
    
    if (bits <= FLINT_BITS)
        halves = mpoly_monomial_halves(q_exp + 0, exp2 + 0, N, mask);
    else
        halves = mpoly_monomial_halves_mp(q_exp + 0, exp2 + 0, N, bits);

    if (!halves)
        goto not_sqrt; /* exponent is not square */

    /* optimisation, compute final exponent */
    if (!check)
    {
        if (bits <= FLINT_BITS)
            halves = mpoly_monomial_halves(exp3, exp2 + (len2 - 1)*N, N, mask);
        else
            halves = mpoly_monomial_halves_mp(exp3, exp2 + (len2 - 1)*N, N, bits);

        if (!halves)
            goto not_sqrt; /* exponent is not square */

        if (bits <= FLINT_BITS)
            mpoly_monomial_add(exp3, exp3, q_exp + 0, N);
        else
            mpoly_monomial_add_mp(exp3, exp3, q_exp + 0, N);
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
            lt_divides = mpoly_monomial_divides(q_exp + q_len*N, exp, q_exp + 0, N, mask);
        else
            lt_divides = mpoly_monomial_divides_mp(q_exp + q_len*N, exp, q_exp + 0, N, bits);

        /* take nodes from heap with exponent matching exp */

        process_coeffs = check || !mpoly_monomial_gt(exp, exp3 + 0, N, cmpmask);

        if (!lt_divides && !check)
        {
            do
            {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
                do
                {
                    *store++ = x->i;
                    *store++ = x->j;
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

                    if (process_coeffs)
                    {
                        if (x->i == -WORD(1))
                            _fmpz_mpoly_add_uiuiui_fmpz(acc_sm, poly2 + x->j);
                        else if (x->i == x->j)
                            _fmpz_mpoly_submul_uiuiui_fmpz(acc_sm, q_coeff[x->i], q_coeff[x->j]);
                        else
                        {
                            /* TODO: write function that performs this operation */
                            _fmpz_mpoly_submul_uiuiui_fmpz(acc_sm, q_coeff[x->i], q_coeff[x->j]);
                            _fmpz_mpoly_submul_uiuiui_fmpz(acc_sm, q_coeff[x->i], q_coeff[x->j]);
                        }
                    }
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

                    if (process_coeffs)
                    {
                        if (x->i == -WORD(1))
                            fmpz_add(acc_lg, acc_lg, poly2 + x->j);
                        else
                        {
                            fmpz_mul(temp, q_coeff + x->i, q_coeff + x->j);
                            if (x->i != x->j)
                                fmpz_mul_2exp(temp, temp, 1);
                            fmpz_sub(acc_lg, acc_lg, temp);
                        }
                    }
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
                /* take next input term */
                if (j + 1 < len2)
                {
                    x = chain + 0;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    if (check || !mpoly_monomial_gt(exp3 + 0, exp2 + x->j*N, N, cmpmask))
                    {
                        mpoly_monomial_set(exp_list[exp_next], exp2 + x->j*N, N);
                        exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                    }
                }
            } else
            {
                /* should we go right */
                if (j < i)
                {
                    x = chain + i;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;

                    if (bits <= FLINT_BITS)
                        mpoly_monomial_add(exp_list[exp_next], q_exp + x->i*N,
                                                              q_exp + x->j*N, N);
                    else
                        mpoly_monomial_add_mp(exp_list[exp_next], q_exp + x->i*N,
                                                                 q_exp + x->j*N, N);
                    if (check || !mpoly_monomial_gt(exp3 + 0, exp_list[exp_next], N, cmpmask))
                    {
                        exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                    }
                }
            }
        }

        /* try to divide accumulated term by leading term */
        if (!check && !lt_divides)
            continue;

        if (process_coeffs)
        {
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

                    if (rr != 0)
                        goto not_sqrt;

                    if (qq == 0)
                        continue;

                    if (qq <= COEFF_MAX)
                    {
                        _fmpz_demote(q_coeff + q_len);
                        q_coeff[q_len] = qq;
                        if (ds != 0)
                            q_coeff[q_len] = -q_coeff[q_len];
                    }
                    else
                    {
                        small = 0;
                        fmpz_set_ui(q_coeff + q_len, qq);
                        if (ds != 0)
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

                fmpz_fdiv_qr(q_coeff + q_len, r, acc_lg, q_coeff + 0);

                if (!fmpz_is_zero(r))
                    goto not_sqrt;
            }
        }

        /* put (q_len, 1) in heap */
        i = q_len;
        x = chain + i;
        x->i = i;
        x->j = 1;
        x->next = NULL;

        if (bits <= FLINT_BITS)
                mpoly_monomial_add(exp_list[exp_next], q_exp + x->i*N,
                                                      q_exp + x->j*N, N);
        else
                mpoly_monomial_add_mp(exp_list[exp_next], q_exp + x->i*N,
                                                         q_exp + x->j*N, N);
        
        exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);

        q_len++;
    }

    /* divide extra factor of 2 back out of leading coefficient */
    fmpz_fdiv_q_2exp(q_coeff + 0, q_coeff + 0, 1);

cleanup:

    fmpz_clear(acc_lg);
    fmpz_clear(r);
    fmpz_clear(temp);

    (*polyq) = q_coeff;
    (*expq) = q_exp;

    TMP_END;

    /* return sqrt poly length, or zero if not a square root */
    return q_len;

not_sqrt:
    for (i = 0; i < q_len; i++)
        _fmpz_demote(q_coeff + i);
    q_len = 0;
    goto cleanup;

exp_overflow:
    for (i = 0; i < q_len; i++)
        _fmpz_demote(q_coeff + i);
    q_len = -WORD(1);
    goto cleanup;

}

int fmpz_mpoly_sqrt_heap(fmpz_mpoly_t q, const fmpz_mpoly_t poly2, 
                          const fmpz_mpoly_ctx_t ctx, int check)
{
    slong exp_bits, N, lenq = 0;
    ulong * exp2 = poly2->exps;
    ulong * cmpmask;
    int free2 = 0;
    fmpz_mpoly_t temp1;
    fmpz_mpoly_struct * tq;

    /* input zero, write out sqrt */
    if (poly2->length == 0)
    {
        fmpz_mpoly_zero(q, ctx);

        return 1;
    }

    /* maximum bits in sqrt exps and input is max for poly2 */
    exp_bits = poly2->bits;

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

    /* take care of aliasing */
    if (q == poly2)
    {
        fmpz_mpoly_init2(temp1, n_sqrt(poly2->length), ctx);
        fmpz_mpoly_fit_bits(temp1, exp_bits, ctx);
        temp1->bits = exp_bits;
        tq = temp1;
    } else
    {
        fmpz_mpoly_fit_length(q, n_sqrt(poly2->length), ctx);
        fmpz_mpoly_fit_bits(q, exp_bits, ctx);
        q->bits = exp_bits;
        tq = q;
    }

    /* do sqrt, check for overflow */
    while ((lenq = _fmpz_mpoly_sqrt_heap(&tq->coeffs, &tq->exps,
                         &tq->alloc, poly2->coeffs, exp2, poly2->length, 
                                      exp_bits, N, cmpmask, check)) == -WORD(1))
    {
        ulong * old_exp2 = exp2;
        slong old_exp_bits = exp_bits;

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

        fmpz_mpoly_fit_bits(tq, exp_bits, ctx);
        tq->bits = exp_bits;
    }

    /* take care of aliasing */
    if (q == poly2)
    {
        fmpz_mpoly_swap(temp1, q, ctx);
        fmpz_mpoly_clear(temp1, ctx);
    } 

    _fmpz_mpoly_set_length(q, lenq, ctx);

    if (free2)
      flint_free(exp2);

    flint_free(cmpmask);

    return lenq != 0;
}
