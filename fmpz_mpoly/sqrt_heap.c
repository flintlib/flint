/*
    Copyright (C) 2016, 2020 William Hart
    Copyright (C) 2020 Daniel Schultz

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


/* try to prove that A is not a square */
static int _is_definitely_not_square(
    int count,
    mp_limb_t * p,
    flint_rand_t state,
    const fmpz * Acoeffs,
    const ulong * Aexps,
    slong Alen,
    flint_bitcnt_t Abits,
    const mpoly_ctx_t mctx)
{
    int success = 0;
    slong i, N = mpoly_words_per_exp(Abits, mctx);
    ulong mask = Abits <= FLINT_BITS ? mpoly_overflow_mask_sp(Abits) : 0;
    mp_limb_t eval, * alphas;
    nmod_t mod;
    ulong * t;
    TMP_INIT;

    FLINT_ASSERT(Alen > 0);

    TMP_START;
    t = (ulong *) TMP_ALLOC(FLINT_MAX(Alen, N)*sizeof(ulong));

    if (count == 1)
    {
        /* check for odd degrees & check total degree too in degree orderings */
        mpoly_monomial_set(t, Aexps + N*0, N);
        if (Abits <= FLINT_BITS)
        {
            for (i = 1; i < Alen; i++)
                mpoly_monomial_max(t, t, Aexps + N*i, Abits, N, mask);

            success = !mpoly_monomial_halves(t, t, N, mask);
        }
        else
        {
            for (i = 1; i < Alen; i++)
                mpoly_monomial_max_mp(t, t, Aexps + N*i, Abits, N);

            success = !mpoly_monomial_halves_mp(t, t, N, Abits);
        }

        if (success)
            goto cleanup;
    }

    /* try at most 3*count evaluations */
    count *= 3;

    alphas = (mp_limb_t *) TMP_ALLOC(mctx->nvars*sizeof(mp_limb_t));

next_p:

    if (*p >= UWORD_MAX_PRIME)
        *p = UWORD(1) << (FLINT_BITS - 2);
    *p = n_nextprime(*p, 1);
    nmod_init(&mod, *p);

    for (i = 0; i < mctx->nvars; i++)
        alphas[i] = n_urandint(state, mod.n);

    _fmpz_vec_get_nmod_vec(t, Acoeffs, Alen, mod);
    eval = _nmod_mpoly_eval_all_ui(t, Aexps, Alen, Abits, alphas, mctx, mod);

    success = n_jacobi_unsigned(eval, mod.n) < 0;

    if (!success && --count >= 0)
        goto next_p;

cleanup:

    TMP_END;

    return success;
}


#define DEBUG 0

/*
   Set polyq to the square root of poly2 and return the length of the square
   root if it exists or zero otherwise. This version of the function assumes
   the exponent vectors all fit in a single word. The exponent vectors are
   assumed to have fields with the given number of bits. Assumes input poly
   is nonzero. Implements "Heap based multivariate square root" by William
   Hart. Square root is from left to right with a
   heap with largest exponent at the head. Output poly is written in order.
*/
slong _fmpz_mpoly_sqrt_heap1(
    fmpz ** polyq, ulong ** expq, slong * allocq,
    const fmpz * poly2, const ulong * exp2, slong len2,
    flint_bitcnt_t bits,
    const mpoly_ctx_t mctx,
    int check)
{
    slong i, j, q_len;
    slong next_loc, heap_len = 1, heap_alloc;
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain_nodes[64];
    mpoly_heap_t ** chain;
    slong exp_alloc;
    slong * store, * store_base;
    mpoly_heap_t * x;
    fmpz * q_coeff = *polyq;
    ulong * q_exp = *expq;
    ulong mask, exp, exp3 = 0;
    ulong maskhi;
    fmpz_t r, acc_lg, temp;
    ulong acc_sm[3]; /* three word accumulation for small coefficients */
    int lt_divides, small;
    slong bits2;
    ulong lc_abs = 0; /* sqrt of lc (positive sign) */
    ulong lc_norm = 0; /* number of bits to shift sqrt(lc) to normalise */
    ulong lc_n = 0; /* sqrt of lc normalised */
    ulong lc_i = 0; /* precomputed inverse of sqrt(lc) */
    flint_rand_t heuristic_state;
    mp_limb_t heuristic_p = UWORD(1) << (FLINT_BITS - 2);
    int heuristic_count = 0;

#if DEBUG
    printf("Small case\n");
#endif

    FLINT_ASSERT(mpoly_words_per_exp(bits, mctx) == 1);
    mpoly_get_cmpmask(&maskhi, 1, bits, mctx);

    flint_randinit(heuristic_state);

    fmpz_init(acc_lg);
    fmpz_init(r);
    fmpz_init(temp);

   /* if intermediate computations a - c_i*c_j likely fit in three words */
   bits2 = _fmpz_vec_max_bits(poly2, len2);
   /* add 1 for sign, 1 for subtraction and 1 for multiplication by 2 */
   small = FLINT_ABS(bits2) + FLINT_BIT_COUNT(len2) + 3 <= 2*FLINT_BITS && 
           FLINT_ABS(bits2) <= 2*(FLINT_BITS - 3);

    /* alloc array of heap nodes which can be chained together */
    next_loc = 2*n_sqrt(len2) + 4;   /* something bigger than heap can ever be */
    heap_alloc = next_loc - 3;
    heap = (mpoly_heap1_s *) flint_malloc((heap_alloc + 1)*sizeof(mpoly_heap1_s));
    chain_nodes[0] = (mpoly_heap_t *) flint_malloc(heap_alloc*sizeof(mpoly_heap_t));
    chain = (mpoly_heap_t **) flint_malloc(heap_alloc*sizeof(mpoly_heap_t*));
    store = store_base = (slong *) flint_malloc(2*heap_alloc*sizeof(mpoly_heap_t *));
    
    for (i = 0; i < heap_alloc; i++)
       chain[i] = chain_nodes[0] + i;

    exp_alloc = 1;

    mask = mpoly_overflow_mask_sp(bits);

    q_len = 0;
   
    /* insert (-1, 1, exp2[1]) into heap */
    if (len2 > 1)
    {
#if DEBUG
       printf("insert (-1, 1)\n");
#endif
       x = chain[0];
       x->i = -WORD(1);
       x->j = 1;
       x->next = NULL;
       HEAP_ASSIGN(heap[1], exp2[1], x);
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
    } 

    if (!mpoly_monomial_halves1(q_exp + 0, exp2[0], mask))
        goto not_sqrt; /* exponent is not square */

    /* optimisation, compute final exponent */
    {
        if (!fmpz_is_square(poly2 + len2 - 1))
            goto not_sqrt;

        if (!mpoly_monomial_halves1(&exp3, exp2[len2 - 1], mask))
            goto not_sqrt; /* exponent is not square */

        exp3 += q_exp[0]; /* overflow not possible */
    }

    while (heap_len > 1)
    {
        exp = heap[1].exp;

        /*
            All input fields are < 2^(bits - 1). If poly2 is a square, all
            output fields of q are < 2^(bits - 2). Since the heap only contains
                (1) exp2[j], or
                (2) q_exp[i] + q_exp[j],
            if something in the heap does overflow, a bad q_exp must have been
            generated, possibly from bad input.
        */
        if (mpoly_monomial_overflows1(exp, mask))
            goto not_sqrt;

        _fmpz_mpoly_fit_length(&q_coeff, &q_exp, allocq, q_len + 1, 1);

        lt_divides = mpoly_monomial_divides1(q_exp + q_len, exp, q_exp[0], mask);

        /* take nodes from heap with exponent matching exp */

        if (!lt_divides && !check)
        {
            do
            {
                x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
                do
                {
#if DEBUG
                    flint_printf("pop1 (%wd, %wd)\n", x->i, x->j);
#endif
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
#if DEBUG
                    flint_printf("pop2 (%wd, %wd)\n", x->i, x->j);
#endif
                    *store++ = x->i;
                    *store++ = x->j;

                    if (x->i == -WORD(1))
                        _fmpz_mpoly_add_uiuiui_fmpz(acc_sm, poly2 + x->j);
                    else if (x->i == x->j)
                        _fmpz_mpoly_submul_uiuiui_fmpz(acc_sm, q_coeff[x->i], q_coeff[x->j]);
                    else
                        _fmpz_mpoly_submul2_uiuiui_fmpz(acc_sm, q_coeff[x->i], q_coeff[x->j]);

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
#if DEBUG
                    flint_printf("pop3 (%wd, %wd)\n", x->i, x->j);
#endif
                    *store++ = x->i;
                    *store++ = x->j;

                    if (x->i == -WORD(1))
                        fmpz_add(acc_lg, acc_lg, poly2 + x->j);
                    else if (x->i == x->j)
                        fmpz_submul(acc_lg, q_coeff + x->i, q_coeff + x->i);
                    else
                    {
                        fmpz_add(temp, q_coeff + x->j, q_coeff + x->j);
                        fmpz_submul(acc_lg, q_coeff + x->i, temp);
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
                    x = chain[0];
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    if (check || !mpoly_monomial_gt1(exp3, exp2[x->j], maskhi))
                    {
#if DEBUG
                        flint_printf("insert1 (%wd, %wd)\n", x->i, x->j);
#endif
                        _mpoly_heap_insert1(heap, exp2[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
                    }
                }
            } else
            {
                /* should we go right */
                if (j < i)
                {
                    x = chain[i];
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;

                    if (check || !mpoly_monomial_gt1(exp3, q_exp[x->i] + q_exp[x->j], maskhi))
                    {
#if DEBUG
                        flint_printf("insert2 (%wd, %wd)\n", x->i, x->j);
#endif
                        _mpoly_heap_insert1(heap, q_exp[x->i] + q_exp[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
                    }
                }
            }
        }

        /* try to divide accumulated term by leading term */
        if (!check && !lt_divides)
            continue;

        if (small)
        {
            ulong d0, d1, ds = acc_sm[2];

            /* d1:d0 = abs(acc_sm[1:0]) assuming ds is sign extension of acc_sm[1] */
            sub_ddmmss(d1, d0, acc_sm[1]^ds, acc_sm[0]^ds, ds, ds);
            
            if ((acc_sm[0] | acc_sm[1] | acc_sm[2]) == 0)
                continue;

            if (!lt_divides)
                goto not_sqrt;

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

            if (!lt_divides)
                goto not_sqrt;

large_lt_divides:

            fmpz_fdiv_qr(q_coeff + q_len, r, acc_lg, q_coeff + 0);

            if (!fmpz_is_zero(r))
            {
                q_len++;
                goto not_sqrt;
            }
        }

        if (q_len >= heap_alloc)
        {
            /* run some tests if the square root is getting long */
            if (q_len > len2 && _is_definitely_not_square(
                            ++heuristic_count, &heuristic_p, heuristic_state,
                                                poly2, exp2, len2, bits, mctx))
            {
                q_len++; /* for demotion */
                goto not_sqrt;
            }

            heap_alloc *= 2;
            heap = (mpoly_heap1_s *) flint_realloc(heap, (heap_alloc + 1)*sizeof(mpoly_heap1_s));
            chain_nodes[exp_alloc] = (mpoly_heap_t *) flint_malloc((heap_alloc/2)*sizeof(mpoly_heap_t));
            chain = (mpoly_heap_t **) flint_realloc(chain, heap_alloc*sizeof(mpoly_heap_t*));
            store = store_base = (slong *) flint_realloc(store_base, 2*heap_alloc*sizeof(mpoly_heap_t *));
            for (i = 0; i < heap_alloc/2; i++)
                chain[i + heap_alloc/2] = chain_nodes[exp_alloc] + i;
            exp_alloc++;
        }

        /* put (q_len, 1) in heap */
        i = q_len;
        x = chain[i];
        x->i = i;
        x->j = 1;
        x->next = NULL;

        if (check || !mpoly_monomial_gt1(exp3, q_exp[i] + q_exp[1], maskhi))
        {
#if DEBUG
           flint_printf("insert3 (%wd, %wd)\n", x->i, x->j);
#endif

           _mpoly_heap_insert1(heap, q_exp[i] + q_exp[1], x,
                                                 &next_loc, &heap_len, maskhi);
        }

        q_len++;
    }

    /* divide extra factor of 2 back out of leading coefficient */
    fmpz_fdiv_q_2exp(q_coeff + 0, q_coeff + 0, 1);

cleanup:

    flint_randclear(heuristic_state);

    fmpz_clear(acc_lg);
    fmpz_clear(r);
    fmpz_clear(temp);

    (*polyq) = q_coeff;
    (*expq) = q_exp;

    flint_free(heap);
    flint_free(chain);
    flint_free(store_base);
    for (i = 0; i < exp_alloc; i++)
        flint_free(chain_nodes[i]);

    /* return sqrt poly length, or zero if not a square root */
    return q_len;

not_sqrt:
    for (i = 0; i < q_len; i++)
        _fmpz_demote(q_coeff + i);
    q_len = 0;
    goto cleanup;
}


slong _fmpz_mpoly_sqrt_heap(
    fmpz ** polyq, ulong ** expq, slong * allocq,
    const fmpz * poly2, const ulong * exp2, slong len2,
    flint_bitcnt_t bits,
    const mpoly_ctx_t mctx,
    int check)
{
    slong N = mpoly_words_per_exp(bits, mctx);
    ulong * cmpmask;
    slong i, j, q_len;
    slong next_loc;
    slong heap_len = 1, heap_alloc;
    int exp_alloc;
    mpoly_heap_s * heap;
    mpoly_heap_t * chain_nodes[64];
    mpoly_heap_t ** chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    fmpz * q_coeff = *polyq;
    ulong * q_exp = *expq;
    ulong * exp, * exp3;
    ulong * exps[64];
    ulong ** exp_list;
    slong exp_next;
    ulong mask;
    fmpz_t r, acc_lg, temp;
    ulong acc_sm[3];
    int lt_divides, small, halves;
    slong bits2;
    ulong lc_abs = 0; /* sqrt of lc (positive sign) */
    ulong lc_norm = 0; /* number of bits to shift sqrt(lc) to normalise */
    ulong lc_n = 0; /* sqrt of lc normalised */
    ulong lc_i = 0; /* precomputed inverse of sqrt(lc) */
    flint_rand_t heuristic_state;
    mp_limb_t heuristic_p = UWORD(1) << (FLINT_BITS - 2);
    int heuristic_count = 0;

    TMP_INIT;

    if (N == 1)
        return _fmpz_mpoly_sqrt_heap1(polyq, expq, allocq,
                                        poly2, exp2, len2, bits, mctx, check);

#if DEBUG
    printf("Large case\n");
#endif

    TMP_START;

    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, bits, mctx);

    flint_randinit(heuristic_state);

    fmpz_init(acc_lg);
    fmpz_init(r);
    fmpz_init(temp);

    /* if intermediate computations a - c_i*c_j likely fit in three words */
    bits2 = _fmpz_vec_max_bits(poly2, len2);
    /* add 1 for sign, 1 for subtraction and 1 for multiplication by 2 */
    small = FLINT_ABS(bits2) + FLINT_BIT_COUNT(len2) + 3 <= 2*FLINT_BITS && 
           FLINT_ABS(bits2) <= 2*(FLINT_BITS - 3);

    /* alloc array of heap nodes which can be chained together */
    next_loc = 2*sqrt(len2) + 4;   /* something bigger than heap can ever be */
    heap_alloc = next_loc - 3;
    heap = (mpoly_heap_s *) flint_malloc((heap_alloc + 1)*sizeof(mpoly_heap_s));
    chain_nodes[0] = (mpoly_heap_t *) flint_malloc(heap_alloc*sizeof(mpoly_heap_t));
    chain = (mpoly_heap_t **) flint_malloc(heap_alloc*sizeof(mpoly_heap_t*));
    store = store_base = (slong *) flint_malloc(2*heap_alloc*sizeof(mpoly_heap_t *));

    for (i = 0; i < heap_alloc; i++)
       chain[i] = chain_nodes[0] + i;

    /* array of exponent vectors, each of "N" words */
    exps[0] = (ulong *) flint_malloc(heap_alloc*N*sizeof(ulong));
    exp_alloc = 1;
    /* list of pointers to available exponent vectors */
    exp_list = (ulong **) flint_malloc(heap_alloc*sizeof(ulong *));
    /* space to save copy of current exponent vector */
    exp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    /* final exponent */
    exp3 = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    /* set up list of available exponent vectors */
    exp_next = 0;
    for (i = 0; i < heap_alloc; i++)
        exp_list[i] = exps[0] + i*N;

    mask = (bits <= FLINT_BITS) ? mpoly_overflow_mask_sp(bits) : 0;

    q_len = 0;
   
    /* insert (-1, 1, exp2[0]) into heap */
    if (len2 > 1)
    {
#if DEBUG
       printf("insert (-1, 1)\n");
#endif
        x = chain[0];
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
    }
    
    if (bits <= FLINT_BITS)
        halves = mpoly_monomial_halves(q_exp + 0, exp2 + 0, N, mask);
    else
        halves = mpoly_monomial_halves_mp(q_exp + 0, exp2 + 0, N, bits);

    if (!halves)
        goto not_sqrt; /* exponent is not square */

    /* optimisation, compute final term */
    {
        if (!fmpz_is_square(poly2 + len2 - 1))
            goto not_sqrt;

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
        _fmpz_mpoly_fit_length(&q_coeff, &q_exp, allocq, q_len + 1, N);

        mpoly_monomial_set(exp, heap[1].exp, N);

        /*
            All input fields are < 2^(bits - 1). If poly2 is a square, all
            output fields of q are < 2^(bits - 2). Since the heap only contains
                (1) exp2[j], or
                (2) q_exp[i] + q_exp[j],
            if something in the heap does overflow, a bad q_exp must have been
            generated, possibly from bad input.
        */
        if (bits <= FLINT_BITS)
        {
            if (mpoly_monomial_overflows(exp, N, mask))
                goto not_sqrt;

            lt_divides = mpoly_monomial_divides(q_exp + q_len*N,
                                                      exp, q_exp + 0, N, mask);
        }
        else
        {
            if (mpoly_monomial_overflows_mp(exp, N, bits))
                goto not_sqrt;

            lt_divides = mpoly_monomial_divides_mp(q_exp + q_len*N,
                                                      exp, q_exp + 0, N, bits);
        }

        /* take nodes from heap with exponent matching exp */

        if (!lt_divides && !check)
        {
            do
            {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
                do
                {
#if DEBUG
                    flint_printf("pop1 (%wd, %wd)\n", x->i, x->j);
#endif
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
#if DEBUG
                    flint_printf("pop2 (%wd, %wd)\n", x->i, x->j);
#endif
                    *store++ = x->i;
                    *store++ = x->j;

                    if (x->i == -WORD(1))
                        _fmpz_mpoly_add_uiuiui_fmpz(acc_sm, poly2 + x->j);
                    else if (x->i == x->j)
                        _fmpz_mpoly_submul_uiuiui_fmpz(acc_sm, q_coeff[x->i], q_coeff[x->j]);
                    else
                        _fmpz_mpoly_submul2_uiuiui_fmpz(acc_sm, q_coeff[x->i], q_coeff[x->j]);

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
#if DEBUG
                    flint_printf("pop2 (%wd, %wd)\n", x->i, x->j);
#endif
                    *store++ = x->i;
                    *store++ = x->j;

                    if (x->i == -WORD(1))
                        fmpz_add(acc_lg, acc_lg, poly2 + x->j);
                    else if (x->i == x->j)
                        fmpz_submul(acc_lg, q_coeff + x->i, q_coeff + x->i);
                    else
                    {
                        fmpz_add(temp, q_coeff + x->j, q_coeff + x->j);
                        fmpz_submul(acc_lg, q_coeff + x->i, temp);
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
                    x = chain[0];
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    if (check || !mpoly_monomial_gt(exp3 + 0, exp2 + x->j*N, N, cmpmask))
                    {
#if DEBUG
                        flint_printf("insert1 (%wd, %wd)\n", x->i, x->j);
#endif
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
                    x = chain[i];
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
#if DEBUG
                        flint_printf("insert2 (%wd, %wd)\n", x->i, x->j);
#endif
                        exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                    }
                }
            }
        }

        /* try to divide accumulated term by leading term */
        if (!check && !lt_divides)
            continue;

        if (small)
        {
            ulong d0, d1, ds = acc_sm[2];

            /* d1:d0 = abs(acc_sm[1:0]) assuming ds is sign extension of acc_sm[1] */
            sub_ddmmss(d1, d0, acc_sm[1]^ds, acc_sm[0]^ds, ds, ds);
            
            if ((acc_sm[0] | acc_sm[1] | acc_sm[2]) == 0)
                continue;

            if (!lt_divides)
                goto not_sqrt;

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

            if (!lt_divides)
                goto not_sqrt;

large_lt_divides:

            fmpz_fdiv_qr(q_coeff + q_len, r, acc_lg, q_coeff + 0);

            if (!fmpz_is_zero(r))
            {
                q_len++;
                goto not_sqrt;
            }
        }

        if (q_len >= heap_alloc)
        {
            /* run some tests if the square root is getting long */
            if (q_len > len2 && _is_definitely_not_square(
                            ++heuristic_count, &heuristic_p, heuristic_state,
                                                poly2, exp2, len2, bits, mctx))
            {
                q_len++; /* for demotion */
                goto not_sqrt;
            }

            heap_alloc *= 2;
            heap = (mpoly_heap_s *) flint_realloc(heap, (heap_alloc + 1)*sizeof(mpoly_heap_s));
            chain_nodes[exp_alloc] = (mpoly_heap_t *) flint_malloc((heap_alloc/2)*sizeof(mpoly_heap_t));
            chain = (mpoly_heap_t **) flint_realloc(chain, heap_alloc*sizeof(mpoly_heap_t*));
            store = store_base = (slong *) flint_realloc(store_base, 2*heap_alloc*sizeof(mpoly_heap_t *));
            exps[exp_alloc] = (ulong *) flint_malloc((heap_alloc/2)*N*sizeof(ulong));
            exp_list = (ulong **) flint_realloc(exp_list, heap_alloc*sizeof(ulong *));
            for (i = 0; i < heap_alloc/2; i++)
            {
               chain[i + heap_alloc/2] = chain_nodes[exp_alloc] + i;
               exp_list[i + heap_alloc/2] = exps[exp_alloc] + i*N;
            }
            exp_alloc++;
        }

        /* put (q_len, 1) in heap */
        i = q_len;
        x = chain[i];
        x->i = i;
        x->j = 1;
        x->next = NULL;

        if (bits <= FLINT_BITS)
            mpoly_monomial_add(exp_list[exp_next], q_exp + x->i*N,
                                                      q_exp + x->j*N, N);
        else
            mpoly_monomial_add_mp(exp_list[exp_next], q_exp + x->i*N,
                                                         q_exp + x->j*N, N);
        
        if (check || !mpoly_monomial_gt(exp3 + 0, exp_list[exp_next], N, cmpmask))
        {
#if DEBUG
            flint_printf("insert3 (%wd, %wd)\n", x->i, x->j);
#endif           
            exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
        }

        q_len++;
    }

    /* divide extra factor of 2 back out of leading coefficient */
    fmpz_fdiv_q_2exp(q_coeff + 0, q_coeff + 0, 1);

cleanup:

    flint_randclear(heuristic_state);

    fmpz_clear(acc_lg);
    fmpz_clear(r);
    fmpz_clear(temp);

    (*polyq) = q_coeff;
    (*expq) = q_exp;

    flint_free(heap);
    flint_free(chain);
    flint_free(store_base);
    flint_free(exp_list);
    for (i = 0; i < exp_alloc; i++)
    {
        flint_free(exps[i]);
        flint_free(chain_nodes[i]);
    }

    TMP_END;

    /* return sqrt poly length, or zero if not a square root */
    return q_len;

not_sqrt:
    for (i = 0; i < q_len; i++)
        _fmpz_demote(q_coeff + i);
    q_len = 0;
    goto cleanup;
}

int fmpz_mpoly_sqrt_heap(fmpz_mpoly_t q, const fmpz_mpoly_t poly2, 
                          const fmpz_mpoly_ctx_t ctx, int check)
{
    slong lenq, lenq_est;
    flint_bitcnt_t exp_bits;
    ulong * exp2 = poly2->exps;
    fmpz_mpoly_t temp1;
    fmpz_mpoly_struct * tq;

    if (poly2->length == 0)
    {
        fmpz_mpoly_zero(q, ctx);
        return 1;
    }

    /* maximum bits in sqrt exps and input is max for poly2 */
    exp_bits = poly2->bits;
    lenq_est = n_sqrt(poly2->length);

    /* take care of aliasing */
    if (q == poly2)
    {
        fmpz_mpoly_init3(temp1, lenq_est, exp_bits, ctx);
        tq = temp1;
    }
    else
    {
        fmpz_mpoly_fit_length_reset_bits(q, lenq_est, exp_bits, ctx);
        tq = q;
    }

    lenq = _fmpz_mpoly_sqrt_heap(&tq->coeffs, &tq->exps, &tq->alloc,
                                 poly2->coeffs, exp2, poly2->length,
                                      exp_bits, ctx->minfo, check);
    /* take care of aliasing */
    if (q == poly2)
    {
        fmpz_mpoly_swap(temp1, q, ctx);
        fmpz_mpoly_clear(temp1, ctx);
    } 

    _fmpz_mpoly_set_length(q, lenq, ctx);

    return lenq != 0;
}

#undef DEBUG
