/*
    Copyright (C) 2016, 2020 William Hart
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef __GNUC__
# define sqrt __builtin_sqrt
#else
# include <math.h>
#endif

#include "gmpcompat.h"
#include "mpn_extras.h"
#include "nmod_mpoly.h"
#include "fmpz_mpoly.h"

/*
    if r is the returned mpz, then x = r + sm, where sm is a signed 3 limb integer
    it must be safe to add +-COEFF_MAX^2*2^FLINT_BITS to the returned sm
    t is temp space since x is not to be modified
*/
static mpz_srcptr _fmpz_mpoly_get_mpz_signed_uiuiui(ulong * sm, fmpz x, mpz_ptr t)
{
    mpz_ptr p;
    slong i, abs_size;
    ulong s;

    if (!COEFF_IS_MPZ(x))
    {
        sm[0] = x;
        sm[1] = FLINT_SIGN_EXT(x);
        sm[2] = FLINT_SIGN_EXT(x);
    }
    else
    {
        p = COEFF_TO_PTR(x);

        sm[0] = 0;
        sm[1] = 0;
        sm[2] = 0;

        s = FLINT_SIGN_EXT(p->_mp_size);
        abs_size = FLINT_ABS(p->_mp_size);

        if (abs_size > 3 || (abs_size == 3 && p->_mp_d[2] >= COEFF_MAX))
            return p;

        for (i = 0; i < abs_size; i++)
            sm[i] = p->_mp_d[i];

        sub_dddmmmsss(sm[2], sm[1], sm[0], s^sm[2], s^sm[1], s^sm[0], s, s, s);
    }

    mpz_set_ui(t, 0);
    return t;
}

/* try to prove that A is not a square */
static int _is_proved_not_square(
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
    mp_limb_t eval, * alphas;
    nmod_t mod;
    ulong * t;
    TMP_INIT;

    FLINT_ASSERT(Alen > 0);

    TMP_START;
    t = (ulong *) TMP_ALLOC(FLINT_MAX(Alen, N)*sizeof(ulong));

    if (count == 1)
    {
        success = mpoly_is_proved_not_square(Aexps, Alen, Abits, N, t);
        if (success)
            goto cleanup;
    }

    /* try at most 3*count evaluations */
    count *= 3;

    alphas = (mp_limb_t *) TMP_ALLOC(mctx->nvars*sizeof(mp_limb_t));

next_p:

    if (*p >= UWORD_MAX_PRIME)
        *p = UWORD(1) << (SMALL_FMPZ_BITCOUNT_MAX);
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


/*
   Set polyq to the square root of A and return the length of the square
   root if it exists or zero otherwise. This version of the function assumes
   the exponent vectors all fit in a single word. The exponent vectors are
   assumed to have fields with the given number of bits. Assumes input poly
   is nonzero. Implements "Heap based multivariate square root" by William
   Hart. Square root is from left to right with a
   heap with largest exponent at the head. Output poly is written in order.
   A never explicitly goes into the heap and is only scanned once.
   TODO: copy this strategy for the "small" case (i.e. no fmpz arithmetic)
         to the other fmpz mpoly mul/div functions.
*/
slong _fmpz_mpoly_sqrt_heap1(
    fmpz ** polyq, ulong ** expq, slong * allocq,
    const fmpz * Acoeffs, const ulong * Aexps, slong Alen,
    flint_bitcnt_t bits,
    const mpoly_ctx_t mctx,
    int check)
{
    slong i, j, Qlen, Ai;
    slong next_loc, heap_len = 1, heap_alloc;
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain_nodes[64];
    mpoly_heap_t ** chain;
    slong exp_alloc;
    slong * store, * store_base;
    mpoly_heap_t * x;
    fmpz * Qcoeffs = *polyq;
    ulong * Qexps = *expq;
    ulong mask, exp, exp3 = 0;
    ulong maskhi;
    mpz_t r, acc, acc2;
    mpz_srcptr acc_lg;
    mpz_ptr t;
    ulong acc_sm[3], acc_sm2[3], pp[3];
    int lt_divides, q_rest_small;
    flint_rand_t heuristic_state;
    mp_limb_t heuristic_p = UWORD(1) << (SMALL_FMPZ_BITCOUNT_MAX);
    int heuristic_count = 0;
    ulong lc_abs = 0; /* 2*sqrt(lc) if it fits in ulong, otherwise 0 */
    ulong lc_norm = 0;
    ulong lc_n = 0;
    ulong lc_i = 0;
    mpz_ptr lc_lg = NULL; /* 2*sqrt(lc) if it is large */

    FLINT_ASSERT(mpoly_words_per_exp(bits, mctx) == 1);
    mpoly_get_cmpmask(&maskhi, 1, bits, mctx);

    flint_randinit(heuristic_state);

    mpz_init(r);
    mpz_init(acc);
    mpz_init(acc2);

    /* alloc array of heap nodes which can be chained together */
    next_loc = 2*n_sqrt(Alen) + 4;   /* something bigger than heap can ever be */
    heap_alloc = next_loc - 3;
    heap = (mpoly_heap1_s *) flint_malloc((heap_alloc + 1)*sizeof(mpoly_heap1_s));
    chain_nodes[0] = (mpoly_heap_t *) flint_malloc(heap_alloc*sizeof(mpoly_heap_t));
    chain = (mpoly_heap_t **) flint_malloc(heap_alloc*sizeof(mpoly_heap_t*));
    store = store_base = (slong *) flint_malloc(2*heap_alloc*sizeof(mpoly_heap_t *));

    for (i = 0; i < heap_alloc; i++)
       chain[i] = chain_nodes[0] + i;

    exp_alloc = 1;

    mask = mpoly_overflow_mask_sp(bits);

    Qlen = 0;

    /* "insert" (-1, 1, Aexps[1]) into "heap" */
    Ai = 1;

    /* compute first term */
    if (!fmpz_is_square(Acoeffs + 0))
        goto not_sqrt;

    _fmpz_mpoly_fit_length(&Qcoeffs, &Qexps, allocq, Qlen + 1, 1);

    fmpz_sqrt(Qcoeffs + 0, Acoeffs + 0);

    Qlen++;

    /* multiply by 2, we revert this at the end */
    fmpz_mul_2exp(Qcoeffs + 0, Qcoeffs + 0, 1);

    /* q_rest_small means Qcoeffs[1] ... Qcoeffs[Qlen-1] are small */
    q_rest_small = 1;

    if (fmpz_abs_fits_ui(Qcoeffs + 0))
    {
        lc_abs = fmpz_get_ui(Qcoeffs + 0);
        lc_norm = flint_clz(lc_abs);
        lc_n = lc_abs << lc_norm;
        lc_i = n_preinvert_limb_prenorm(lc_n);
    }
    else
    {
        lc_lg = COEFF_TO_PTR(Qcoeffs[0]);
    }

    if (!mpoly_monomial_halves1(Qexps + 0, Aexps[0], mask))
        goto not_sqrt; /* exponent is not square */

    /* optimisation, compute final exponent */
    {
        if (!fmpz_is_square(Acoeffs + Alen - 1))
            goto not_sqrt;

        if (!mpoly_monomial_halves1(&exp3, Aexps[Alen - 1], mask))
            goto not_sqrt; /* exponent is not square */

        exp3 += Qexps[0]; /* overflow not possible */
    }

    while (heap_len > 1 || Ai < Alen)
    {
        _fmpz_mpoly_fit_length(&Qcoeffs, &Qexps, allocq, Qlen + 1, 1);

        if (heap_len > 1 && Ai < Alen && Aexps[Ai] == heap[1].exp)
        {
            /* take from both A and heap */
            exp = Aexps[Ai];
            acc_lg = _fmpz_mpoly_get_mpz_signed_uiuiui(acc_sm, Acoeffs[Ai], acc);
            Ai++;
        }
        else if (heap_len > 1 && (Ai >= Alen ||
                           mpoly_monomial_gt1(heap[1].exp, Aexps[Ai], maskhi)))
        {
            /* take only from heap */
            exp = heap[1].exp;
            mpz_set_ui(acc, 0);
            acc_lg = acc;
            acc_sm[2] = acc_sm[1] = acc_sm[0] = 0;

            if (mpoly_monomial_overflows1(exp, mask))
                goto not_sqrt;
        }
        else
        {
            FLINT_ASSERT(Ai < Alen);

            /* take only from A */
            exp = Aexps[Ai];
            acc_lg = _fmpz_mpoly_get_mpz_signed_uiuiui(acc_sm, Acoeffs[Ai], acc);
            Ai++;

            if (!check && mpoly_monomial_gt1(exp3, exp, maskhi))
                break;

            lt_divides = mpoly_monomial_divides1(Qexps + Qlen, exp, Qexps[0], mask);

            goto skip_heap;
        }

        lt_divides = mpoly_monomial_divides1(Qexps + Qlen, exp, Qexps[0], mask);

        /* take nodes from heap with exponent matching exp */

        if (!lt_divides && !check)
        {
            do {
                x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
                do {
                    *store++ = x->i;
                    *store++ = x->j;
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);

            mpz_set_ui(acc, 0);
            acc_lg = acc;
        }
        else if (q_rest_small)
        {
            /* optimization: small coeff arithmetic */

            acc_sm2[2] = acc_sm2[1] = acc_sm2[0] = 0;
            do {
                x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
                do {
                    *store++ = x->i;
                    *store++ = x->j;

                    smul_ppmm(pp[1], pp[0], Qcoeffs[x->i], Qcoeffs[x->j]);
                    pp[2] = FLINT_SIGN_EXT(pp[1]);

                    if (x->i != x->j)
                        sub_dddmmmsss(acc_sm2[2], acc_sm2[1], acc_sm2[0],
                                      acc_sm2[2], acc_sm2[1], acc_sm2[0],
                                      pp[2], pp[1], pp[0]);
                    else
                        sub_dddmmmsss(acc_sm[2], acc_sm[1], acc_sm[0],
                                      acc_sm[2], acc_sm[1], acc_sm[0],
                                      pp[2], pp[1], pp[0]);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);

            add_sssaaaaaa(acc_sm[2], acc_sm[1], acc_sm[0],
                          acc_sm[2], acc_sm[1], acc_sm[0],
                          acc_sm2[2], acc_sm2[1], acc_sm2[0]);
            add_sssaaaaaa(acc_sm[2], acc_sm[1], acc_sm[0],
                          acc_sm[2], acc_sm[1], acc_sm[0],
                          acc_sm2[2], acc_sm2[1], acc_sm2[0]);

            if (mpz_sgn(acc_lg) != 0)
            {
                flint_mpz_add_signed_uiuiui(acc, acc_lg,
                                              acc_sm[2], acc_sm[1], acc_sm[0]);
                acc_lg = acc;
                acc_sm[2] = acc_sm[1] = acc_sm[0] = 0;
            }
        }
        else
        {
            acc_sm2[2] = acc_sm2[1] = acc_sm2[0] = 0;

            /* total is always acc + acc_sm + 2*(acc2 + acc_sm2) */
            mpz_tdiv_q_2exp(acc2, acc_lg, 1);
            mpz_tdiv_r_2exp(acc, acc_lg, 1);

            do {
                x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
                do {
                    fmpz Qi, Qj;

                    *store++ = x->i;
                    *store++ = x->j;

                    Qi = Qcoeffs[x->i];
                    Qj = Qcoeffs[x->j];
                    t = (x->i != x->j) ? acc2 : acc;

                    if (!COEFF_IS_MPZ(Qi) && !COEFF_IS_MPZ(Qj))
                    {
                        smul_ppmm(pp[1], pp[0], Qi, Qj);
                        pp[2] = FLINT_SIGN_EXT(pp[1]);

                        if (x->i != x->j)
                            sub_dddmmmsss(acc_sm2[2], acc_sm2[1], acc_sm2[0],
                                          acc_sm2[2], acc_sm2[1], acc_sm2[0],
                                          pp[2], pp[1], pp[0]);
                        else
                            sub_dddmmmsss(acc_sm[2], acc_sm[1], acc_sm[0],
                                          acc_sm[2], acc_sm[1], acc_sm[0],
                                          pp[2], pp[1], pp[0]);
                    }
                    else if (!COEFF_IS_MPZ(Qi) && COEFF_IS_MPZ(Qj))
                    {
                        if (Qi < WORD(0))
                            flint_mpz_addmul_ui(t, COEFF_TO_PTR(Qj), -Qi);
                        else
                            flint_mpz_submul_ui(t, COEFF_TO_PTR(Qj), Qi);
                    }
                    else if (COEFF_IS_MPZ(Qi) && !COEFF_IS_MPZ(Qj))
                    {
                        if (Qj < WORD(0))
                            flint_mpz_addmul_ui(t, COEFF_TO_PTR(Qi), -Qj);
                        else
                            flint_mpz_submul_ui(t, COEFF_TO_PTR(Qi), Qj);
                    }
                    else
                    {
                        mpz_submul(t, COEFF_TO_PTR(Qi), COEFF_TO_PTR(Qj));
                    }

                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);

            add_sssaaaaaa(acc_sm[2], acc_sm[1], acc_sm[0],
                          acc_sm[2], acc_sm[1], acc_sm[0],
                          acc_sm2[2], acc_sm2[1], acc_sm2[0]);
            add_sssaaaaaa(acc_sm[2], acc_sm[1], acc_sm[0],
                          acc_sm[2], acc_sm[1], acc_sm[0],
                          acc_sm2[2], acc_sm2[1], acc_sm2[0]);

            flint_mpz_add_signed_uiuiui(acc, acc, acc_sm[2], acc_sm[1], acc_sm[0]);
            mpz_addmul_ui(acc, acc2, 2);
            acc_lg = acc;
            acc_sm[2] = acc_sm[1] = acc_sm[0] = 0;
        }

        /* process nodes taken from the heap */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            /* should we go right */
            if (j < i)
            {
                x = chain[i];
                x->i = i;
                x->j = j + 1;
                x->next = NULL;

                if (check || !mpoly_monomial_gt1(exp3, Qexps[x->i] + Qexps[x->j], maskhi))
                {
                    _mpoly_heap_insert1(heap, Qexps[x->i] + Qexps[x->j], x,
                                             &next_loc, &heap_len, maskhi);
                }
            }
        }

    skip_heap:

        /* try to divide accumulated term by leading term */
        if (!check && !lt_divides)
            continue;

        if (mpz_sgn(acc_lg) == 0)
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
                nhi = MPN_LEFT_SHIFT_HI(d1, d0, lc_norm);
                nlo = d0 << lc_norm;
                udiv_qrnnd_preinv(qq, rr, nhi, nlo, lc_n, lc_i);

                if (rr != 0)
                    goto not_sqrt;

                if (qq == 0)
                    continue;

                if (qq <= COEFF_MAX)
                {
                    _fmpz_demote(Qcoeffs + Qlen);
                    Qcoeffs[Qlen] = qq;
                    if (ds != 0)
                        Qcoeffs[Qlen] = -Qcoeffs[Qlen];
                }
                else
                {
                    q_rest_small = 0;
                    if (ds == 0)
                        fmpz_set_ui(Qcoeffs + Qlen, qq);
                    else
                        fmpz_neg_ui(Qcoeffs + Qlen, qq);
                }
            }
            else
            {
                flint_mpz_add_signed_uiuiui(acc, acc_lg, acc_sm[2], acc_sm[1], acc_sm[0]);
                goto large_lt_divides;
            }
        }
        else
        {
            flint_mpz_add_signed_uiuiui(acc, acc_lg, acc_sm[2], acc_sm[1], acc_sm[0]);

            if (mpz_sgn(acc) == 0)
                continue;

            if (!lt_divides)
                goto not_sqrt;

        large_lt_divides:

            t = _fmpz_promote(Qcoeffs + Qlen);
            if (lc_abs > 0)
                flint_mpz_fdiv_qr_ui(t, r, acc, lc_abs);
            else
                mpz_fdiv_qr(t, r, acc, lc_lg);

            _fmpz_demote_val(Qcoeffs + Qlen);
            q_rest_small = q_rest_small && !COEFF_IS_MPZ(Qcoeffs[Qlen]);

            if (mpz_sgn(r) != 0)
            {
                Qlen++;
                goto not_sqrt;
            }
        }

        if (Qlen >= heap_alloc)
        {
            /* run some tests if the square root is getting long */
            if (Qlen > Alen && _is_proved_not_square(
                            ++heuristic_count, &heuristic_p, heuristic_state,
                                                Acoeffs, Aexps, Alen, bits, mctx))
            {
                Qlen++; /* for demotion */
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

        /* put (Qlen, 1) in heap */
        i = Qlen;
        x = chain[i];
        x->i = i;
        x->j = 1;
        x->next = NULL;

        if (check || !mpoly_monomial_gt1(exp3, Qexps[i] + Qexps[1], maskhi))
        {
           _mpoly_heap_insert1(heap, Qexps[i] + Qexps[1], x,
                                                 &next_loc, &heap_len, maskhi);
        }

        Qlen++;
    }

    /* divide extra factor of 2 back out of leading coefficient */
    fmpz_fdiv_q_2exp(Qcoeffs + 0, Qcoeffs + 0, 1);

cleanup:

    flint_randclear(heuristic_state);

    mpz_clear(r);
    mpz_clear(acc);
    mpz_clear(acc2);

    (*polyq) = Qcoeffs;
    (*expq) = Qexps;

    flint_free(heap);
    flint_free(chain);
    flint_free(store_base);
    for (i = 0; i < exp_alloc; i++)
        flint_free(chain_nodes[i]);

    /* return sqrt poly length, or zero if not a square root */
    return Qlen;

not_sqrt:
    for (i = 0; i < Qlen; i++)
        _fmpz_demote(Qcoeffs + i);
    Qlen = 0;
    goto cleanup;
}


slong _fmpz_mpoly_sqrt_heap(
    fmpz ** polyq, ulong ** expq, slong * allocq,
    const fmpz * Acoeffs, const ulong * Aexps, slong Alen,
    flint_bitcnt_t bits,
    const mpoly_ctx_t mctx,
    int check)
{
    slong N = mpoly_words_per_exp(bits, mctx);
    ulong * cmpmask;
    slong i, j, Qlen, Ai;
    slong next_loc;
    slong heap_len = 1, heap_alloc;
    int exp_alloc;
    mpoly_heap_s * heap;
    mpoly_heap_t * chain_nodes[64];
    mpoly_heap_t ** chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    fmpz * Qcoeffs = *polyq;
    ulong * Qexps = *expq;
    ulong * exp, * exp3;
    ulong * exps[64];
    ulong ** exp_list;
    slong exp_next;
    ulong mask;
    mpz_t r, acc, acc2;
    mpz_srcptr acc_lg;
    mpz_ptr t;
    ulong acc_sm[3], acc_sm2[3], pp[3];
    int halves, use_heap, lt_divides, q_rest_small;
    flint_rand_t heuristic_state;
    mp_limb_t heuristic_p = UWORD(1) << (SMALL_FMPZ_BITCOUNT_MAX);
    int heuristic_count = 0;
    ulong lc_abs = 0; /* 2*sqrt(lc) if it fits in ulong, otherwise 0 */
    ulong lc_norm = 0;
    ulong lc_n = 0;
    ulong lc_i = 0;
    mpz_ptr lc_lg = NULL; /* 2*sqrt(lc) if it is large */
    TMP_INIT;

    if (N == 1)
        return _fmpz_mpoly_sqrt_heap1(polyq, expq, allocq,
                                        Acoeffs, Aexps, Alen, bits, mctx, check);

    TMP_START;

    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, bits, mctx);

    flint_randinit(heuristic_state);

    mpz_init(r);
    mpz_init(acc);
    mpz_init(acc2);

    /* alloc array of heap nodes which can be chained together */
    next_loc = 2*sqrt(Alen) + 4;   /* something bigger than heap can ever be */
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

    Qlen = 0;

    /* "insert" (-1, 1, Aexps[0]) into "heap" */
    Ai = 1;

    /* compute first term */
    if (!fmpz_is_square(Acoeffs + 0))
        goto not_sqrt;

    _fmpz_mpoly_fit_length(&Qcoeffs, &Qexps, allocq, Qlen + 1, 1);

    fmpz_sqrt(Qcoeffs + 0, Acoeffs + 0);
    Qlen++;

    /* multiply by 2, we revert this at the end */
    fmpz_mul_2exp(Qcoeffs + 0, Qcoeffs + 0, 1);

    /* q_rest_small means Qcoeffs[1] ... Qcoeffs[Qlen-1] are small */
    q_rest_small = 1;

    if (fmpz_abs_fits_ui(Qcoeffs + 0))
    {
        lc_abs = fmpz_get_ui(Qcoeffs + 0);
        lc_norm = flint_clz(lc_abs);
        lc_n = lc_abs << lc_norm;
        lc_i = n_preinvert_limb_prenorm(lc_n);
    }
    else
    {
        lc_lg = COEFF_TO_PTR(Qcoeffs[0]);
    }

    if (bits <= FLINT_BITS)
        halves = mpoly_monomial_halves(Qexps + 0, Aexps + 0, N, mask);
    else
        halves = mpoly_monomial_halves_mp(Qexps + 0, Aexps + 0, N, bits);

    if (!halves)
        goto not_sqrt; /* exponent is not square */

    /* optimisation, compute final term */
    {
        if (!fmpz_is_square(Acoeffs + Alen - 1))
            goto not_sqrt;

        if (bits <= FLINT_BITS)
            halves = mpoly_monomial_halves(exp3, Aexps + (Alen - 1)*N, N, mask);
        else
            halves = mpoly_monomial_halves_mp(exp3, Aexps + (Alen - 1)*N, N, bits);

        if (!halves)
            goto not_sqrt; /* exponent is not square */

        if (bits <= FLINT_BITS)
            mpoly_monomial_add(exp3, exp3, Qexps + 0, N);
        else
            mpoly_monomial_add_mp(exp3, exp3, Qexps + 0, N);
    }

    while (heap_len > 1 || Ai < Alen)
    {
        _fmpz_mpoly_fit_length(&Qcoeffs, &Qexps, allocq, Qlen + 1, N);

        if (heap_len > 1 && Ai < Alen &&
                            mpoly_monomial_equal(Aexps + N*Ai, heap[1].exp, N))
        {
            /* take from both A and heap */
            mpoly_monomial_set(exp, Aexps + N*Ai, N);
            acc_lg = _fmpz_mpoly_get_mpz_signed_uiuiui(acc_sm, Acoeffs[Ai], acc);
            Ai++;
            use_heap = 1;
        }
        else if (heap_len > 1 && (Ai >= Alen ||
                     mpoly_monomial_lt(Aexps + N*Ai, heap[1].exp, N, cmpmask)))
        {
            /* take only from heap */
            mpoly_monomial_set(exp, heap[1].exp, N);
            mpz_set_ui(acc, 0);
            acc_lg = acc;
            acc_sm[2] = acc_sm[1] = acc_sm[0] = 0;

            if (bits <= FLINT_BITS ? mpoly_monomial_overflows(exp, N, mask)
                                   : mpoly_monomial_overflows_mp(exp, N, bits))
                goto not_sqrt;

            use_heap = 1;
        }
        else
        {
            FLINT_ASSERT(Ai < Alen);

            /* take only from A */
            mpoly_monomial_set(exp, Aexps + N*Ai, N);
            acc_lg = _fmpz_mpoly_get_mpz_signed_uiuiui(acc_sm, Acoeffs[Ai], acc);
            Ai++;

            if (!check && mpoly_monomial_gt(exp3, exp, N, cmpmask))
                break;

            use_heap = 0;
        }

        if (bits <= FLINT_BITS)
            lt_divides = mpoly_monomial_divides(Qexps + N*Qlen,
                                                      exp, Qexps + 0, N, mask);
        else
            lt_divides = mpoly_monomial_divides_mp(Qexps + N*Qlen,
                                                      exp, Qexps + 0, N, bits);

        if (!use_heap)
            goto skip_heap;

        /* take nodes from heap with exponent matching exp */

        if (!lt_divides && !check)
        {
            do {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
                do {
                    *store++ = x->i;
                    *store++ = x->j;
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

            mpz_set_ui(acc, 0);
            acc_lg = acc;
        }
        else if (q_rest_small)
        {
            /* optimization: small coeff arithmetic */

            acc_sm2[2] = acc_sm2[1] = acc_sm2[0] = 0;
            do {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
                do {
                    *store++ = x->i;
                    *store++ = x->j;

                    smul_ppmm(pp[1], pp[0], Qcoeffs[x->i], Qcoeffs[x->j]);
                    pp[2] = FLINT_SIGN_EXT(pp[1]);

                    if (x->i != x->j)
                        sub_dddmmmsss(acc_sm2[2], acc_sm2[1], acc_sm2[0],
                                      acc_sm2[2], acc_sm2[1], acc_sm2[0],
                                      pp[2], pp[1], pp[0]);
                    else
                        sub_dddmmmsss(acc_sm[2], acc_sm[1], acc_sm[0],
                                      acc_sm[2], acc_sm[1], acc_sm[0],
                                      pp[2], pp[1], pp[0]);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

            add_sssaaaaaa(acc_sm[2], acc_sm[1], acc_sm[0],
                          acc_sm[2], acc_sm[1], acc_sm[0],
                          acc_sm2[2], acc_sm2[1], acc_sm2[0]);
            add_sssaaaaaa(acc_sm[2], acc_sm[1], acc_sm[0],
                          acc_sm[2], acc_sm[1], acc_sm[0],
                          acc_sm2[2], acc_sm2[1], acc_sm2[0]);

            if (mpz_sgn(acc_lg) != 0)
            {
                flint_mpz_add_signed_uiuiui(acc, acc_lg,
                                              acc_sm[2], acc_sm[1], acc_sm[0]);
                acc_lg = acc;
                acc_sm[2] = acc_sm[1] = acc_sm[0] = 0;
            }
        }
        else
        {
            acc_sm2[2] = acc_sm2[1] = acc_sm2[0] = 0;

            /* total is always acc + acc_sm + 2*(acc2 + acc_sm2) */
            mpz_tdiv_q_2exp(acc2, acc_lg, 1);
            mpz_tdiv_r_2exp(acc, acc_lg, 1);

            do {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
                do {
                    fmpz Qi, Qj;

                    *store++ = x->i;
                    *store++ = x->j;

                    Qi = Qcoeffs[x->i];
                    Qj = Qcoeffs[x->j];
                    t = (x->i != x->j) ? acc2 : acc;

                    if (!COEFF_IS_MPZ(Qi) && !COEFF_IS_MPZ(Qj))
                    {
                        smul_ppmm(pp[1], pp[0], Qi, Qj);
                        pp[2] = FLINT_SIGN_EXT(pp[1]);

                        if (x->i != x->j)
                            sub_dddmmmsss(acc_sm2[2], acc_sm2[1], acc_sm2[0],
                                          acc_sm2[2], acc_sm2[1], acc_sm2[0],
                                          pp[2], pp[1], pp[0]);
                        else
                            sub_dddmmmsss(acc_sm[2], acc_sm[1], acc_sm[0],
                                          acc_sm[2], acc_sm[1], acc_sm[0],
                                          pp[2], pp[1], pp[0]);
                    }
                    else if (!COEFF_IS_MPZ(Qi) && COEFF_IS_MPZ(Qj))
                    {
                        if (Qi < WORD(0))
                            flint_mpz_addmul_ui(t, COEFF_TO_PTR(Qj), -Qi);
                        else
                            flint_mpz_submul_ui(t, COEFF_TO_PTR(Qj), Qi);
                    }
                    else if (COEFF_IS_MPZ(Qi) && !COEFF_IS_MPZ(Qj))
                    {
                        if (Qj < WORD(0))
                            flint_mpz_addmul_ui(t, COEFF_TO_PTR(Qi), -Qj);
                        else
                            flint_mpz_submul_ui(t, COEFF_TO_PTR(Qi), Qj);
                    }
                    else
                    {
                        mpz_submul(t, COEFF_TO_PTR(Qi), COEFF_TO_PTR(Qj));
                    }
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

            add_sssaaaaaa(acc_sm[2], acc_sm[1], acc_sm[0],
                          acc_sm[2], acc_sm[1], acc_sm[0],
                          acc_sm2[2], acc_sm2[1], acc_sm2[0]);
            add_sssaaaaaa(acc_sm[2], acc_sm[1], acc_sm[0],
                          acc_sm[2], acc_sm[1], acc_sm[0],
                          acc_sm2[2], acc_sm2[1], acc_sm2[0]);

            flint_mpz_add_signed_uiuiui(acc, acc, acc_sm[2], acc_sm[1], acc_sm[0]);
            mpz_addmul_ui(acc, acc2, 2);
            acc_lg = acc;
            acc_sm[2] = acc_sm[1] = acc_sm[0] = 0;
        }

        /* process nodes taken from the heap */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            /* should we go right */
            if (j < i)
            {
                x = chain[i];
                x->i = i;
                x->j = j + 1;
                x->next = NULL;

                if (bits <= FLINT_BITS)
                    mpoly_monomial_add(exp_list[exp_next], Qexps + x->i*N,
                                                          Qexps + x->j*N, N);
                else
                    mpoly_monomial_add_mp(exp_list[exp_next], Qexps + x->i*N,
                                                             Qexps + x->j*N, N);
                if (check || !mpoly_monomial_gt(exp3 + 0, exp_list[exp_next], N, cmpmask))
                {
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                         &next_loc, &heap_len, N, cmpmask);
                }
            }
        }

    skip_heap:

        /* try to divide accumulated term by leading term */
        if (!check && !lt_divides)
            continue;

        if (mpz_sgn(acc_lg) == 0)
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
                nhi = MPN_LEFT_SHIFT_HI(d1, d0, lc_norm);
                nlo = d0 << lc_norm;
                udiv_qrnnd_preinv(qq, rr, nhi, nlo, lc_n, lc_i);

                if (rr != 0)
                    goto not_sqrt;

                if (qq == 0)
                    continue;

                if (qq <= COEFF_MAX)
                {
                    _fmpz_demote(Qcoeffs + Qlen);
                    Qcoeffs[Qlen] = qq;
                    if (ds != 0)
                        Qcoeffs[Qlen] = -Qcoeffs[Qlen];
                }
                else
                {
                    q_rest_small = 0;
                    if (ds == 0)
                        fmpz_set_ui(Qcoeffs + Qlen, qq);
                    else
                        fmpz_neg_ui(Qcoeffs + Qlen, qq);
                }
            }
            else
            {
                flint_mpz_add_signed_uiuiui(acc, acc_lg, acc_sm[2], acc_sm[1], acc_sm[0]);
                goto large_lt_divides;
            }
        }
        else
        {
            flint_mpz_add_signed_uiuiui(acc, acc_lg, acc_sm[2], acc_sm[1], acc_sm[0]);

            if (mpz_sgn(acc) == 0)
                continue;

            if (!lt_divides)
                goto not_sqrt;

        large_lt_divides:

            t = _fmpz_promote(Qcoeffs + Qlen);
            if (lc_abs > 0)
                flint_mpz_fdiv_qr_ui(t, r, acc, lc_abs);
            else
                mpz_fdiv_qr(t, r, acc, lc_lg);

            _fmpz_demote_val(Qcoeffs + Qlen);
            q_rest_small = q_rest_small && !COEFF_IS_MPZ(Qcoeffs[Qlen]);

            if (mpz_sgn(r) != 0)
            {
                Qlen++;
                goto not_sqrt;
            }
        }

        if (Qlen >= heap_alloc)
        {
            /* run some tests if the square root is getting long */
            if (Qlen > Alen && _is_proved_not_square(
                            ++heuristic_count, &heuristic_p, heuristic_state,
                                                Acoeffs, Aexps, Alen, bits, mctx))
            {
                Qlen++; /* for demotion */
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

        /* put (Qlen, 1) in heap */
        i = Qlen;
        x = chain[i];
        x->i = i;
        x->j = 1;
        x->next = NULL;

        if (bits <= FLINT_BITS)
            mpoly_monomial_add(exp_list[exp_next], Qexps + x->i*N,
                                                      Qexps + x->j*N, N);
        else
            mpoly_monomial_add_mp(exp_list[exp_next], Qexps + x->i*N,
                                                         Qexps + x->j*N, N);

        if (check || !mpoly_monomial_gt(exp3 + 0, exp_list[exp_next], N, cmpmask))
        {
            exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
        }

        Qlen++;
    }

    /* divide extra factor of 2 back out of leading coefficient */
    fmpz_fdiv_q_2exp(Qcoeffs + 0, Qcoeffs + 0, 1);

cleanup:

    flint_randclear(heuristic_state);

    mpz_clear(r);
    mpz_clear(acc);
    mpz_clear(acc2);

    (*polyq) = Qcoeffs;
    (*expq) = Qexps;

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
    return Qlen;

not_sqrt:
    for (i = 0; i < Qlen; i++)
        _fmpz_demote(Qcoeffs + i);
    Qlen = 0;
    goto cleanup;
}

int fmpz_mpoly_sqrt_heap(fmpz_mpoly_t Q, const fmpz_mpoly_t A,
                          const fmpz_mpoly_ctx_t ctx, int check)
{
    slong lenq, lenq_est;
    flint_bitcnt_t exp_bits;
    fmpz_mpoly_t T;
    fmpz_mpoly_struct * q;

    if (fmpz_mpoly_is_zero(A, ctx))
    {
        fmpz_mpoly_zero(Q, ctx);
        return 1;
    }

    /* square root fits in A->bits if it exists */
    exp_bits = A->bits;

    /* rought lower estimate on length of square root */
    lenq_est = n_sqrt(A->length);

    if (Q == A)
    {
        fmpz_mpoly_init3(T, lenq_est, exp_bits, ctx);
        q = T;
    }
    else
    {
        fmpz_mpoly_fit_length_reset_bits(Q, lenq_est, exp_bits, ctx);
        q = Q;
    }

    lenq = _fmpz_mpoly_sqrt_heap(&q->coeffs, &q->exps, &q->alloc,
                                 A->coeffs, A->exps, A->length,
                                      exp_bits, ctx->minfo, check);
    if (Q == A)
    {
        fmpz_mpoly_swap(Q, T, ctx);
        fmpz_mpoly_clear(T, ctx);
    }

    _fmpz_mpoly_set_length(Q, lenq, ctx);

    return lenq != 0;
}
