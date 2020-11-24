/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"
#include "fq_zech_mpoly.h"


static int _is_proved_not_square_medprime(
    int count,
    flint_rand_t state,
    const mp_limb_t * Acoeffs,
    const ulong * Aexps,
    slong Alen,
    flint_bitcnt_t Abits,
    const mpoly_ctx_t mctx,
    nmod_t mod)
{
    int success = 0;
    slong i;
    fq_zech_struct eval[1], * t, * alphas, ** alpha_ptrs;
    fq_zech_ctx_t fqctx;
    fmpz_t p;
    slong edeg, max_degree = n_flog(1000000, mod.n);
    TMP_INIT;

    FLINT_ASSERT(Alen > 0);

    edeg = (max_degree + count - 2)/2;
    edeg = FLINT_MAX(2, edeg);
    if (edeg > max_degree)
        return 0;

    fmpz_init_set_ui(p, mod.n);
    fq_zech_ctx_init(fqctx, p, edeg, "#");
    fq_zech_init(eval, fqctx);

    TMP_START;

    alphas = (fq_zech_struct *) TMP_ALLOC(mctx->nvars*sizeof(fq_zech_struct));
    alpha_ptrs = (fq_zech_struct **) TMP_ALLOC(mctx->nvars*sizeof(fq_zech_struct *));
    for (i = 0; i < mctx->nvars; i++)
    {
        alpha_ptrs[i] = alphas + i;
        fq_zech_init(alphas + i, fqctx);
    }

    t = (fq_zech_struct *) TMP_ALLOC(Alen*sizeof(fq_zech_struct));
    for (i = 0; i < Alen; i++)
    {
        fq_zech_init(t + i, fqctx);
        fq_zech_set_ui(t + i, Acoeffs[i], fqctx);
    }

    /* try at most 3*count evaluations */
    count *= 3;

next_p:

    for (i = 0; i < mctx->nvars; i++)
        fq_zech_rand(alphas + i, state, fqctx);

    _fq_zech_mpoly_eval_all_fq_zech(eval, t, Aexps, Alen, Abits,
                                                      alpha_ptrs, mctx, fqctx);

    success = !fq_zech_is_square(eval, fqctx);

    if (!success && --count >= 0)
        goto next_p;

    fmpz_clear(p);
    fq_zech_clear(eval, fqctx);
    fq_zech_ctx_clear(fqctx);

    TMP_END;

    return success;
}


/* try to prove that A is not a square */
static int _is_proved_not_square(
    int count,
    flint_rand_t state,
    const mp_limb_t * Acoeffs,
    const ulong * Aexps,
    slong Alen,
    flint_bitcnt_t Abits,
    const mpoly_ctx_t mctx,
    nmod_t mod)
{
    int tries_left, success = 0;
    slong i, N = mpoly_words_per_exp(Abits, mctx);
    mp_limb_t eval, * alphas;
    ulong * t;
    TMP_INIT;

    FLINT_ASSERT(Alen > 0);

    TMP_START;
    t = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    if (count == 1)
    {
        success = mpoly_is_proved_not_square(Aexps, Alen, Abits, N, t);
        if (success)
            goto cleanup;
    }

    tries_left = 3*count;

    alphas = (mp_limb_t *) TMP_ALLOC(mctx->nvars*sizeof(mp_limb_t));

next_p:

    for (i = 0; i < mctx->nvars; i++)
        alphas[i] = n_urandint(state, mod.n);

    eval = _nmod_mpoly_eval_all_ui(Acoeffs, Aexps, Alen, Abits, alphas, mctx, mod);

    success = n_jacobi_unsigned(eval, mod.n) < 0;

    if (!success && --tries_left >= 0)
        goto next_p;

cleanup:

    TMP_END;

    if (!success)
        success = _is_proved_not_square_medprime(count, state,
                                       Acoeffs, Aexps, Alen, Abits, mctx, mod);
    return success;
}


static int _nmod_mpoly_sqrt_heap1(
    nmod_mpoly_t Q,
    const mp_limb_t * Acoeffs,
    const ulong * Aexps,
    slong Alen,
    flint_bitcnt_t bits,
    const mpoly_ctx_t mctx,
    nmod_t mod)
{
    slong i, j, Qlen, Ai;
    slong next_loc, heap_len = 1, heap_alloc;
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain_nodes[64];
    mpoly_heap_t ** chain;
    slong exp_alloc;
    slong * store, * store_base;
    mpoly_heap_t * x;
    mp_limb_t * Qcoeffs = Q->coeffs;
    ulong * Qexps = Q->exps;
    ulong mask, exp, exp3 = 0;
    ulong maskhi;
    ulong pp1, pp0, acc2, acc1, acc0, tmp2, tmp1, tmp0, lc_minus_inv;
    int lt_divides;
    flint_rand_t heuristic_state;
    int heuristic_count = 0;

    FLINT_ASSERT(mpoly_words_per_exp(bits, mctx) == 1);
    mpoly_get_cmpmask(&maskhi, 1, bits, mctx);

    flint_randinit(heuristic_state);

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

    /* "insert" (-1, 1, Aexps[1]) into "heap" */
    Ai = 1;

    /* compute first term */
    Qlen = 0;
    _nmod_mpoly_fit_length(&Qcoeffs, &Q->coeffs_alloc,
                           &Qexps, &Q->exps_alloc, 1, Qlen + 1);

    Qcoeffs[0] = n_sqrtmod(Acoeffs[0], mod.n);
    if (Qcoeffs[0] == 0)
        goto not_sqrt;

    Qlen = 1;

    /* precompute leading cofficient info */
    lc_minus_inv = mod.n - nmod_inv(nmod_add(Qcoeffs[0], Qcoeffs[0], mod), mod);

    if (!mpoly_monomial_halves1(Qexps + 0, Aexps[0], mask))
        goto not_sqrt; /* exponent is not square */

    /* optimisation, compute final exponent */
    {
        if (0 == n_sqrtmod(Acoeffs[Alen - 1], mod.n))
            goto not_sqrt;

        if (!mpoly_monomial_halves1(&exp3, Aexps[Alen - 1], mask))
            goto not_sqrt; /* exponent is not square */

        exp3 += Qexps[0]; /* overflow not possible */
    }

    while (heap_len > 1 || Ai < Alen)
    {
        _nmod_mpoly_fit_length(&Qcoeffs, &Q->coeffs_alloc,
                               &Qexps, &Q->exps_alloc, 1, Qlen + 1);

        acc2 = acc1 = acc0 = 0;

        if (heap_len > 1 && Ai < Alen && Aexps[Ai] == heap[1].exp)
        {
            /* take from both A and heap */
            exp = Aexps[Ai];
            acc0 = mod.n - Acoeffs[Ai];
            Ai++;
        }
        else if (heap_len > 1 && (Ai >= Alen ||
                           mpoly_monomial_gt1(heap[1].exp, Aexps[Ai], maskhi)))
        {
            /* take only from heap */
            exp = heap[1].exp;
            if (mpoly_monomial_overflows1(exp, mask))
                goto not_sqrt;
        }
        else
        {
            FLINT_ASSERT(Ai < Alen);

            /* take only from A */
            exp = Aexps[Ai];
            acc0 = mod.n - Acoeffs[Ai];
            Ai++;

            goto skip_heap;
        }

        /* total is always acc + 2*tmp */
        tmp2 = tmp1 = tmp0 = 0;

        do {
            x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
            do {
                *store++ = x->i;
                *store++ = x->j;

                umul_ppmm(pp1, pp0, Qcoeffs[x->i], Qcoeffs[x->j]);
                if (x->i != x->j)
                    add_sssaaaaaa(tmp2, tmp1, tmp0, tmp2, tmp1, tmp0,
                                                   UWORD(0), pp1, pp0);
                else
                    add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0,
                                                   UWORD(0), pp1, pp0);
            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && heap[1].exp == exp);

        add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, tmp2, tmp1, tmp0);
        add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, tmp2, tmp1, tmp0);
        NMOD_RED3(acc0, acc2, acc1, acc0, mod);

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

                _mpoly_heap_insert1(heap, Qexps[x->i] + Qexps[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
            }
        }

        if (acc0 == 0)
            continue;

    skip_heap:

        lt_divides = mpoly_monomial_divides1(Qexps + Qlen, exp, Qexps[0], mask);
        if (!lt_divides)
            goto not_sqrt;

        Qcoeffs[Qlen] = nmod_mul(acc0, lc_minus_inv, mod);

        if (Qlen >= heap_alloc)
        {
            /* run some tests if the square root is getting long */
            if (Qlen > Alen && _is_proved_not_square(
                                        ++heuristic_count, heuristic_state,
                                           Acoeffs, Aexps, Alen, bits, mctx, mod))
            {
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

       _mpoly_heap_insert1(heap, Qexps[i] + Qexps[1], x,
                                                 &next_loc, &heap_len, maskhi);

        Qlen++;
    }

cleanup:

    flint_randclear(heuristic_state);

    Q->coeffs = Qcoeffs;
    Q->exps = Qexps;
    Q->length = Qlen;

    flint_free(heap);
    flint_free(chain);
    flint_free(store_base);
    for (i = 0; i < exp_alloc; i++)
        flint_free(chain_nodes[i]);

    return Qlen > 0;

not_sqrt:

    Qlen = 0;
    goto cleanup;
}

static int _nmod_mpoly_sqrt_heap(
    nmod_mpoly_t Q,
    const mp_limb_t * Acoeffs,
    const ulong * Aexps,
    slong Alen,
    flint_bitcnt_t bits,
    const mpoly_ctx_t mctx,
    nmod_t mod)
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
    mp_limb_t * Qcoeffs = Q->coeffs;
    ulong * Qexps = Q->exps;
    ulong * exp, * exp3;
    ulong * exps[64];
    ulong ** exp_list;
    slong exp_next;
    ulong mask;
    ulong pp1, pp0, acc2, acc1, acc0, tmp2, tmp1, tmp0, lc_minus_inv;
    int lt_divides, halves;
    flint_rand_t heuristic_state;
    int heuristic_count = 0;
    TMP_INIT;

    if (N == 1)
        return _nmod_mpoly_sqrt_heap1(Q, Acoeffs, Aexps, Alen, bits, mctx, mod);

    TMP_START;

    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, bits, mctx);

    flint_randinit(heuristic_state);

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

    /* "insert" (-1, 1, Aexps[0]) into "heap" */
    Ai = 1;

    /* compute first term */
    Qlen = 0;
    _nmod_mpoly_fit_length(&Qcoeffs, &Q->coeffs_alloc,
                           &Qexps, &Q->exps_alloc, N, Qlen + 1);

    Qcoeffs[0] = n_sqrtmod(Acoeffs[0], mod.n);
    if (Qcoeffs[0] == 0)
        goto not_sqrt;

    Qlen = 1;

    /* precompute leading cofficient info */
    lc_minus_inv = mod.n - nmod_inv(nmod_add(Qcoeffs[0], Qcoeffs[0], mod), mod);

    if (bits <= FLINT_BITS)
        halves = mpoly_monomial_halves(Qexps + 0, Aexps + 0, N, mask);
    else
        halves = mpoly_monomial_halves_mp(Qexps + 0, Aexps + 0, N, bits);

    if (!halves)
        goto not_sqrt; /* exponent is not square */

    /* optimisation, compute final term */
    {
        if (0 == n_sqrtmod(Acoeffs[Alen - 1], mod.n))
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
        _nmod_mpoly_fit_length(&Qcoeffs, &Q->coeffs_alloc,
                               &Qexps, &Q->exps_alloc, N, Qlen + 1);

        acc2 = acc1 = acc0 = 0;

        if (heap_len > 1 && Ai < Alen &&
            mpoly_monomial_equal(Aexps + N*Ai, heap[1].exp, N))
        {
            /* take from both A and heap */
            mpoly_monomial_set(exp, Aexps + N*Ai, N);
            acc0 = mod.n - Acoeffs[Ai];
            Ai++;
        }
        else if (heap_len > 1 && (Ai >= Alen || mpoly_monomial_lt(
                                       Aexps + N*Ai, heap[1].exp, N, cmpmask)))
        {
            /* take only from heap */
            mpoly_monomial_set(exp, heap[1].exp, N);

            if (bits <= FLINT_BITS ? mpoly_monomial_overflows(exp, N, mask)
                                   : mpoly_monomial_overflows_mp(exp, N, bits))
                goto not_sqrt;
        }
        else
        {
            FLINT_ASSERT(Ai < Alen);

            /* take only from A */
            mpoly_monomial_set(exp, Aexps + N*Ai, N);
            acc0 = mod.n - Acoeffs[Ai];
            Ai++;

            goto skip_heap;
        }

        /* total is always -acc - 2*tmp */
        tmp2 = tmp1 = tmp0 = 0;
        do {
            exp_list[--exp_next] = heap[1].exp;
            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
            do {
                *store++ = x->i;
                *store++ = x->j;

                FLINT_ASSERT(x->i > 0 && x->j > 0);
                umul_ppmm(pp1, pp0, Qcoeffs[x->i], Qcoeffs[x->j]);
                if (x->i != x->j)
                    add_sssaaaaaa(tmp2, tmp1, tmp0, tmp2, tmp1, tmp0,
                                                       UWORD(0), pp1, pp0);
                else
                    add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0,
                                                       UWORD(0), pp1, pp0);
            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

        add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, tmp2, tmp1, tmp0);
        add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, tmp2, tmp1, tmp0);
        NMOD_RED3(acc0, acc2, acc1, acc0, mod);

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
                    mpoly_monomial_add(exp_list[exp_next], Qexps + N*x->i,
                                                            Qexps + N*x->j, N);
                else
                    mpoly_monomial_add_mp(exp_list[exp_next], Qexps + N*x->i,
                                                            Qexps + N*x->j, N);

                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
            }
        }

        if (acc0 == 0)
            continue;

    skip_heap:

        if (bits <= FLINT_BITS)
            lt_divides = mpoly_monomial_divides(Qexps + N*Qlen,
                                                exp, Qexps + N*0, N, mask);
        else
            lt_divides = mpoly_monomial_divides_mp(Qexps + N*Qlen,
                                                exp, Qexps + N*0, N, bits);
        if (!lt_divides)
            goto not_sqrt;

        Qcoeffs[Qlen] = nmod_mul(acc0, lc_minus_inv, mod);

        if (Qlen >= heap_alloc)
        {
            /* run some tests if the square root is getting long */
            if (Qlen > Alen && _is_proved_not_square(
                                          ++heuristic_count, heuristic_state,
                                        Acoeffs, Aexps, Alen, bits, mctx, mod))
            {
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

        exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                         &next_loc, &heap_len, N, cmpmask);

        Qlen++;
    }

cleanup:

    flint_randclear(heuristic_state);

    Q->coeffs = Qcoeffs;
    Q->exps = Qexps;
    Q->length = Qlen;

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

    return Qlen > 0;

not_sqrt:
    Qlen = 0;
    goto cleanup;
}

int nmod_mpoly_sqrt_heap(nmod_mpoly_t Q, const nmod_mpoly_t A,
                                                    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong lenq_est;

    if ((ctx->mod.n % 2) == 0)
    {
        flint_bitcnt_t bits = A->bits;
        mp_limb_t * Aexps = A->exps;
        slong Alen = A->length;
        slong i, N = mpoly_words_per_exp(bits, ctx->minfo);
        ulong mask = (bits <= FLINT_BITS) ? mpoly_overflow_mask_sp(bits) : 0;

        if (ctx->mod.n != 2)
            flint_throw(FLINT_IMPINV, "nmod_mpoly_sqrt_heap: "
                        "cannot compute sqrt modulo %wd*%wd", 2, ctx->mod.n/2);

        if (Q != A)
        {
            nmod_mpoly_fit_length_reset_bits(Q, Alen, bits, ctx);
            for (i = 0; i < Alen; i++)
                Q->coeffs[i] = 1;
        }

        for (i = 0; i < Alen; i++)
        {
            if (bits <= FLINT_BITS ?
                !mpoly_monomial_halves(Q->exps + N*i, Aexps + N*i, N, mask) :
                !mpoly_monomial_halves_mp(Q->exps + N*i, Aexps + N*i, N, bits))
            {
                Q->length = 0;
                return 0;
            }
        }

        Q->length = Alen;
        return 1;
    }

    if (nmod_mpoly_is_zero(A, ctx))
    {
        nmod_mpoly_zero(Q, ctx);
        return 1;
    }

    lenq_est = n_sqrt(A->length);

    if (Q == A)
    {
        nmod_mpoly_t T;
        nmod_mpoly_init3(T, lenq_est, A->bits, ctx);
        success = _nmod_mpoly_sqrt_heap(T, A->coeffs, A->exps, A->length,
                                                A->bits, ctx->minfo, ctx->mod);
        nmod_mpoly_swap(Q, T, ctx);
        nmod_mpoly_clear(T, ctx);
    }
    else
    {
        nmod_mpoly_fit_length_reset_bits(Q, lenq_est, A->bits, ctx);
        success = _nmod_mpoly_sqrt_heap(Q, A->coeffs, A->exps, A->length,
                                                A->bits, ctx->minfo, ctx->mod);
    }

    return success;
}

