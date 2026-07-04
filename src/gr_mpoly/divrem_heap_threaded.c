/*
    Copyright (C) 2017 William Hart
    Copyright (C) 2018 Daniel Schultz
    Copyright (C) 2019 Daniel Schultz
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdint.h>
#include <string.h>
#include "thread_pool.h"
#include "thread_support.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mpoly.h"
#include "mpoly.h"
#include "gr_generic.h"
#include "gr_mpoly.h"
#include "threaded_divmod.h"

/*
    Multithreaded gr_mpoly_divrem_heap / gr_mpoly_div / gr_mpoly_div_weak /
    gr_mpoly_divrem_weak, sharing the chunked producer/consumer pipeline
    (gr_mpoly_ts, _gr_mpoly_mulsub_stripe*, the chunk/base structs,
    chunk_mulsub, trychunk, worker_loop) defined in divides_heap_threaded.c
    and declared in threaded_divmod.h.

    The only new pieces here are:

        - _gr_mpoly_divrem_stripe1 / _gr_mpoly_divrem_stripe: the bounded
          ("stripe") analogue of _gr_mpoly_divrem_heap1 / _gr_mpoly_divrem_heap
          (src/gr_mpoly/divrem_heap.c), i.e. exactly the transformation that
          already turns _gr_mpoly_divides_heap1/heap into
          _gr_mpoly_divides_stripe1/stripe in divides_heap_threaded.c: add an
          S->emin lower bound below which heap nodes are never (re)inserted.
          Unlike the divides stripe leaves, these never abort on a
          non-reducible term -- it is routed to a local remainder output
          instead -- so they always fully drain their input.

        - _gr_mpoly_divrem_heap_threaded_pool and the four public entry
          points. This mirrors _gr_mpoly_divides_heap_threaded_pool /
          gr_mpoly_divides_heap_threaded, with two differences forced by
          "does not require exactness": mpoly_divides_select_exps's failure
          return (which only detects "provably not exact") is treated as
          "could not find a clean partition, fall back to the serial
          engine" rather than GR_DOMAIN, and the initial burst of leading
          terms routes non-reducible terms to R instead of failing.
*/

/* gather product node (Bcoeff[i], Qcoeff[j]) into dot_a/dot_b; record
   dividend node. Identical in shape to STRIPE_DIVIDES_FETCH_NODE in
   divides_heap_threaded.c (kept as a separate, tiny macro rather than
   shared, exactly as the serial divides_heap.c / divrem_heap.c already
   each keep their own copy of this pattern). */
#define STRIPE_DIVREM_FETCH_NODE(SET_SHALLOW) \
    do { \
        store[store_len++] = x->i; \
        store[store_len++] = x->j; \
        if (x->i == -UWORD(1)) \
        { \
            have_dividend = 1; \
            dividend_j = x->j; \
        } \
        else \
        { \
            hind[x->i] |= WORD(1); \
            SET_SHALLOW(dot_a, dot_len, Bcoeff, x->i); \
            SET_SHALLOW(dot_b, dot_len, Qcoeff, x->j); \
            dot_len++; \
        } \
    } while (0)

/* append acc (or rem) as a new remainder term at monomial exp/expN */
#define STRIPE_ADD_W_TERM1(val) \
    do { \
        _gr_mpoly_fit_length(W_coeff, W_alloc, W_exp, W_exps_alloc, 1, Wlen + 1, ctx); \
        status |= gr_set(GR_ENTRY(*W_coeff, Wlen, sz), val, cctx); \
        (*W_exp)[Wlen] = exp; \
        Wlen++; \
    } while (0)

#define STRIPE_ADD_W_TERMN(val) \
    do { \
        _gr_mpoly_fit_length(W_coeff, W_alloc, W_exp, W_exps_alloc, N, Wlen + 1, ctx); \
        status |= gr_set(GR_ENTRY(*W_coeff, Wlen, sz), val, cctx); \
        mpoly_monomial_set(*W_exp + Wlen*N, exp, N); \
        Wlen++; \
    } while (0)

/* exact-coefficient path: same as STRIPE_DIVIDES_COEFF in
   divides_heap_threaded.c, except for the extra "unchecked" branch used by
   gr_mpoly_divexact_heap_threaded (S->unchecked implies S->nonfield == 0,
   i.e. the two are mutually exclusive -- see the base struct comment). */
#define STRIPE_DIVREM_QCOEFF(acc) \
    (lc_is_one  ? gr_set(GR_ENTRY(Qcoeff, Qlen, sz), acc, cctx) \
     : lc_is_unit ? gr_mul(GR_ENTRY(Qcoeff, Qlen, sz), acc, lc_inv, cctx) \
     : S->unchecked ? gr_divexact(GR_ENTRY(Qcoeff, Qlen, sz), acc, Bcoeff, cctx) \
                    : gr_div(GR_ENTRY(Qcoeff, Qlen, sz), acc, Bcoeff, cctx))

/*
    Decide the fate of the accumulated coefficient acc at monomial exp:
    commit a quotient term (bumping Qlen after the caller inserts the
    follow-up heap node) or route val to the remainder via ADD_W_TERM.
    Mirrors DIVREM_COEFF_STEP in gr_mpoly/divrem_heap.c, restricted to a
    single term (the heap bookkeeping around it is identical for a
    committed or a non-committed term, so it lives outside this macro,
    exactly as in _gr_mpoly_divides_stripe1/stripe).
*/
#define STRIPE_DIVREM_COEFF_STEP(ADD_W_TERM) \
{ \
    if (!lt_divides) \
    { \
        ADD_W_TERM(acc); \
        commit = 0; \
    } \
    else if (S->nonfield) \
    { \
        cstatus = gr_euclidean_divrem(GR_ENTRY(Qcoeff, Qlen, sz), rem, acc, Bcoeff, cctx); \
        if (cstatus != GR_SUCCESS) { status |= cstatus; goto unable; } \
        /* an unknown zero-status is treated as "not zero" (i.e. kept in \
           W): we only need the euclidean division itself to succeed, not \
           to know for certain whether its remainder happens to vanish. */ \
        if (gr_is_zero(rem, cctx) != T_TRUE) \
            ADD_W_TERM(rem); \
        commit = (gr_is_zero(GR_ENTRY(Qcoeff, Qlen, sz), cctx) != T_TRUE); \
    } \
    else \
    { \
        cstatus = STRIPE_DIVREM_QCOEFF(acc); \
        if (cstatus != GR_SUCCESS) \
        { \
            status |= cstatus; goto unable; \
        } \
        else if (gr_is_zero(GR_ENTRY(Qcoeff, Qlen, sz), cctx) == T_TRUE) \
        { \
            ADD_W_TERM(acc); \
            commit = 0; \
        } \
        else \
        { \
            commit = 1; \
        } \
    } \
}


/*
    Bounded ("stripe") division-with-remainder, single-word exponent
    version, assuming Alen > 0. This is _gr_mpoly_divrem_heap1
    (gr_mpoly/divrem_heap.c) restricted to exponents >= S->emin, in exactly
    the same way _gr_mpoly_divides_stripe1 restricts _gr_mpoly_divides_heap1.

    Never aborts merely because a term is not divisible by lm(B) (that is
    precisely what routes it to W instead of Q). For the non-weak variant
    (nonfield == 0), an exactly-reducible-by-monomial term whose
    coefficient does not divide exactly *is* a hard abort, exactly as in
    the serial gr_mpoly_divrem_heap kernel (see the doc comment at the top
    of gr_mpoly/divrem_heap.c): *res_status is then GR_DOMAIN. Otherwise
    *res_status is GR_SUCCESS, or GR_UNABLE if some ring operation could
    not decide.
*/
slong _gr_mpoly_divrem_stripe1(
    gr_ptr * Q_coeff, ulong ** Q_exp, slong * Q_alloc, slong * Q_exps_alloc,
    gr_ptr * W_coeff, ulong ** W_exp, slong * W_alloc, slong * W_exps_alloc,
    slong * Wlen_out,
    gr_srcptr Acoeff, const ulong * Aexp, slong Alen,
    gr_srcptr Bcoeff, const ulong * Bexp, slong Blen,
    const _gr_mpoly_stripe_t S, int * res_status)
{
    gr_mpoly_ctx_struct * ctx = S->ctx;
    gr_ctx_struct * cctx = S->cctx;
    slong sz = S->sz;
    int have_fast_dot = S->have_fast_dot;
    int lc_is_one = S->lc_is_one;
    int lc_is_unit = S->lc_is_unit;
    gr_srcptr lc_inv = S->lc_inv;
    flint_bitcnt_t bits = S->bits;
    ulong emin = S->emin[0];
    ulong cmpmask = S->cmpmask[0];
    ulong texp;
    int lt_divides, commit;
    slong i, j, s;
    slong next_loc, heap_len;
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base, store_len;
    mpoly_heap_t * x;
    slong Qlen, Wlen;
    slong Qalloc = *Q_alloc;
    slong Qexps_alloc = *Q_exps_alloc;
    gr_ptr Qcoeff = *Q_coeff;
    ulong * Qexp = *Q_exp;
    ulong exp;
    ulong mask;
    slong * hind;
    gr_ptr acc, pp, rem, dot_a, dot_b;
    slong dot_len;
    int have_dividend, cstatus;
    slong dividend_j = 0;
    int status = GR_SUCCESS;

    FLINT_ASSERT(Alen > 0);
    FLINT_ASSERT(Blen > 0);
    FLINT_ASSERT(S->N == 1);

    next_loc = Blen + 4;

    i = 0;
    hind = (slong *)(S->big_mem + i);
    i += Blen*sizeof(slong);
    store = store_base = (slong *) (S->big_mem + i);
    i += 2*Blen*sizeof(slong);
    heap = (mpoly_heap1_s *)(S->big_mem + i);
    i += (Blen + 1)*sizeof(mpoly_heap1_s);
    chain = (mpoly_heap_t *)(S->big_mem + i);
    i += Blen*sizeof(mpoly_heap_t);
    FLINT_ASSERT(i <= S->big_mem_alloc);

    GR_TMP_INIT3(acc, pp, rem, cctx);
    dot_a = flint_malloc(2 * Blen * sz);
    dot_b = GR_ENTRY(dot_a, Blen, sz);

    for (i = 0; i < Blen; i++)
        hind[i] = 1;

    mask = mpoly_overflow_mask_sp(bits);

    Qlen = 0;
    Wlen = 0;
    s = Blen;

    heap_len = 2;
    x = chain + 0;
    x->i = -WORD(1);
    x->j = 0;
    x->next = NULL;
    heap[1].next = x;
    HEAP_ASSIGN(heap[1], Aexp[0], x);

    FLINT_ASSERT(mpoly_monomial_cmp1(Aexp[0], emin, cmpmask) >= 0);

    while (heap_len > 1)
    {
        exp = heap[1].exp;

        /* An exponent window is chosen (by the caller, via
           mpoly_divides_select_exps) from the actual monomials appearing
           in A and B, so this cannot overflow in practice; treat it like
           any other "cannot proceed" condition rather than introducing a
           GR_DOMAIN-shaped abort that the divrem family has no use for. */
        if (mpoly_monomial_overflows1(exp, mask))
        {
            status |= GR_UNABLE;
            goto unable;
        }

        FLINT_ASSERT(mpoly_monomial_cmp1(exp, emin, cmpmask) >= 0);

        _gr_mpoly_fit_length(&Qcoeff, &Qalloc, &Qexp, &Qexps_alloc, 1, Qlen + 1, ctx);

        lt_divides = mpoly_monomial_divides1(Qexp + Qlen, exp, Bexp[0], mask);

        dot_len = 0;
        have_dividend = 0;
        store_len = 0;

        if (have_fast_dot)
        {
            do
            {
                x = _mpoly_heap_pop1(heap, &heap_len, cmpmask);
                do
                {
                    if (sz == 1)       { STRIPE_DIVREM_FETCH_NODE(SET_SHALLOW1); }
                    else if (sz == 2)  { STRIPE_DIVREM_FETCH_NODE(SET_SHALLOW2); }
                    else if (sz == 4)  { STRIPE_DIVREM_FETCH_NODE(SET_SHALLOW4); }
                    else if (sz == 8)  { STRIPE_DIVREM_FETCH_NODE(SET_SHALLOW8); }
                    else if (sz == 16) { STRIPE_DIVREM_FETCH_NODE(SET_SHALLOW16); }
                    else               { STRIPE_DIVREM_FETCH_NODE(SET_SHALLOW_GENERIC); }
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);

            status |= _gr_vec_dot(acc,
                        have_dividend ? GR_ENTRY(Acoeff, dividend_j, sz) : NULL,
                        1, dot_a, dot_b, dot_len, cctx);
        }
        else
        {
            status |= gr_zero(acc, cctx);
            do
            {
                x = _mpoly_heap_pop1(heap, &heap_len, cmpmask);
                do
                {
                    store[store_len++] = x->i;
                    store[store_len++] = x->j;
                    if (x->i == -UWORD(1))
                    {
                        status |= gr_add(acc, acc, GR_ENTRY(Acoeff, x->j, sz), cctx);
                    }
                    else
                    {
                        hind[x->i] |= WORD(1);
                        status |= gr_mul(pp, GR_ENTRY(Bcoeff, x->i, sz), GR_ENTRY(Qcoeff, x->j, sz), cctx);
                        status |= gr_sub(acc, acc, pp, cctx);
                    }
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);
        }

        /* process nodes taken from the heap */
        while (store_len > 0)
        {
            j = store_base[--store_len];
            i = store_base[--store_len];

            if (i == -WORD(1))
            {
                if (j + 1 < Alen)
                {
                    x = chain + 0;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;

                    FLINT_ASSERT(mpoly_monomial_cmp1(Aexp[x->j], emin, cmpmask) >= 0);

                    _mpoly_heap_insert1(heap, Aexp[x->j], x,
                                                &next_loc, &heap_len, cmpmask);
                }
            }
            else
            {
                /* should we go right? */
                if ((i + 1 < Blen) && (hind[i + 1] == 2*j + 1))
                {
                    x = chain + i + 1;
                    x->i = i + 1;
                    x->j = j;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;

                    texp = Bexp[x->i] + Qexp[x->j];
                    if (mpoly_monomial_cmp1(texp, emin, cmpmask) >= 0)
                        _mpoly_heap_insert1(heap, texp, x, &next_loc, &heap_len, cmpmask);
                    else
                        hind[x->i] |= 1;
                }
                /* should we go up? */
                if (j + 1 == Qlen)
                {
                    s++;
                }
                else if (((hind[i] & 1) == 1) &&
                         ((i == 1) || (hind[i - 1] >= 2*(j + 2) + 1)))
                {
                    x = chain + i;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;

                    texp = Bexp[x->i] + Qexp[x->j];
                    if (mpoly_monomial_cmp1(texp, emin, cmpmask) >= 0)
                        _mpoly_heap_insert1(heap, texp, x, &next_loc, &heap_len, cmpmask);
                    else
                        hind[x->i] |= 1;
                }
            }
        }

        /* an unknown zero-status is treated as "not zero" -- see the
           discussion at gr_mpoly_divexact for why this is safe: we only
           need the coefficient division below to succeed, not to know for
           certain that acc is nonzero. */
        if (gr_is_zero(acc, cctx) == T_TRUE)
            continue;

        commit = 0;
        STRIPE_DIVREM_COEFF_STEP(STRIPE_ADD_W_TERM1);

        if (commit)
        {
            if (s > 1)
            {
                i = 1;
                x = chain + i;
                x->i = i;
                x->j = Qlen;
                x->next = NULL;
                hind[x->i] = 2*(x->j + 1) + 0;

                texp = Bexp[x->i] + Qexp[x->j];
                if (mpoly_monomial_cmp1(texp, emin, cmpmask) >= 0)
                    _mpoly_heap_insert1(heap, texp, x, &next_loc, &heap_len, cmpmask);
                else
                    hind[x->i] |= 1;
            }
            s = 1;
            Qlen++;
        }
    }

    *res_status = status;
    goto cleanup;

unable:
    Qlen = 0;
    Wlen = 0;
    *res_status = status;

cleanup:

    GR_TMP_CLEAR3(acc, pp, rem, cctx);
    flint_free(dot_a);

    *Q_coeff = Qcoeff;
    *Q_exp = Qexp;
    *Q_alloc = Qalloc;
    *Q_exps_alloc = Qexps_alloc;
    *Wlen_out = Wlen;

    return Qlen;
}

/* Multi-word exponent version of the above. */
slong _gr_mpoly_divrem_stripe(
    gr_ptr * Q_coeff, ulong ** Q_exp, slong * Q_alloc, slong * Q_exps_alloc,
    gr_ptr * W_coeff, ulong ** W_exp, slong * W_alloc, slong * W_exps_alloc,
    slong * Wlen_out,
    gr_srcptr Acoeff, const ulong * Aexp, slong Alen,
    gr_srcptr Bcoeff, const ulong * Bexp, slong Blen,
    const _gr_mpoly_stripe_t S, int * res_status)
{
    gr_mpoly_ctx_struct * ctx = S->ctx;
    gr_ctx_struct * cctx = S->cctx;
    slong sz = S->sz;
    int have_fast_dot = S->have_fast_dot;
    int lc_is_one = S->lc_is_one;
    int lc_is_unit = S->lc_is_unit;
    gr_srcptr lc_inv = S->lc_inv;
    flint_bitcnt_t bits = S->bits;
    slong N = S->N;
    int lt_divides, commit;
    slong i, j, s;
    slong next_loc, heap_len;
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base, store_len;
    mpoly_heap_t * x;
    slong Qlen, Wlen;
    slong Qalloc = *Q_alloc;
    slong Qexps_alloc = *Q_exps_alloc;
    gr_ptr Qcoeff = *Q_coeff;
    ulong * Qexp = *Q_exp;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    ulong mask;
    slong * hind;
    gr_ptr acc, pp, rem, dot_a, dot_b;
    slong dot_len;
    int have_dividend, cstatus;
    slong dividend_j = 0;
    int status = GR_SUCCESS;

    FLINT_ASSERT(Alen > 0);
    FLINT_ASSERT(Blen > 0);

    next_loc = Blen + 4;

    i = 0;
    hind = (slong *) (S->big_mem + i);
    i += Blen*sizeof(slong);
    store = store_base = (slong *) (S->big_mem + i);
    i += 2*Blen*sizeof(slong);
    heap = (mpoly_heap_s *)(S->big_mem + i);
    i += (Blen + 1)*sizeof(mpoly_heap_s);
    chain = (mpoly_heap_t *)(S->big_mem + i);
    i += Blen*sizeof(mpoly_heap_t);
    exps = (ulong *)(S->big_mem + i);
    i +=  Blen*N*sizeof(ulong);
    exp_list = (ulong **)(S->big_mem + i);
    i +=  Blen*sizeof(ulong *);
    exp = (ulong *)(S->big_mem + i);
    i +=  N*sizeof(ulong);
    FLINT_ASSERT(i <= S->big_mem_alloc);

    GR_TMP_INIT3(acc, pp, rem, cctx);
    dot_a = flint_malloc(2 * Blen * sz);
    dot_b = GR_ENTRY(dot_a, Blen, sz);

    exp_next = 0;
    for (i = 0; i < Blen; i++)
        exp_list[i] = exps + i*N;

    for (i = 0; i < Blen; i++)
        hind[i] = 1;

    mask = bits <= FLINT_BITS ? mpoly_overflow_mask_sp(bits) : 0;

    Qlen = 0;
    Wlen = 0;
    s = Blen;

    heap_len = 2;
    x = chain + 0;
    x->i = -WORD(1);
    x->j = 0;
    x->next = NULL;
    heap[1].next = x;
    heap[1].exp = exp_list[exp_next++];

    FLINT_ASSERT(mpoly_monomial_cmp(Aexp + N*0, S->emin, N, S->cmpmask) >= 0);

    mpoly_monomial_set(heap[1].exp, Aexp + N*0, N);

    while (heap_len > 1)
    {
        _gr_mpoly_fit_length(&Qcoeff, &Qalloc, &Qexp, &Qexps_alloc, N, Qlen + 1, ctx);

        mpoly_monomial_set(exp, heap[1].exp, N);

        if (bits <= FLINT_BITS)
        {
            if (mpoly_monomial_overflows(exp, N, mask))
            {
                status |= GR_UNABLE;
                goto unable;
            }
            lt_divides = mpoly_monomial_divides(Qexp + N*Qlen, exp, Bexp + N*0, N, mask);
        }
        else
        {
            if (mpoly_monomial_overflows_mp(exp, N, bits))
            {
                status |= GR_UNABLE;
                goto unable;
            }
            lt_divides = mpoly_monomial_divides_mp(Qexp + N*Qlen, exp, Bexp + N*0, N, bits);
        }

        FLINT_ASSERT(mpoly_monomial_cmp(exp, S->emin, N, S->cmpmask) >= 0);

        dot_len = 0;
        have_dividend = 0;
        store_len = 0;

        if (have_fast_dot)
        {
            do
            {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, S->cmpmask);
                do
                {
                    if (sz == 1)       { STRIPE_DIVREM_FETCH_NODE(SET_SHALLOW1); }
                    else if (sz == 2)  { STRIPE_DIVREM_FETCH_NODE(SET_SHALLOW2); }
                    else if (sz == 4)  { STRIPE_DIVREM_FETCH_NODE(SET_SHALLOW4); }
                    else if (sz == 8)  { STRIPE_DIVREM_FETCH_NODE(SET_SHALLOW8); }
                    else if (sz == 16) { STRIPE_DIVREM_FETCH_NODE(SET_SHALLOW16); }
                    else               { STRIPE_DIVREM_FETCH_NODE(SET_SHALLOW_GENERIC); }
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

            status |= _gr_vec_dot(acc,
                        have_dividend ? GR_ENTRY(Acoeff, dividend_j, sz) : NULL,
                        1, dot_a, dot_b, dot_len, cctx);
        }
        else
        {
            status |= gr_zero(acc, cctx);
            do
            {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, S->cmpmask);
                do
                {
                    store[store_len++] = x->i;
                    store[store_len++] = x->j;
                    if (x->i == -UWORD(1))
                    {
                        status |= gr_add(acc, acc, GR_ENTRY(Acoeff, x->j, sz), cctx);
                    }
                    else
                    {
                        hind[x->i] |= WORD(1);
                        status |= gr_mul(pp, GR_ENTRY(Bcoeff, x->i, sz), GR_ENTRY(Qcoeff, x->j, sz), cctx);
                        status |= gr_sub(acc, acc, pp, cctx);
                    }
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));
        }

        /* process nodes taken from the heap */
        while (store_len > 0)
        {
            j = store_base[--store_len];
            i = store_base[--store_len];

            if (i == -WORD(1))
            {
                if (j + 1 < Alen)
                {
                    x = chain + 0;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    mpoly_monomial_set(exp_list[exp_next], Aexp + x->j*N, N);

                    FLINT_ASSERT(mpoly_monomial_cmp(exp_list[exp_next], S->emin, N, S->cmpmask) >= 0);

                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                          &next_loc, &heap_len, N, S->cmpmask);
                }
            }
            else
            {
                /* should we go right? */
                if ((i + 1 < Blen) && (hind[i + 1] == 2*j + 1))
                {
                    x = chain + i + 1;
                    x->i = i + 1;
                    x->j = j;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;

                    mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i,
                                                              Qexp + N*x->j, N);

                    if (mpoly_monomial_cmp(exp_list[exp_next], S->emin, N, S->cmpmask) >= 0)
                        exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                          &next_loc, &heap_len, N, S->cmpmask);
                    else
                        hind[x->i] |= 1;
                }
                /* should we go up? */
                if (j + 1 == Qlen)
                {
                    s++;
                }
                else if (((hind[i] & 1) == 1) &&
                         ((i == 1) || (hind[i - 1] >= 2*(j + 2) + 1)))
                {
                    x = chain + i;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;

                    mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i,
                                                              Qexp + N*x->j, N);

                    if (mpoly_monomial_cmp(exp_list[exp_next], S->emin, N, S->cmpmask) >= 0)
                        exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                          &next_loc, &heap_len, N, S->cmpmask);
                    else
                        hind[x->i] |= 1;
                }
            }
        }

        /* an unknown zero-status is treated as "not zero" -- see the
           discussion at gr_mpoly_divexact for why this is safe: we only
           need the coefficient division below to succeed, not to know for
           certain that acc is nonzero. */
        if (gr_is_zero(acc, cctx) == T_TRUE)
            continue;

        commit = 0;
        STRIPE_DIVREM_COEFF_STEP(STRIPE_ADD_W_TERMN);

        if (commit)
        {
            if (s > 1)
            {
                i = 1;
                x = chain + i;
                x->i = i;
                x->j = Qlen;
                x->next = NULL;
                hind[x->i] = 2*(x->j + 1) + 0;

                mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i, Qexp + N*x->j, N);

                if (mpoly_monomial_cmp(exp_list[exp_next], S->emin, N, S->cmpmask) >= 0)
                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                          &next_loc, &heap_len, N, S->cmpmask);
                else
                    hind[x->i] |= 1;
            }
            s = 1;
            Qlen++;
        }
    }

    *res_status = status;
    goto cleanup;

unable:
    Qlen = 0;
    Wlen = 0;
    *res_status = status;

cleanup:

    GR_TMP_CLEAR3(acc, pp, rem, cctx);
    flint_free(dot_a);

    *Q_coeff = Qcoeff;
    *Q_exp = Qexp;
    *Q_alloc = Qalloc;
    *Q_exps_alloc = Qexps_alloc;
    *Wlen_out = Wlen;

    return Qlen;
}


static int
_gr_mpoly_divrem_serial(gr_mpoly_t Q, gr_mpoly_t R,
    const gr_mpoly_t A, const gr_mpoly_t B, int nonfield, int unchecked,
    gr_mpoly_ctx_t ctx)
{
    /*
        Call the guaranteed-serial kernel directly, NOT the public
        gr_mpoly_div/gr_mpoly_divrem/... dispatchers: those may route large
        instances straight back into this threaded engine, which would
        recurse indefinitely whenever this function is reached as a
        fallback from within the threaded engine itself (A, B unchanged).
    */
    return _gr_mpoly_divrem_mp(Q, R, A, B, nonfield, unchecked, ctx);
}

/*
    R may be NULL, in which case the remainder is not written at all (the
    gr_mpoly_div(_weak)_heap_threaded case) rather than merely discarded --
    this matches how _gr_mpoly_divrem_mp handles R == NULL in the serial
    engine, and lets a chunk skip producing/appending remainder terms
    entirely rather than computing and immediately throwing them away.
*/
static int _gr_mpoly_divrem_heap_threaded_pool(
    gr_mpoly_t Q,
    gr_mpoly_t R,
    const gr_mpoly_t A,
    const gr_mpoly_t B,
    gr_mpoly_ctx_t ctx,
    int nonfield,
    int unchecked,
    const thread_pool_handle * handles,
    slong num_handles)
{
    mpoly_ctx_struct * mctx = GR_MPOLY_MCTX(ctx);
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    int status;
    fmpz_mpoly_ctx_t zctx;
    fmpz_mpoly_t S;
    slong i, N;
    flint_bitcnt_t exp_bits;
    ulong * cmpmask;
    ulong * Aexp, * Bexp;
    int freeAexp, freeBexp;
    worker_arg_struct * worker_args;
    gr_ptr lc_inv;
    divides_heap_base_t H;
    int lc_is_one, lc_is_unit;
    TMP_INIT;

    if (B->length < 2 || A->length < 2)
        return _gr_mpoly_divrem_serial(Q, R, A, B, nonfield, unchecked, ctx);

    TMP_START;

    GR_TMP_INIT(lc_inv, cctx);

    exp_bits = MPOLY_MIN_BITS;
    exp_bits = FLINT_MAX(exp_bits, A->bits);
    exp_bits = FLINT_MAX(exp_bits, B->bits);
    exp_bits = mpoly_fix_bits(exp_bits, mctx);

    N = mpoly_words_per_exp(exp_bits, mctx);
    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, exp_bits, mctx);

    /* ensure input exponents packed to same size as output exponents */
    Aexp = A->exps;
    freeAexp = 0;
    if (exp_bits > A->bits)
    {
        freeAexp = 1;
        Aexp = (ulong *) flint_malloc(N*A->length*sizeof(ulong));
        mpoly_repack_monomials(Aexp, exp_bits, A->exps, A->bits,
                                                        A->length, mctx);
    }

    Bexp = B->exps;
    freeBexp = 0;
    if (exp_bits > B->bits)
    {
        freeBexp = 1;
        Bexp = (ulong *) flint_malloc(N*B->length*sizeof(ulong));
        mpoly_repack_monomials(Bexp, exp_bits, B->exps, B->bits,
                                                    B->length, mctx);
    }

    fmpz_mpoly_ctx_init(zctx, mctx->nvars, mctx->ord);
    fmpz_mpoly_init(S, zctx);

    lc_is_one = lc_is_unit = 0;
    if (!nonfield)
    {
        lc_is_one = (gr_is_one(B->coeffs, cctx) == T_TRUE);
        lc_is_unit = lc_is_one || (gr_inv(lc_inv, B->coeffs, cctx) == GR_SUCCESS);
    }

    /*
        mpoly_divides_select_exps's failure return only ever means "this
        division is provably not exact" -- irrelevant for div/divrem, which
        have no notion of an "impossible" division. Rather than try to
        patch up a partition in that case, just fall back to the serial
        engine (this can only happen on genuinely pathological / huge
        exponent spreads).
    */
    if (mpoly_divides_select_exps(S, zctx, num_handles,
                                   Aexp, A->length, Bexp, B->length, exp_bits))
    {
        status = _gr_mpoly_divrem_serial(Q, R, A, B, nonfield, unchecked, ctx);
        goto cleanup1;
    }

    divides_heap_base_init(H);

    H->polyA->coeffs = A->coeffs;
    H->polyA->exps = Aexp;
    H->polyA->bits = exp_bits;
    H->polyA->length = A->length;
    H->polyA->coeffs_alloc = A->coeffs_alloc;
    H->polyA->exps_alloc = N*A->length;

    H->polyB->coeffs = B->coeffs;
    H->polyB->exps = Bexp;
    H->polyB->bits = exp_bits;
    H->polyB->length = B->length;
    H->polyB->coeffs_alloc = B->coeffs_alloc;
    H->polyB->exps_alloc = N*B->length;

    H->ctx = ctx;
    H->cctx = cctx;
    H->bits = exp_bits;
    H->N = N;
    H->cmpmask = cmpmask;
    H->lc_is_one = lc_is_one;
    H->lc_is_unit = lc_is_unit;
    H->lc_inv = lc_inv;
    H->require_exact = 0;
    H->want_remainder = (R != NULL);
    H->nonfield = nonfield;
    H->unchecked = unchecked;
    H->failed = 0;
    H->have_domain = 0;
    H->have_unable = 0;

    for (i = 0; i + 1 < S->length; i++)
    {
        divides_heap_chunk_struct * L;
        L = (divides_heap_chunk_struct *) flint_malloc(
                                            sizeof(divides_heap_chunk_struct));
        L->ma = 0;
        L->mq = 0;
        L->emax = S->exps + N*i;
        L->emin = S->exps + N*(i + 1);
        L->upperclosed = 0;
        L->startidx = B->length;
        L->endidx = B->length;
        L->producer = 0;
        L->Cinited = 0;
        L->lock = -2;
        divides_heap_base_add_chunk(H, L);
    }

    H->head->upperclosed = 1;
    H->head->producer = 1;
    H->cur = H->head;

    /*
        Unlike gr_mpoly_divides_heap_threaded, no "generate an initial
        burst of leading terms" pre-pass is needed here: that pre-pass is
        purely a latency optimisation there (so that H->polyQ has *some*
        content before the thread pool is woken), never a correctness
        requirement -- a chunk that is not waiting on any higher quotient
        terms (the head chunk, here marked as producer from the start) can
        always be fully resolved on first touch, regardless of how many
        terms of A it covers, via the ordinary stripe finalisation path in
        trychunk. So we simply seed both thread-safe arrays empty.
    */
    gr_mpoly_ts_init(H->polyQ, NULL, NULL, 0, H->bits, H->N, cctx);
    if (H->want_remainder)
        gr_mpoly_ts_init(H->polyR, NULL, NULL, 0, H->bits, H->N, cctx);

    /* start the workers */

#if FLINT_USES_PTHREAD
    pthread_mutex_init(&H->mutex, NULL);
#endif

    worker_args = (worker_arg_struct *) flint_malloc((num_handles + 1)
                                                        *sizeof(worker_arg_t));

    for (i = 0; i < num_handles; i++)
    {
        (worker_args + i)->H = H;
        thread_pool_wake(global_thread_pool, handles[i], 0,
                                                 worker_loop, worker_args + i);
    }
    (worker_args + num_handles)->H = H;
    worker_loop(worker_args + num_handles);
    for (i = 0; i < num_handles; i++)
        thread_pool_wait(global_thread_pool, handles[i]);

    flint_free(worker_args);

#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&H->mutex);
#endif

    status = divides_heap_base_clear(Q, R, H);

cleanup1:

    fmpz_mpoly_clear(S, zctx);
    fmpz_mpoly_ctx_clear(zctx);

    if (freeAexp)
        flint_free(Aexp);

    if (freeBexp)
        flint_free(Bexp);

    GR_TMP_CLEAR(lc_inv, cctx);

    TMP_END;

    return status;
}


static int _gr_mpoly_divrem_heap_threaded_dispatch(
    gr_mpoly_t Q,
    gr_mpoly_t R,
    const gr_mpoly_t A,
    const gr_mpoly_t B,
    gr_mpoly_ctx_t ctx,
    int nonfield,
    int unchecked)
{
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    thread_pool_handle * handles;
    slong num_handles;
    slong thread_limit;
    int status;

    if (B->length == 0)
        return GR_DOMAIN;

    if (A->length == 0)
    {
        status = gr_mpoly_zero(Q, ctx);
        if (R != NULL)
            status |= gr_mpoly_zero(R, ctx);
        return status;
    }

    if (gr_is_zero(B->coeffs, cctx) != T_FALSE)
        return GR_UNABLE;

    /* fall back to the single-threaded algorithm for small inputs, or when
       the coefficient ring does not allow concurrent operations */
    if (B->length < 2 || A->length < 2 || gr_ctx_is_threadsafe(cctx) != T_TRUE)
        return _gr_mpoly_divrem_serial(Q, R, A, B, nonfield, unchecked, ctx);

    thread_limit = A->length/32;

    num_handles = flint_request_threads(&handles, thread_limit);

    status = _gr_mpoly_divrem_heap_threaded_pool(Q, R, A, B, ctx,
                                    nonfield, unchecked, handles, num_handles);

    flint_give_back_threads(handles, num_handles);

    return status;
}


int gr_mpoly_divrem_heap_threaded(
    gr_mpoly_t Q, gr_mpoly_t R,
    const gr_mpoly_t A, const gr_mpoly_t B,
    gr_mpoly_ctx_t ctx)
{
    return _gr_mpoly_divrem_heap_threaded_dispatch(Q, R, A, B, ctx, 0, 0);
}


int gr_mpoly_divrem_weak_heap_threaded(
    gr_mpoly_t Q, gr_mpoly_t R,
    const gr_mpoly_t A, const gr_mpoly_t B,
    gr_mpoly_ctx_t ctx)
{
    return _gr_mpoly_divrem_heap_threaded_dispatch(Q, R, A, B, ctx, 1, 0);
}


int gr_mpoly_div_heap_threaded(
    gr_mpoly_t Q,
    const gr_mpoly_t A, const gr_mpoly_t B,
    gr_mpoly_ctx_t ctx)
{
    return _gr_mpoly_divrem_heap_threaded_dispatch(Q, NULL, A, B, ctx, 0, 0);
}


int gr_mpoly_div_weak_heap_threaded(
    gr_mpoly_t Q,
    const gr_mpoly_t A, const gr_mpoly_t B,
    gr_mpoly_ctx_t ctx)
{
    return _gr_mpoly_divrem_heap_threaded_dispatch(Q, NULL, A, B, ctx, 1, 0);
}


int gr_mpoly_divexact_heap_threaded(
    gr_mpoly_t Q,
    const gr_mpoly_t A, const gr_mpoly_t B,
    gr_mpoly_ctx_t ctx)
{
    return _gr_mpoly_divrem_heap_threaded_dispatch(Q, NULL, A, B, ctx, 0, 1);
}
