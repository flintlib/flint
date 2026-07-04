/*
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
    GR port of fmpz_mpoly_divides_heap_threaded (equivalently
    nmod_mpoly_divides_heap_threaded).

    The dividend A is split, by monomial, into a number of chunks (using the
    same purely exponent-based selection routine, mpoly_divides_select_exps,
    used by the fmpz/nmod versions).  Each chunk is responsible for producing
    the quotient terms whose monomial lies in its exponent range.  A chunk can
    only start producing terms once all quotient terms of every higher chunk
    are known, since those may contribute (via multiplication by B) to the
    remainder in the current chunk; this dependency is realised by a singly
    linked list of chunks together with a lock-protected "producer" flag,
    exactly as in the fmpz/nmod code.

    As in gr_mpoly_divides_heap and gr_mpoly_mul_heap_threaded, the two
    small/large coefficient code paths used by the fmpz version (packing the
    accumulator into two or three words) correspond here to the two
    accumulation strategies of the generic GR heap arithmetic: a fast path
    that gathers the operands of a diagonal into two shallow vectors and
    accumulates with a single _gr_vec_dot call (used when the coefficient
    ring overloads VEC_DOT), and a generic path that accumulates term by term
    using a single preallocated product temporary.

    Exact division of coefficients follows gr_poly_divexact/gr_mpoly_divides_heap:
    if the leading coefficient of B is a unit, we invert it once (in the
    calling thread, before any workers are started) and multiply; otherwise
    we fall back to gr_div, which over an integral domain such as Z returns
    GR_DOMAIN exactly when the division has a nonzero remainder.

    Status handling: ring operations can return GR_UNABLE (undecidable, e.g.
    over an inexact ring) in addition to GR_SUCCESS/GR_DOMAIN.  A definite
    GR_DOMAIN (the division is provably not exact) takes priority and causes
    all workers to stop as soon as they notice; GR_UNABLE is recorded but by
    itself also causes the computation to stop early, since none of the
    partial results can be trusted at that point (the same as with a
    definite failure).  The two conditions are tracked separately in the
    base structure so that the correct final status can be reported.
*/

/* shallow copy helpers (SET_SHALLOW*), MULSUB_FETCH_NODE, and the
   gr_mpoly_ts_* / divides_heap_chunk_struct / divides_heap_base_struct /
   _gr_mpoly_stripe_struct / worker_arg_struct types now live in
   threaded_divmod.h, shared with divrem_heap_threaded.c. */

/* gather product node (Bcoeff[i], Qcoeff[j]) into dot_a/dot_b; record
   dividend node (used by the divides stripe functions) */
#define STRIPE_DIVIDES_FETCH_NODE(SET_SHALLOW) \
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

/*
    Divide the accumulated value by the leading coefficient of B, storing the
    quotient coefficient in GR_ENTRY(Qcoeff, Qlen, sz).  Mirrors the
    DIVIDES_COEFF macro of gr_mpoly/divides_heap.c.
*/
#define STRIPE_DIVIDES_COEFF(acc) \
    (lc_is_one  ? gr_set(GR_ENTRY(Qcoeff, Qlen, sz), acc, cctx) \
     : lc_is_unit ? gr_mul(GR_ENTRY(Qcoeff, Qlen, sz), acc, lc_inv, cctx) \
                  : gr_div(GR_ENTRY(Qcoeff, Qlen, sz), acc, Bcoeff, cctx))

/*
    gr_mpoly_ts_struct/gr_mpoly_ts_t is declared in threaded_divmod.h; the
    functions below (used by both this file and divrem_heap_threaded.c) are
    defined here with external linkage.
*/

/* Bcoeff is emptied (its elements are moved into A) */
void gr_mpoly_ts_init(gr_mpoly_ts_t A,
        gr_ptr Bcoeff, ulong * Bexp, slong Blen,
        flint_bitcnt_t bits, slong N, gr_ctx_t cctx)
{
    slong i;
    slong sz = cctx->sizeof_elem;
    flint_bitcnt_t idx = FLINT_BIT_COUNT(Blen);
    idx = (idx <= 8) ? 0 : idx - 8;
    for (i = 0; i < FLINT_BITS; i++)
    {
        A->exp_array[i] = NULL;
        A->coeff_array[i] = NULL;
    }
    A->bits = bits;
    A->N = N;
    A->idx = idx;
    A->alloc = WORD(256) << idx;
    A->exps = A->exp_array[idx]
            = (ulong *) flint_malloc(N*A->alloc*sizeof(ulong));
    A->coeffs = A->coeff_array[idx]
              = (gr_ptr) flint_malloc(A->alloc*sz);
    _gr_vec_init(A->coeffs, A->alloc, cctx);
    A->length = Blen;
    for (i = 0; i < Blen; i++)
    {
        gr_swap(GR_ENTRY(A->coeffs, i, sz), GR_ENTRY(Bcoeff, i, sz), cctx);
        mpoly_monomial_set(A->exps + N*i, Bexp + N*i, N);
    }
}

void gr_mpoly_ts_clear(gr_mpoly_ts_t A, gr_ctx_t cctx)
{
    slong i;

    for (i = 0; i < FLINT_BITS; i++)
    {
        if (A->exp_array[i] != NULL)
        {
            slong level_alloc = UWORD(256) << i;
            FLINT_ASSERT(A->coeff_array[i] != NULL);
            _gr_vec_clear(A->coeff_array[i], level_alloc, cctx);
            flint_free(A->coeff_array[i]);
            flint_free(A->exp_array[i]);
        }
    }
}

void gr_mpoly_ts_clear_poly(gr_mpoly_t Q, gr_mpoly_ts_t A, gr_ctx_t cctx)
{
    if (Q->coeffs_alloc != 0)
    {
        _gr_vec_clear(Q->coeffs, Q->coeffs_alloc, cctx);
        flint_free(Q->coeffs);
    }
    if (Q->exps_alloc != 0)
        flint_free(Q->exps);

    Q->exps = A->exps;
    Q->coeffs = A->coeffs;
    Q->bits = A->bits;
    Q->coeffs_alloc = A->alloc;
    Q->exps_alloc = A->N * A->alloc;
    Q->length = A->length;

    A->coeff_array[A->idx] = NULL;
    A->exp_array[A->idx] = NULL;
    gr_mpoly_ts_clear(A, cctx);
}

/* put B on the end of A -- Bcoeff is emptied */
void gr_mpoly_ts_append(gr_mpoly_ts_t A,
                gr_ptr Bcoeff, ulong * Bexps, slong Blen, slong N, gr_ctx_t cctx)
{
/* TODO: this needs barriers on non-x86 */

    slong i;
    slong sz = cctx->sizeof_elem;
    ulong * oldexps = A->exps;
    gr_ptr oldcoeffs = A->coeffs;
    slong oldlength = A->length;
    slong newlength = A->length + Blen;

    if (newlength <= A->alloc)
    {
        /* write new terms first */
        for (i = 0; i < Blen; i++)
        {
            gr_swap(GR_ENTRY(oldcoeffs, oldlength + i, sz),
                    GR_ENTRY(Bcoeff, i, sz), cctx);
            mpoly_monomial_set(oldexps + N*(oldlength + i), Bexps + N*i, N);
        }
    }
    else
    {
        slong newalloc;
        ulong * newexps;
        gr_ptr newcoeffs;
        flint_bitcnt_t newidx;
        newidx = FLINT_BIT_COUNT(newlength - 1);
        newidx = (newidx > 8) ? newidx - 8 : 0;
        FLINT_ASSERT(newidx > A->idx);

        newalloc = UWORD(256) << newidx;
        FLINT_ASSERT(newlength <= newalloc);
        newexps = A->exp_array[newidx]
                = (ulong *) flint_malloc(N*newalloc*sizeof(ulong));
        newcoeffs = A->coeff_array[newidx]
                  = (gr_ptr) flint_malloc(newalloc*sz);
        _gr_vec_init(newcoeffs, newalloc, cctx);

        for (i = 0; i < oldlength; i++)
        {
            /* NOTE: fmpz_mpoly_ts_append uses a shallow set here.
               This should be possible in the GR setting too but would
               require some care regarding how we init/clear
               the coefficients. */
            GR_MUST_SUCCEED(gr_set(GR_ENTRY(newcoeffs, i, sz), GR_ENTRY(oldcoeffs, i, sz), cctx));
            mpoly_monomial_set(newexps + N*i, oldexps + N*i, N);
        }
        for (i = 0; i < Blen; i++)
        {
            gr_swap(GR_ENTRY(newcoeffs, oldlength + i, sz),
                    GR_ENTRY(Bcoeff, i, sz), cctx);
            mpoly_monomial_set(newexps + N*(oldlength + i), Bexps + N*i, N);
        }

        A->alloc = newalloc;
        A->exps = newexps;
        A->coeffs = newcoeffs;
        A->idx = newidx;

        /* do not free oldcoeffs/oldexps as other threads may be using them */
    }

    /* update length at the very end */
    A->length = newlength;
}


/*
    divides_heap_chunk_struct / divides_heap_base_struct / worker_arg_struct
    and _gr_mpoly_stripe_struct are declared in threaded_divmod.h.
*/

void divides_heap_base_init(divides_heap_base_t H)
{
    H->head = NULL;
    H->tail = NULL;
    H->cur = NULL;
    H->ctx = NULL;
    H->length = 0;
    H->N = 0;
    H->bits = 0;
    H->cmpmask = NULL;
}

void divides_heap_chunk_clear(divides_heap_chunk_t L, divides_heap_base_t H)
{
    if (L->Cinited)
        gr_mpoly_clear(L->polyC, H->ctx);
}

int divides_heap_base_clear(gr_mpoly_t Q, gr_mpoly_t R, divides_heap_base_t H)
{
    divides_heap_chunk_struct * L = H->head;
    int status;

    FLINT_ASSERT((R != NULL) == H->want_remainder);

    while (L != NULL)
    {
        divides_heap_chunk_struct * nextL = L->next;
        divides_heap_chunk_clear(L, H);
        flint_free(L);
        L = nextL;
    }
    H->head = NULL;
    H->tail = NULL;
    H->cur = NULL;
    H->length = 0;
    H->N = 0;
    H->bits = 0;
    H->cmpmask = NULL;

    if (H->have_domain)
        status = GR_DOMAIN;
    else if (H->have_unable || H->overflowed)
        status = GR_UNABLE;
    else
        status = GR_SUCCESS;

    if (status != GR_SUCCESS)
    {
        GR_IGNORE(gr_mpoly_zero(Q, H->ctx));
        gr_mpoly_ts_clear(H->polyQ, H->cctx);
        if (H->want_remainder)
        {
            GR_IGNORE(gr_mpoly_zero(R, H->ctx));
            gr_mpoly_ts_clear(H->polyR, H->cctx);
        }
    }
    else
    {
        gr_mpoly_ts_clear_poly(Q, H->polyQ, H->cctx);
        if (H->want_remainder)
            gr_mpoly_ts_clear_poly(R, H->polyR, H->cctx);
    }

    return status;
}

void divides_heap_base_add_chunk(divides_heap_base_t H, divides_heap_chunk_t L)
{
    L->next = NULL;

    if (H->tail == NULL)
    {
        FLINT_ASSERT(H->head == NULL);
        H->tail = L;
        H->head = L;
    }
    else
    {
        divides_heap_chunk_struct * tail = H->tail;
        FLINT_ASSERT(tail->next == NULL);
        tail->next = L;
        H->tail = L;
    }
    H->length++;
}


/*
    A = D - (a stripe of B * C), single-word exponent version.

    B is the small "increment" -- the newly available slice of the growing
    quotient -- and is swept afresh on every call (S->big_mem is sized
    relative to Blen).  C is the actual, fixed divisor; S->startidx/S->endidx
    are persistent (monotonically shrinking) indices into C, valid provided
    B keeps decreasing (i.e. this chunk's exponent window keeps narrowing)
    while C stays the same across calls.  This (initially confusing) choice
    of roles mirrors fmpz_mpoly_divides_heap_threaded / nmod_mpoly_divides_heap_threaded,
    where it is an important optimisation: the per-call heap setup is O(Blen)
    (the small increment) rather than O(Clen) (the possibly much larger B).

    saveD: 0 means the coefficients of D may be moved out of D (D is a
           scratch remainder), 1 means D must be left undisturbed (D may be
           a slice of the original dividend, shared with other chunks).

    Returns the length of A and writes GR_SUCCESS/GR_UNABLE to *res_status
    (this routine, unlike the plain division, cannot detect that a division
    is impossible -- it is simply forming a difference of products).
*/
slong _gr_mpoly_mulsub_stripe1(
    gr_ptr * A_coeff, ulong ** A_exp, slong * A_alloc, slong * A_exps_alloc,
    gr_srcptr Dcoeff, const ulong * Dexp, slong Dlen, int saveD,
    gr_srcptr Bcoeff, const ulong * Bexp, slong Blen,
    gr_srcptr Ccoeff, const ulong * Cexp, slong Clen,
    const _gr_mpoly_stripe_t S, int * res_status, int * overflowed)
{
    gr_mpoly_ctx_struct * ctx = S->ctx;
    gr_ctx_struct * cctx = S->cctx;
    slong sz = S->sz;
    int have_fast_dot = S->have_fast_dot;
    int upperclosed;
    slong startidx, endidx;
    ulong prev_startidx;
    ulong maskhi = S->cmpmask[0];
    ulong mask = mpoly_overflow_mask_sp(S->bits);
    ulong emax, emin;
    slong i, j;
    slong next_loc = Blen + 4;
    slong heap_len = 1;
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base, store_len;
    mpoly_heap_t * x;
    slong Di;
    slong Alen;
    slong Aalloc = *A_alloc;
    slong Aexps_alloc = *A_exps_alloc;
    gr_ptr Acoeff = *A_coeff;
    ulong * Aexp = *A_exp;
    ulong exp;
    slong * ends;
    ulong texp;
    slong * hind;
    gr_ptr pp, dot_a, dot_b;
    slong dot_len;
    int status = GR_SUCCESS;

    *overflowed = 0;

    FLINT_ASSERT(S->N == 1);

    i = 0;
    hind = (slong *)(S->big_mem + i);
    i += Blen*sizeof(slong);
    ends = (slong *)(S->big_mem + i);
    i += Blen*sizeof(slong);
    store = store_base = (slong *) (S->big_mem + i);
    i += 2*Blen*sizeof(slong);
    heap = (mpoly_heap1_s *)(S->big_mem + i);
    i += (Blen + 1)*sizeof(mpoly_heap1_s);
    chain = (mpoly_heap_t *)(S->big_mem + i);
    i += Blen*sizeof(mpoly_heap_t);
    FLINT_ASSERT(i <= S->big_mem_alloc);

    GR_TMP_INIT(pp, cctx);
    dot_a = flint_malloc(2 * Blen * sz);
    dot_b = GR_ENTRY(dot_a, Blen, sz);

    startidx = *S->startidx;
    endidx = *S->endidx;
    upperclosed = S->upperclosed;
    emax = S->emax[0];
    emin = S->emin[0];

    /* put all the starting nodes on the heap */
    prev_startidx = -UWORD(1);
    for (i = 0; i < Blen; i++)
    {
        if (startidx < Clen)
        {
            texp = Bexp[i] + Cexp[startidx];
            FLINT_ASSERT(mpoly_monomial_cmp1(emax, texp, maskhi)
                                                               > -upperclosed);
        }
        while (startidx > 0)
        {
            texp = Bexp[i] + Cexp[startidx - 1];
            if (mpoly_monomial_cmp1(emax, texp, maskhi) <= -upperclosed)
                break;
            startidx--;
        }

        if (endidx < Clen)
        {
            texp = Bexp[i] + Cexp[endidx];
            FLINT_ASSERT(mpoly_monomial_cmp1(emin, texp, maskhi) > 0);
        }
        while (endidx > 0)
        {
            texp = Bexp[i] + Cexp[endidx - 1];
            if (mpoly_monomial_cmp1(emin, texp, maskhi) <= 0)
                break;
            endidx--;
        }

        ends[i] = endidx;
        hind[i] = 2*startidx + 1;

        if ((startidx < endidx) && (((ulong) startidx) < prev_startidx))
        {
            x = chain + i;
            x->i = i;
            x->j = startidx;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;
            _mpoly_heap_insert1(heap, Bexp[x->i] + Cexp[x->j], x,
                                      &next_loc, &heap_len, maskhi);
        }

        prev_startidx = startidx;
    }

    *S->startidx = startidx;
    *S->endidx = endidx;

    Alen = 0;
    Di = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        if (mpoly_monomial_overflows1(exp, mask))
        {
            *overflowed = 1;
            goto cleanup;
        }

        while (Di < Dlen && mpoly_monomial_gt1(Dexp[Di], exp, maskhi))
        {
            _gr_mpoly_fit_length(&Acoeff, &Aalloc, &Aexp, &Aexps_alloc, 1, Alen + 1, ctx);
            Aexp[Alen] = Dexp[Di];
            if (saveD)
                status |= gr_set(GR_ENTRY(Acoeff, Alen, sz), GR_ENTRY(Dcoeff, Di, sz), cctx);
            else
                gr_swap(GR_ENTRY(Acoeff, Alen, sz), GR_ENTRY(Dcoeff, Di, sz), cctx);
            Alen++;
            Di++;
        }

        _gr_mpoly_fit_length(&Acoeff, &Aalloc, &Aexp, &Aexps_alloc, 1, Alen + 1, ctx);

        Aexp[Alen] = exp;

        store_len = 0;

        if (have_fast_dot)
        {
            int have_D;

            dot_len = 0;
            do
            {
                x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
                do
                {
                    if (sz == 1)       { MULSUB_FETCH_NODE(SET_SHALLOW1); }
                    else if (sz == 2)  { MULSUB_FETCH_NODE(SET_SHALLOW2); }
                    else if (sz == 4)  { MULSUB_FETCH_NODE(SET_SHALLOW4); }
                    else if (sz == 8)  { MULSUB_FETCH_NODE(SET_SHALLOW8); }
                    else if (sz == 16) { MULSUB_FETCH_NODE(SET_SHALLOW16); }
                    else               { MULSUB_FETCH_NODE(SET_SHALLOW_GENERIC); }
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);

            have_D = (Di < Dlen && Dexp[Di] == exp);
            status |= _gr_vec_dot(GR_ENTRY(Acoeff, Alen, sz),
                        have_D ? GR_ENTRY(Dcoeff, Di, sz) : NULL, 1,
                        dot_a, dot_b, dot_len, cctx);
            if (have_D)
                Di++;
        }
        else
        {
            if (Di < Dlen && Dexp[Di] == exp)
            {
                status |= gr_set(GR_ENTRY(Acoeff, Alen, sz), GR_ENTRY(Dcoeff, Di, sz), cctx);
                Di++;
            }
            else
            {
                status |= gr_zero(GR_ENTRY(Acoeff, Alen, sz), cctx);
            }

            do
            {
                x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
                do
                {
                    hind[x->i] |= WORD(1);
                    store[store_len++] = x->i;
                    store[store_len++] = x->j;
                    status |= gr_mul(pp, GR_ENTRY(Bcoeff, x->i, sz), GR_ENTRY(Ccoeff, x->j, sz), cctx);
                    status |= gr_sub(GR_ENTRY(Acoeff, Alen, sz), GR_ENTRY(Acoeff, Alen, sz), pp, cctx);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);
        }

        /* process nodes taken from the heap */
        while (store_len > 0)
        {
            j = store_base[--store_len];
            i = store_base[--store_len];

            /* should we go right? */
            if ((i + 1 < Blen) && (j + 0 < ends[i + 1]) && (hind[i + 1] == 2*j + 1))
            {
                x = chain + i + 1;
                x->i = i + 1;
                x->j = j;
                x->next = NULL;
                hind[x->i] = 2*(x->j + 1) + 0;
                _mpoly_heap_insert1(heap, Bexp[x->i] + Cexp[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
            }

            /* should we go up? */
            if ((j + 1 < ends[i + 0]) && ((hind[i] & 1) == 1) &&
                ((i == 0) || (hind[i - 1] >= 2*(j + 2) + 1)))
            {
                x = chain + i;
                x->i = i;
                x->j = j + 1;
                x->next = NULL;
                hind[x->i] = 2*(x->j + 1) + 0;
                _mpoly_heap_insert1(heap, Bexp[x->i] + Cexp[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
            }
        }

        /* Note: the zero status of this coefficient does not affect divisibility
           testing at this point. */
        if (gr_is_zero(GR_ENTRY(Acoeff, Alen, sz), cctx) != T_TRUE)
            Alen++;
    }

    FLINT_ASSERT(Di <= Dlen);
    _gr_mpoly_fit_length(&Acoeff, &Aalloc, &Aexp, &Aexps_alloc, 1, Alen + Dlen - Di, ctx);
    if (saveD)
        status |= _gr_vec_set(GR_ENTRY(Acoeff, Alen, sz), GR_ENTRY(Dcoeff, Di, sz), Dlen - Di, cctx);
    else
        _gr_vec_swap(GR_ENTRY(Acoeff, Alen, sz), (gr_ptr) GR_ENTRY(Dcoeff, Di, sz), Dlen - Di, cctx);
    mpoly_copy_monomials(Aexp + Alen, Dexp + Di, Dlen - Di, 1);
    Alen += Dlen - Di;

cleanup:

    GR_TMP_CLEAR(pp, cctx);
    flint_free(dot_a);

    *A_coeff = Acoeff;
    *A_exp = Aexp;
    *A_alloc = Aalloc;
    *A_exps_alloc = Aexps_alloc;
    *res_status = status;
    return Alen;
}


/* Multi-word exponent version of the above. */
slong _gr_mpoly_mulsub_stripe(
    gr_ptr * A_coeff, ulong ** A_exp, slong * A_alloc, slong * A_exps_alloc,
    gr_srcptr Dcoeff, const ulong * Dexp, slong Dlen, int saveD,
    gr_srcptr Bcoeff, const ulong * Bexp, slong Blen,
    gr_srcptr Ccoeff, const ulong * Cexp, slong Clen,
    const _gr_mpoly_stripe_t S, int * res_status, int * overflowed)
{
    gr_mpoly_ctx_struct * ctx = S->ctx;
    gr_ctx_struct * cctx = S->cctx;
    slong sz = S->sz;
    int have_fast_dot = S->have_fast_dot;
    int upperclosed;
    slong startidx, endidx;
    ulong prev_startidx;
    ulong * emax = S->emax;
    ulong * emin = S->emin;
    slong N = S->N;
    flint_bitcnt_t bits = S->bits;
    ulong mask = bits <= FLINT_BITS ? mpoly_overflow_mask_sp(bits) : 0;
    slong i, j;
    slong next_loc = Blen + 4;
    slong heap_len = 1;
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base, store_len;
    mpoly_heap_t * x;
    slong Di;
    slong Alen;
    slong Aalloc = *A_alloc;
    slong Aexps_alloc = *A_exps_alloc;
    gr_ptr Acoeff = *A_coeff;
    ulong * Aexp = *A_exp;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    slong * ends;
    ulong * texp;
    slong * hind;
    gr_ptr pp, dot_a, dot_b;
    slong dot_len;
    int status = GR_SUCCESS;

    *overflowed = 0;

    i = 0;
    hind = (slong *)(S->big_mem + i);
    i += Blen*sizeof(slong);
    ends = (slong *)(S->big_mem + i);
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
    texp = (ulong *)(S->big_mem + i);
    i +=  N*sizeof(ulong);
    FLINT_ASSERT(i <= S->big_mem_alloc);

    GR_TMP_INIT(pp, cctx);
    dot_a = flint_malloc(2 * Blen * sz);
    dot_b = GR_ENTRY(dot_a, Blen, sz);

    exp_next = 0;

    startidx = *S->startidx;
    endidx = *S->endidx;
    upperclosed = S->upperclosed;

    for (i = 0; i < Blen; i++)
        exp_list[i] = exps + i*N;

    /* put all the starting nodes on the heap */
    prev_startidx = -UWORD(1);
    for (i = 0; i < Blen; i++)
    {
        if (startidx < Clen)
        {
            mpoly_monomial_add_mp(texp, Bexp + N*i, Cexp + N*startidx, N);
            FLINT_ASSERT(mpoly_monomial_cmp(emax, texp, N, S->cmpmask)
                                                               > -upperclosed);
        }
        while (startidx > 0)
        {
            mpoly_monomial_add_mp(texp, Bexp + N*i, Cexp + N*(startidx - 1), N);
            if (mpoly_monomial_cmp(emax, texp, N, S->cmpmask) <= -upperclosed)
                break;
            startidx--;
        }

        if (endidx < Clen)
        {
            mpoly_monomial_add_mp(texp, Bexp + N*i, Cexp + N*endidx, N);
            FLINT_ASSERT(mpoly_monomial_cmp(emin, texp, N, S->cmpmask) > 0);
        }
        while (endidx > 0)
        {
            mpoly_monomial_add_mp(texp, Bexp + N*i, Cexp + N*(endidx - 1), N);
            if (mpoly_monomial_cmp(emin, texp, N, S->cmpmask) <= 0)
                break;
            endidx--;
        }

        ends[i] = endidx;
        hind[i] = 2*startidx + 1;

        if ((startidx < endidx) && (((ulong) startidx) < prev_startidx))
        {
            x = chain + i;
            x->i = i;
            x->j = startidx;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;

            mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i,
                                                      Cexp + N*x->j, N);

            exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                          &next_loc, &heap_len, N, S->cmpmask);
        }

        prev_startidx = startidx;
    }

    *S->startidx = startidx;
    *S->endidx = endidx;

    Alen = 0;
    Di = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        if (bits <= FLINT_BITS
                ? mpoly_monomial_overflows(exp, N, mask)
                : mpoly_monomial_overflows_mp(exp, N, bits))
        {
            *overflowed = 1;
            goto cleanup;
        }

        while (Di < Dlen && mpoly_monomial_gt(Dexp + N*Di, exp, N, S->cmpmask))
        {
            _gr_mpoly_fit_length(&Acoeff, &Aalloc, &Aexp, &Aexps_alloc, N, Alen + 1, ctx);
            mpoly_monomial_set(Aexp + N*Alen, Dexp + N*Di, N);
            if (saveD)
                status |= gr_set(GR_ENTRY(Acoeff, Alen, sz), GR_ENTRY(Dcoeff, Di, sz), cctx);
            else
                gr_swap(GR_ENTRY(Acoeff, Alen, sz), GR_ENTRY(Dcoeff, Di, sz), cctx);
            Alen++;
            Di++;
        }

        _gr_mpoly_fit_length(&Acoeff, &Aalloc, &Aexp, &Aexps_alloc, N, Alen + 1, ctx);

        mpoly_monomial_set(Aexp + N*Alen, exp, N);

        store_len = 0;

        if (have_fast_dot)
        {
            int have_D;

            dot_len = 0;
            do
            {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, S->cmpmask);
                do
                {
                    if (sz == 1)       { MULSUB_FETCH_NODE(SET_SHALLOW1); }
                    else if (sz == 2)  { MULSUB_FETCH_NODE(SET_SHALLOW2); }
                    else if (sz == 4)  { MULSUB_FETCH_NODE(SET_SHALLOW4); }
                    else if (sz == 8)  { MULSUB_FETCH_NODE(SET_SHALLOW8); }
                    else if (sz == 16) { MULSUB_FETCH_NODE(SET_SHALLOW16); }
                    else               { MULSUB_FETCH_NODE(SET_SHALLOW_GENERIC); }
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

            have_D = (Di < Dlen && mpoly_monomial_equal(Dexp + N*Di, exp, N));
            status |= _gr_vec_dot(GR_ENTRY(Acoeff, Alen, sz),
                        have_D ? GR_ENTRY(Dcoeff, Di, sz) : NULL, 1,
                        dot_a, dot_b, dot_len, cctx);
            if (have_D)
                Di++;
        }
        else
        {
            if (Di < Dlen && mpoly_monomial_equal(Dexp + N*Di, exp, N))
            {
                status |= gr_set(GR_ENTRY(Acoeff, Alen, sz), GR_ENTRY(Dcoeff, Di, sz), cctx);
                Di++;
            }
            else
            {
                status |= gr_zero(GR_ENTRY(Acoeff, Alen, sz), cctx);
            }

            do
            {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, S->cmpmask);
                do
                {
                    hind[x->i] |= WORD(1);
                    store[store_len++] = x->i;
                    store[store_len++] = x->j;
                    status |= gr_mul(pp, GR_ENTRY(Bcoeff, x->i, sz), GR_ENTRY(Ccoeff, x->j, sz), cctx);
                    status |= gr_sub(GR_ENTRY(Acoeff, Alen, sz), GR_ENTRY(Acoeff, Alen, sz), pp, cctx);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));
        }

        /* process nodes taken from the heap */
        while (store_len > 0)
        {
            j = store_base[--store_len];
            i = store_base[--store_len];

            /* should we go right? */
            if ((i + 1 < Blen) && (j + 0 < ends[i + 1]) && (hind[i + 1] == 2*j + 1))
            {
                x = chain + i + 1;
                x->i = i + 1;
                x->j = j;
                x->next = NULL;
                hind[x->i] = 2*(x->j + 1) + 0;

                mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i,
                                                          Cexp + N*x->j, N);

                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                          &next_loc, &heap_len, N, S->cmpmask);
            }

            /* should we go up? */
            if ((j + 1 < ends[i + 0]) && ((hind[i] & 1) == 1) &&
                ((i == 0) || (hind[i - 1] >= 2*(j + 2) + 1)))
            {
                x = chain + i;
                x->i = i;
                x->j = j + 1;
                x->next = NULL;
                hind[x->i] = 2*(x->j + 1) + 0;

                mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i,
                                                          Cexp + N*x->j, N);

                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                          &next_loc, &heap_len, N, S->cmpmask);
            }
        }

        /* an unknown zero-status is kept as a (possibly redundant) term
           rather than causing an abort -- this mirrors the "acc" checks
           fixed for divrem/divexact (see the discussion at
           gr_mpoly_divexact), and is unconditionally safe here regardless
           of that discussion's open questions for gr_mpoly_divides itself:
           mulsub_stripe never makes any exactness determination, it only
           chooses a (possibly non-maximally-reduced) representation of the
           difference D - stripe(B*C). This was already keeping the term in
           the "unknown" case; the only change is to stop also tainting the
           overall status with GR_UNABLE merely for being unable to prove
           it away. */
        if (gr_is_zero(GR_ENTRY(Acoeff, Alen, sz), cctx) != T_TRUE)
            Alen++;
    }

    FLINT_ASSERT(Di <= Dlen);
    _gr_mpoly_fit_length(&Acoeff, &Aalloc, &Aexp, &Aexps_alloc, N, Alen + Dlen - Di, ctx);
    if (saveD)
        status |= _gr_vec_set(GR_ENTRY(Acoeff, Alen, sz), GR_ENTRY(Dcoeff, Di, sz), Dlen - Di, cctx);
    else
        _gr_vec_swap(GR_ENTRY(Acoeff, Alen, sz), (gr_ptr) GR_ENTRY(Dcoeff, Di, sz), Dlen - Di, cctx);
    mpoly_copy_monomials(Aexp + N*Alen, Dexp + N*Di, Dlen - Di, N);
    Alen += Dlen - Di;

cleanup:

    GR_TMP_CLEAR(pp, cctx);
    flint_free(dot_a);

    *A_coeff = Acoeff;
    *A_exp = Aexp;
    *A_alloc = Aalloc;
    *A_exps_alloc = Aexps_alloc;
    *res_status = status;
    return Alen;
}


/*
    Q = stripe of R/B (single-word exponent version), assuming R != 0.

    New quotient terms/heap nodes are only produced while their exponent
    stays within [S->emin, infinity) (the fmpz/nmod analogue restricts the
    heap to a chunk's exponent window in the same way): nodes whose exponent
    would fall below S->emin are never inserted (see the "should we go
    right/up" logic below), matching the fact that lower quotient terms are
    the responsibility of a different chunk.

    Returns Qlen and writes GR_SUCCESS / GR_DOMAIN / GR_UNABLE to
    *res_status.  Qlen is meaningful only when *res_status == GR_SUCCESS.
*/
static slong _gr_mpoly_divides_stripe1(
    gr_ptr * Q_coeff, ulong ** Q_exp, slong * Q_alloc, slong * Q_exps_alloc,
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
    int lt_divides;
    slong i, j, s;
    slong next_loc, heap_len;
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base, store_len;
    mpoly_heap_t * x;
    slong Qlen;
    slong Qalloc = *Q_alloc;
    slong Qexps_alloc = *Q_exps_alloc;
    gr_ptr Qcoeff = *Q_coeff;
    ulong * Qexp = *Q_exp;
    ulong exp;
    ulong mask;
    slong * hind;
    gr_ptr acc, pp, dot_a, dot_b;
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

    GR_TMP_INIT2(acc, pp, cctx);
    dot_a = flint_malloc(2 * Blen * sz);
    dot_b = GR_ENTRY(dot_a, Blen, sz);

    for (i = 0; i < Blen; i++)
        hind[i] = 1;

    mask = mpoly_overflow_mask_sp(bits);

    Qlen = 0;
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

        if (mpoly_monomial_overflows1(exp, mask))
            goto not_exact_division;

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
                    if (sz == 1)       { STRIPE_DIVIDES_FETCH_NODE(SET_SHALLOW1); }
                    else if (sz == 2)  { STRIPE_DIVIDES_FETCH_NODE(SET_SHALLOW2); }
                    else if (sz == 4)  { STRIPE_DIVIDES_FETCH_NODE(SET_SHALLOW4); }
                    else if (sz == 8)  { STRIPE_DIVIDES_FETCH_NODE(SET_SHALLOW8); }
                    else if (sz == 16) { STRIPE_DIVIDES_FETCH_NODE(SET_SHALLOW16); }
                    else               { STRIPE_DIVIDES_FETCH_NODE(SET_SHALLOW_GENERIC); }
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

        switch (gr_is_zero(acc, cctx))
        {
            case T_TRUE:
                continue;
            case T_FALSE:
                break;
            default:
                status |= GR_UNABLE;
                goto unable;
        }

        cstatus = STRIPE_DIVIDES_COEFF(acc);
        if (cstatus == GR_DOMAIN)
            goto not_exact_division;
        if (cstatus != GR_SUCCESS)
        {
            status |= cstatus;
            goto unable;
        }

        if (!lt_divides)
            goto not_exact_division;

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

    *res_status = status;
    goto cleanup;

not_exact_division:
    Qlen = 0;
    /* preserve any GR_UNABLE accumulated while computing this term's acc:
       if an earlier operation could not decide, this term's apparent
       non-exactness may just be an artifact of that, not a genuine proof --
       see the analogous fix in gr_mpoly_divides_heap. */
    status |= GR_DOMAIN;
    *res_status = status;
    goto cleanup;

unable:
    Qlen = 0;
    *res_status = status;

cleanup:

    GR_TMP_CLEAR2(acc, pp, cctx);
    flint_free(dot_a);

    *Q_coeff = Qcoeff;
    *Q_exp = Qexp;
    *Q_alloc = Qalloc;
    *Q_exps_alloc = Qexps_alloc;

    return Qlen;
}

/* Multi-word exponent version of the above. */
static slong _gr_mpoly_divides_stripe(
    gr_ptr * Q_coeff, ulong ** Q_exp, slong * Q_alloc, slong * Q_exps_alloc,
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
    int lt_divides;
    slong i, j, s;
    slong next_loc, heap_len;
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base, store_len;
    mpoly_heap_t * x;
    slong Qlen;
    slong Qalloc = *Q_alloc;
    slong Qexps_alloc = *Q_exps_alloc;
    gr_ptr Qcoeff = *Q_coeff;
    ulong * Qexp = *Q_exp;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    ulong mask;
    slong * hind;
    gr_ptr acc, pp, dot_a, dot_b;
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

    GR_TMP_INIT2(acc, pp, cctx);
    dot_a = flint_malloc(2 * Blen * sz);
    dot_b = GR_ENTRY(dot_a, Blen, sz);

    exp_next = 0;
    for (i = 0; i < Blen; i++)
        exp_list[i] = exps + i*N;

    for (i = 0; i < Blen; i++)
        hind[i] = 1;

    mask = bits <= FLINT_BITS ? mpoly_overflow_mask_sp(bits) : 0;

    Qlen = 0;
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
                goto not_exact_division;
            lt_divides = mpoly_monomial_divides(Qexp + N*Qlen, exp, Bexp + N*0, N, mask);
        }
        else
        {
            if (mpoly_monomial_overflows_mp(exp, N, bits))
                goto not_exact_division;
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
                    if (sz == 1)       { STRIPE_DIVIDES_FETCH_NODE(SET_SHALLOW1); }
                    else if (sz == 2)  { STRIPE_DIVIDES_FETCH_NODE(SET_SHALLOW2); }
                    else if (sz == 4)  { STRIPE_DIVIDES_FETCH_NODE(SET_SHALLOW4); }
                    else if (sz == 8)  { STRIPE_DIVIDES_FETCH_NODE(SET_SHALLOW8); }
                    else if (sz == 16) { STRIPE_DIVIDES_FETCH_NODE(SET_SHALLOW16); }
                    else               { STRIPE_DIVIDES_FETCH_NODE(SET_SHALLOW_GENERIC); }
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

        switch (gr_is_zero(acc, cctx))
        {
            case T_TRUE:
                continue;
            case T_FALSE:
                break;
            default:
                status |= GR_UNABLE;
                goto unable;
        }

        cstatus = STRIPE_DIVIDES_COEFF(acc);
        if (cstatus == GR_DOMAIN)
            goto not_exact_division;
        if (cstatus != GR_SUCCESS)
        {
            status |= cstatus;
            goto unable;
        }

        if (!lt_divides)
            goto not_exact_division;

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

    *res_status = status;
    goto cleanup;

not_exact_division:
    Qlen = 0;
    /* preserve any GR_UNABLE accumulated while computing this term's acc:
       if an earlier operation could not decide, this term's apparent
       non-exactness may just be an artifact of that, not a genuine proof --
       see the analogous fix in gr_mpoly_divides_heap. */
    status |= GR_DOMAIN;
    *res_status = status;
    goto cleanup;

unable:
    Qlen = 0;
    *res_status = status;

cleanup:

    GR_TMP_CLEAR2(acc, pp, cctx);
    flint_free(dot_a);

    *Q_coeff = Qcoeff;
    *Q_exp = Qexp;
    *Q_alloc = Qalloc;
    *Q_exps_alloc = Qexps_alloc;

    return Qlen;
}


slong chunk_find_exp(ulong * exp, slong a,
                     const ulong * Aexp, slong Alen, slong N, const ulong * cmpmask)
{
    slong b = Alen;

try_again:
    FLINT_ASSERT(b >= a);

    FLINT_ASSERT(a > 0);
    FLINT_ASSERT(mpoly_monomial_cmp(Aexp + N*(a - 1), exp, N, cmpmask) >= 0);
    FLINT_ASSERT(b >= Alen
                  ||  mpoly_monomial_cmp(Aexp + N*b, exp, N, cmpmask) < 0);

    if (b - a < 5)
    {
        slong i = a;
        while (i < b
                && mpoly_monomial_cmp(Aexp + N*i, exp, N, cmpmask) >= 0)
        {
            i++;
        }
        return i;
    }
    else
    {
        slong c = a + (b - a)/2;
        if (mpoly_monomial_cmp(Aexp + N*c, exp, N, cmpmask) < 0)
            b = c;
        else
            a = c;
        goto try_again;
    }
}

void stripe_fit_length(_gr_mpoly_stripe_struct * S, slong new_len)
{
    slong N = S->N;
    slong new_alloc;
    new_alloc = 0;
    if (N == 1)
    {
        new_alloc += new_len*sizeof(slong);
        new_alloc += new_len*sizeof(slong);
        new_alloc += 2*new_len*sizeof(slong);
        new_alloc += (new_len + 1)*sizeof(mpoly_heap1_s);
        new_alloc += new_len*sizeof(mpoly_heap_t);
    }
    else
    {
        new_alloc += new_len*sizeof(slong);
        new_alloc += new_len*sizeof(slong);
        new_alloc += 2*new_len*sizeof(slong);
        new_alloc += (new_len + 1)*sizeof(mpoly_heap_s);
        new_alloc += new_len*sizeof(mpoly_heap_t);
        new_alloc += new_len*N*sizeof(ulong);
        new_alloc += new_len*sizeof(ulong *);
        new_alloc += N*sizeof(ulong);
    }

    if (S->big_mem_alloc >= new_alloc)
        return;

    new_alloc = FLINT_MAX(new_alloc, S->big_mem_alloc + S->big_mem_alloc/4);
    S->big_mem_alloc = new_alloc;

    if (S->big_mem != NULL)
        S->big_mem = (char *) flint_realloc(S->big_mem, new_alloc);
    else
        S->big_mem = (char *) flint_malloc(new_alloc);
}


/*
    Set L->polyC = D - (a stripe of B*Q_new), where Q_new is the newly
    available part of the quotient (indices [L->mq, q_prev_length) of
    H->polyQ), and D is the previous L->polyC (if L->Cinited) or, the first
    time this chunk is touched, the slice of the original dividend A lying
    in this chunk's exponent window.
*/
void chunk_mulsub(worker_arg_t W, divides_heap_chunk_t L, slong q_prev_length)
{
    divides_heap_base_struct * H = W->H;
    slong N = H->N;
    gr_ctx_struct * cctx = H->cctx;
    gr_mpoly_struct * C = L->polyC;
    const gr_mpoly_struct * B = H->polyB;
    const gr_mpoly_struct * A = H->polyA;
    gr_mpoly_ts_struct * Q = H->polyQ;
    gr_mpoly_struct * T1 = W->polyT1;
    _gr_mpoly_stripe_struct * S = W->S;
    int status;
    int overflowed;

    S->startidx = &L->startidx;
    S->endidx = &L->endidx;
    S->emin = L->emin;
    S->emax = L->emax;
    S->upperclosed = L->upperclosed;
    FLINT_ASSERT(S->N == N);
    stripe_fit_length(S, q_prev_length - L->mq);

    if (L->Cinited)
    {
        if (N == 1)
        {
            T1->length = _gr_mpoly_mulsub_stripe1(
                    &T1->coeffs, &T1->exps, &T1->coeffs_alloc, &T1->exps_alloc,
                    C->coeffs, C->exps, C->length, 1,
                    GR_ENTRY(Q->coeffs, L->mq, cctx->sizeof_elem), Q->exps + N*L->mq, q_prev_length - L->mq,
                    B->coeffs, B->exps, B->length,
                    S, &status, &overflowed);
        }
        else
        {
            T1->length = _gr_mpoly_mulsub_stripe(
                    &T1->coeffs, &T1->exps, &T1->coeffs_alloc, &T1->exps_alloc,
                    C->coeffs, C->exps, C->length, 1,
                    GR_ENTRY(Q->coeffs, L->mq, cctx->sizeof_elem), Q->exps + N*L->mq, q_prev_length - L->mq,
                    B->coeffs, B->exps, B->length,
                    S, &status, &overflowed);
        }
        if (!overflowed)
            gr_mpoly_swap(C, T1, H->ctx);
    }
    else
    {
        slong startidx, stopidx;
        if (L->upperclosed)
        {
            startidx = 0;
            stopidx = chunk_find_exp(L->emin, 1, H->polyA->exps, H->polyA->length, H->N, H->cmpmask);
        }
        else
        {
            startidx = chunk_find_exp(L->emax, 1, H->polyA->exps, H->polyA->length, H->N, H->cmpmask);
            stopidx = chunk_find_exp(L->emin, startidx, H->polyA->exps, H->polyA->length, H->N, H->cmpmask);
        }

        L->Cinited = 1;
        gr_mpoly_init3(C, 16 + stopidx - startidx, H->bits, H->ctx); /* any is OK */

        if (N == 1)
        {
            C->length = _gr_mpoly_mulsub_stripe1(
                    &C->coeffs, &C->exps, &C->coeffs_alloc, &C->exps_alloc,
                    GR_ENTRY(A->coeffs, startidx, cctx->sizeof_elem), A->exps + N*startidx, stopidx - startidx, 1,
                    GR_ENTRY(Q->coeffs, L->mq, cctx->sizeof_elem), Q->exps + N*L->mq, q_prev_length - L->mq,
                    B->coeffs, B->exps, B->length,
                    S, &status, &overflowed);
        }
        else
        {
            C->length = _gr_mpoly_mulsub_stripe(
                    &C->coeffs, &C->exps, &C->coeffs_alloc, &C->exps_alloc,
                    GR_ENTRY(A->coeffs, startidx, cctx->sizeof_elem), A->exps + N*startidx, stopidx - startidx, 1,
                    GR_ENTRY(Q->coeffs, L->mq, cctx->sizeof_elem), Q->exps + N*L->mq, q_prev_length - L->mq,
                    B->coeffs, B->exps, B->length,
                    S, &status, &overflowed);
        }
    }

    if (overflowed)
    {
#if FLINT_USES_PTHREAD
        pthread_mutex_lock(&H->mutex);
#endif
        H->overflowed = 1;
        H->failed = 1;
#if FLINT_USES_PTHREAD
        pthread_mutex_unlock(&H->mutex);
#endif
        return;
    }

    if (status != GR_SUCCESS)
    {
#if FLINT_USES_PTHREAD
        pthread_mutex_lock(&H->mutex);
#endif
        H->have_unable = 1;
#if FLINT_USES_PTHREAD
        pthread_mutex_unlock(&H->mutex);
#endif
    }

    L->mq = q_prev_length;
}

void trychunk(worker_arg_t W, divides_heap_chunk_t L)
{
    divides_heap_base_struct * H = W->H;
    slong N = H->N;
    gr_ctx_struct * cctx = H->cctx;
    gr_mpoly_struct * C = L->polyC;
    slong q_prev_length;
    const gr_mpoly_struct * B = H->polyB;
    const gr_mpoly_struct * A = H->polyA;
    gr_mpoly_ts_struct * Q = H->polyQ;
    gr_mpoly_struct * T2 = W->polyT2;

    /* return if this section has already finished processing */
    if (L->mq < 0)
        return;

    /* process more quotient terms if available */
    q_prev_length = Q->length;
    if (q_prev_length > L->mq)
    {
        if (L->producer == 0 && q_prev_length - L->mq < 20)
            return;

        chunk_mulsub(W, L, q_prev_length);
    }

    if (L->producer == 1)
    {
        divides_heap_chunk_struct * next;
        gr_srcptr Rcoeff;
        ulong * Rexp;
        slong Rlen;

        /* process the remaining quotient terms */
        q_prev_length = Q->length;
        if (q_prev_length > L->mq)
            chunk_mulsub(W, L, q_prev_length);

        /* find location of remaining terms */
        if (L->Cinited)
        {
            Rlen = C->length;
            Rexp = C->exps;
            Rcoeff = C->coeffs;
        }
        else
        {
            slong startidx, stopidx;
            if (L->upperclosed)
            {
                startidx = 0;
                stopidx = chunk_find_exp(L->emin, 1, H->polyA->exps, H->polyA->length, H->N, H->cmpmask);
            }
            else
            {
                startidx = chunk_find_exp(L->emax, 1, H->polyA->exps, H->polyA->length, H->N, H->cmpmask);
                stopidx = chunk_find_exp(L->emin, startidx, H->polyA->exps, H->polyA->length, H->N, H->cmpmask);
            }
            Rlen = stopidx - startidx;
            Rcoeff = GR_ENTRY(A->coeffs, startidx, cctx->sizeof_elem);
            Rexp = A->exps + N*startidx;
        }

        /* if we have remaining terms, add to quotient (and remainder) */
        if (Rlen > 0)
        {
            int rstatus;
            _gr_mpoly_stripe_struct * S = W->S;
            S->startidx = &L->startidx;
            S->endidx = &L->endidx;
            S->emin = L->emin;
            S->emax = L->emax;
            S->upperclosed = L->upperclosed;

            if (H->require_exact)
            {
                if (N == 1)
                {
                    T2->length = _gr_mpoly_divides_stripe1(
                                        &T2->coeffs, &T2->exps, &T2->coeffs_alloc, &T2->exps_alloc,
                                           Rcoeff, Rexp, Rlen,
                                           B->coeffs, B->exps, B->length, S, &rstatus);
                }
                else
                {
                    T2->length = _gr_mpoly_divides_stripe(
                                        &T2->coeffs, &T2->exps, &T2->coeffs_alloc, &T2->exps_alloc,
                                           Rcoeff, Rexp, Rlen,
                                           B->coeffs, B->exps, B->length, S, &rstatus);
                }

                if (rstatus != GR_SUCCESS)
                {
#if FLINT_USES_PTHREAD
                    pthread_mutex_lock(&H->mutex);
#endif
                    if (rstatus == GR_DOMAIN)
                        H->have_domain = 1;
                    else
                        H->have_unable = 1;
                    H->failed = 1;
#if FLINT_USES_PTHREAD
                    pthread_mutex_unlock(&H->mutex);
#endif
                    return;
                }
                else
                {
                    gr_mpoly_ts_append(H->polyQ, T2->coeffs, T2->exps, T2->length, N, cctx);
                }
            }
            else
            {
                /* divrem family: never a hard failure on non-reducible
                   terms -- they are routed to the remainder instead. */
                gr_mpoly_struct * T3 = W->polyT3;
                slong Wlen;
                int overflowed;

                if (N == 1)
                {
                    T2->length = _gr_mpoly_divrem_stripe1(
                                        &T2->coeffs, &T2->exps, &T2->coeffs_alloc, &T2->exps_alloc,
                                        &T3->coeffs, &T3->exps, &T3->coeffs_alloc, &T3->exps_alloc, &Wlen,
                                           Rcoeff, Rexp, Rlen,
                                           B->coeffs, B->exps, B->length, S, &rstatus, &overflowed);
                }
                else
                {
                    T2->length = _gr_mpoly_divrem_stripe(
                                        &T2->coeffs, &T2->exps, &T2->coeffs_alloc, &T2->exps_alloc,
                                        &T3->coeffs, &T3->exps, &T3->coeffs_alloc, &T3->exps_alloc, &Wlen,
                                           Rcoeff, Rexp, Rlen,
                                           B->coeffs, B->exps, B->length, S, &rstatus, &overflowed);
                }

                if (overflowed)
                {
#if FLINT_USES_PTHREAD
                    pthread_mutex_lock(&H->mutex);
#endif
                    H->overflowed = 1;
                    H->failed = 1;
#if FLINT_USES_PTHREAD
                    pthread_mutex_unlock(&H->mutex);
#endif
                    return;
                }

                if (rstatus != GR_SUCCESS)
                {
                    /* GR_DOMAIN: an exactly-reducible-by-monomial term's
                       coefficient did not divide exactly (nonfield == 0
                       only -- see _gr_mpoly_divrem_stripe1/stripe).
                       Otherwise GR_UNABLE. */
#if FLINT_USES_PTHREAD
                    pthread_mutex_lock(&H->mutex);
#endif
                    if (rstatus == GR_DOMAIN)
                        H->have_domain = 1;
                    else
                        H->have_unable = 1;
                    H->failed = 1;
#if FLINT_USES_PTHREAD
                    pthread_mutex_unlock(&H->mutex);
#endif
                    return;
                }

                if (T2->length > 0)
                    gr_mpoly_ts_append(H->polyQ, T2->coeffs, T2->exps, T2->length, N, cctx);
                if (H->want_remainder && Wlen > 0)
                    gr_mpoly_ts_append(H->polyR, T3->coeffs, T3->exps, Wlen, N, cctx);
            }
        }

        next = L->next;
        H->length--;
        H->cur = next;

        if (next != NULL)
            next->producer = 1;

        L->producer = 0;
        L->mq = -1;
    }

    return;
}


void worker_loop(void * varg)
{
    worker_arg_struct * W = (worker_arg_struct *) varg;
    divides_heap_base_struct * H = W->H;
    _gr_mpoly_stripe_struct * S = W->S;
    const gr_mpoly_struct * B = H->polyB;
    gr_mpoly_struct * T1 = W->polyT1;
    gr_mpoly_struct * T2 = W->polyT2;
    gr_mpoly_struct * T3 = W->polyT3;
    slong N = H->N;
    slong Blen = B->length;

    /* initialize stripe working memory */
    S->N = N;
    S->bits = H->bits;
    S->ctx = H->ctx;
    S->cctx = H->cctx;
    S->sz = H->cctx->sizeof_elem;
    S->cmpmask = H->cmpmask;
    S->have_fast_dot = (GR_VEC_DOT_OP(H->cctx, VEC_DOT) != (gr_method_vec_dot_op) gr_generic_vec_dot);
    S->lc_is_one = H->lc_is_one;
    S->lc_is_unit = H->lc_is_unit;
    S->lc_inv = H->lc_inv;
    S->nonfield = H->nonfield;
    S->unchecked = H->unchecked;
    S->big_mem_alloc = 0;
    S->big_mem = NULL;

    stripe_fit_length(S, Blen);

    gr_mpoly_init3(T1, 16, H->bits, H->ctx);
    gr_mpoly_init3(T2, 16, H->bits, H->ctx);
    gr_mpoly_init3(T3, 16, H->bits, H->ctx);

    while (!H->failed)
    {
        divides_heap_chunk_struct * L;
        L = H->cur;

        if (L == NULL)
            break;

        while (L != NULL)
        {
#if FLINT_USES_PTHREAD
            pthread_mutex_lock(&H->mutex);
#endif
            if (L->lock != -1)
            {
                L->lock = -1;
#if FLINT_USES_PTHREAD
                pthread_mutex_unlock(&H->mutex);
#endif
                trychunk(W, L);
#if FLINT_USES_PTHREAD
                pthread_mutex_lock(&H->mutex);
#endif
                L->lock = 0;
#if FLINT_USES_PTHREAD
                pthread_mutex_unlock(&H->mutex);
#endif
                break;
            }
            else
            {
#if FLINT_USES_PTHREAD
                pthread_mutex_unlock(&H->mutex);
#endif
            }

            L = L->next;
        }
    }

    gr_mpoly_clear(T1, H->ctx);
    gr_mpoly_clear(T2, H->ctx);
    gr_mpoly_clear(T3, H->ctx);
    flint_free(S->big_mem);

    return;
}


/*
    return the overall GR status (GR_SUCCESS / GR_DOMAIN / GR_UNABLE).
    The leading coefficient of B need not be a unit; gr_div is used as a
    fallback exact-division test (see gr_mpoly_divides_heap / DIVIDES_COEFF).
*/
static int _gr_mpoly_divides_heap_threaded_pool(
    gr_mpoly_t Q,
    const gr_mpoly_t A,
    const gr_mpoly_t B,
    gr_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    mpoly_ctx_struct * mctx = GR_MPOLY_MCTX(ctx);
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    slong sz = cctx->sizeof_elem;
    ulong mask;
    int status;
    fmpz_mpoly_ctx_t zctx;
    fmpz_mpoly_t S;
    slong i, k, N;
    flint_bitcnt_t exp_bits;
    ulong * cmpmask;
    ulong * Aexp, * Bexp;
    int freeAexp, freeBexp;
    worker_arg_struct * worker_args;
    gr_ptr qcoeff, lc_inv;
    ulong * texps, * qexps;
    divides_heap_base_t H;
    int lc_is_one, lc_is_unit;
    int cstatus;
    TMP_INIT;

    if (B->length < 2 || A->length < 2)
        return gr_mpoly_divides_heap(Q, A, B, ctx);

    TMP_START;

    GR_TMP_INIT2(qcoeff, lc_inv, cctx);

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

    /* precompute 1/lc(B) if it is a unit, and the first quotient coeff */
    lc_is_one = (gr_is_one(B->coeffs, cctx) == T_TRUE);
    lc_is_unit = lc_is_one || (gr_inv(lc_inv, B->coeffs, cctx) == GR_SUCCESS);

    cstatus = lc_is_one    ? gr_set(qcoeff, A->coeffs, cctx)
              : lc_is_unit ? gr_mul(qcoeff, A->coeffs, lc_inv, cctx)
                           : gr_div(qcoeff, A->coeffs, B->coeffs, cctx);
    if (cstatus == GR_DOMAIN)
    {
        status = GR_DOMAIN;
        GR_IGNORE(gr_mpoly_zero(Q, ctx));
        goto cleanup1;
    }
    if (cstatus != GR_SUCCESS)
    {
        status = GR_UNABLE;
        GR_IGNORE(gr_mpoly_zero(Q, ctx));
        goto cleanup1;
    }

    if (mpoly_divides_select_exps(S, zctx, num_handles,
                                   Aexp, A->length, Bexp, B->length, exp_bits))
    {
        status = GR_DOMAIN;
        GR_IGNORE(gr_mpoly_zero(Q, ctx));
        goto cleanup1;
    }

    /*
        At this point A and B both have at least two terms, the leading
        coefficient of A divides that of B, and the exponent selection did
        not give an easy exit.
    */
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
    H->require_exact = 1;
    H->want_remainder = 0;
    H->nonfield = 0;
    H->unchecked = 0;
    H->failed = 0;
    H->have_domain = 0;
    H->have_unable = 0;
    H->overflowed = 0;

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

    /* generate at least the first quotient term */

    texps = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    qexps = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    mpoly_monomial_sub_mp(qexps + N*0, Aexp + N*0, Bexp + N*0, N);

    gr_mpoly_ts_init(H->polyQ, qcoeff, qexps, 1, H->bits, H->N, cctx);

    mpoly_monomial_add_mp(texps, qexps + N*0, Bexp + N*1, N);

    mask = 0;
    for (i = 0; (flint_bitcnt_t) i < FLINT_BITS/exp_bits; i++)
        mask = (mask << exp_bits) + (UWORD(1) << (exp_bits - 1));

    k = 1;
    while (k < A->length && mpoly_monomial_gt(Aexp + N*k, texps, N, cmpmask))
    {
        int lt_divides;
        if (exp_bits <= FLINT_BITS)
            lt_divides = mpoly_monomial_divides(qexps, Aexp + N*k,
                                                      Bexp + N*0, N, mask);
        else
            lt_divides = mpoly_monomial_divides_mp(qexps, Aexp + N*k,
                                                      Bexp + N*0, N, exp_bits);
        if (!lt_divides)
        {
            H->failed = 1;
            H->have_domain = 1;
            break;
        }

        cstatus = lc_is_one    ? gr_set(qcoeff, GR_ENTRY(A->coeffs, k, sz), cctx)
                  : lc_is_unit ? gr_mul(qcoeff, GR_ENTRY(A->coeffs, k, sz), lc_inv, cctx)
                               : gr_div(qcoeff, GR_ENTRY(A->coeffs, k, sz), B->coeffs, cctx);
        if (cstatus == GR_DOMAIN)
        {
            H->failed = 1;
            H->have_domain = 1;
            break;
        }
        if (cstatus != GR_SUCCESS)
        {
            H->failed = 1;
            H->have_unable = 1;
            break;
        }

        gr_mpoly_ts_append(H->polyQ, qcoeff, qexps, 1, H->N, cctx);
        k++;
    }

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

    status = divides_heap_base_clear(Q, NULL, H);

cleanup1:

    fmpz_mpoly_clear(S, zctx);
    fmpz_mpoly_ctx_clear(zctx);

    if (freeAexp)
        flint_free(Aexp);

    if (freeBexp)
        flint_free(Bexp);

    GR_TMP_CLEAR2(qcoeff, lc_inv, cctx);

    TMP_END;

    return status;
}


int gr_mpoly_divides_heap_threaded(
    gr_mpoly_t Q,
    const gr_mpoly_t A,
    const gr_mpoly_t B,
    gr_mpoly_ctx_t ctx)
{
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    thread_pool_handle * handles;
    slong num_handles;
    slong thread_limit;
    int status;

    if (B->length == 0)
        return GR_DOMAIN;     /* division by zero: not divisible */

    if (A->length == 0)
        return gr_mpoly_zero(Q, ctx);

    /* The leading monomials must be known */
    if (gr_is_zero(B->coeffs, cctx) != T_FALSE ||
        gr_is_zero(A->coeffs, cctx) != T_FALSE)
        return GR_UNABLE;

    /* fall back to the single-threaded algorithm for small inputs, or when
       the coefficient ring does not allow concurrent operations */
    if (B->length < 2 || A->length < 2 || gr_ctx_is_threadsafe(cctx) != T_TRUE)
        return gr_mpoly_divides_heap(Q, A, B, ctx);

    thread_limit = A->length/32;

    num_handles = flint_request_threads(&handles, thread_limit);

    status = _gr_mpoly_divides_heap_threaded_pool(Q, A, B, ctx,
                                                         handles, num_handles);

    flint_give_back_threads(handles, num_handles);

    return status;
}
