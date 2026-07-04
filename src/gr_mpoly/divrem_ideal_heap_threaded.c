/*
    Copyright (C) 2018 Daniel Schultz
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
#include "fmpz_mpoly.h"
#include "mpoly.h"
#include "gr_generic.h"
#include "gr_mpoly.h"
#include "threaded_divmod.h"

/*
    Multithreaded gr_mpoly_divrem_ideal / gr_mpoly_divrem_ideal_weak: divide A
    by a list of divisors B[0..len-1], producing quotients Q[0..len-1] and a
    remainder R with A = sum_w Q[w]*B[w] + R, where a leading term is reduced
    by the first divisor whose leading monomial divides it (and, field-like,
    whose leading coefficient divides the term coefficient exactly); a term
    that no divisor can reduce is moved to R. See gr_mpoly/divrem_ideal.c for
    the serial algorithm and semantics (in particular: unlike plain divrem,
    an inexact field-like coefficient division is not a hard GR_DOMAIN error
    here -- it just means *that* divisor cannot reduce the term, so the term
    falls through to the next candidate divisor and ultimately to R if none
    work; this routine therefore never returns GR_DOMAIN).

    This shares the chunked producer/consumer pipeline infrastructure
    (gr_mpoly_ts, _gr_mpoly_mulsub_stripe1/stripe, stripe_fit_length,
    chunk_find_exp) from threaded_divmod.h / divides_heap_threaded.c, but
    uses its own chunk/base/worker-arg types (defined locally, not shared):
    a chunk's running partial remainder must be reduced against *every*
    divisor's newly available quotient increment (so chunk_mulsub below
    simply calls the existing, unmodified _gr_mpoly_mulsub_stripe1/stripe
    once per divisor, chaining the output of one call into the next), and a
    chunk's finalisation step must run a genuinely new bounded multi-chain
    heap (_gr_mpoly_divrem_ideal_stripe1/stripe below, a direct descendant
    of _gr_mpoly_divrem_ideal_mp1/mp in divrem_ideal.c with the same
    "restrict to exponents >= emin" transformation already used to derive
    the single-divisor stripe leaves from their heap1/heap counterparts).
*/

/* gather one heap node tagged with which divisor (p) it belongs to; the
   single dividend node (if any, always p == -1) becomes the dot initial
   value. Mirrors IDEAL_FETCH_NODE in divrem_ideal.c. */
#define IDEAL_STRIPE_FETCH_NODE(SET_SHALLOW) \
    do { \
        store[store_len++] = x->i; \
        store[store_len++] = x->j; \
        store[store_len++] = x->p; \
        if (x->i == -UWORD(1)) \
        { \
            have_dividend = 1; \
            dividend_j = x->j; \
        } \
        else \
        { \
            hinds[x->p][x->i] |= WORD(1); \
            SET_SHALLOW(dot_a, dot_len, poly3[x->p]->coeffs, x->i); \
            SET_SHALLOW(dot_b, dot_len, Qcoeff[x->p], x->j); \
            dot_len++; \
        } \
    } while (0)

/* field-like exact quotient coefficient = val / lc(B[w]) */
#define IDEAL_STRIPE_QCOEFF(dst, w, val) \
    (lc_is_one[w]  ? gr_set(dst, val, cctx) \
     : lc_is_unit[w] ? gr_mul(dst, val, GR_ENTRY(lc_inv, w, sz), cctx) \
                     : gr_div(dst, val, GR_ENTRY(poly3[w]->coeffs, 0, sz), cctx))


/*
    A chunk holds an exponent range on the dividend, exactly like
    divides_heap_chunk_struct, except that the persistent per-B-term
    "startidx"/"endidx" state used by chunk_mulsub (see the comment there)
    now needs one entry per divisor, since a chunk is reduced against
    *every* divisor's growing quotient.
*/
typedef struct _ideal_chunk_struct
{
    gr_mpoly_t polyC;
    struct _ideal_chunk_struct * next;
    ulong * emin;
    ulong * emax;
    slong * startidx;  /* array of length H->len */
    slong * endidx;    /* array of length H->len */
    slong * mq;        /* array of length H->len: consumed-so-far per divisor */
    int upperclosed;
    volatile int lock;
    volatile int producer;
    volatile int done;
    int Cinited;
} ideal_chunk_struct;

typedef ideal_chunk_struct ideal_chunk_t[1];

/*
    The base struct: a linked list of chunks, one thread-safe growable
    quotient array per divisor, one thread-safe growable remainder array,
    and the divisor data (shallow, repacked to a common bit width) and
    precomputed per-divisor leading-coefficient data.
*/
typedef struct
{
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif
    ideal_chunk_struct * head;
    ideal_chunk_struct * tail;
    ideal_chunk_struct * volatile cur;
    gr_mpoly_t polyA;
    gr_mpoly_struct * const * polyB; /* array of length `len`, shallow */
    gr_mpoly_ts_struct * polyQ;      /* array of length `len` */
    gr_mpoly_ts_t polyR;
    gr_mpoly_ctx_struct * ctx;
    gr_ctx_struct * cctx;
    slong len;                        /* number of divisors */
    slong length;                     /* number of chunks (unrelated to len) */
    slong N;
    flint_bitcnt_t bits;
    ulong * cmpmask;
    int * lc_is_one;    /* array of length `len` */
    int * lc_is_unit;   /* array of length `len` */
    gr_ptr lc_inv;      /* vector of length `len` */
    int nonfield;
    volatile int failed;
    volatile int have_unable;
} ideal_base_struct;

typedef ideal_base_struct ideal_base_t[1];

/*
    The worker struct has the (shared, unmodified) stripe scratch memory
    used for chunk_mulsub, plus one scratch poly per divisor (for a newly
    finalised quotient increment) and one for a newly finalised remainder
    increment.
*/
typedef struct
{
    ideal_base_struct * H;
    _gr_mpoly_stripe_t S;
    gr_mpoly_t polyT1;   /* chunk_mulsub scratch, reused per divisor */
    gr_mpoly_struct * polyT2; /* array of length H->len */
    gr_mpoly_t polyT3;   /* remainder increment scratch */
} ideal_worker_arg_struct;

typedef ideal_worker_arg_struct ideal_worker_arg_t[1];


static void ideal_base_init(ideal_base_t H)
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

static void ideal_chunk_clear(ideal_chunk_t L, ideal_base_t H)
{
    if (L->Cinited)
        gr_mpoly_clear(L->polyC, H->ctx);
}

/* Finalise H into (Q[0..len-1], R). Always clears/frees H's chunk list and
   thread-safe arrays. */
static int ideal_base_clear(gr_mpoly_struct * const * Q, gr_mpoly_t R, ideal_base_t H)
{
    ideal_chunk_struct * L = H->head;
    slong w;
    int status;

    while (L != NULL)
    {
        ideal_chunk_struct * nextL = L->next;
        ideal_chunk_clear(L, H);
        flint_free(L->startidx);
        flint_free(L);
        L = nextL;
    }
    H->head = NULL;
    H->tail = NULL;
    H->cur = NULL;
    H->length = 0;

    status = H->have_unable ? GR_UNABLE : GR_SUCCESS;

    if (status != GR_SUCCESS)
    {
        for (w = 0; w < H->len; w++)
        {
            GR_IGNORE(gr_mpoly_zero(Q[w], H->ctx));
            gr_mpoly_ts_clear(H->polyQ + w, H->cctx);
        }
        GR_IGNORE(gr_mpoly_zero(R, H->ctx));
        gr_mpoly_ts_clear(H->polyR, H->cctx);
    }
    else
    {
        for (w = 0; w < H->len; w++)
            gr_mpoly_ts_clear_poly(Q[w], H->polyQ + w, H->cctx);
        gr_mpoly_ts_clear_poly(R, H->polyR, H->cctx);
    }

    return status;
}

static void ideal_base_add_chunk(ideal_base_t H, ideal_chunk_t L)
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
        ideal_chunk_struct * tail = H->tail;
        FLINT_ASSERT(tail->next == NULL);
        tail->next = L;
        H->tail = L;
    }
    H->length++;
}


/* Bundles the (coeffs, exps, alloc, exps_alloc, length) quintuple for one
   growable output array (a quotient increment for one divisor, or the
   remainder increment), avoiding a long parallel-array parameter list. */
typedef struct
{
    gr_ptr coeffs;
    ulong * exps;
    slong alloc;
    slong exps_alloc;
    slong len;
} ideal_stripe_out_struct;

typedef ideal_stripe_out_struct ideal_stripe_out_t[1];


/*
    Bounded ("stripe") multi-divisor division-with-remainder, single-word
    exponent version: a direct descendant of _gr_mpoly_divrem_ideal_mp1
    (divrem_ideal.c), restricted to exponents >= emin in exactly the way
    the single-divisor stripe leaves restrict _gr_mpoly_divrem_heap1, and
    producing fresh local output increments (Qout[0..len-1], Wout) instead
    of writing into shared/growing state.

    Never a hard failure due to "no divisor reduces this term" (that is
    precisely what routes it to Wout instead); *res_status only ever
    reports GR_SUCCESS or GR_UNABLE (this routine, like
    _gr_mpoly_divrem_ideal_mp1, never returns GR_DOMAIN -- an inexact
    field-like coefficient division for one candidate divisor just means
    that divisor cannot reduce the term, not that the whole computation is
    impossible).
*/
static void _gr_mpoly_divrem_ideal_stripe1(
    ideal_stripe_out_t * Qout, ideal_stripe_out_t Wout,
    gr_srcptr Acoeff, const ulong * Aexp, slong Alen,
    gr_mpoly_struct * const * poly3, slong len,
    ulong emin, ulong cmpmask, flint_bitcnt_t bits,
    int nonfield, int * lc_is_one, int * lc_is_unit, gr_srcptr lc_inv,
    gr_mpoly_ctx_t ctx, int * res_status)
{
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    slong sz = cctx->sizeof_elem;
    int have_fast_dot = (GR_VEC_DOT_OP(cctx, VEC_DOT) != (gr_method_vec_dot_op) gr_generic_vec_dot);
    slong i, j, p, w, len3 = 0;
    slong next_loc, heap_len = 2;
    mpoly_heap1_s * heap;
    mpoly_nheap_t ** chains, * chains_ptr;
    slong ** hinds, * hinds_ptr;
    mpoly_nheap_t * x;
    ulong texp;
    ulong mask;
    slong * store, * store_base, store_len;
    slong * k, * s;
    gr_ptr acc, pp, rem, dot_a, dot_b;
    slong dot_len;
    int have_dividend, cstatus;
    slong dividend_j = 0;
    gr_ptr * Qcoeff;
    gr_ptr r_coeff;
    ulong * r_exp;
    slong r_len;
    int status = GR_SUCCESS;
    TMP_INIT;

    TMP_START;

    GR_TMP_INIT3(acc, pp, rem, cctx);

    Qcoeff = (gr_ptr *) TMP_ALLOC(len*sizeof(gr_ptr));
    for (w = 0; w < len; w++)
        Qcoeff[w] = Qout[w]->coeffs;
    r_coeff = Wout->coeffs;
    r_exp = Wout->exps;
    r_len = 0;

    chains = (mpoly_nheap_t **) TMP_ALLOC(len*sizeof(mpoly_nheap_t *));
    hinds = (slong **) TMP_ALLOC(len*sizeof(slong *));
    for (w = 0; w < len; w++)
        len3 += poly3[w]->length;
    chains_ptr = (mpoly_nheap_t *) TMP_ALLOC(len3*sizeof(mpoly_nheap_t));
    hinds_ptr = (slong *) TMP_ALLOC(len3*sizeof(slong));

    len3 = 0;
    for (w = 0; w < len; w++)
    {
        chains[w] = chains_ptr + len3;
        hinds[w] = hinds_ptr + len3;
        len3 += poly3[w]->length;
        for (i = 0; i < poly3[w]->length; i++)
            hinds[w][i] = 1;
    }

    next_loc = len3 + 4;
    heap = (mpoly_heap1_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap1_s));
    store = store_base = (slong *) TMP_ALLOC(3*len3*sizeof(slong));

    dot_a = flint_malloc(2 * len3 * sz);
    dot_b = GR_ENTRY(dot_a, len3, sz);

    k = (slong *) TMP_ALLOC(len*sizeof(slong));
    s = (slong *) TMP_ALLOC(len*sizeof(slong));
    for (w = 0; w < len; w++)
    {
        k[w] = -WORD(1);
        s[w] = poly3[w]->length;
    }

    mask = mpoly_overflow_mask_sp(bits);

    x = chains[0] + 0;
    x->i = -UWORD(1);
    x->j = 0;
    x->p = -WORD(1);
    x->next = NULL;
    HEAP_ASSIGN(heap[1], Aexp[0], x);

    FLINT_ASSERT(mpoly_monomial_cmp1(Aexp[0], emin, cmpmask) >= 0);

    while (heap_len > 1)
    {
        ulong exp = heap[1].exp;

        if (mpoly_monomial_overflows1(exp, mask))
        {
            status |= GR_UNABLE;
            goto unable;
        }

        FLINT_ASSERT(mpoly_monomial_cmp1(exp, emin, cmpmask) >= 0);

        store_len = 0;
        dot_len = 0;
        have_dividend = 0;

        if (have_fast_dot)
        {
            do
            {
                x = _mpoly_heap_pop1(heap, &heap_len, cmpmask);
                do
                {
                    if (sz == 1)       { IDEAL_STRIPE_FETCH_NODE(SET_SHALLOW1); }
                    else if (sz == 2)  { IDEAL_STRIPE_FETCH_NODE(SET_SHALLOW2); }
                    else if (sz == 4)  { IDEAL_STRIPE_FETCH_NODE(SET_SHALLOW4); }
                    else if (sz == 8)  { IDEAL_STRIPE_FETCH_NODE(SET_SHALLOW8); }
                    else if (sz == 16) { IDEAL_STRIPE_FETCH_NODE(SET_SHALLOW16); }
                    else               { IDEAL_STRIPE_FETCH_NODE(SET_SHALLOW_GENERIC); }
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
                    store[store_len++] = x->p;
                    if (x->i == -UWORD(1))
                        status |= gr_add(acc, acc, GR_ENTRY(Acoeff, x->j, sz), cctx);
                    else
                    {
                        hinds[x->p][x->i] |= WORD(1);
                        status |= gr_mul(pp, GR_ENTRY(poly3[x->p]->coeffs, x->i, sz),
                                             GR_ENTRY(Qcoeff[x->p], x->j, sz), cctx);
                        status |= gr_sub(acc, acc, pp, cctx);
                    }
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);
        }

        /* process popped nodes */
        while (store_len > 0)
        {
            p = store_base[--store_len];
            j = store_base[--store_len];
            i = store_base[--store_len];

            if (i == -WORD(1))
            {
                if (j + 1 < Alen)
                {
                    x = chains[0] + 0;
                    x->i = -UWORD(1);
                    x->j = j + 1;
                    x->p = p;
                    x->next = NULL;

                    FLINT_ASSERT(mpoly_monomial_cmp1(Aexp[x->j], emin, cmpmask) >= 0);

                    _mpoly_heap_insert1(heap, Aexp[x->j], x,
                                                &next_loc, &heap_len, cmpmask);
                }
            }
            else
            {
                if ((i + 1 < poly3[p]->length) && (hinds[p][i + 1] == 2*j + 1))
                {
                    x = chains[p] + i + 1;
                    x->i = i + 1;
                    x->j = j;
                    x->p = p;
                    x->next = NULL;
                    hinds[p][x->i] = 2*(x->j + 1) + 0;

                    texp = poly3[p]->exps[x->i] + Qout[p]->exps[x->j];
                    if (mpoly_monomial_cmp1(texp, emin, cmpmask) >= 0)
                        _mpoly_heap_insert1(heap, texp, x, &next_loc, &heap_len, cmpmask);
                    else
                        hinds[p][x->i] |= 1;
                }
                if (j == k[p])
                {
                    s[p]++;
                }
                else if (((hinds[p][i] & 1) == 1) &&
                         ((i == 1) || (hinds[p][i - 1] >= 2*(j + 2) + 1)))
                {
                    x = chains[p] + i;
                    x->i = i;
                    x->j = j + 1;
                    x->p = p;
                    x->next = NULL;
                    hinds[p][x->i] = 2*(x->j + 1) + 0;

                    texp = poly3[p]->exps[x->i] + Qout[p]->exps[x->j];
                    if (mpoly_monomial_cmp1(texp, emin, cmpmask) >= 0)
                        _mpoly_heap_insert1(heap, texp, x, &next_loc, &heap_len, cmpmask);
                    else
                        hinds[p][x->i] |= 1;
                }
            }
        }

        /* an unknown zero-status is treated as "not zero" -- see the
           discussion at gr_mpoly_divexact for why this is safe: we only
           need the coefficient division below to succeed, not to know for
           certain that acc is nonzero. */
        if (gr_is_zero(acc, cctx) == T_TRUE)
            continue;

        /* try to reduce the term acc*x^exp by the divisors in order */
        {
            int div_flag = 0;

            for (w = 0; w < len; w++)
            {
                int d1, d2;

                d1 = mpoly_monomial_divides1(&texp, exp, poly3[w]->exps[0], mask);

                if (!d1)
                    continue;

                _gr_mpoly_fit_length(&Qcoeff[w], &Qout[w]->alloc,
                                     &Qout[w]->exps, &Qout[w]->exps_alloc, 1, k[w] + 2, ctx);

                if (nonfield)
                {
                    cstatus = gr_euclidean_divrem(GR_ENTRY(Qcoeff[w], k[w] + 1, sz),
                                     rem, acc, GR_ENTRY(poly3[w]->coeffs, 0, sz), cctx);
                    if (cstatus != GR_SUCCESS) { status |= cstatus; goto unable; }
                    gr_swap(acc, rem, cctx);
                    d2 = (gr_is_zero(GR_ENTRY(Qcoeff[w], k[w] + 1, sz), cctx) != T_TRUE);
                }
                else
                {
                    cstatus = IDEAL_STRIPE_QCOEFF(GR_ENTRY(Qcoeff[w], k[w] + 1, sz), w, acc);
                    if (cstatus == GR_DOMAIN)
                        continue;
                    if (cstatus != GR_SUCCESS) { status |= cstatus; goto unable; }
                    d2 = (gr_is_zero(GR_ENTRY(Qcoeff[w], k[w] + 1, sz), cctx) != T_TRUE);
                }

                if (d2)
                {
                    k[w]++;
                    Qout[w]->exps[k[w]] = texp;

                    if (s[w] > 1)
                    {
                        x = chains[w] + 1;
                        x->i = 1;
                        x->j = k[w];
                        x->p = w;
                        x->next = NULL;
                        hinds[w][x->i] = 2*(x->j + 1) + 0;

                        texp = poly3[w]->exps[1] + Qout[w]->exps[k[w]];
                        if (mpoly_monomial_cmp1(texp, emin, cmpmask) >= 0)
                            _mpoly_heap_insert1(heap, texp, x, &next_loc, &heap_len, cmpmask);
                        else
                            hinds[w][x->i] |= 1;
                    }
                    s[w] = 1;
                }

                if (nonfield)
                {
                    if (gr_is_zero(acc, cctx) == T_TRUE)
                    {
                        div_flag = 1;
                        break;
                    }
                }
                else
                {
                    div_flag = 1;
                    break;
                }
            }

            if (!div_flag)
            {
                _gr_mpoly_fit_length(&r_coeff, &Wout->alloc, &r_exp, &Wout->exps_alloc, 1, r_len + 1, ctx);
                status |= gr_set(GR_ENTRY(r_coeff, r_len, sz), acc, cctx);
                r_exp[r_len] = exp;
                r_len++;
            }
        }
    }

    *res_status = status;
    goto cleanup;

unable:
    for (w = 0; w < len; w++)
        k[w] = -1;
    r_len = 0;
    *res_status = status;

cleanup:

    GR_TMP_CLEAR3(acc, pp, rem, cctx);
    flint_free(dot_a);

    for (w = 0; w < len; w++)
    {
        Qout[w]->coeffs = Qcoeff[w];
        Qout[w]->len = k[w] + 1;
    }
    Wout->coeffs = r_coeff;
    Wout->exps = r_exp;
    Wout->len = r_len;

    TMP_END;
}


/* Multi-word exponent version of the above. */
static void _gr_mpoly_divrem_ideal_stripe(
    ideal_stripe_out_t * Qout, ideal_stripe_out_t Wout,
    gr_srcptr Acoeff, const ulong * Aexp, slong Alen,
    gr_mpoly_struct * const * poly3, slong len,
    const ulong * emin, const ulong * cmpmask, slong N, flint_bitcnt_t bits,
    int nonfield, int * lc_is_one, int * lc_is_unit, gr_srcptr lc_inv,
    gr_mpoly_ctx_t ctx, int * res_status)
{
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    slong sz = cctx->sizeof_elem;
    int have_fast_dot = (GR_VEC_DOT_OP(cctx, VEC_DOT) != (gr_method_vec_dot_op) gr_generic_vec_dot);
    slong i, j, p, w, len3 = 0;
    slong next_loc, heap_len = 2;
    mpoly_heap_s * heap;
    mpoly_nheap_t ** chains, * chains_ptr;
    slong ** hinds, * hinds_ptr;
    mpoly_nheap_t * x;
    ulong * exp, * exps, * texp;
    ulong ** exp_list;
    slong exp_next;
    ulong mask;
    slong * store, * store_base, store_len;
    slong * k, * s;
    gr_ptr acc, pp, rem, dot_a, dot_b;
    slong dot_len;
    int have_dividend, cstatus;
    slong dividend_j = 0;
    gr_ptr * Qcoeff;
    gr_ptr r_coeff;
    ulong * r_exp;
    slong r_len;
    int status = GR_SUCCESS;
    TMP_INIT;

    TMP_START;

    GR_TMP_INIT3(acc, pp, rem, cctx);

    Qcoeff = (gr_ptr *) TMP_ALLOC(len*sizeof(gr_ptr));
    for (w = 0; w < len; w++)
        Qcoeff[w] = Qout[w]->coeffs;
    r_coeff = Wout->coeffs;
    r_exp = Wout->exps;
    r_len = 0;

    chains = (mpoly_nheap_t **) TMP_ALLOC(len*sizeof(mpoly_nheap_t *));
    hinds = (slong **) TMP_ALLOC(len*sizeof(slong *));
    for (w = 0; w < len; w++)
        len3 += poly3[w]->length;
    chains_ptr = (mpoly_nheap_t *) TMP_ALLOC(len3*sizeof(mpoly_nheap_t));
    hinds_ptr = (slong *) TMP_ALLOC(len3*sizeof(slong));

    len3 = 0;
    for (w = 0; w < len; w++)
    {
        chains[w] = chains_ptr + len3;
        hinds[w] = hinds_ptr + len3;
        len3 += poly3[w]->length;
        for (i = 0; i < poly3[w]->length; i++)
            hinds[w][i] = 1;
    }

    next_loc = len3 + 4;
    heap = (mpoly_heap_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap_s));
    store = store_base = (slong *) TMP_ALLOC(3*len3*sizeof(slong));

    dot_a = flint_malloc(2 * len3 * sz);
    dot_b = GR_ENTRY(dot_a, len3, sz);

    exps = (ulong *) TMP_ALLOC(len3*N*sizeof(ulong));
    exp_list = (ulong **) TMP_ALLOC(len3*sizeof(ulong *));
    exp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    texp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    exp_next = 0;
    for (i = 0; i < len3; i++)
        exp_list[i] = exps + i*N;

    k = (slong *) TMP_ALLOC(len*sizeof(slong));
    s = (slong *) TMP_ALLOC(len*sizeof(slong));
    for (w = 0; w < len; w++)
    {
        k[w] = -WORD(1);
        s[w] = poly3[w]->length;
    }

    mask = bits <= FLINT_BITS ? mpoly_overflow_mask_sp(bits) : 0;

    x = chains[0] + 0;
    x->i = -UWORD(1);
    x->j = 0;
    x->p = -WORD(1);
    x->next = NULL;
    heap[1].next = x;
    heap[1].exp = exp_list[exp_next++];
    mpoly_monomial_set(heap[1].exp, Aexp, N);

    FLINT_ASSERT(mpoly_monomial_cmp(Aexp, emin, N, cmpmask) >= 0);

    while (heap_len > 1)
    {
        mpoly_monomial_set(exp, heap[1].exp, N);

        if (bits <= FLINT_BITS)
        {
            if (mpoly_monomial_overflows(exp, N, mask))
            {
                status |= GR_UNABLE;
                goto unable;
            }
        }
        else
        {
            if (mpoly_monomial_overflows_mp(exp, N, bits))
            {
                status |= GR_UNABLE;
                goto unable;
            }
        }

        FLINT_ASSERT(mpoly_monomial_cmp(exp, emin, N, cmpmask) >= 0);

        store_len = 0;
        dot_len = 0;
        have_dividend = 0;

        if (have_fast_dot)
        {
            do
            {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
                do
                {
                    if (sz == 1)       { IDEAL_STRIPE_FETCH_NODE(SET_SHALLOW1); }
                    else if (sz == 2)  { IDEAL_STRIPE_FETCH_NODE(SET_SHALLOW2); }
                    else if (sz == 4)  { IDEAL_STRIPE_FETCH_NODE(SET_SHALLOW4); }
                    else if (sz == 8)  { IDEAL_STRIPE_FETCH_NODE(SET_SHALLOW8); }
                    else if (sz == 16) { IDEAL_STRIPE_FETCH_NODE(SET_SHALLOW16); }
                    else               { IDEAL_STRIPE_FETCH_NODE(SET_SHALLOW_GENERIC); }
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
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
                do
                {
                    store[store_len++] = x->i;
                    store[store_len++] = x->j;
                    store[store_len++] = x->p;
                    if (x->i == -UWORD(1))
                        status |= gr_add(acc, acc, GR_ENTRY(Acoeff, x->j, sz), cctx);
                    else
                    {
                        hinds[x->p][x->i] |= WORD(1);
                        status |= gr_mul(pp, GR_ENTRY(poly3[x->p]->coeffs, x->i, sz),
                                             GR_ENTRY(Qcoeff[x->p], x->j, sz), cctx);
                        status |= gr_sub(acc, acc, pp, cctx);
                    }
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));
        }

        /* process popped nodes */
        while (store_len > 0)
        {
            p = store_base[--store_len];
            j = store_base[--store_len];
            i = store_base[--store_len];

            if (i == -WORD(1))
            {
                if (j + 1 < Alen)
                {
                    x = chains[0] + 0;
                    x->i = -UWORD(1);
                    x->j = j + 1;
                    x->p = p;
                    x->next = NULL;
                    mpoly_monomial_set(exp_list[exp_next], Aexp + x->j*N, N);

                    FLINT_ASSERT(mpoly_monomial_cmp(exp_list[exp_next], emin, N, cmpmask) >= 0);

                    exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
                }
            }
            else
            {
                if ((i + 1 < poly3[p]->length) && (hinds[p][i + 1] == 2*j + 1))
                {
                    x = chains[p] + i + 1;
                    x->i = i + 1;
                    x->j = j;
                    x->p = p;
                    x->next = NULL;
                    hinds[p][x->i] = 2*(x->j + 1) + 0;

                    mpoly_monomial_add_mp(exp_list[exp_next], poly3[p]->exps + x->i*N,
                                                        Qout[p]->exps + x->j*N, N);

                    if (mpoly_monomial_cmp(exp_list[exp_next], emin, N, cmpmask) >= 0)
                        exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                                 &next_loc, &heap_len, N, cmpmask);
                    else
                        hinds[p][x->i] |= 1;
                }
                if (j == k[p])
                {
                    s[p]++;
                }
                else if (((hinds[p][i] & 1) == 1) &&
                         ((i == 1) || (hinds[p][i - 1] >= 2*(j + 2) + 1)))
                {
                    x = chains[p] + i;
                    x->i = i;
                    x->j = j + 1;
                    x->p = p;
                    x->next = NULL;
                    hinds[p][x->i] = 2*(x->j + 1) + 0;

                    mpoly_monomial_add_mp(exp_list[exp_next], poly3[p]->exps + x->i*N,
                                                        Qout[p]->exps + x->j*N, N);

                    if (mpoly_monomial_cmp(exp_list[exp_next], emin, N, cmpmask) >= 0)
                        exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                                 &next_loc, &heap_len, N, cmpmask);
                    else
                        hinds[p][x->i] |= 1;
                }
            }
        }

        /* an unknown zero-status is treated as "not zero" -- see the
           discussion at gr_mpoly_divexact for why this is safe: we only
           need the coefficient division below to succeed, not to know for
           certain that acc is nonzero. */
        if (gr_is_zero(acc, cctx) == T_TRUE)
            continue;

        /* try to reduce the term acc*x^exp by the divisors in order */
        {
            int div_flag = 0;

            for (w = 0; w < len; w++)
            {
                int d1, d2;

                if (bits <= FLINT_BITS)
                    d1 = mpoly_monomial_divides(texp, exp, poly3[w]->exps, N, mask);
                else
                    d1 = mpoly_monomial_divides_mp(texp, exp, poly3[w]->exps, N, bits);

                if (!d1)
                    continue;

                _gr_mpoly_fit_length(&Qcoeff[w], &Qout[w]->alloc,
                                     &Qout[w]->exps, &Qout[w]->exps_alloc, N, k[w] + 2, ctx);

                if (nonfield)
                {
                    cstatus = gr_euclidean_divrem(GR_ENTRY(Qcoeff[w], k[w] + 1, sz),
                                     rem, acc, GR_ENTRY(poly3[w]->coeffs, 0, sz), cctx);
                    if (cstatus != GR_SUCCESS) { status |= cstatus; goto unable; }
                    gr_swap(acc, rem, cctx);
                    d2 = (gr_is_zero(GR_ENTRY(Qcoeff[w], k[w] + 1, sz), cctx) != T_TRUE);
                }
                else
                {
                    cstatus = IDEAL_STRIPE_QCOEFF(GR_ENTRY(Qcoeff[w], k[w] + 1, sz), w, acc);
                    if (cstatus == GR_DOMAIN)
                        continue;
                    if (cstatus != GR_SUCCESS) { status |= cstatus; goto unable; }
                    d2 = (gr_is_zero(GR_ENTRY(Qcoeff[w], k[w] + 1, sz), cctx) != T_TRUE);
                }

                if (d2)
                {
                    k[w]++;
                    mpoly_monomial_set(Qout[w]->exps + k[w]*N, texp, N);

                    if (s[w] > 1)
                    {
                        x = chains[w] + 1;
                        x->i = 1;
                        x->j = k[w];
                        x->p = w;
                        x->next = NULL;
                        hinds[w][x->i] = 2*(x->j + 1) + 0;

                        mpoly_monomial_add_mp(exp_list[exp_next], poly3[w]->exps + N,
                                                            Qout[w]->exps + k[w]*N, N);

                        if (mpoly_monomial_cmp(exp_list[exp_next], emin, N, cmpmask) >= 0)
                            exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                                     &next_loc, &heap_len, N, cmpmask);
                        else
                            hinds[w][x->i] |= 1;
                    }
                    s[w] = 1;
                }

                if (nonfield)
                {
                    if (gr_is_zero(acc, cctx) == T_TRUE)
                    {
                        div_flag = 1;
                        break;
                    }
                }
                else
                {
                    div_flag = 1;
                    break;
                }
            }

            if (!div_flag)
            {
                _gr_mpoly_fit_length(&r_coeff, &Wout->alloc, &r_exp, &Wout->exps_alloc, N, r_len + 1, ctx);
                status |= gr_set(GR_ENTRY(r_coeff, r_len, sz), acc, cctx);
                mpoly_monomial_set(r_exp + r_len*N, exp, N);
                r_len++;
            }
        }
    }

    *res_status = status;
    goto cleanup;

unable:
    for (w = 0; w < len; w++)
        k[w] = -1;
    r_len = 0;
    *res_status = status;

cleanup:

    GR_TMP_CLEAR3(acc, pp, rem, cctx);
    flint_free(dot_a);

    for (w = 0; w < len; w++)
    {
        Qout[w]->coeffs = Qcoeff[w];
        Qout[w]->len = k[w] + 1;
    }
    Wout->coeffs = r_coeff;
    Wout->exps = r_exp;
    Wout->len = r_len;

    TMP_END;
}


static slong ideal_chunk_find_exp(ulong * exp, slong a, const ideal_base_t H)
{
    return chunk_find_exp(exp, a, H->polyA->exps, H->polyA->length, H->N, H->cmpmask);
}

/*
    Reduce L->polyC by every divisor's newly available quotient increment,
    one divisor at a time, chaining the output of one call into the input
    of the next. Each call is exactly the existing, unmodified
    _gr_mpoly_mulsub_stripe1/stripe (see divides_heap_threaded.c): this
    routine has no notion of "ideal" beyond looping over the divisor list.
*/
static void ideal_chunk_mulsub(ideal_worker_arg_t W, ideal_chunk_t L, const slong * q_prev_length)
{
    ideal_base_struct * H = W->H;
    slong N = H->N;
    gr_ctx_struct * cctx = H->cctx;
    slong sz = cctx->sizeof_elem;
    const gr_mpoly_struct * A = H->polyA;
    gr_mpoly_struct * C;
    gr_mpoly_struct * T1 = W->polyT1;
    _gr_mpoly_stripe_struct * S = W->S;
    slong w;

    S->emin = L->emin;
    S->emax = L->emax;
    S->upperclosed = L->upperclosed;
    FLINT_ASSERT(S->N == N);

    if (!L->Cinited)
    {
        /*
            First touch: L->polyC has never been initialised. Seed it with
            the slice of the original dividend A lying in this chunk's
            exponent window, exactly like the single-divisor chunk_mulsub
            (divides_heap_threaded.c) does for its own first-touch case.
            (For simplicity we always do this as a plain copy rather than
            trying to fold it into the first divisor's mulsub call as the
            single-divisor version does -- this chunk-lifetime-one-off copy
            is not worth the extra complexity of special-casing "whichever
            divisor happens to have a pending increment first".)
        */
        slong startidx, stopidx, Alen;
        int status;

        if (L->upperclosed)
        {
            startidx = 0;
            stopidx = ideal_chunk_find_exp(L->emin, 1, H);
        }
        else
        {
            startidx = ideal_chunk_find_exp(L->emax, 1, H);
            stopidx = ideal_chunk_find_exp(L->emin, startidx, H);
        }
        Alen = stopidx - startidx;

        L->Cinited = 1;
        gr_mpoly_init3(L->polyC, 16 + Alen, H->bits, H->ctx);
        C = L->polyC;

        gr_mpoly_fit_length(C, Alen, H->ctx);
        status = _gr_vec_set(C->coeffs, GR_ENTRY(A->coeffs, startidx, sz), Alen, cctx);
        mpoly_copy_monomials(C->exps, A->exps + N*startidx, Alen, N);
        C->length = Alen;

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
    }

    C = L->polyC;

    for (w = 0; w < H->len; w++)
    {
        gr_mpoly_struct * B = H->polyB[w];
        gr_mpoly_ts_struct * Q = H->polyQ + w;
        slong mq = L->mq[w];
        slong new_length = q_prev_length[w];
        int status;

        if (new_length <= mq)
            continue;

        stripe_fit_length(S, new_length - mq);
        S->startidx = &L->startidx[w];
        S->endidx = &L->endidx[w];

        if (N == 1)
        {
            T1->length = _gr_mpoly_mulsub_stripe1(
                    &T1->coeffs, &T1->exps, &T1->coeffs_alloc, &T1->exps_alloc,
                    C->coeffs, C->exps, C->length, 1,
                    GR_ENTRY(Q->coeffs, mq, sz), Q->exps + N*mq, new_length - mq,
                    B->coeffs, B->exps, B->length,
                    S, &status);
        }
        else
        {
            T1->length = _gr_mpoly_mulsub_stripe(
                    &T1->coeffs, &T1->exps, &T1->coeffs_alloc, &T1->exps_alloc,
                    C->coeffs, C->exps, C->length, 1,
                    GR_ENTRY(Q->coeffs, mq, sz), Q->exps + N*mq, new_length - mq,
                    B->coeffs, B->exps, B->length,
                    S, &status);
        }

        gr_mpoly_swap(C, T1, H->ctx);

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

        L->mq[w] = new_length;
    }
}

static void ideal_trychunk(ideal_worker_arg_t W, ideal_chunk_t L)
{
    ideal_base_struct * H = W->H;
    slong N = H->N;
    gr_ctx_struct * cctx = H->cctx;
    slong sz = cctx->sizeof_elem;
    gr_mpoly_struct * C = L->polyC;
    slong * q_prev_length;
    const gr_mpoly_struct * A = H->polyA;
    slong w;

    if (L->done)
        return;

    q_prev_length = (slong *) flint_malloc(H->len*sizeof(slong));

    for (w = 0; w < H->len; w++)
        q_prev_length[w] = (H->polyQ + w)->length;

    for (w = 0; w < H->len; w++)
        if (q_prev_length[w] > L->mq[w])
            break;

    if (w < H->len)
    {
        if (L->producer == 0)
        {
            /* only bother with a small trickle if we are not the producer */
            slong total = 0;
            for (w = 0; w < H->len; w++)
                total += q_prev_length[w] - L->mq[w];
            if (total < 20)
            {
                flint_free(q_prev_length);
                return;
            }
        }

        ideal_chunk_mulsub(W, L, q_prev_length);
    }

    if (L->producer == 1)
    {
        ideal_chunk_struct * next;
        gr_srcptr Rcoeff;
        ulong * Rexp;
        slong Rlen;

        /* process any further quotient terms that trickled in */
        for (w = 0; w < H->len; w++)
            q_prev_length[w] = (H->polyQ + w)->length;
        for (w = 0; w < H->len; w++)
            if (q_prev_length[w] > L->mq[w])
                break;
        if (w < H->len)
            ideal_chunk_mulsub(W, L, q_prev_length);

        flint_free(q_prev_length);

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
                stopidx = ideal_chunk_find_exp(L->emin, 1, H);
            }
            else
            {
                startidx = ideal_chunk_find_exp(L->emax, 1, H);
                stopidx = ideal_chunk_find_exp(L->emin, startidx, H);
            }
            Rlen = stopidx - startidx;
            Rcoeff = GR_ENTRY(A->coeffs, startidx, sz);
            Rexp = A->exps + N*startidx;
        }

        if (Rlen > 0)
        {
            gr_mpoly_struct * T2 = W->polyT2;
            gr_mpoly_struct * T3 = W->polyT3;
            ideal_stripe_out_t * Qout;
            ideal_stripe_out_t Wout;
            int rstatus;

            Qout = (ideal_stripe_out_t *) flint_malloc(H->len*sizeof(ideal_stripe_out_t));
            for (w = 0; w < H->len; w++)
            {
                Qout[w]->coeffs = T2[w].coeffs;
                Qout[w]->exps = T2[w].exps;
                Qout[w]->alloc = T2[w].coeffs_alloc;
                Qout[w]->exps_alloc = T2[w].exps_alloc;
                Qout[w]->len = 0;
            }
            Wout->coeffs = T3->coeffs;
            Wout->exps = T3->exps;
            Wout->alloc = T3->coeffs_alloc;
            Wout->exps_alloc = T3->exps_alloc;
            Wout->len = 0;

            if (N == 1)
                _gr_mpoly_divrem_ideal_stripe1(Qout, Wout,
                        Rcoeff, Rexp, Rlen, H->polyB, H->len,
                        L->emin[0], H->cmpmask[0], H->bits,
                        H->nonfield, H->lc_is_one, H->lc_is_unit, H->lc_inv,
                        H->ctx, &rstatus);
            else
                _gr_mpoly_divrem_ideal_stripe(Qout, Wout,
                        Rcoeff, Rexp, Rlen, H->polyB, H->len,
                        L->emin, H->cmpmask, N, H->bits,
                        H->nonfield, H->lc_is_one, H->lc_is_unit, H->lc_inv,
                        H->ctx, &rstatus);

            for (w = 0; w < H->len; w++)
            {
                T2[w].coeffs = Qout[w]->coeffs;
                T2[w].exps = Qout[w]->exps;
                T2[w].coeffs_alloc = Qout[w]->alloc;
                T2[w].exps_alloc = Qout[w]->exps_alloc;
                T2[w].length = Qout[w]->len;
            }
            T3->coeffs = Wout->coeffs;
            T3->exps = Wout->exps;
            T3->coeffs_alloc = Wout->alloc;
            T3->exps_alloc = Wout->exps_alloc;
            T3->length = Wout->len;

            flint_free(Qout);

            if (rstatus != GR_SUCCESS)
            {
#if FLINT_USES_PTHREAD
                pthread_mutex_lock(&H->mutex);
#endif
                H->have_unable = 1;
                H->failed = 1;
#if FLINT_USES_PTHREAD
                pthread_mutex_unlock(&H->mutex);
#endif
                return;
            }

            for (w = 0; w < H->len; w++)
                if (T2[w].length > 0)
                    gr_mpoly_ts_append(H->polyQ + w, T2[w].coeffs, T2[w].exps, T2[w].length, N, cctx);
            if (T3->length > 0)
                gr_mpoly_ts_append(H->polyR, T3->coeffs, T3->exps, T3->length, N, cctx);
        }

        next = L->next;
        H->length--;
        H->cur = next;

        if (next != NULL)
            next->producer = 1;

        L->producer = 0;
        L->done = 1;
    }
    else
    {
        flint_free(q_prev_length);
    }
}

static void ideal_worker_loop(void * varg)
{
    ideal_worker_arg_struct * W = (ideal_worker_arg_struct *) varg;
    ideal_base_struct * H = W->H;
    _gr_mpoly_stripe_struct * S = W->S;
    gr_mpoly_struct * T1 = W->polyT1;
    gr_mpoly_struct * T2 = W->polyT2;
    slong N = H->N;
    slong w, maxBlen = 0;

    S->N = N;
    S->bits = H->bits;
    S->ctx = H->ctx;
    S->cctx = H->cctx;
    S->sz = H->cctx->sizeof_elem;
    S->cmpmask = H->cmpmask;
    S->have_fast_dot = (GR_VEC_DOT_OP(H->cctx, VEC_DOT) != (gr_method_vec_dot_op) gr_generic_vec_dot);
    S->lc_is_one = 0;
    S->lc_is_unit = 0;
    S->lc_inv = NULL;
    S->nonfield = H->nonfield;
    S->unchecked = 0;
    S->big_mem_alloc = 0;
    S->big_mem = NULL;

    for (w = 0; w < H->len; w++)
        maxBlen = FLINT_MAX(maxBlen, H->polyB[w]->length);
    stripe_fit_length(S, maxBlen);

    gr_mpoly_init3(T1, 16, H->bits, H->ctx);
    for (w = 0; w < H->len; w++)
        gr_mpoly_init3(T2 + w, 16, H->bits, H->ctx);
    gr_mpoly_init3(W->polyT3, 16, H->bits, H->ctx);

    while (!H->failed)
    {
        ideal_chunk_struct * L;
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
                ideal_trychunk(W, L);
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
    for (w = 0; w < H->len; w++)
        gr_mpoly_clear(T2 + w, H->ctx);
    gr_mpoly_clear(W->polyT3, H->ctx);
    flint_free(S->big_mem);
}


static int
_gr_mpoly_divrem_ideal_serial(gr_mpoly_struct * const * Q, gr_mpoly_t R,
    const gr_mpoly_t A, gr_mpoly_struct * const * B, slong len,
    int nonfield, gr_mpoly_ctx_t ctx)
{
    /*
        Call the guaranteed-serial kernel directly, NOT the public
        gr_mpoly_divrem_ideal/_weak dispatchers: those may route large
        instances straight back into this threaded engine, which would
        recurse indefinitely whenever this function is reached as a
        fallback from within the threaded engine itself (A, B unchanged).
    */
    return _gr_mpoly_divrem_ideal((gr_mpoly_struct **) Q, R, A, B, len, nonfield, ctx);
}


static int _gr_mpoly_divrem_ideal_heap_threaded_pool(
    gr_mpoly_struct * const * Q,
    gr_mpoly_t R,
    const gr_mpoly_t A,
    gr_mpoly_struct * const * B,
    slong len,
    gr_mpoly_ctx_t ctx,
    int nonfield,
    const thread_pool_handle * handles,
    slong num_handles)
{
    mpoly_ctx_struct * mctx = GR_MPOLY_MCTX(ctx);
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    slong sz = cctx->sizeof_elem;
    int status;
    fmpz_mpoly_ctx_t zctx;
    fmpz_mpoly_t S;
    slong i, w, N;
    flint_bitcnt_t exp_bits;
    ulong * cmpmask;
    ulong * Aexp;
    ulong ** Bexp;
    int freeAexp, * freeBexp;
    ideal_worker_arg_struct * worker_args;
    ideal_base_t H;
    gr_mpoly_struct * polyB_storage;
    gr_mpoly_struct ** polyB_ptrs;
    TMP_INIT;

    if (A->length < 2 || B[0]->length < 2)
        return _gr_mpoly_divrem_ideal_serial(Q, R, A, B, len, nonfield, ctx);

    TMP_START;

    exp_bits = MPOLY_MIN_BITS;
    exp_bits = FLINT_MAX(exp_bits, A->bits);
    for (w = 0; w < len; w++)
        exp_bits = FLINT_MAX(exp_bits, B[w]->bits);
    exp_bits = mpoly_fix_bits(exp_bits, mctx);

    N = mpoly_words_per_exp(exp_bits, mctx);
    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, exp_bits, mctx);

    Aexp = A->exps;
    freeAexp = 0;
    if (exp_bits > A->bits)
    {
        freeAexp = 1;
        Aexp = (ulong *) flint_malloc(N*A->length*sizeof(ulong));
        mpoly_repack_monomials(Aexp, exp_bits, A->exps, A->bits, A->length, mctx);
    }

    Bexp = (ulong **) TMP_ALLOC(len*sizeof(ulong *));
    freeBexp = (int *) TMP_ALLOC(len*sizeof(int));
    for (w = 0; w < len; w++)
    {
        Bexp[w] = B[w]->exps;
        freeBexp[w] = 0;
        if (exp_bits > B[w]->bits)
        {
            freeBexp[w] = 1;
            Bexp[w] = (ulong *) flint_malloc(N*B[w]->length*sizeof(ulong));
            mpoly_repack_monomials(Bexp[w], exp_bits, B[w]->exps, B[w]->bits, B[w]->length, mctx);
        }
    }

    polyB_storage = (gr_mpoly_struct *) flint_malloc(len*sizeof(gr_mpoly_struct));
    polyB_ptrs = (gr_mpoly_struct **) flint_malloc(len*sizeof(gr_mpoly_struct *));
    for (w = 0; w < len; w++)
    {
        polyB_ptrs[w] = polyB_storage + w;
        polyB_storage[w].coeffs = B[w]->coeffs;
        polyB_storage[w].exps = Bexp[w];
        polyB_storage[w].bits = exp_bits;
        polyB_storage[w].length = B[w]->length;
        polyB_storage[w].coeffs_alloc = B[w]->coeffs_alloc;
        polyB_storage[w].exps_alloc = N*B[w]->length;
    }

    fmpz_mpoly_ctx_init(zctx, mctx->nvars, mctx->ord);
    fmpz_mpoly_init(S, zctx);

    /*
        mpoly_divides_select_exps only needs *a* representative divisor to
        pick reasonable chunk boundaries -- these are purely a load
        balancing hint (unlike the plain-divide case, correctness never
        depends on them), so B[0] is a fine, simple choice regardless of
        how many divisors there are.
    */
    if (mpoly_divides_select_exps(S, zctx, num_handles,
                                   Aexp, A->length, Bexp[0], B[0]->length, exp_bits))
    {
        /* select_exps's failure notion ("provably not exact") does not
           apply to divrem_ideal (which always succeeds via the
           remainder), but a degenerate/pathological exponent spread is
           not worth trying to patch up -- fall back to the serial engine. */
        status = _gr_mpoly_divrem_ideal_serial(Q, R, A, B, len, nonfield, ctx);
        goto cleanup1;
    }

    ideal_base_init(H);

    H->polyA->coeffs = A->coeffs;
    H->polyA->exps = Aexp;
    H->polyA->bits = exp_bits;
    H->polyA->length = A->length;
    H->polyA->coeffs_alloc = A->coeffs_alloc;
    H->polyA->exps_alloc = N*A->length;

    H->polyB = polyB_ptrs;
    H->len = len;
    H->ctx = ctx;
    H->cctx = cctx;
    H->bits = exp_bits;
    H->N = N;
    H->cmpmask = cmpmask;
    H->nonfield = nonfield;
    H->failed = 0;
    H->have_unable = 0;

    H->polyQ = (gr_mpoly_ts_struct *) flint_malloc(len*sizeof(gr_mpoly_ts_struct));
    for (w = 0; w < len; w++)
        gr_mpoly_ts_init(H->polyQ + w, NULL, NULL, 0, exp_bits, N, cctx);
    gr_mpoly_ts_init(H->polyR, NULL, NULL, 0, exp_bits, N, cctx);

    H->lc_inv = flint_malloc(len*sz);
    _gr_vec_init(H->lc_inv, len, cctx);
    H->lc_is_one = (int *) flint_malloc(len*sizeof(int));
    H->lc_is_unit = (int *) flint_malloc(len*sizeof(int));
    for (w = 0; w < len; w++)
    {
        H->lc_is_one[w] = H->lc_is_unit[w] = 0;
        if (!nonfield)
        {
            H->lc_is_one[w] = (gr_is_one(GR_ENTRY(B[w]->coeffs, 0, sz), cctx) == T_TRUE);
            H->lc_is_unit[w] = H->lc_is_one[w] ||
                (gr_inv(GR_ENTRY(H->lc_inv, w, sz), GR_ENTRY(B[w]->coeffs, 0, sz), cctx) == GR_SUCCESS);
        }
    }

    for (i = 0; i + 1 < S->length; i++)
    {
        ideal_chunk_struct * L;
        L = (ideal_chunk_struct *) flint_malloc(sizeof(ideal_chunk_struct));
        L->startidx = (slong *) flint_malloc(3*len*sizeof(slong));
        L->endidx = L->startidx + len;
        L->mq = L->startidx + 2*len;
        for (w = 0; w < len; w++)
        {
            L->startidx[w] = B[w]->length;
            L->endidx[w] = B[w]->length;
            L->mq[w] = 0;
        }
        L->emax = S->exps + N*i;
        L->emin = S->exps + N*(i + 1);
        L->upperclosed = 0;
        L->producer = 0;
        L->Cinited = 0;
        L->done = 0;
        L->lock = -2;
        ideal_base_add_chunk(H, L);
    }

    H->head->upperclosed = 1;
    H->head->producer = 1;
    H->cur = H->head;

    /*
        No "initial burst" pre-pass here, for the same reason
        gr_mpoly_divrem_heap_threaded doesn't need one (see the long
        comment there): the head chunk, already marked as producer, can be
        fully resolved from scratch on first touch regardless of how many
        terms of A it covers.
    */

#if FLINT_USES_PTHREAD
    pthread_mutex_init(&H->mutex, NULL);
#endif

    worker_args = (ideal_worker_arg_struct *) flint_malloc((num_handles + 1)
                                                        *sizeof(ideal_worker_arg_struct));
    for (i = 0; i < num_handles; i++)
    {
        (worker_args + i)->H = H;
        (worker_args + i)->polyT2 = (gr_mpoly_struct *) flint_malloc(len*sizeof(gr_mpoly_struct));
        thread_pool_wake(global_thread_pool, handles[i], 0,
                                                 ideal_worker_loop, worker_args + i);
    }
    (worker_args + num_handles)->H = H;
    (worker_args + num_handles)->polyT2 = (gr_mpoly_struct *) flint_malloc(len*sizeof(gr_mpoly_struct));
    ideal_worker_loop(worker_args + num_handles);
    for (i = 0; i < num_handles; i++)
        thread_pool_wait(global_thread_pool, handles[i]);

    for (i = 0; i <= num_handles; i++)
        flint_free((worker_args + i)->polyT2);
    flint_free(worker_args);

#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&H->mutex);
#endif

    status = ideal_base_clear(Q, R, H);

    flint_free(H->polyQ);
    _gr_vec_clear(H->lc_inv, len, cctx);
    flint_free(H->lc_inv);
    flint_free(H->lc_is_one);
    flint_free(H->lc_is_unit);

cleanup1:

    fmpz_mpoly_clear(S, zctx);
    fmpz_mpoly_ctx_clear(zctx);

    flint_free(polyB_storage);
    flint_free(polyB_ptrs);

    if (freeAexp)
        flint_free(Aexp);
    for (w = 0; w < len; w++)
        if (freeBexp[w])
            flint_free(Bexp[w]);

    TMP_END;

    return status;
}


static int _gr_mpoly_divrem_ideal_heap_threaded_dispatch(
    gr_mpoly_struct * const * Q,
    gr_mpoly_t R,
    const gr_mpoly_t A,
    gr_mpoly_struct * const * B,
    slong len,
    gr_mpoly_ctx_t ctx,
    int nonfield)
{
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    thread_pool_handle * handles;
    slong num_handles;
    slong thread_limit;
    slong w;
    int status;

    for (w = 0; w < len; w++)
    {
        if (B[w]->length == 0)
            return GR_DOMAIN;
        if (gr_is_zero(B[w]->coeffs, cctx) != T_FALSE)
            return GR_UNABLE;
    }

    if (A->length == 0)
    {
        status = GR_SUCCESS;
        for (w = 0; w < len; w++)
            status |= gr_mpoly_zero(Q[w], ctx);
        status |= gr_mpoly_zero(R, ctx);
        return status;
    }

    /* fall back to the single-threaded algorithm for small inputs, or when
       the coefficient ring does not allow concurrent operations */
    if (A->length < 2 || B[0]->length < 2 || gr_ctx_is_threadsafe(cctx) != T_TRUE)
        return _gr_mpoly_divrem_ideal_serial(Q, R, A, B, len, nonfield, ctx);

    thread_limit = A->length/32;

    num_handles = flint_request_threads(&handles, thread_limit);

    status = _gr_mpoly_divrem_ideal_heap_threaded_pool(Q, R, A, B, len, ctx,
                                                nonfield, handles, num_handles);

    flint_give_back_threads(handles, num_handles);

    return status;
}


static int
_gr_mpoly_divrem_ideal_vec_threaded(gr_mpoly_vec_t Q, gr_mpoly_t R,
    const gr_mpoly_t A, const gr_mpoly_vec_t B, int nonfield, gr_mpoly_ctx_t ctx)
{
    slong w, len = B->length;
    gr_mpoly_struct ** Qptr;
    gr_mpoly_struct ** Bptr;
    int status;

    if (len == 0)
    {
        gr_mpoly_vec_set_length(Q, 0, ctx);
        return gr_mpoly_set(R, A, ctx);
    }

    gr_mpoly_vec_set_length(Q, len, ctx);

    Qptr = (gr_mpoly_struct **) flint_malloc(len*sizeof(gr_mpoly_struct *));
    Bptr = (gr_mpoly_struct **) flint_malloc(len*sizeof(gr_mpoly_struct *));
    for (w = 0; w < len; w++)
    {
        Qptr[w] = Q->entries + w;
        Bptr[w] = B->entries + w;
    }

    status = _gr_mpoly_divrem_ideal_heap_threaded_dispatch(Qptr, R, A, Bptr, len, ctx, nonfield);

    flint_free(Qptr);
    flint_free(Bptr);

    return status;
}


int gr_mpoly_divrem_ideal_heap_threaded(
    gr_mpoly_vec_t Q, gr_mpoly_t R,
    const gr_mpoly_t A, const gr_mpoly_vec_t B,
    gr_mpoly_ctx_t ctx)
{
    return _gr_mpoly_divrem_ideal_vec_threaded(Q, R, A, B, 0, ctx);
}


int gr_mpoly_divrem_ideal_weak_heap_threaded(
    gr_mpoly_vec_t Q, gr_mpoly_t R,
    const gr_mpoly_t A, const gr_mpoly_vec_t B,
    gr_mpoly_ctx_t ctx)
{
    return _gr_mpoly_divrem_ideal_vec_threaded(Q, R, A, B, 1, ctx);
}
