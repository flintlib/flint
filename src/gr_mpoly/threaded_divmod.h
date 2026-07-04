/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef GR_MPOLY_THREADED_DIVMOD_H
#define GR_MPOLY_THREADED_DIVMOD_H

#include "thread_pool.h"
#include "thread_support.h"
#include "gr_mpoly.h"

/*
    Infrastructure shared by gr_mpoly_divides_heap_threaded (in
    divides_heap_threaded.c) and gr_mpoly_divrem_heap_threaded /
    gr_mpoly_div_heap_threaded / gr_mpoly_div_weak_heap_threaded /
    gr_mpoly_divrem_weak_heap_threaded (in divrem_heap_threaded.c).

    Both algorithms partition the dividend A into exponent-range chunks and
    resolve them top-down through the same producer/consumer chunk pipeline,
    accumulating quotient terms (and, for the divrem family, remainder terms)
    into thread-safe growable arrays.  Only the *leaf* computation performed
    once a chunk becomes the producer differs:

        - gr_mpoly_divides_heap_threaded requires every surviving term of a
          finalized chunk to become a quotient term, and aborts (GR_DOMAIN)
          otherwise; see _gr_mpoly_divides_stripe1/_gr_mpoly_divides_stripe
          in divides_heap_threaded.c.

        - the divrem family instead always succeeds (barring GR_UNABLE /
          exponent overflow): terms that are not divisible by lm(B), or
          whose coefficient does not divide exactly (unless the "weak"
          Euclidean mode is requested), are routed to the remainder instead
          of aborting; see _gr_mpoly_divrem_stripe1/_gr_mpoly_divrem_stripe
          in divrem_heap_threaded.c.

    Everything else -- the thread-safe quotient/remainder arrays, the
    B-times-(quotient increment) subtraction used to keep each chunk's
    partial remainder up to date, the chunk struct and its lock-protected
    producer handoff, and the worker loop -- does not care which flavour of
    leaf is used and is defined once, in divides_heap_threaded.c, with
    linkage visible here so that divrem_heap_threaded.c can build its own
    top-level dispatch on top of it.
*/

/* shallow copy helpers (see mul_johnson.c, gr_mpoly/divides_heap.c) */
#define SET_SHALLOW_GENERIC(a,i,b,j) memcpy(GR_ENTRY(a, i, sz), GR_ENTRY(b, j, sz), sz)
#define SET_SHALLOW1(a,i,b,j) ((uint8_t *) a)[i] = ((uint8_t *) b)[j]
#define SET_SHALLOW2(a,i,b,j) ((uint16_t *) a)[i] = ((uint16_t *) b)[j]
#define SET_SHALLOW4(a,i,b,j) ((uint32_t *) a)[i] = ((uint32_t *) b)[j]

#if FLINT_BITS == 64
#define SET_SHALLOW8(a,i,b,j) ((uint64_t *) a)[i] = ((uint64_t *) b)[j]
#define SET_SHALLOW16(a,i,b,j) ((uint64_t *) a)[2 * (i)] = ((uint64_t *) b)[2 * (j)]; ((uint64_t *) a)[2 * (i) + 1] = ((uint64_t *) b)[2 * (j) + 1]
#else
#define SET_SHALLOW8(a,i,b,j) gr_set_shallow(GR_ENTRY(a, i, sz), GR_ENTRY(b, j, sz), cctx)
#define SET_SHALLOW16(a,i,b,j) gr_set_shallow(GR_ENTRY(a, i, sz), GR_ENTRY(b, j, sz), cctx)
#endif

/* gather product node (Bcoeff[i], Ccoeff[j]) into dot_a/dot_b (used by the
   mulsub stripe functions, where -- unlike the plain division heaps --
   there is no dividend node to special-case) */
#define MULSUB_FETCH_NODE(SET_SHALLOW) \
    do { \
        hind[x->i] |= WORD(1); \
        store[store_len++] = x->i; \
        store[store_len++] = x->j; \
        SET_SHALLOW(dot_a, dot_len, Bcoeff, x->i); \
        SET_SHALLOW(dot_b, dot_len, Ccoeff, x->j); \
        dot_len++; \
    } while (0)

/*
    A thread safe growable gr_mpoly.  Only three mutating operations are
    supported:
        - init from an array of terms
        - append an array of terms (the only operation performed while other
          threads may be concurrently reading)
        - clear out contents to a normal gr_mpoly_t

    Coefficients are moved by value (gr_swap) rather than copied, and gr
    elements are assumed to be relocatable by memcpy (the same assumption
    already made by the "shallow copy" helpers used throughout gr_mpoly),
    which lets us grow the coefficient array with a plain flint_malloc-based
    scheme just like the fmpz/nmod originals.

    IMPORTANT: when the backing array has to grow, *already published*
    entries must be migrated into the new buffer with a non-mutating copy
    (gr_set), never a swap/move -- other threads may still be concurrently
    reading through a stale pointer to the old buffer (which is
    deliberately never freed), so the old buffer's contents must not be
    disturbed. Only the newly incoming terms (which come from a private,
    unshared scratch buffer) may be moved by swap.
*/
typedef struct _gr_mpoly_ts_struct
{
    gr_ptr volatile coeffs; /* this is coeff_array[idx] */
    ulong * volatile exps;  /* this is exp_array[idx] */
    volatile slong length;
    slong alloc;
    slong N;
    flint_bitcnt_t bits;
    flint_bitcnt_t idx;
    ulong * exp_array[FLINT_BITS];
    gr_ptr coeff_array[FLINT_BITS];
} gr_mpoly_ts_struct;

typedef gr_mpoly_ts_struct gr_mpoly_ts_t[1];

void gr_mpoly_ts_init(gr_mpoly_ts_t A,
        gr_ptr Bcoeff, ulong * Bexp, slong Blen,
        flint_bitcnt_t bits, slong N, gr_ctx_t cctx);

void gr_mpoly_ts_clear(gr_mpoly_ts_t A, gr_ctx_t cctx);

void gr_mpoly_ts_clear_poly(gr_mpoly_t Q, gr_mpoly_ts_t A, gr_ctx_t cctx);

void gr_mpoly_ts_append(gr_mpoly_ts_t A,
                gr_ptr Bcoeff, ulong * Bexps, slong Blen, slong N, gr_ctx_t cctx);

/*
    A chunk holds an exponent range on the dividend.  Shared verbatim by
    both algorithms: chunk_mulsub (which reduces L->polyC by newly available
    quotient terms) has no notion of "divides" vs "divrem".
*/
typedef struct _divides_heap_chunk_struct
{
    gr_mpoly_t polyC;
    struct _divides_heap_chunk_struct * next;
    ulong * emin;
    ulong * emax;
    slong startidx;
    slong endidx;
    int upperclosed;
    volatile int lock;
    volatile int producer;
    volatile slong ma;
    volatile slong mq;
    int Cinited;
} divides_heap_chunk_struct;

typedef divides_heap_chunk_struct divides_heap_chunk_t[1];

/*
    The base struct includes a linked list of chunks together with the
    thread safe growable quotient (and, for the divrem family, remainder).

    require_exact / want_remainder / nonfield select which algorithm
    trychunk's producer-finalisation step runs:

        require_exact = 1  ->  gr_mpoly_divides_heap_threaded semantics:
                                a chunk with leftover, non-quotient terms is
                                a hard failure (have_domain).  want_remainder
                                and nonfield are unused; polyR is untouched.

        require_exact = 0  ->  divrem family: leftover terms are routed to
                                polyR (if want_remainder) or simply discarded
                                (if !want_remainder, i.e. gr_mpoly_div(_weak)
                                _heap_threaded).  nonfield selects exact
                                (gr_div) vs Euclidean (gr_euclidean_divrem)
                                coefficient handling, exactly like the
                                serial gr_mpoly_divrem_heap kernel.  unchecked
                                (only meaningful when nonfield == 0) instead
                                selects gr_divexact for a non-unit lc(B),
                                i.e. gr_mpoly_divexact_heap_threaded; leftover
                                terms are always discarded in that case
                                (want_remainder == 0), matching the fact that
                                gr_mpoly_divexact has no remainder-producing
                                counterpart.
*/
typedef struct
{
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif
    divides_heap_chunk_struct * head;
    divides_heap_chunk_struct * tail;
    divides_heap_chunk_struct * volatile cur;
    gr_mpoly_t polyA;
    gr_mpoly_t polyB;
    gr_mpoly_ts_t polyQ;
    gr_mpoly_ts_t polyR;       /* only used/initialised if want_remainder */
    gr_mpoly_ctx_struct * ctx;
    gr_ctx_struct * cctx;
    slong length;
    slong N;
    flint_bitcnt_t bits;
    ulong * cmpmask;
    int lc_is_one;
    int lc_is_unit;
    gr_srcptr lc_inv;    /* 1/lc(B), valid iff lc_is_unit */
    int require_exact;   /* 1: gr_mpoly_divides semantics; 0: div/divrem */
    int want_remainder;  /* 0: gr_mpoly_div(_weak); 1: gr_mpoly_divrem(_weak) */
    int nonfield;         /* 1: Euclidean ("weak") per-term coefficient divrem */
    int unchecked;        /* 1: gr_mpoly_divexact -- gr_divexact instead of
                              gr_div for a non-unit lc(B); mutually exclusive
                              with nonfield (both require_exact == 0) */
    volatile int failed;      /* stop: either not exact, or arithmetic unable */
    volatile int have_domain; /* division is provably not exact (require_exact only) */
    volatile int have_unable; /* some ring operation returned GR_UNABLE */
    volatile int overflowed;  /* an internal exponent computation overflowed
                                  the current bit width (see the comment on
                                  _gr_mpoly_mulsub_stripe1 above): the whole
                                  computation must be redone at a wider bit
                                  width, it cannot be salvaged in place.
                                  Setting this also sets failed, to stop all
                                  workers immediately, exactly like
                                  have_domain/have_unable. Whether this is
                                  acted on (retried) or just folded into an
                                  UNABLE verdict depends on the caller: the
                                  require_exact (gr_mpoly_divides) pool
                                  function does the latter (unchanged
                                  behaviour, now merely made safe rather than
                                  silently undetected -- see the long comment
                                  at _gr_mpoly_mulsub_stripe1), while the
                                  divrem/divrem_ideal pool functions retry. */
} divides_heap_base_struct;

typedef divides_heap_base_struct divides_heap_base_t[1];

void divides_heap_base_init(divides_heap_base_t H);
void divides_heap_chunk_clear(divides_heap_chunk_t L, divides_heap_base_t H);

/*
    Finalise H into (Q, R): R may be NULL (matches H->want_remainder == 0).
    Returns the overall GR status.  Always clears/frees H's chunk list and
    thread-safe arrays.
*/
int divides_heap_base_clear(gr_mpoly_t Q, gr_mpoly_t R, divides_heap_base_t H);

void divides_heap_base_add_chunk(divides_heap_base_t H, divides_heap_chunk_t L);

/*
    Binary search for the first index >= a (in the caller-supplied,
    monomial-decreasing array Aexp of length Alen) whose exponent no longer
    exceeds exp. Generic over which gr_mpoly's exponent array is searched,
    so it does not need to know about divides_heap_base_t.
*/
slong chunk_find_exp(ulong * exp, slong a,
                     const ulong * Aexp, slong Alen, slong N, const ulong * cmpmask);

/*
    Per-worker scratch memory (analogous to fmpz_mpoly_stripe_t / the
    _gr_mpoly_stripe_struct of gr_mpoly_mul_heap_threaded, extended with the
    exponent-window and leading-coefficient data needed for division).
*/
typedef struct
{
    char * big_mem;
    slong big_mem_alloc;
    slong N;
    flint_bitcnt_t bits;
    const ulong * cmpmask;
    slong * startidx;
    slong * endidx;
    ulong * emin;
    ulong * emax;
    int upperclosed;
    gr_mpoly_ctx_struct * ctx;
    gr_ctx_struct * cctx;
    slong sz;
    int have_fast_dot;
    int lc_is_one;
    int lc_is_unit;
    gr_srcptr lc_inv;
    int unchecked;  /* only consulted by the divrem-family stripe leaves,
                       and only meaningful together with nonfield == 0 */
    int nonfield;   /* only consulted by the divrem-family stripe leaves */
} _gr_mpoly_stripe_struct;

typedef _gr_mpoly_stripe_struct _gr_mpoly_stripe_t[1];

void stripe_fit_length(_gr_mpoly_stripe_struct * S, slong new_len);

/* A = D - (a stripe of B * C); see the long comment at its definition in
   divides_heap_threaded.c. Shared verbatim: has no notion of "divides" vs
   "divrem", it just forms a difference of products.

   *overflowed is set to 1 (and *res_status/the returned length should be
   ignored -- they are not meaningful) if an exponent computed internally
   (e.g. as a sum of a term of B and a term of C) does not fit in the
   current bit width. This can happen even when every term of D, B and C
   individually fits, since two individually-valid packed exponents can
   overflow a shared field when added; the caller must repack to a wider
   bit width and redo the whole computation -- there is no way to recover
   in place. See the retry loops in gr_mpoly_divrem_heap_threaded /
   gr_mpoly_divrem_ideal_heap_threaded. */
slong _gr_mpoly_mulsub_stripe1(
    gr_ptr * A_coeff, ulong ** A_exp, slong * A_alloc, slong * A_exps_alloc,
    gr_srcptr Dcoeff, const ulong * Dexp, slong Dlen, int saveD,
    gr_srcptr Bcoeff, const ulong * Bexp, slong Blen,
    gr_srcptr Ccoeff, const ulong * Cexp, slong Clen,
    const _gr_mpoly_stripe_t S, int * res_status, int * overflowed);

slong _gr_mpoly_mulsub_stripe(
    gr_ptr * A_coeff, ulong ** A_exp, slong * A_alloc, slong * A_exps_alloc,
    gr_srcptr Dcoeff, const ulong * Dexp, slong Dlen, int saveD,
    gr_srcptr Bcoeff, const ulong * Bexp, slong Blen,
    gr_srcptr Ccoeff, const ulong * Cexp, slong Clen,
    const _gr_mpoly_stripe_t S, int * res_status, int * overflowed);

/*
    Bounded ("stripe") division-with-remainder leaves, defined in
    divrem_heap_threaded.c, analogous to _gr_mpoly_divides_stripe1/
    _gr_mpoly_divides_stripe in divides_heap_threaded.c but never aborting
    on a non-reducible term: it is routed to the local remainder output
    (Wcoeff/Wexp/Wlen) instead. *res_status only ever reports GR_SUCCESS or
    GR_UNABLE (never GR_DOMAIN -- that concept does not apply once
    require_exact == 0).

    *overflowed has the same meaning as for _gr_mpoly_mulsub_stripe1 above:
    an internal exponent computation overflowed the current bit width, and
    *res_status / the returned Qlen / *Wlen_out are not meaningful.
*/
slong _gr_mpoly_divrem_stripe1(
    gr_ptr * Q_coeff, ulong ** Q_exp, slong * Q_alloc, slong * Q_exps_alloc,
    gr_ptr * W_coeff, ulong ** W_exp, slong * W_alloc, slong * W_exps_alloc,
    slong * Wlen_out,
    gr_srcptr Acoeff, const ulong * Aexp, slong Alen,
    gr_srcptr Bcoeff, const ulong * Bexp, slong Blen,
    const _gr_mpoly_stripe_t S, int * res_status, int * overflowed);

slong _gr_mpoly_divrem_stripe(
    gr_ptr * Q_coeff, ulong ** Q_exp, slong * Q_alloc, slong * Q_exps_alloc,
    gr_ptr * W_coeff, ulong ** W_exp, slong * W_alloc, slong * W_exps_alloc,
    slong * Wlen_out,
    gr_srcptr Acoeff, const ulong * Aexp, slong Alen,
    gr_srcptr Bcoeff, const ulong * Bexp, slong Blen,
    const _gr_mpoly_stripe_t S, int * res_status, int * overflowed);

/*
    The worker struct has the stripe scratch memory and three polys for
    workspace: T1/T2 as before (mulsub output / newly finalised quotient
    increment), T3 an additional scratch poly used to hold a newly
    finalised *remainder* increment when H->want_remainder.
*/
typedef struct _worker_arg_struct
{
    divides_heap_base_struct * H;
    _gr_mpoly_stripe_t S;
    gr_mpoly_t polyT1;
    gr_mpoly_t polyT2;
    gr_mpoly_t polyT3;
} worker_arg_struct;

typedef worker_arg_struct worker_arg_t[1];

/* the chunk pipeline itself: defined in divides_heap_threaded.c, mode-aware
   via H->require_exact / H->want_remainder / H->nonfield */
void trychunk(worker_arg_t W, divides_heap_chunk_t L);
void worker_loop(void * varg);
void chunk_mulsub(worker_arg_t W, divides_heap_chunk_t L, slong q_prev_length);

/*
    The guaranteed-serial kernel for gr_mpoly_div(_weak)/gr_mpoly_divrem(_weak),
    defined in divrem_heap.c (R == NULL for the div/div_weak case). Declared
    here so that the threaded engine in divrem_heap_threaded.c can fall back
    to it directly, rather than through the public dispatchers (gr_mpoly_div,
    gr_mpoly_divrem, ...), which may route large instances straight back into
    the threaded engine.
*/
int _gr_mpoly_divrem_mp(
    gr_mpoly_t Q, gr_mpoly_t R,
    const gr_mpoly_t A, const gr_mpoly_t B, int nonfield, int unchecked,
    gr_mpoly_ctx_t ctx);

/*
    The guaranteed-serial divrem_ideal kernel, defined in divrem_ideal.c.
    Declared here so that divrem_ideal_heap_threaded.c can fall back to it
    directly, for the same reason _gr_mpoly_divrem_mp is exposed above.
*/
int _gr_mpoly_divrem_ideal(
    gr_mpoly_struct ** Q, gr_mpoly_t R,
    const gr_mpoly_t A, gr_mpoly_struct * const * B, slong len,
    int nonfield, gr_mpoly_ctx_t ctx);

#endif
