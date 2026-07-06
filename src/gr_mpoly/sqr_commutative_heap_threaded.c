/*
    Copyright (C) 2017-2019 Daniel Schultz
    Copyright (C) 2022, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdint.h>
#include <string.h>
#include "longlong.h"
#include "thread_pool.h"
#include "thread_support.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "mpoly.h"
#include "gr_generic.h"
#include "gr_mpoly.h"

/*
    Multithreaded squaring of a sparse multivariate polynomial, assuming a
    commutative coefficient ring.  It is modeled on gr_mpoly_mul_heap_threaded
    (which parallelizes gr_mpoly_mul_heap) and reuses that routine's
    partitioning, heap traversal and join machinery verbatim.  The only change
    is the per-monomial coefficient accumulation.

    Writing f = sum_i a_i X^{e_i} (with the e_i distinct),

        f^2 = sum_i a_i^2 X^{2 e_i}  +  2 sum_{i<j} a_i a_j X^{e_i + e_j}.

    The product f*f is partitioned by output monomial into disjoint, already
    sorted chunks exactly as in the multiplication (each chunk is a Johnson
    heap over the len x len grid restricted to a column range [start, end) per
    row).  Because identical output monomials always fall in a single chunk,
    both orientations (i,j) and (j,i) of every cross pair are enumerated inside
    the same diagonal batch.  We therefore accumulate only the terms with
    i <= j: the off-diagonal products (i < j) are summed and doubled, and the
    diagonal square a_d^2 is added whenever a diagonal pair (d,d) shares the
    monomial.  The terms with i > j are still popped and cascaded through the
    heap (so the traversal is byte-for-byte identical to the multiplication),
    but their coefficient product is skipped.  This performs only about
    len*(len+1)/2 coefficient multiplications, half of the len^2 of a plain
    squaring, while the partitioning and correctness follow immediately from
    the multiplication together with a_i a_j = a_j a_i.
*/

/* Helpers to place operands in two temporary shallow arrays so that we can
   accumulate the off-diagonal part with a single _gr_vec_dot call.  (Same as
   in mul_heap_threaded.c.) */

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

/* Fast-dot fetch: pop the whole chain of nodes sharing the current exponent.
   Every node updates hind and Q (to drive the heap cascade exactly as in the
   multiplication), but only i < j nodes go into the dot vectors and only the
   unique i == j node is recorded as the diagonal.  i > j nodes are skipped. */

#define SQR_FETCH_TERMS_INNER(SET_SHALLOW) \
    do { \
        slong _ii = x->i, _jj = x->j; \
        hind[_ii] |= WORD(1); \
        Q[Q_len++] = _ii; \
        Q[Q_len++] = _jj; \
        if (_ii < _jj) \
        { \
            SET_SHALLOW(dot_a, dot_len, coeff2, _ii); \
            SET_SHALLOW(dot_b, dot_len, coeff2, _jj); \
            dot_len++; \
        } \
        else if (_ii == _jj) \
        { \
            diag_i = _ii; \
            have_diag = 1; \
        } \
    } while (0); \
    while ((x = x->next) != NULL) \
    { \
        slong _ii = x->i, _jj = x->j; \
        hind[_ii] |= WORD(1); \
        Q[Q_len++] = _ii; \
        Q[Q_len++] = _jj; \
        if (_ii < _jj) \
        { \
            SET_SHALLOW(dot_a, dot_len, coeff2, _ii); \
            SET_SHALLOW(dot_b, dot_len, coeff2, _jj); \
            dot_len++; \
        } \
        else if (_ii == _jj) \
        { \
            diag_i = _ii; \
            have_diag = 1; \
        } \
    }

#define SQR_FETCH_TERMS1(SET_SHALLOW) \
    do { \
        while (heap_len > 1 && heap[1].exp == exp) \
        { \
            x = _mpoly_heap_pop1(heap, &heap_len, maskhi); \
            SQR_FETCH_TERMS_INNER(SET_SHALLOW) \
        } \
    } while (0)

#define SQR_FETCH_TERMS(SET_SHALLOW) \
    do \
    { \
        exp_list[--exp_next] = heap[1].exp; \
        x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask); \
        SQR_FETCH_TERMS_INNER(SET_SHALLOW) \
    } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N))

/* Combine the accumulated off-diagonal dot product D with the optional
   diagonal square:  result = 2*D + a_diag^2. */

#define SQR_FINALIZE_DOT \
    do { \
        gr_ptr _out = GR_ENTRY(p1, len1, sz); \
        if (dot_len > 0) \
        { \
            status |= _gr_vec_dot(_out, NULL, 0, dot_a, dot_b, dot_len, cctx); \
            status |= gr_mul_two(_out, _out, cctx); \
            if (have_diag) \
            { \
                status |= gr_sqr(t, GR_ENTRY(coeff2, diag_i, sz), cctx); \
                status |= gr_add(_out, _out, t, cctx); \
            } \
        } \
        else \
        { \
            /* pure diagonal monomial 2*e_diag */ \
            status |= gr_sqr(_out, GR_ENTRY(coeff2, diag_i, sz), cctx); \
        } \
    } while (0)

/* Generic accumulation (no fast dot product): off-diagonal products (i < j)
   are summed with gr_mul / gr_addmul, doubled once, and the diagonal square is
   added.  i > j nodes update hind/Q but are otherwise skipped. */

#define SQR_DOT_TERMS_GENERIC \
    do { \
        do { \
            slong _ii = x->i, _jj = x->j; \
            hind[_ii] |= WORD(1); \
            Q[Q_len++] = _ii; \
            Q[Q_len++] = _jj; \
            if (_ii < _jj) \
            { \
                if (first_off) \
                { \
                    status |= gr_mul(GR_ENTRY(p1, len1, sz), GR_ENTRY(coeff2, _ii, sz), GR_ENTRY(coeff2, _jj, sz), cctx); \
                    first_off = 0; \
                } \
                else \
                { \
                    status |= gr_addmul(GR_ENTRY(p1, len1, sz), GR_ENTRY(coeff2, _ii, sz), GR_ENTRY(coeff2, _jj, sz), cctx); \
                } \
            } \
            else if (_ii == _jj) \
            { \
                status |= gr_sqr(t, GR_ENTRY(coeff2, _ii, sz), cctx); \
                have_diag = 1; \
            } \
        } while ((x = x->next) != NULL); \
    } while (0)

#define SQR_FINALIZE_GENERIC \
    do { \
        gr_ptr _out = GR_ENTRY(p1, len1, sz); \
        if (!first_off) \
        { \
            status |= gr_mul_two(_out, _out, cctx); \
            if (have_diag) \
                status |= gr_add(_out, _out, t, cctx); \
        } \
        else \
        { \
            /* only a diagonal term contributed */ \
            status |= gr_set(_out, t, cctx); \
        } \
    } while (0)


/* Per-worker scratch + ring/order parameters shared by the part functions. */
typedef struct
{
    gr_mpoly_ctx_struct * ctx;
    gr_ctx_struct * cctx;
    slong sz;
    int have_fast_dot;
    flint_bitcnt_t bits;
    slong N;
    const ulong * cmpmask;
    char * big_mem;
    slong big_mem_alloc;
}
_gr_mpoly_stripe_struct;

typedef _gr_mpoly_stripe_struct _gr_mpoly_stripe_t[1];


/* Compute the squaring terms whose output monomial lies in the range covered
   by [start, end) of the (single) polynomial coeff2/exp2/len2, single-word
   exponent version.  Both operand roles read from coeff2/exp2 since f is being
   squared.  Returns the number of terms produced and writes the accumulated gr
   status to *res_status. */
static slong _gr_mpoly_sqr_commutative_heap_part1(
    gr_ptr * A_coeff, ulong ** A_exp, slong * A_alloc, slong * A_exps_alloc,
    gr_srcptr coeff2, const ulong * exp2, slong len2,
    slong * start, slong * end, slong * hind,
    const _gr_mpoly_stripe_t S, int * res_status)
{
    gr_mpoly_ctx_struct * ctx = S->ctx;
    gr_ctx_struct * cctx = S->cctx;
    slong sz = S->sz;
    int have_fast_dot = S->have_fast_dot;
    const ulong maskhi = S->cmpmask[0];
    slong i, j;
    ulong exp;
    mpoly_heap_t * x;
    slong next_loc;
    slong heap_len;
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * Q;
    slong Q_len;
    slong len1;
    gr_ptr p1 = *A_coeff;
    ulong * e1 = *A_exp;
    gr_ptr t, dot_a, dot_b;
    slong dot_len;
    int first_off;
    int have_diag;
    slong diag_i;
    int status = GR_SUCCESS;

    i = 0;
    Q = (slong *) (S->big_mem + i);
    i += 2*len2*sizeof(slong);
    heap = (mpoly_heap1_s *) (S->big_mem + i);
    i += (len2 + 1)*sizeof(mpoly_heap1_s);
    chain = (mpoly_heap_t *) (S->big_mem + i);
    i += len2*sizeof(mpoly_heap_t);
    FLINT_ASSERT(i <= S->big_mem_alloc);

    GR_TMP_INIT(t, cctx);
    dot_a = flint_malloc(2 * len2 * sz);
    dot_b = GR_ENTRY(dot_a, len2, sz);

    /* put all the starting nodes on the heap */
    Q_len = 0;
    heap_len = 1; /* heap zero index unused */
    next_loc = len2 + 4;   /* something bigger than heap can ever be */
    for (i = 0; i < len2; i++)
        hind[i] = 2*start[i] + 1;
    for (i = 0; i < len2; i++)
    {
        if ((start[i] < end[i]) && ((i == 0) || (start[i] < start[i - 1])))
        {
            x = chain + i;
            x->i = i;
            x->j = start[i];
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;
            _mpoly_heap_insert1(heap, exp2[x->i] + exp2[x->j], x,
                                                &next_loc, &heap_len, maskhi);
        }
    }

    len1 = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        _gr_mpoly_fit_length(&p1, A_alloc, &e1, A_exps_alloc, 1, len1 + 1, ctx);

        e1[len1] = exp;

        have_diag = 0;
        diag_i = 0;

        if (!have_fast_dot)
        {
            first_off = 1;

            while (heap_len > 1 && heap[1].exp == exp)
            {
                x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
                SQR_DOT_TERMS_GENERIC;
            }

            SQR_FINALIZE_GENERIC;
        }
        else
        {
            dot_len = 0;

            if (sz == 1)
            {
                SQR_FETCH_TERMS1(SET_SHALLOW1);
            }
            else if (sz == 2)
            {
                SQR_FETCH_TERMS1(SET_SHALLOW2);
            }
            else if (sz == 4)
            {
                SQR_FETCH_TERMS1(SET_SHALLOW4);
            }
            else if (sz == 8)
            {
                SQR_FETCH_TERMS1(SET_SHALLOW8);
            }
            else if (sz == 16)
            {
                SQR_FETCH_TERMS1(SET_SHALLOW16);
            }
            else
            {
                SQR_FETCH_TERMS1(SET_SHALLOW_GENERIC);
            }

            SQR_FINALIZE_DOT;
        }

        len1 += (gr_is_zero(GR_ENTRY(p1, len1, sz), cctx) != T_TRUE);

        /* for each node temporarily stored */
        while (Q_len > 0)
        {
            j = Q[--Q_len];
            i = Q[--Q_len];

            /* should we go right? */
            if ((i + 1 < len2) && (j + 0 < end[i + 1]) && (hind[i + 1] == 2*j + 1))
            {
                x = chain + i + 1;
                x->i = i + 1;
                x->j = j;
                x->next = NULL;
                hind[x->i] = 2*(x->j + 1) + 0;
                _mpoly_heap_insert1(heap, exp2[x->i] + exp2[x->j], x,
                                                &next_loc, &heap_len, maskhi);
            }

            /* should we go up? */
            if ((j + 1 < end[i + 0]) && ((hind[i] & 1) == 1) &&
                ((i == 0) || (hind[i - 1] >= 2*(j + 2) + 1)))
            {
                x = chain + i;
                x->i = i;
                x->j = j + 1;
                x->next = NULL;
                hind[x->i] = 2*(x->j + 1) + 0;
                _mpoly_heap_insert1(heap, exp2[x->i] + exp2[x->j], x,
                                                &next_loc, &heap_len, maskhi);
            }
        }
    }

    flint_free(dot_a);
    GR_TMP_CLEAR(t, cctx);

    *A_coeff = p1;
    *A_exp = e1;
    *res_status = status;
    return len1;
}


/* Multi-word exponent version of the above. */
static slong _gr_mpoly_sqr_commutative_heap_part(
    gr_ptr * A_coeff, ulong ** A_exp, slong * A_alloc, slong * A_exps_alloc,
    gr_srcptr coeff2, const ulong * exp2, slong len2,
    slong * start, slong * end, slong * hind,
    const _gr_mpoly_stripe_t S, int * res_status)
{
    gr_mpoly_ctx_struct * ctx = S->ctx;
    gr_ctx_struct * cctx = S->cctx;
    slong sz = S->sz;
    int have_fast_dot = S->have_fast_dot;
    flint_bitcnt_t bits = S->bits;
    slong N = S->N;
    const ulong * cmpmask = S->cmpmask;
    slong i, j;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    mpoly_heap_t * x;
    slong next_loc;
    slong heap_len;
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * Q;
    slong Q_len;
    slong len1;
    gr_ptr p1 = *A_coeff;
    ulong * e1 = *A_exp;
    gr_ptr t, dot_a, dot_b;
    slong dot_len;
    int first_off;
    int have_diag;
    slong diag_i;
    int status = GR_SUCCESS;

    i = 0;
    Q = (slong *) (S->big_mem + i);
    i += 2*len2*sizeof(slong);
    exp_list = (ulong **) (S->big_mem + i);
    i += len2*sizeof(ulong *);
    exps = (ulong *) (S->big_mem + i);
    i += len2*N*sizeof(ulong);
    heap = (mpoly_heap_s *) (S->big_mem + i);
    i += (len2 + 1)*sizeof(mpoly_heap_s);
    chain = (mpoly_heap_t *) (S->big_mem + i);
    i += len2*sizeof(mpoly_heap_t);
    FLINT_ASSERT(i <= S->big_mem_alloc);

    GR_TMP_INIT(t, cctx);
    dot_a = flint_malloc(2 * len2 * sz);
    dot_b = GR_ENTRY(dot_a, len2, sz);

    /* put all the starting nodes on the heap */
    Q_len = 0;
    heap_len = 1; /* heap zero index unused */
    next_loc = len2 + 4;   /* something bigger than heap can ever be */
    exp_next = 0;
    for (i = 0; i < len2; i++)
        exp_list[i] = exps + i*N;
    for (i = 0; i < len2; i++)
        hind[i] = 2*start[i] + 1;
    for (i = 0; i < len2; i++)
    {
        if ((start[i] < end[i]) && ((i == 0) || (start[i] < start[i - 1])))
        {
            x = chain + i;
            x->i = i;
            x->j = start[i];
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;
            mpoly_monomial_add_any_bits(exp_list[exp_next], exp2 + x->i*N, exp2 + x->j*N, N, bits);

            exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
        }
    }

    len1 = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        _gr_mpoly_fit_length(&p1, A_alloc, &e1, A_exps_alloc, N, len1 + 1, ctx);

        mpoly_monomial_set(e1 + len1*N, exp, N);

        have_diag = 0;
        diag_i = 0;

        if (!have_fast_dot)
        {
            first_off = 1;

            do
            {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
                SQR_DOT_TERMS_GENERIC;
            }
            while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

            SQR_FINALIZE_GENERIC;
        }
        else
        {
            dot_len = 0;

            if (sz == 1)
            {
                SQR_FETCH_TERMS(SET_SHALLOW1);
            }
            else if (sz == 2)
            {
                SQR_FETCH_TERMS(SET_SHALLOW2);
            }
            else if (sz == 4)
            {
                SQR_FETCH_TERMS(SET_SHALLOW4);
            }
            else if (sz == 8)
            {
                SQR_FETCH_TERMS(SET_SHALLOW8);
            }
            else if (sz == 16)
            {
                SQR_FETCH_TERMS(SET_SHALLOW16);
            }
            else
            {
                SQR_FETCH_TERMS(SET_SHALLOW_GENERIC);
            }

            SQR_FINALIZE_DOT;
        }

        len1 += (gr_is_zero(GR_ENTRY(p1, len1, sz), cctx) != T_TRUE);

        /* for each node temporarily stored */
        while (Q_len > 0)
        {
            j = Q[--Q_len];
            i = Q[--Q_len];

            /* should we go right? */
            if ((i + 1 < len2) && (j + 0 < end[i + 1]) && (hind[i + 1] == 2*j + 1))
            {
                x = chain + i + 1;
                x->i = i + 1;
                x->j = j;
                x->next = NULL;
                hind[x->i] = 2*(x->j + 1) + 0;

                mpoly_monomial_add_any_bits(exp_list[exp_next], exp2 + x->i*N, exp2 + x->j*N, N, bits);

                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
            }

            /* should we go up? */
            if ((j + 1 < end[i + 0]) && ((hind[i] & 1) == 1) &&
                ((i == 0) || (hind[i - 1] >= 2*(j + 2) + 1)))
            {
                x = chain + i;
                x->i = i;
                x->j = j + 1;
                x->next = NULL;
                hind[x->i] = 2*(x->j + 1) + 0;

                mpoly_monomial_add_any_bits(exp_list[exp_next], exp2 + x->i*N, exp2 + x->j*N, N, bits);

                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
            }
        }
    }

    flint_free(dot_a);
    GR_TMP_CLEAR(t, cctx);

    *A_coeff = p1;
    *A_exp = e1;
    *res_status = status;
    return len1;
}


/*
    The workers calculate squaring terms from 4*n divisions, where n is the
    number of threads.
*/

typedef struct
{
    volatile int idx;
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif
    slong nthreads;
    slong ndivs;
    gr_mpoly_ctx_struct * ctx;
    gr_ptr Acoeff;
    ulong * Aexp;
    gr_srcptr coeff2;
    const ulong * exp2;
    slong len2;
    slong N;
    flint_bitcnt_t bits;
    const ulong * cmpmask;
    int have_fast_dot;
}
_base_struct;

typedef _base_struct _base_t[1];

typedef struct
{
    slong lower;
    slong upper;
    slong thread_idx;
    slong Aoffset;
    slong Alen;
    slong Aalloc;
    slong Aexps_alloc;
    ulong * Aexp;
    gr_ptr Acoeff;
    int status;
}
_div_struct;

typedef struct
{
    _gr_mpoly_stripe_t S;
    slong idx;
    _base_struct * base;
    _div_struct * divs;
    ulong * exp;
}
_worker_arg_struct;


#define SWAP_PTRS(xx, yy) \
   do { \
      tt = xx; \
      xx = yy; \
      yy = tt; \
   } while (0)

static void _gr_mpoly_sqr_commutative_heap_threaded_worker(void * varg)
{
    _worker_arg_struct * arg = (_worker_arg_struct *) varg;
    _gr_mpoly_stripe_struct * S = arg->S;
    _div_struct * divs = arg->divs;
    _base_struct * base = arg->base;
    slong len2 = base->len2;
    slong N = base->N;
    slong i, j;
    ulong * exp;
    slong score;
    slong * start, * end, * t1, * t2, * t3, * t4, * tt;

    exp = (ulong *) flint_malloc(N*sizeof(ulong));
    t1 = (slong *) flint_malloc(len2*sizeof(slong));
    t2 = (slong *) flint_malloc(len2*sizeof(slong));
    t3 = (slong *) flint_malloc(len2*sizeof(slong));
    t4 = (slong *) flint_malloc(len2*sizeof(slong));

    S->ctx = base->ctx;
    S->cctx = GR_MPOLY_CCTX(base->ctx);
    S->sz = S->cctx->sizeof_elem;
    S->have_fast_dot = base->have_fast_dot;
    S->bits = base->bits;
    S->N = N;
    S->cmpmask = base->cmpmask;

    S->big_mem_alloc = 0;
    if (N == 1)
    {
        S->big_mem_alloc += 2*len2*sizeof(slong);
        S->big_mem_alloc += (len2 + 1)*sizeof(mpoly_heap1_s);
        S->big_mem_alloc += len2*sizeof(mpoly_heap_t);
    }
    else
    {
        S->big_mem_alloc += 2*len2*sizeof(slong);
        S->big_mem_alloc += len2*sizeof(ulong *);
        S->big_mem_alloc += len2*N*sizeof(ulong);
        S->big_mem_alloc += (len2 + 1)*sizeof(mpoly_heap_s);
        S->big_mem_alloc += len2*sizeof(mpoly_heap_t);
    }
    S->big_mem = (char *) flint_malloc(S->big_mem_alloc);

    /* get index to start working on */
    if (arg->idx + 1 < base->nthreads)
    {
#if FLINT_USES_PTHREAD
        pthread_mutex_lock(&base->mutex);
#endif
        i = base->idx - 1;
        base->idx = i;
#if FLINT_USES_PTHREAD
        pthread_mutex_unlock(&base->mutex);
#endif
    }
    else
    {
        i = base->ndivs - 1;
    }

    while (i >= 0)
    {
        FLINT_ASSERT(divs[i].thread_idx == -WORD(1));
        divs[i].thread_idx = arg->idx;

        /* calculate start */
        if (i + 1 < base->ndivs)
        {
            mpoly_search_monomials(
                &start, exp, &score, t1, t2, t3,
                            divs[i].lower, divs[i].lower,
                            base->exp2, base->len2, base->exp2, base->len2,
                                          base->N, base->cmpmask);
            if (start == t2)
            {
                SWAP_PTRS(t1, t2);
            }
            else if (start == t3)
            {
                SWAP_PTRS(t1, t3);
            }
        }
        else
        {
            start = t1;
            for (j = 0; j < base->len2; j++)
                start[j] = 0;
        }

        /* calculate end */
        if (i > 0)
        {
            mpoly_search_monomials(
                &end, exp, &score, t2, t3, t4,
                            divs[i - 1].lower, divs[i - 1].lower,
                            base->exp2, base->len2, base->exp2, base->len2,
                                          base->N, base->cmpmask);
            if (end == t3)
            {
                SWAP_PTRS(t2, t3);
            }
            else if (end == t4)
            {
                SWAP_PTRS(t2, t4);
            }
        }
        else
        {
            end = t2;
            for (j = 0; j < base->len2; j++)
                end[j] = base->len2;
        }
        /* t3 and t4 are free for workspace at this point */

        /* join code assumes all divisions have been allocated */
        _gr_mpoly_fit_length(&divs[i].Acoeff, &divs[i].Aalloc,
                             &divs[i].Aexp, &divs[i].Aexps_alloc, N, 256, base->ctx);

        /* calculate products in [start, end) */
        if (N == 1)
        {
            divs[i].Alen = _gr_mpoly_sqr_commutative_heap_part1(
                         &divs[i].Acoeff, &divs[i].Aexp,
                         &divs[i].Aalloc, &divs[i].Aexps_alloc,
                              base->coeff2, base->exp2, base->len2,
                                          start, end, t3, S, &divs[i].status);
        }
        else
        {
            divs[i].Alen = _gr_mpoly_sqr_commutative_heap_part(
                         &divs[i].Acoeff, &divs[i].Aexp,
                         &divs[i].Aalloc, &divs[i].Aexps_alloc,
                              base->coeff2, base->exp2, base->len2,
                                          start, end, t3, S, &divs[i].status);
        }

        /* get next index to work on */
#if FLINT_USES_PTHREAD
        pthread_mutex_lock(&base->mutex);
#endif
        i = base->idx - 1;
        base->idx = i;
#if FLINT_USES_PTHREAD
        pthread_mutex_unlock(&base->mutex);
#endif
    }

    flint_free(S->big_mem);
    flint_free(t4);
    flint_free(t3);
    flint_free(t2);
    flint_free(t1);
    flint_free(exp);
}

static void _join_worker(void * varg)
{
    _worker_arg_struct * arg = (_worker_arg_struct *) varg;
    _div_struct * divs = arg->divs;
    _base_struct * base = arg->base;
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(base->ctx);
    slong sz = cctx->sizeof_elem;
    slong N = base->N;
    slong i;

    for (i = base->ndivs - 2; i >= 0; i--)
    {
        FLINT_ASSERT(divs[i].thread_idx != -WORD(1));

        if (divs[i].thread_idx != arg->idx)
            continue;

        FLINT_ASSERT(divs[i].Acoeff != NULL);
        FLINT_ASSERT(divs[i].Aexp != NULL);

        /* Move this division's coefficients into the output buffer.  The
           destination slots have been initialised (to zero) by the main
           thread, so we exchange rather than overwrite to avoid leaking and
           leave clearable elements behind in the worker buffer. */
        _gr_vec_swap(GR_ENTRY(base->Acoeff, divs[i].Aoffset, sz),
                     divs[i].Acoeff, divs[i].Alen, cctx);

        memcpy(base->Aexp + N*divs[i].Aoffset, divs[i].Aexp,
                                                 N*divs[i].Alen*sizeof(ulong));

        _gr_vec_clear(divs[i].Acoeff, divs[i].Aalloc, cctx);
        flint_free(divs[i].Acoeff);
        flint_free(divs[i].Aexp);
    }
}

static int _gr_mpoly_sqr_commutative_heap_threaded(
    gr_mpoly_t A,
    gr_srcptr coeff2, const ulong * exp2, slong len2,
    flint_bitcnt_t bits,
    slong N,
    const ulong * cmpmask,
    gr_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    slong i;
    slong BClen, hi;
    _base_t base;
    _div_struct * divs;
    _worker_arg_struct * args;
    slong Aalloc, Aexps_alloc;
    slong Alen;
    gr_ptr Acoeff;
    ulong * Aexp;
    int status = GR_SUCCESS;

    base->nthreads = num_handles + 1;
    base->ndivs = base->nthreads*4;  /* number of divisions */
    base->ctx = ctx;
    base->coeff2 = coeff2;
    base->exp2 = exp2;
    base->len2 = len2;
    base->bits = bits;
    base->N = N;
    base->cmpmask = cmpmask;
    base->idx = base->ndivs - 1;    /* decremented by worker threads */
    base->have_fast_dot = (GR_VEC_DOT_OP(cctx, VEC_DOT) != (gr_method_vec_dot_op) gr_generic_vec_dot);

    /* product of lengths cannot overflow: checked by caller */
    umul_ppmm(hi, BClen, len2, len2);
    FLINT_ASSERT(hi == 0 && BClen >= 0);
    (void) hi;

    divs = (_div_struct *) flint_malloc(base->ndivs*sizeof(_div_struct));
    args = (_worker_arg_struct *) flint_malloc(base->nthreads
                                                  *sizeof(_worker_arg_struct));

    /* allocate space and set the boundary for each division */
    for (i = base->ndivs - 1; i >= 0; i--)
    {
        double d = (double)(i + 1) / (double)(base->ndivs);

        /* divisions decrease in size so that no worker finishes too early */
        divs[i].lower = (d * d) * BClen;
        divs[i].lower = FLINT_MIN(divs[i].lower, BClen);
        divs[i].lower = FLINT_MAX(divs[i].lower, WORD(0));
        divs[i].upper = divs[i].lower;
        divs[i].Aoffset = -WORD(1);
        divs[i].thread_idx = -WORD(1);
        divs[i].status = GR_SUCCESS;

        divs[i].Alen = 0;
        if (i == base->ndivs - 1)
        {
            /* highest division writes to original poly */
            divs[i].Aalloc = A->coeffs_alloc;
            divs[i].Aexps_alloc = A->exps_alloc;
            divs[i].Aexp = A->exps;
            divs[i].Acoeff = A->coeffs;
        }
        else
        {
            /* lower divisions write to a new worker poly */
            divs[i].Aalloc = 0;
            divs[i].Aexps_alloc = 0;
            divs[i].Aexp = NULL;
            divs[i].Acoeff = NULL;
        }
    }

    /* compute each chunk in parallel */
#if FLINT_USES_PTHREAD
    pthread_mutex_init(&base->mutex, NULL);
#endif
    for (i = 0; i < num_handles; i++)
    {
        args[i].idx = i;
        args[i].base = base;
        args[i].divs = divs;
        thread_pool_wake(global_thread_pool, handles[i], 0,
                               _gr_mpoly_sqr_commutative_heap_threaded_worker, &args[i]);
    }
    i = num_handles;
    args[i].idx = i;
    args[i].base = base;
    args[i].divs = divs;
    _gr_mpoly_sqr_commutative_heap_threaded_worker(&args[i]);
    for (i = 0; i < num_handles; i++)
    {
        thread_pool_wait(global_thread_pool, handles[i]);
    }

    /* collect status from all divisions */
    for (i = 0; i < base->ndivs; i++)
        status |= divs[i].status;

    /* calculate and allocate space for final answer */
    i = base->ndivs - 1;
    Alen = divs[i].Alen;
    Acoeff = divs[i].Acoeff;
    Aexp = divs[i].Aexp;
    Aalloc = divs[i].Aalloc;
    Aexps_alloc = divs[i].Aexps_alloc;
    for (i = base->ndivs - 2; i >= 0; i--)
    {
        divs[i].Aoffset = Alen;
        Alen += divs[i].Alen;
    }

    /* grow the output buffer; new slots are initialised to zero so the join
       can safely exchange the worker coefficients into place */
    _gr_mpoly_fit_length(&Acoeff, &Aalloc, &Aexp, &Aexps_alloc, N, Alen, ctx);

    base->Acoeff = Acoeff;
    base->Aexp = Aexp;

    /* join answers */
    for (i = 0; i < num_handles; i++)
    {
        thread_pool_wake(global_thread_pool, handles[i], 0, _join_worker, &args[i]);
    }
    _join_worker(&args[num_handles]);
    for (i = 0; i < num_handles; i++)
    {
        thread_pool_wait(global_thread_pool, handles[i]);
    }

#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&base->mutex);
#endif

    flint_free(args);
    flint_free(divs);

    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->coeffs_alloc = Aalloc;
    A->exps_alloc = Aexps_alloc;
    A->length = Alen;

    return status;
}


/* maxBfields gets clobbered */
static int _gr_mpoly_sqr_commutative_heap_threaded_maxfields(
    gr_mpoly_t A,
    const gr_mpoly_t B, fmpz * maxBfields,
    gr_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    mpoly_ctx_struct * mctx = GR_MPOLY_MCTX(ctx);
    slong N;
    flint_bitcnt_t exp_bits;
    ulong * cmpmask;
    ulong * Bexp;
    int freeBexp;
    slong hi, BClen;
    int status = GR_SUCCESS;
    TMP_INIT;

    TMP_START;

    /* the exponents of f^2 are twice those of f */
    _fmpz_vec_add(maxBfields, maxBfields, maxBfields, mctx->nfields);

    exp_bits = _fmpz_vec_max_bits(maxBfields, mctx->nfields);
    exp_bits = FLINT_MAX(MPOLY_MIN_BITS, exp_bits + 1);
    exp_bits = FLINT_MAX(exp_bits, B->bits);
    exp_bits = mpoly_fix_bits(exp_bits, mctx);

    N = mpoly_words_per_exp(exp_bits, mctx);
    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, exp_bits, mctx);

    /* ensure input exponents are packed into same sized fields as output */
    freeBexp = 0;
    Bexp = B->exps;
    if (exp_bits > B->bits)
    {
        freeBexp = 1;
        Bexp = (ulong *) flint_malloc(N*B->length*sizeof(ulong));
        mpoly_repack_monomials(Bexp, exp_bits, B->exps, B->bits,
                                                        B->length, mctx);
    }

    /* the threaded algorithm requires the product of lengths to fit a word;
       fall back to the serial algorithm otherwise */
    umul_ppmm(hi, BClen, B->length, B->length);

    /* deal with aliasing and do squaring */
    if (A == B)
    {
        gr_mpoly_t T;
        gr_mpoly_init3(T, 0, exp_bits, ctx);

        if (hi != 0 || BClen < 0)
            status = gr_mpoly_sqr_commutative_heap(T, B, ctx);
        else
            status = _gr_mpoly_sqr_commutative_heap_threaded(T,
                                B->coeffs, Bexp, B->length,
                                exp_bits, N, cmpmask, ctx, handles, num_handles);

        gr_mpoly_swap(T, A, ctx);
        gr_mpoly_clear(T, ctx);
    }
    else
    {
        gr_mpoly_fit_length_reset_bits(A, B->length + B->length, exp_bits, ctx);

        if (hi != 0 || BClen < 0)
            status = gr_mpoly_sqr_commutative_heap(A, B, ctx);
        else
            status = _gr_mpoly_sqr_commutative_heap_threaded(A,
                                B->coeffs, Bexp, B->length,
                                exp_bits, N, cmpmask, ctx, handles, num_handles);
    }

    if (freeBexp)
        flint_free(Bexp);

    TMP_END;

    return status;
}


int gr_mpoly_sqr_commutative_heap_threaded(
    gr_mpoly_t A,
    const gr_mpoly_t B,
    gr_mpoly_ctx_t ctx)
{
    mpoly_ctx_struct * mctx = GR_MPOLY_MCTX(ctx);
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    slong i;
    fmpz * maxBfields;
    thread_pool_handle * handles;
    slong num_handles;
    slong thread_limit;
    int status;
    TMP_INIT;

    if (B->length == 0)
        return gr_mpoly_zero(A, ctx);

    /* The squaring identity f^2 = sum a_i^2 X^{2e_i} + 2 sum_{i<j} a_i a_j ...
       assumes a commutative (or approximately commutative) coefficient ring;
       this is guaranteed a priori by the caller. */

    /* the coefficient ring must allow concurrent operations */
    if (gr_ctx_is_threadsafe(cctx) != T_TRUE)
        return gr_mpoly_sqr_commutative_heap(A, B, ctx);

    thread_limit = 1 + (B->length * B->length) / 10000;
    thread_limit = FLINT_MIN(thread_limit, B->length / 2);

    TMP_START;

    maxBfields = (fmpz *) TMP_ALLOC(mctx->nfields*sizeof(fmpz));
    for (i = 0; i < mctx->nfields; i++)
        fmpz_init(maxBfields + i);
    mpoly_max_fields_fmpz(maxBfields, B->exps, B->length, B->bits, mctx);

    num_handles = flint_request_threads(&handles, thread_limit);

    status = _gr_mpoly_sqr_commutative_heap_threaded_maxfields(A, B, maxBfields,
                                                    ctx, handles, num_handles);

    flint_give_back_threads(handles, num_handles);

    for (i = 0; i < mctx->nfields; i++)
        fmpz_clear(maxBfields + i);

    TMP_END;

    return status;
}
