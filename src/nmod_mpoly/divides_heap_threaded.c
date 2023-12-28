/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

#if FLINT_KNOW_STRONG_ORDER

#include "thread_support.h"
#include "fmpz_mpoly.h"

#define PROFILE_THIS 0

#if PROFILE_THIS

#include "profiler.h"

typedef struct _vec_slong_struct
{
    slong * array;
    slong alloc;
    slong length;
} vec_slong_struct;

typedef vec_slong_struct vec_slong_t[1];

static void vec_slong_init(vec_slong_t v)
{
    v->length = 0;
    v->alloc = 16;
    v->array = (slong *) flint_malloc(v->alloc*sizeof(slong));
}

static void vec_slong_clear(vec_slong_t v)
{
    flint_free(v->array);
}

static void vec_slong_push_back(vec_slong_t v, slong a)
{
    v->length++;
    if (v->length > v->alloc)
    {
        v->alloc = FLINT_MAX(v->length, 2*v->alloc);
        v->array = (slong *) flint_realloc(v->array, v->alloc*sizeof(slong));
    }
    v->array[v->length - 1] = a;
}

static void vec_slong_print(const vec_slong_t v)
{
    slong i;
    flint_printf("[");
    for (i = 0; i < v->length; i++)
    {
        flint_printf("%wd",v->array[i]);
        if (i + 1 < v->length)
        {
            flint_printf(",",v->array[i]);
        }
    }
    flint_printf("]");
}
#endif



/*
    a thread safe mpoly supports three mutating operations
    - init from an array of terms
    - append an array of terms
    - clear out contents to a normal mpoly
*/
typedef struct _nmod_mpoly_ts_struct
{
    mp_limb_t * volatile coeffs; /* this is coeff_array[idx] */
    ulong * volatile exps;       /* this is exp_array[idx] */
    volatile slong length;
    slong alloc;
    flint_bitcnt_t bits;
    flint_bitcnt_t idx;
    mp_limb_t * exp_array[FLINT_BITS];
    ulong * coeff_array[FLINT_BITS];
} nmod_mpoly_ts_struct;

typedef nmod_mpoly_ts_struct nmod_mpoly_ts_t[1];

static void nmod_mpoly_ts_init(nmod_mpoly_ts_t A,
                              mp_limb_t * Bcoeff, ulong * Bexp, slong Blen,
                                                    flint_bitcnt_t bits, slong N)
{
    slong i;
    flint_bitcnt_t idx = FLINT_BIT_COUNT(Blen);
    idx = (idx <= 8) ? 0 : idx - 8;
    for (i = 0; i < FLINT_BITS; i++)
    {
        A->exp_array[i] = NULL;
        A->coeff_array[i] = NULL;
    }
    A->bits = bits;
    A->idx = idx;
    A->alloc = WORD(256) << idx;
    A->exps = A->exp_array[idx]
            = (ulong *) flint_malloc(N*A->alloc*sizeof(ulong));
    A->coeffs = A->coeff_array[idx]
              = (mp_limb_t *) flint_malloc(A->alloc*sizeof(mp_limb_t));
    A->length = Blen;
    for (i = 0; i < Blen; i++)
    {
        A->coeffs[i] = Bcoeff[i];
        mpoly_monomial_set(A->exps + N*i, Bexp + N*i, N);
    }
}

static void nmod_mpoly_ts_clear(nmod_mpoly_ts_t A)
{
    slong i;
    for (i = 0; i < FLINT_BITS; i++)
    {
        if (A->exp_array[i] != NULL)
        {
            FLINT_ASSERT(A->coeff_array[i] != NULL);
            flint_free(A->coeff_array[i]);
            flint_free(A->exp_array[i]);
        }
    }
}


/* put B on the end of A */
static void nmod_mpoly_ts_append(nmod_mpoly_ts_t A,
                        mp_limb_t * Bcoeff, ulong * Bexps, slong Blen, slong N)
{
/* TODO: this needs barriers on non-x86 */

    slong i;
    ulong * oldexps = A->exps;
    ulong * oldcoeffs = A->coeffs;
    slong oldlength = A->length;
    slong newlength = A->length + Blen;

    if (newlength <= A->alloc)
    {
        /* write new terms first */
        for (i = 0; i < Blen; i++)
        {
            oldcoeffs[oldlength + i] = Bcoeff[i];
            mpoly_monomial_set(oldexps + N*(oldlength + i), Bexps + N*i, N);
        }
    }
    else
    {
        slong newalloc;
        ulong * newexps;
        mp_limb_t * newcoeffs;
        flint_bitcnt_t newidx;
        newidx = FLINT_BIT_COUNT(newlength - 1);
        newidx = (newidx > 8) ? newidx - 8 : 0;
        FLINT_ASSERT(newidx > A->idx);

        newalloc = UWORD(256) << newidx;
        FLINT_ASSERT(newlength <= newalloc);
        newexps = A->exp_array[newidx]
                = (ulong *) flint_malloc(N*newalloc*sizeof(ulong));
        newcoeffs = A->coeff_array[newidx]
                  = (mp_limb_t *) flint_malloc(newalloc*sizeof(mp_limb_t));

        for (i = 0; i < oldlength; i++)
        {
            newcoeffs[i] = oldcoeffs[i];
            mpoly_monomial_set(newexps + N*i, oldexps + N*i, N);
        }
        for (i = 0; i < Blen; i++)
        {
            newcoeffs[oldlength + i] = Bcoeff[i];
            mpoly_monomial_set(newexps + N*(oldlength + i), Bexps + N*i, N);
        }

        A->alloc = newalloc;
        A->exps = newexps;
        A->coeffs = newcoeffs;
        A->idx = newidx;

        /* do not free oldcoeff/exps as other threads may be using them */
    }

    /* update length at the very end */
    A->length = newlength;
}


/*
    a chunk holds an exponent range on the dividend
*/
typedef struct _divides_heap_chunk_struct
{
    nmod_mpoly_t polyC;
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
#if PROFILE_THIS
    slong idx;
#endif
} divides_heap_chunk_struct;

typedef divides_heap_chunk_struct divides_heap_chunk_t[1];

/*
    the base struct includes a linked list of chunks
*/
typedef struct
{
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif
    divides_heap_chunk_struct * head;
    divides_heap_chunk_struct * tail;
    divides_heap_chunk_struct * volatile cur;
    nmod_mpoly_t polyA;
    nmod_mpoly_t polyB;
    nmod_mpoly_ts_t polyQ;
    const nmod_mpoly_ctx_struct * ctx;
    slong length;
    slong N;
    flint_bitcnt_t bits;
    mp_limb_t lc_inv;
    ulong * cmpmask;
    int failed;
#if PROFILE_THIS
    timeit_t timer;
#endif
} divides_heap_base_struct;

typedef divides_heap_base_struct divides_heap_base_t[1];

/*
    the worker struct has a big chunk of memory in the stripe_t
    and two polys for work space
*/
typedef struct _worker_arg_struct
{
    divides_heap_base_struct * H;
    nmod_mpoly_stripe_t S;
    nmod_mpoly_t polyT1;
    nmod_mpoly_t polyT2;
#if PROFILE_THIS
    vec_slong_t time_data;
#endif
} worker_arg_struct;

typedef worker_arg_struct worker_arg_t[1];


static void divides_heap_base_init(divides_heap_base_t H)
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

static void divides_heap_chunk_clear(divides_heap_chunk_t L, divides_heap_base_t H)
{
    if (L->Cinited)
    {
        nmod_mpoly_clear(L->polyC, H->ctx);
    }
}
static int divides_heap_base_clear(nmod_mpoly_t Q, divides_heap_base_t H)
{
    divides_heap_chunk_struct * L = H->head;
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

    if (H->failed)
    {
        nmod_mpoly_zero(Q, H->ctx);
        nmod_mpoly_ts_clear(H->polyQ);
        return 0;
    }
    else
    {
        nmod_mpoly_ts_struct * A = H->polyQ;
        slong N = mpoly_words_per_exp(A->bits, H->ctx->minfo);

        if (Q->exps)
            flint_free(Q->exps);
        if (Q->coeffs)
            flint_free(Q->coeffs);

        Q->exps = A->exps;
        Q->coeffs = A->coeffs;
        Q->bits = A->bits;
        Q->length = A->length;
        Q->coeffs_alloc = A->alloc;
        Q->exps_alloc = N*A->alloc;

        A->coeff_array[A->idx] = NULL;
        A->exp_array[A->idx] = NULL;
        nmod_mpoly_ts_clear(A);
        return 1;
    }
}

static void divides_heap_base_add_chunk(divides_heap_base_t H, divides_heap_chunk_t L)
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
    A = D - (a stripe of B * C)
    S->startidx and S->endidx are assumed to be correct
        that is, we expect and successive calls to keep
            B decreasing
            C the same
*/
static void _nmod_mpoly_mulsub_stripe1(
    nmod_mpoly_t A,
    const mp_limb_t * Dcoeff, const ulong * Dexp, slong Dlen,
    const mp_limb_t * Bcoeff, const ulong * Bexp, slong Blen,
    const mp_limb_t * Ccoeff, const ulong * Cexp, slong Clen,
    const nmod_mpoly_stripe_t S)
{
    int upperclosed;
    slong startidx, endidx;
    ulong prev_startidx;
    ulong maskhi = S->cmpmask[0];
    ulong emax = S->emax[0];
    ulong emin = S->emin[0];
    slong i, j;
    slong next_loc = Blen + 4;   /* something bigger than heap can ever be */
    slong heap_len = 1; /* heap zero index unused */
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    slong Di;
    slong Alen;
    mp_limb_t * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    ulong acc0, acc1, acc2, pp0, pp1;
    ulong exp;
    slong * ends;
    ulong texp;
    slong * hind;

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
            FLINT_ASSERT(mpoly_monomial_cmp1(emax, texp, maskhi) > -upperclosed);
        }
        while (startidx > 0)
        {
            texp = Bexp[i] + Cexp[startidx - 1];
            if (mpoly_monomial_cmp1(emax, texp, maskhi) <= -upperclosed)
            {
                break;
            }
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
            {
                break;
            }
            endidx--;
        }

        ends[i] = endidx;

        hind[i] = 2*startidx + 1;

        if (  (startidx < endidx)
           && (((ulong)startidx) < prev_startidx)
           )
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

    /* set the indices for the next time mul is called */
    *S->startidx = startidx;
    *S->endidx = endidx;

    Alen = 0;
    Di = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        while (Di < Dlen && mpoly_monomial_gt1(Dexp[Di], exp, maskhi))
        {
            _nmod_mpoly_fit_length(&Acoeff, &A->coeffs_alloc,
                                   &Aexp, &A->exps_alloc, 1, Alen + 1);
            Acoeff[Alen] = Dcoeff[Di];
            Aexp[Alen] = Dexp[Di];
            Alen++;
            Di++;
        }

        _nmod_mpoly_fit_length(&Acoeff, &A->coeffs_alloc,
                               &Aexp, &A->exps_alloc, 1, Alen + 1);

        Aexp[Alen] = exp;

        acc0 = acc1 = acc2 = 0;
        if (Di < Dlen && Dexp[Di] == exp)
        {
            acc0 = S->mod.n - Dcoeff[Di];
            Di++;
        }

        do
        {
            x = _mpoly_heap_pop1(heap, &heap_len, maskhi);

            hind[x->i] |= WORD(1);
            *store++ = x->i;
            *store++ = x->j;
            umul_ppmm(pp1, pp0, Bcoeff[x->i], Ccoeff[x->j]);
            add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);

            while ((x = x->next) != NULL)
            {
                hind[x->i] |= WORD(1);
                *store++ = x->i;
                *store++ = x->j;
                umul_ppmm(pp1, pp0, Bcoeff[x->i], Ccoeff[x->j]);
                add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);
            }
        } while (heap_len > 1 && heap[1].exp == exp);

        NMOD_RED3(Acoeff[Alen], acc2, acc1, acc0, S->mod);
        if (Acoeff[Alen] != 0)
        {
            Acoeff[Alen] = S->mod.n - Acoeff[Alen];
            Alen++;
        }

        /* process nodes taken from the heap */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            /* should we go right? */
            if (  (i + 1 < Blen)
               && (j + 0 < ends[i + 1])
               && (hind[i + 1] == 2*j + 1)
               )
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
            if (  (j + 1 < ends[i + 0])
               && ((hind[i] & 1) == 1)
               && (  (i == 0)
                  || (hind[i - 1] >= 2*(j + 2) + 1)
                  )
               )
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
    }

    FLINT_ASSERT(Di <= Dlen);
    if (Di < Dlen)
    {
        _nmod_mpoly_fit_length(&Acoeff, &A->coeffs_alloc,
                               &Aexp, &A->exps_alloc, 1, Alen + Dlen - Di);
        flint_mpn_copyi(Acoeff + Alen, Dcoeff + Di, Dlen - Di);
        mpoly_copy_monomials(Aexp + 1*Alen, Dexp + 1*Di, Dlen - Di, 1);
        Alen += Dlen - Di;
    }

    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->length = Alen;
}


static void _nmod_mpoly_mulsub_stripe(
    nmod_mpoly_t A,
    const mp_limb_t * Dcoeff, const ulong * Dexp, slong Dlen,
    const mp_limb_t * Bcoeff, const ulong * Bexp, slong Blen,
    const mp_limb_t * Ccoeff, const ulong * Cexp, slong Clen,
    const nmod_mpoly_stripe_t S)
{
    int upperclosed;
    slong startidx, endidx;
    ulong prev_startidx;
    ulong * emax = S->emax;
    ulong * emin = S->emin;
    slong N = S->N;
    slong i, j;
    slong next_loc = Blen + 4;   /* something bigger than heap can ever be */
    slong heap_len = 1; /* heap zero index unused */
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    slong Di;
    slong Alen;
    mp_limb_t * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    ulong acc0, acc1, acc2, pp0, pp1;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    slong * ends;
    ulong * texp;
    slong * hind;

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
            FLINT_ASSERT(mpoly_monomial_cmp(emax, texp, N, S->cmpmask) > -upperclosed);
        }
        while (startidx > 0)
        {
            mpoly_monomial_add_mp(texp, Bexp + N*i, Cexp + N*(startidx - 1), N);
            if (mpoly_monomial_cmp(emax, texp, N, S->cmpmask) <= -upperclosed)
            {
                break;
            }
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
            {
                break;
            }
            endidx--;
        }

        ends[i] = endidx;

        hind[i] = 2*startidx + 1;

        if (  (startidx < endidx)
           && (((ulong)startidx) < prev_startidx)
           )
        {
            x = chain + i;
            x->i = i;
            x->j = startidx;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;

            mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i, Cexp + N*x->j, N);

            if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                      &next_loc, &heap_len, N, S->cmpmask))
               exp_next--;
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

        while (Di < Dlen && mpoly_monomial_gt(Dexp + N*Di, exp, N, S->cmpmask))
        {
            _nmod_mpoly_fit_length(&Acoeff, &A->coeffs_alloc,
                                   &Aexp, &A->exps_alloc, N, Alen + 1);
            mpoly_monomial_set(Aexp + N*Alen, Dexp + N*Di, N);
            Acoeff[Alen] = Dcoeff[Di];
            Alen++;
            Di++;
        }

        _nmod_mpoly_fit_length(&Acoeff, &A->coeffs_alloc,
                               &Aexp, &A->exps_alloc, N, Alen + 1);

        mpoly_monomial_set(Aexp + N*Alen, exp, N);

        acc0 = acc1 = acc2 = 0;
        if (Di < Dlen && mpoly_monomial_equal(Dexp + N*Di, exp, N))
        {
            acc0 = S->mod.n - Dcoeff[Di];
            Di++;
        }

        do
        {
            exp_list[--exp_next] = heap[1].exp;

            x = _mpoly_heap_pop(heap, &heap_len, N, S->cmpmask);

            hind[x->i] |= WORD(1);
            *store++ = x->i;
            *store++ = x->j;
            umul_ppmm(pp1, pp0, Bcoeff[x->i], Ccoeff[x->j]);
            add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);

            while ((x = x->next) != NULL)
            {
                hind[x->i] |= WORD(1);
                *store++ = x->i;
                *store++ = x->j;
                umul_ppmm(pp1, pp0, Bcoeff[x->i], Ccoeff[x->j]);
                add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);
            }
        } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

        NMOD_RED3(Acoeff[Alen], acc2, acc1, acc0, S->mod);
        if (Acoeff[Alen] != 0)
        {
            Acoeff[Alen] = S->mod.n - Acoeff[Alen];
            Alen++;
        }

        /* process nodes taken from the heap */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            /* should we go right? */
            if (  (i + 1 < Blen)
               && (j + 0 < ends[i + 1])
               && (hind[i + 1] == 2*j + 1)
               )
            {
                x = chain + i + 1;
                x->i = i + 1;
                x->j = j;
                x->next = NULL;

                hind[x->i] = 2*(x->j + 1) + 0;

                mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i, Cexp + N*x->j, N);

                if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                      &next_loc, &heap_len, N, S->cmpmask))
                    exp_next--;
            }

            /* should we go up? */
            if (  (j + 1 < ends[i + 0])
               && ((hind[i] & 1) == 1)
               && (  (i == 0)
                  || (hind[i - 1] >= 2*(j + 2) + 1)
                  )
               )
            {
                x = chain + i;
                x->i = i;
                x->j = j + 1;
                x->next = NULL;

                hind[x->i] = 2*(x->j + 1) + 0;

                mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i, Cexp + N*x->j, N);

                if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                      &next_loc, &heap_len, N, S->cmpmask))
                    exp_next--;
            }
        }
    }

    FLINT_ASSERT(Di <= Dlen);
    if (Di < Dlen)
    {
        _nmod_mpoly_fit_length(&Acoeff, &A->coeffs_alloc,
                               &Aexp, &A->exps_alloc, N, Alen + Dlen - Di);
        flint_mpn_copyi(Acoeff + Alen, Dcoeff + Di, Dlen - Di);
        mpoly_copy_monomials(Aexp + N*Alen, Dexp + N*Di, Dlen - Di, N);
        Alen += Dlen - Di;
    }

    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->length = Alen;
}

/*
    Q = stripe of A/B (assume A != 0)
    return Qlen = 0 if exact division is impossible
*/
static int _nmod_mpoly_divides_stripe1(
    nmod_mpoly_t Q,
    const mp_limb_t * Acoeff, const ulong * Aexp, slong Alen,
    const mp_limb_t * Bcoeff, const ulong * Bexp, slong Blen,
    const nmod_mpoly_stripe_t S)
{
    flint_bitcnt_t bits = S->bits;
    ulong emin = S->emin[0];
    ulong cmpmask = S->cmpmask[0];
    ulong texp;
    int lt_divides;
    slong i, j, s;
    slong next_loc, heap_len;
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    slong Qlen;
    mp_limb_t * Qcoeff = Q->coeffs;
    ulong * Qexp = Q->exps;
    ulong exp;
    mp_limb_t acc0, acc1, acc2, pp1, pp0;
    ulong mask;
    slong * hind;

    FLINT_ASSERT(Alen > 0);
    FLINT_ASSERT(Blen > 0);
    FLINT_ASSERT(S->N == 1);

    next_loc = Blen + 4;   /* something bigger than heap can ever be */

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

    for (i = 0; i < Blen; i++)
        hind[i] = 1;

    mask = mpoly_overflow_mask_sp(bits);

    Qlen = WORD(0);

    /* s is the number of terms * (latest quotient) we should put into heap */
    s = Blen;

    /* insert (-1, 0, exp2[0]) into heap */
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

        _nmod_mpoly_fit_length(&Qcoeff, &Q->coeffs_alloc,
                               &Qexp, &Q->exps_alloc, 1, Alen + 1);

        lt_divides = mpoly_monomial_divides1(Qexp + Qlen, exp, Bexp[0], mask);

        acc0 = acc1 = acc2 = 0;
        do
        {
            x = _mpoly_heap_pop1(heap, &heap_len, cmpmask);
            do
            {
                *store++ = x->i;
                *store++ = x->j;

                if (x->i == -WORD(1))
                {
                    add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0,
                                    WORD(0), WORD(0), S->mod.n - Acoeff[x->j]);
                }
                else
                {
                    hind[x->i] |= WORD(1);
                    umul_ppmm(pp1, pp0, Bcoeff[x->i], Qcoeff[x->j]);
                    add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);
                }
            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && heap[1].exp == exp);

        NMOD_RED3(Qcoeff[Qlen], acc2, acc1, acc0, S->mod);

        /* process nodes taken from the heap */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            if (i == -WORD(1))
            {
                /* take next dividend term */
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
                /* should we go up */
                if (  (i + 1 < Blen)
                   && (hind[i + 1] == 2*j + 1)
                   )
                {
                    x = chain + i + 1;
                    x->i = i + 1;
                    x->j = j;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;

                    texp = Bexp[x->i] + Qexp[x->j];
                    if (mpoly_monomial_cmp1(texp, emin, cmpmask) >= 0)
                    {
                        _mpoly_heap_insert1(heap, texp, x,
                                                &next_loc, &heap_len, cmpmask);
                    }
                    else
                    {
                        hind[x->i] |= 1;
                    }
                }
                /* should we go up? */
                if (j + 1 == Qlen)
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

                    texp = Bexp[x->i] + Qexp[x->j];
                    if (mpoly_monomial_cmp1(texp, emin, cmpmask) >= 0)
                    {
                        _mpoly_heap_insert1(heap, texp, x,
                                                &next_loc, &heap_len, cmpmask);
                    }
                    else
                    {
                        hind[x->i] |= 1;
                    }
                }
            }
        }

        Qcoeff[Qlen] = nmod_mul(Qcoeff[Qlen], S->lc_minus_inv, S->mod);
        if (Qcoeff[Qlen] == 0)
        {
            continue;
        }

        if (!lt_divides)
        {
            goto not_exact_division;
        }

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
            {
                _mpoly_heap_insert1(heap, texp, x,
                                         &next_loc, &heap_len, cmpmask);
            }
            else
            {
                hind[x->i] |= 1;
            }
        }
        s = 1;
        Qlen++;
    }

    Q->coeffs = Qcoeff;
    Q->exps = Qexp;
    Q->length = Qlen;

    return 1;

not_exact_division:

    Q->coeffs = Qcoeff;
    Q->exps = Qexp;
    Q->length = 0;

    return 0;
}

static int _nmod_mpoly_divides_stripe(
    nmod_mpoly_t Q,
    const mp_limb_t * Acoeff, const ulong * Aexp, slong Alen,
    const mp_limb_t * Bcoeff, const ulong * Bexp, slong Blen,
    const nmod_mpoly_stripe_t S)
{
    flint_bitcnt_t bits = S->bits;
    slong N = S->N;
    int lt_divides;
    slong i, j, s;
    slong next_loc, heap_len;
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    slong Qlen;
    mp_limb_t * Qcoeff = Q->coeffs;
    ulong * Qexp = Q->exps;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    mp_limb_t acc0, acc1, acc2, pp1, pp0;
    ulong mask;
    slong * hind;

    FLINT_ASSERT(Alen > 0);
    FLINT_ASSERT(Blen > 0);

    next_loc = Blen + 4;   /* something bigger than heap can ever be */

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

    exp_next = 0;
    for (i = 0; i < Blen; i++)
        exp_list[i] = exps + i*N;

    for (i = 0; i < Blen; i++)
        hind[i] = 1;

    mask = bits <= FLINT_BITS ? mpoly_overflow_mask_sp(bits) : 0;

    Qlen = WORD(0);

    /* s is the number of terms * (latest quotient) we should put into heap */
    s = Blen;

    /* insert (-1, 0, exp2[0]) into heap */
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
        _nmod_mpoly_fit_length(&Qcoeff, &Q->coeffs_alloc,
                               &Qexp, &Q->exps_alloc, N, Qlen + 1);

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

        acc0 = acc1 = acc2 = 0;
        do {
            exp_list[--exp_next] = heap[1].exp;
            x = _mpoly_heap_pop(heap, &heap_len, N, S->cmpmask);
            do {
                *store++ = x->i;
                *store++ = x->j;

                if (x->i == -WORD(1))
                {
                    add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0,
                                 WORD(0), WORD(0), S->mod.n - Acoeff[x->j]);
                }
                else
                {
                    hind[x->i] |= WORD(1);
                    umul_ppmm(pp1, pp0, Bcoeff[x->i], Qcoeff[x->j]);
                    add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);
                }
            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

        NMOD_RED3(Qcoeff[Qlen], acc2, acc1, acc0, S->mod);

        /* process nodes taken from the heap */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            if (i == -WORD(1))
            {
                /* take next dividend term */
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
                /* should we go up */
                if (  (i + 1 < Blen)
                   && (hind[i + 1] == 2*j + 1)
                   )
                {
                    x = chain + i + 1;
                    x->i = i + 1;
                    x->j = j;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;

                    mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i,
                                                              Qexp + N*x->j, N);

                    if (mpoly_monomial_cmp(exp_list[exp_next], S->emin, N, S->cmpmask) >= 0)
                    {
                        exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                          &next_loc, &heap_len, N, S->cmpmask);
                    }
                    else
                    {
                        hind[x->i] |= 1;
                    }
                }
                /* should we go up? */
                if (j + 1 == Qlen)
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

                    mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i,
                                                              Qexp + N*x->j, N);

                    if (mpoly_monomial_cmp(exp_list[exp_next], S->emin, N, S->cmpmask) >= 0)
                    {
                        exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                          &next_loc, &heap_len, N, S->cmpmask);
                    }
                    else
                    {
                        hind[x->i] |= 1;
                    }
                }
            }
        }

        Qcoeff[Qlen] = nmod_mul(Qcoeff[Qlen], S->lc_minus_inv, S->mod);
        if (Qcoeff[Qlen] == 0)
        {
            continue;
        }

        if (!lt_divides)
        {
            goto not_exact_division;
        }

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
            {

                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                          &next_loc, &heap_len, N, S->cmpmask);
            }
            else
            {
                hind[x->i] |= 1;
            }
        }
        s = 1;
        Qlen++;
    }

    Q->coeffs = Qcoeff;
    Q->exps = Qexp;
    Q->length = Qlen;

    return 1;

not_exact_division:

    Q->coeffs = Qcoeff;
    Q->exps = Qexp;
    Q->length = 0;

    return 0;
}


static slong chunk_find_exp(ulong * exp, slong a, const divides_heap_base_t H)
{
    slong N = H->N;
    slong b = H->polyA->length;
    const ulong * Aexp = H->polyA->exps;

try_again:
    FLINT_ASSERT(b >= a);

    FLINT_ASSERT(a > 0);
    FLINT_ASSERT(mpoly_monomial_cmp(Aexp + N*(a - 1), exp, N, H->cmpmask) >= 0);
    FLINT_ASSERT(b >= H->polyA->length
                  ||  mpoly_monomial_cmp(Aexp + N*b, exp, N, H->cmpmask) < 0);

    if (b - a < 5)
    {
        slong i = a;
        while (i < b
                && mpoly_monomial_cmp(Aexp + N*i, exp, N, H->cmpmask) >= 0)
        {
            i++;
        }
        return i;
    }
    else
    {
        slong c = a + (b - a)/2;
        if (mpoly_monomial_cmp(Aexp + N*c, exp, N, H->cmpmask) < 0)
        {
            b = c;
        }
        else
        {
            a = c;
        }
        goto try_again;
    }
}

static void stripe_fit_length(nmod_mpoly_stripe_struct * S, slong new_len)
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
    {
        return;
    }

    new_alloc = FLINT_MAX(new_alloc, S->big_mem_alloc + S->big_mem_alloc/4);
    S->big_mem_alloc = new_alloc;

    if (S->big_mem != NULL)
    {
        S->big_mem = (char *) flint_realloc(S->big_mem, new_alloc);
    }
    else
    {
        S->big_mem = (char *) flint_malloc(new_alloc);
    }

}


static void chunk_mulsub(worker_arg_t W, divides_heap_chunk_t L, slong q_prev_length)
{
    divides_heap_base_struct * H = W->H;
    slong N = H->N;
    nmod_mpoly_struct * C = L->polyC;
    const nmod_mpoly_struct * B = H->polyB;
    const nmod_mpoly_struct * A = H->polyA;
    nmod_mpoly_ts_struct * Q = H->polyQ;
    nmod_mpoly_struct * T1 = W->polyT1;
    nmod_mpoly_stripe_struct * S = W->S;

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
            _nmod_mpoly_mulsub_stripe1(T1,
                    C->coeffs, C->exps, C->length,
                    Q->coeffs + L->mq, Q->exps + N*L->mq, q_prev_length - L->mq,
                    B->coeffs, B->exps, B->length, S);
        }
        else
        {
            _nmod_mpoly_mulsub_stripe(T1,
                    C->coeffs, C->exps, C->length,
                    Q->coeffs + L->mq, Q->exps + N*L->mq, q_prev_length - L->mq,
                    B->coeffs, B->exps, B->length, S);
        }
        nmod_mpoly_swap(C, T1, H->ctx);
    }
    else
    {
        slong startidx, stopidx;
        if (L->upperclosed)
        {
            startidx = 0;
            stopidx = chunk_find_exp(L->emin, 1, H);
        }
        else
        {
            startidx = chunk_find_exp(L->emax, 1, H);
            stopidx = chunk_find_exp(L->emin, startidx, H);
        }

        L->Cinited = 1;
        nmod_mpoly_init3(C, 16 + stopidx - startidx, H->bits, H->ctx); /*any is OK*/

        if (N == 1)
        {
            _nmod_mpoly_mulsub_stripe1(C,
                    A->coeffs + startidx, A->exps + N*startidx, stopidx - startidx,
                    Q->coeffs + L->mq, Q->exps + N*L->mq, q_prev_length - L->mq,
                    B->coeffs, B->exps, B->length, S);
        }
        else
        {
            _nmod_mpoly_mulsub_stripe(C,
                    A->coeffs + startidx, A->exps + N*startidx, stopidx - startidx,
                    Q->coeffs + L->mq, Q->exps + N*L->mq, q_prev_length - L->mq,
                    B->coeffs, B->exps, B->length, S);
        }
    }

    L->mq = q_prev_length;
}

static void trychunk(worker_arg_t W, divides_heap_chunk_t L)
{
    divides_heap_base_struct * H = W->H;
    slong i;
    slong N = H->N;
    nmod_mpoly_struct * C = L->polyC;
    slong q_prev_length;
    ulong mask;
    const nmod_mpoly_struct * B = H->polyB;
    const nmod_mpoly_struct * A = H->polyA;
    nmod_mpoly_ts_struct * Q = H->polyQ;
    nmod_mpoly_struct * T2 = W->polyT2;

    mask = 0;
    for (i = 0; i < FLINT_BITS/H->bits; i++)
        mask = (mask << H->bits) + (UWORD(1) << (H->bits - 1));

    /* return if this section has already finished processing */
    if (L->mq < 0)
    {
        return;
    }

    /* process more quotient terms if available */
    q_prev_length = Q->length;
    if (q_prev_length > L->mq)
    {
        if (L->producer == 0 && q_prev_length - L->mq < 20)
            return;

#if PROFILE_THIS
        vec_slong_push_back(W->time_data, 4*L->idx + 0);
        vec_slong_push_back(W->time_data, timeit_query_wall(H->timer));
#endif
        chunk_mulsub(W, L, q_prev_length);
#if PROFILE_THIS
        vec_slong_push_back(W->time_data, 4*L->idx + 1);
        vec_slong_push_back(W->time_data, timeit_query_wall(H->timer));
#endif
    }

    if (L->producer == 1)
    {
        divides_heap_chunk_struct * next;
        mp_limb_t * Rcoeff;
        ulong * Rexp;
        slong Rlen;

#if PROFILE_THIS
        vec_slong_push_back(W->time_data, 4*L->idx + 2);
        vec_slong_push_back(W->time_data, timeit_query_wall(H->timer));
#endif

        /* process the remaining quotient terms */
        q_prev_length = Q->length;
        if (q_prev_length > L->mq)
        {
            chunk_mulsub(W, L, q_prev_length);
        }

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
                stopidx = chunk_find_exp(L->emin, 1, H);
            }
            else
            {
                startidx = chunk_find_exp(L->emax, 1, H);
                stopidx = chunk_find_exp(L->emin, startidx, H);
            }
            Rlen = stopidx - startidx;
            Rcoeff = A->coeffs + startidx;
            Rexp = A->exps + N*startidx;
        }

        /* if we have remaining terms, add to quotient  */
        if (Rlen > 0)
        {
            int divides;
            nmod_mpoly_stripe_struct * S = W->S;
            S->startidx = &L->startidx;
            S->endidx = &L->endidx;
            S->emin = L->emin;
            S->emax = L->emax;
            S->upperclosed = L->upperclosed;
            if (N == 1)
            {
                divides = _nmod_mpoly_divides_stripe1(T2, Rcoeff, Rexp, Rlen,
                                             B->coeffs, B->exps, B->length, S);
            }
            else
            {
                divides = _nmod_mpoly_divides_stripe(T2, Rcoeff, Rexp, Rlen,
                                            B->coeffs, B->exps, B->length,  S);
            }

            if (!divides)
            {
#if PROFILE_THIS
                vec_slong_push_back(W->time_data, 4*L->idx + 3);
                vec_slong_push_back(W->time_data, timeit_query_wall(H->timer));
#endif
                H->failed = 1;
                return;
            }
            else
            {
                nmod_mpoly_ts_append(H->polyQ, T2->coeffs, T2->exps, T2->length, N);
            }
        }

#if PROFILE_THIS
        vec_slong_push_back(W->time_data, 4*L->idx + 3);
        vec_slong_push_back(W->time_data, timeit_query_wall(H->timer));
#endif

        next = L->next;
        H->length--;
        H->cur = next;

        if (next != NULL)
        {
            next->producer = 1;
        }

        L->producer = 0;
        L->mq = -1;
    }

    return;
}


static void worker_loop(void * varg)
{
    worker_arg_struct * W = (worker_arg_struct *) varg;
    divides_heap_base_struct * H = W->H;
    nmod_mpoly_stripe_struct * S = W->S;
    const nmod_mpoly_struct * B = H->polyB;
    nmod_mpoly_struct * T1 = W->polyT1;
    nmod_mpoly_struct * T2 = W->polyT2;
    slong N = H->N;
    slong Blen = B->length;

    /* initialize stripe working memory */
    S->N = N;
    S->bits = H->bits;
    S->ctx = H->ctx;
    S->cmpmask = H->cmpmask;
    S->big_mem_alloc = 0;
    S->big_mem = NULL;
    S->mod = H->ctx->mod;
    S->lc_minus_inv = S->mod.n - H->lc_inv;

    stripe_fit_length(S, Blen);

    nmod_mpoly_init3(T1, 16, H->bits, H->ctx);
    nmod_mpoly_init3(T2, 16, H->bits, H->ctx);

    while (!H->failed)
    {
        divides_heap_chunk_struct * L;
        L = H->cur;

        if (L == NULL)
        {
            break;
        }
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

    nmod_mpoly_clear(T1, H->ctx);
    nmod_mpoly_clear(T2, H->ctx);
    flint_free(S->big_mem);

    return;
}


/*
    return 1 if quotient is exact.
    The leading coefficient of B should be invertible.
*/
int _nmod_mpoly_divides_heap_threaded_pool(
    nmod_mpoly_t Q,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    ulong mask;
    int divides;
    fmpz_mpoly_ctx_t zctx;
    fmpz_mpoly_t S;
    slong i, k, N;
    flint_bitcnt_t exp_bits;
    ulong * cmpmask;
    ulong * Aexp, * Bexp;
    int freeAexp, freeBexp;
    worker_arg_struct * worker_args;
    mp_limb_t qcoeff;
    ulong * texps, * qexps;
    divides_heap_base_t H;
#if PROFILE_THIS
    slong idx = 0;
#endif
    TMP_INIT;

    if (B->length < 2 || A->length < 2)
    {
        return nmod_mpoly_divides_monagan_pearce(Q, A, B, ctx);
    }

    TMP_START;

    exp_bits = MPOLY_MIN_BITS;
    exp_bits = FLINT_MAX(exp_bits, A->bits);
    exp_bits = FLINT_MAX(exp_bits, B->bits);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);

    N = mpoly_words_per_exp(exp_bits, ctx->minfo);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, exp_bits, ctx->minfo);

    /* ensure input exponents packed to same size as output exponents */
    Aexp = A->exps;
    freeAexp = 0;
    if (exp_bits > A->bits)
    {
        freeAexp = 1;
        Aexp = (ulong *) flint_malloc(N*A->length*sizeof(ulong));
        mpoly_repack_monomials(Aexp, exp_bits, A->exps, A->bits,
                                                        A->length, ctx->minfo);
    }

    Bexp = B->exps;
    freeBexp = 0;
    if (exp_bits > B->bits)
    {
        freeBexp = 1;
        Bexp = (ulong *) flint_malloc(N*B->length*sizeof(ulong));
        mpoly_repack_monomials(Bexp, exp_bits, B->exps, B->bits,
                                                    B->length, ctx->minfo);
    }

    fmpz_mpoly_ctx_init(zctx, ctx->minfo->nvars, ctx->minfo->ord);
    fmpz_mpoly_init(S, zctx);

    if (mpoly_divides_select_exps(S, zctx, num_handles,
                                   Aexp, A->length, Bexp, B->length, exp_bits))
    {
        divides = 0;
        nmod_mpoly_zero(Q, ctx);
        goto cleanup1;
    }

    /*
        At this point A and B both have at least two terms and the exponent
        selection did not give an easy exit. Since we are possibly holding
        threads, we should try to not throw, so the leading coefficient should
        already have been checked for invertibility.
    */
    divides_heap_base_init(H);
    qcoeff = n_gcdinv(&H->lc_inv, B->coeffs[0], ctx->mod.n);
    FLINT_ASSERT(qcoeff == 1); /* gcd should be one */
    H->polyA->coeffs = A->coeffs;
    H->polyA->exps = Aexp;
    H->polyA->bits = exp_bits;
    H->polyA->length = A->length;
    H->polyA->coeffs_alloc = A->coeffs_alloc;
    H->polyA->exps_alloc = A->exps_alloc;

    H->polyB->coeffs = B->coeffs;
    H->polyB->exps = Bexp;
    H->polyB->bits = exp_bits;
    H->polyB->length = B->length;
    H->polyB->coeffs_alloc = B->coeffs_alloc;
    H->polyB->exps_alloc = B->coeffs_alloc;

    H->ctx = ctx;
    H->bits = exp_bits;
    H->N = N;
    H->cmpmask = cmpmask;
    H->failed = 0;

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
#if PROFILE_THIS
        L->idx = idx++;
#endif
        divides_heap_base_add_chunk(H, L);
    }

    H->head->upperclosed = 1;
    H->head->producer = 1;
    H->cur = H->head;

    /* generate at least the first quotient terms */

    texps = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    qexps = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    mpoly_monomial_sub_mp(qexps + N*0, Aexp + N*0, Bexp + N*0, N);
    qcoeff = nmod_mul(H->lc_inv, A->coeffs[0], ctx->mod);

    nmod_mpoly_ts_init(H->polyQ, &qcoeff, qexps, 1, H->bits, H->N);

    mpoly_monomial_add_mp(texps, qexps + N*0, Bexp + N*1, N);

    mask = 0;
    for (i = 0; i < FLINT_BITS/exp_bits; i++)
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
            break;
        }
        qcoeff = nmod_mul(H->lc_inv, A->coeffs[k], ctx->mod);
        nmod_mpoly_ts_append(H->polyQ, &qcoeff, qexps, 1, H->N);
        k++;
    }

    /* start the workers */

#if FLINT_USES_PTHREAD
    pthread_mutex_init(&H->mutex, NULL);
#endif

    worker_args = (worker_arg_struct *) flint_malloc((num_handles + 1)
                                                        *sizeof(worker_arg_t));

#if PROFILE_THIS
    for (i = 0; i < num_handles + 1; i++)
    {
        vec_slong_init((worker_args + i)->time_data);
    }
    timeit_start(H->timer);
#endif

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

#if PROFILE_THIS
    timeit_stop(H->timer);
    flint_printf("data = [");
    for (i = 0; i < num_handles + 1; i++)
    {
        flint_printf("[%wd,", i);
        vec_slong_print((worker_args + i)->time_data);
        flint_printf("],\n");
        vec_slong_clear((worker_args + i)->time_data);
    }
    flint_printf("%wd]\n", H->timer->wall);
#endif

    flint_free(worker_args);

#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&H->mutex);
#endif

    divides = divides_heap_base_clear(Q, H);

cleanup1:

    fmpz_mpoly_clear(S, zctx);
    fmpz_mpoly_ctx_clear(zctx);

    if (freeAexp)
        flint_free(Aexp);

    if (freeBexp)
        flint_free(Bexp);

    TMP_END;

    return divides;
}


int nmod_mpoly_divides_heap_threaded(
    nmod_mpoly_t Q,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx)
{
    thread_pool_handle * handles;
    slong num_handles;
    int divides;
    slong thread_limit = A->length/32;

    if (B->length == 0)
    {
        if (A->length == 0 || nmod_mpoly_ctx_modulus(ctx) == 1)
        {
            nmod_mpoly_set(Q, A, ctx);
            return 1;
        }
        else
        {
            flint_throw(FLINT_DIVZERO, "nmod_mpoly_divides_heap_threaded: divide by zero");
        }
    }

    if (B->length < 2 || A->length < 2)
    {
        if (A->length == 0)
        {
            nmod_mpoly_zero(Q, ctx);
            return 1;
        }

        return nmod_mpoly_divides_monagan_pearce(Q, A, B, ctx);
    }

    if (1 != n_gcd(B->coeffs[0], ctx->mod.n))
    {
        flint_throw(FLINT_IMPINV, "nmod_mpoly_divides_heap_threaded: Cannot invert leading coefficient");
    }

    num_handles = flint_request_threads(&handles, thread_limit);

    divides = _nmod_mpoly_divides_heap_threaded_pool(Q, A, B, ctx,
                                                         handles, num_handles);

    flint_give_back_threads(handles, num_handles);

    return divides;
}
#else
typedef int this_file_is_empty;
#endif
