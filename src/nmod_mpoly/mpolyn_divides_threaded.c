/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

#if FLINT_KNOW_STRONG_ORDER

#include "thread_pool.h"
#include "fmpz_mpoly.h"

typedef struct _nmod_mpolyn_stripe_struct
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
    const nmod_mpoly_ctx_struct * ctx;
} nmod_mpolyn_stripe_struct;

typedef nmod_mpolyn_stripe_struct nmod_mpolyn_stripe_t[1];

static void stripe_fit_length(nmod_mpolyn_stripe_struct * S, slong new_len)
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


/*
    a thread safe mpolyn supports three mutating operations
    - init from an array of terms
    - append an array of terms
    - clear out contents to a normal mpoly
*/
typedef struct _nmod_mpolyn_ts_struct
{
    n_poly_struct * volatile coeffs; /* this is coeff_array[idx] */
    ulong * volatile exps;       /* this is exp_array[idx] */
    volatile slong length;
    slong alloc;
    flint_bitcnt_t bits;
    flint_bitcnt_t idx;
    ulong * exp_array[FLINT_BITS];
    n_poly_struct * coeff_array[FLINT_BITS];
} nmod_mpolyn_ts_struct;

typedef nmod_mpolyn_ts_struct nmod_mpolyn_ts_t[1];

/* Bcoeff is changed */
static void nmod_mpolyn_ts_init(nmod_mpolyn_ts_t A,
                        n_poly_struct * Bcoeff, ulong * Bexp, slong Blen,
                      flint_bitcnt_t bits, slong N, const nmod_mpoly_ctx_t ctx)
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
              = (n_poly_struct *) flint_malloc(A->alloc*sizeof(n_poly_struct));

    for (i = 0; i < A->alloc; i++)
    {
        n_poly_init(A->coeffs + i);
    }

    A->length = Blen;
    for (i = 0; i < Blen; i++)
    {
        n_poly_swap(A->coeffs + i, Bcoeff + i);
        mpoly_monomial_set(A->exps + N*i, Bexp + N*i, N);
    }
}

static void nmod_mpolyn_ts_clear(nmod_mpolyn_ts_t A)
{
    slong i;

    for (i = 0; i < A->length; i++)
    {
        n_poly_clear(A->coeffs + i);
    }

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

static void nmod_mpolyn_ts_clear_poly(nmod_mpolyn_t Q, nmod_mpolyn_ts_t A)
{
    if (Q->alloc != 0)
    {
        slong i;

        FLINT_ASSERT(Q->exps != NULL);
        FLINT_ASSERT(Q->coeffs != NULL);

        for (i = 0; i < Q->alloc; i++)
        {
            n_poly_clear(Q->coeffs + i);
        }

        flint_free(Q->exps);
        flint_free(Q->coeffs);
    }

    Q->exps = A->exps;
    Q->coeffs = A->coeffs;
    Q->bits = A->bits;
    Q->alloc = A->alloc;
    Q->length = A->length;

    A->length = 0;
    A->coeff_array[A->idx] = NULL;
    A->exp_array[A->idx] = NULL;
    nmod_mpolyn_ts_clear(A);
}


/* put B on the end of A - Bcoeff is changed*/
static void nmod_mpolyn_ts_append(nmod_mpolyn_ts_t A,
                       n_poly_struct * Bcoeff, ulong * Bexps, slong Blen,
                                           slong N, const nmod_mpoly_ctx_t ctx)
{
/* TODO: this needs barriers on non-x86 */

    slong i;
    ulong * oldexps = A->exps;
    n_poly_struct * oldcoeffs = A->coeffs;
    slong oldlength = A->length;
    slong newlength = A->length + Blen;

    if (newlength <= A->alloc)
    {
        /* write new terms first */
        for (i = 0; i < Blen; i++)
        {
            n_poly_swap(oldcoeffs + oldlength + i, Bcoeff + i);
            mpoly_monomial_set(oldexps + N*(oldlength + i), Bexps + N*i, N);
        }
    }
    else
    {
        slong newalloc;
        ulong * newexps;
        n_poly_struct * newcoeffs;
        flint_bitcnt_t newidx;
        newidx = FLINT_BIT_COUNT(newlength - 1);
        newidx = (newidx > 8) ? newidx - 8 : 0;
        FLINT_ASSERT(newidx > A->idx);

        newalloc = UWORD(256) << newidx;
        FLINT_ASSERT(newlength <= newalloc);
        newexps = A->exp_array[newidx]
                = (ulong *) flint_malloc(N*newalloc*sizeof(ulong));
        newcoeffs = A->coeff_array[newidx]
                  = (n_poly_struct *) flint_malloc(newalloc*sizeof(n_poly_struct));

        for (i = 0; i < newalloc; i++)
        {
            n_poly_init(newcoeffs + i);
        }

        for (i = 0; i < oldlength; i++)
        {
            newcoeffs[i] = oldcoeffs[i]; /* just copy the bits */
            mpoly_monomial_set(newexps + N*i, oldexps + N*i, N);
        }
        for (i = 0; i < Blen; i++)
        {
            n_poly_swap(newcoeffs + oldlength + i, Bcoeff + i);
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
    nmod_mpolyn_t polyC;
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
    nmod_mpolyn_t polyA;
    nmod_mpolyn_t polyB;
    nmod_mpolyn_ts_t polyQ;
    const nmod_mpoly_ctx_struct * ctx;
    slong length;
    slong N;
    flint_bitcnt_t bits;
    ulong * cmpmask;
    int failed;
} divides_heap_base_struct;

typedef divides_heap_base_struct divides_heap_base_t[1];

/*
    the worker struct has a big chunk of memory in the stripe_t
    and two polys for work space
*/
typedef struct _worker_arg_struct
{
    divides_heap_base_struct * H;
    nmod_mpolyn_stripe_t S;
    nmod_mpolyn_t polyT1;
    nmod_mpolyn_t polyT2;
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
        nmod_mpolyn_clear(L->polyC, H->ctx);
    }
}
static int divides_heap_base_clear(nmod_mpolyn_t Q, divides_heap_base_t H)
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
    H->ctx = NULL;
    H->length = 0;
    H->N = 0;
    H->bits = 0;
    H->cmpmask = NULL;

    if (H->failed)
    {
        nmod_mpolyn_zero(Q, H->ctx);
        nmod_mpolyn_ts_clear(H->polyQ);
        return 0;
    }
    else
    {
        nmod_mpolyn_ts_clear_poly(Q, H->polyQ);
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


static void _nmod_mpolyn_fit_length(n_poly_struct ** coeffs,
                            ulong ** exps, slong * alloc, slong length,
                                    slong N, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong old_alloc = *alloc;
    slong new_alloc = FLINT_MAX(length, 2*old_alloc);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            *exps = (ulong *) flint_malloc(new_alloc*N*sizeof(ulong));
            *coeffs = (n_poly_struct *) flint_malloc(new_alloc*sizeof(n_poly_struct));
        } else
        {
            *exps = (ulong *) flint_realloc(*exps, new_alloc*N*sizeof(ulong));
            *coeffs = (n_poly_struct *) flint_realloc(*coeffs, new_alloc*sizeof(n_poly_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            n_poly_init(*coeffs + i);
        }
        *alloc = new_alloc;
    }
}


static slong _nmod_mpolyn_mulsub_stripe1(
    n_poly_struct ** A_coeff, ulong ** A_exp, slong * A_alloc,
    const n_poly_struct * Dcoeff, const ulong * Dexp, slong Dlen, int saveD,
    const n_poly_struct * Bcoeff, const ulong * Bexp, slong Blen,
    const n_poly_struct * Ccoeff, const ulong * Cexp, slong Clen,
    const nmod_mpolyn_stripe_t S)
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
    slong Aalloc = *A_alloc;
    n_poly_struct * Acoeff = *A_coeff;
    ulong * Aexp = *A_exp;
    ulong exp;
    slong * ends;
    ulong texp;
    slong * hind;
    n_poly_t pp;

    FLINT_ASSERT(S->N == 1);
    FLINT_ASSERT(S->bits <= FLINT_BITS);

    n_poly_init(pp);

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
            FLINT_ASSERT(mpoly_monomial_cmp1(emax, texp, maskhi)
                                                               > -upperclosed);
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

    FLINT_ASSERT(ends[0] >= startidx);

    Alen = 0;
    Di = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        while (Di < Dlen && mpoly_monomial_gt1(Dexp[Di], exp, maskhi))
        {
            _nmod_mpolyn_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, 1, S->ctx);
            Aexp[Alen] = Dexp[Di];
            if (saveD)
                n_poly_set(Acoeff + Alen, Dcoeff + Di);
            else
                n_poly_swap(Acoeff + Alen, (n_poly_struct *)(Dcoeff + Di));
            Alen++;
            Di++;
        }

        _nmod_mpolyn_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, 1, S->ctx);

        Aexp[Alen] = exp;

        if (Di < Dlen && Dexp[Di] == exp)
        {
            n_poly_set(Acoeff + Alen, Dcoeff + Di);
            Di++;
        }
        else
        {
            n_poly_zero(Acoeff + Alen);
        }

        do
        {
            x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
            do {
                hind[x->i] |= WORD(1);
                *store++ = x->i;
                *store++ = x->j;
                n_poly_mod_mul(pp, Bcoeff + x->i, Ccoeff + x->j, S->ctx->mod);
                n_poly_mod_sub(Acoeff + Alen, Acoeff + Alen, pp, S->ctx->mod);
            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && heap[1].exp == exp);

        Alen += !n_poly_is_zero(Acoeff + Alen);

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
    _nmod_mpolyn_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + Dlen - Di, 1, S->ctx);
    if (saveD)
    {
        for (i = 0; i < Dlen - Di; i++)
            n_poly_set(Acoeff + Alen + i, Dcoeff + Di + i);
    }
    else
    {
        for (i = 0; i < Dlen - Di; i++)
            n_poly_swap(Acoeff + Alen + i, (n_poly_struct *)(Dcoeff + Di + i));
    }

    mpoly_copy_monomials(Aexp + 1*Alen, Dexp + 1*Di, Dlen - Di, 1);
    Alen += Dlen - Di;

    *A_coeff = Acoeff;
    *A_exp = Aexp;
    *A_alloc = Aalloc;

    n_poly_clear(pp);

    return Alen;
}


static slong _nmod_mpolyn_mulsub_stripe(
    n_poly_struct ** A_coeff, ulong ** A_exp, slong * A_alloc,
    const n_poly_struct * Dcoeff, const ulong * Dexp, slong Dlen, int saveD,
    const n_poly_struct * Bcoeff, const ulong * Bexp, slong Blen,
    const n_poly_struct * Ccoeff, const ulong * Cexp, slong Clen,
    const nmod_mpolyn_stripe_t S)
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
    slong Aalloc = *A_alloc;
    n_poly_struct * Acoeff = *A_coeff;
    ulong * Aexp = *A_exp;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    slong * ends;
    ulong * texp;
    slong * hind;
    n_poly_t pp;

    FLINT_ASSERT(S->bits <= FLINT_BITS);

    n_poly_init(pp);

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
            mpoly_monomial_add(texp, Bexp + N*i, Cexp + N*startidx, N);
            FLINT_ASSERT(mpoly_monomial_cmp(emax, texp, N, S->cmpmask)
                                                               > -upperclosed);
        }
        while (startidx > 0)
        {
            mpoly_monomial_add(texp, Bexp + N*i, Cexp + N*(startidx - 1), N);
            if (mpoly_monomial_cmp(emax, texp, N, S->cmpmask) <= -upperclosed)
            {
                break;
            }
            startidx--;
        }

        if (endidx < Clen)
        {
            mpoly_monomial_add(texp, Bexp + N*i, Cexp + N*endidx, N);
            FLINT_ASSERT(mpoly_monomial_cmp(emin, texp, N, S->cmpmask) > 0);
        }
        while (endidx > 0)
        {
            mpoly_monomial_add(texp, Bexp + N*i, Cexp + N*(endidx - 1), N);
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

            mpoly_monomial_add(exp_list[exp_next], Bexp + N*x->i,
                                                   Cexp + N*x->j, N);

            exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                          &next_loc, &heap_len, N, S->cmpmask);
        }

        prev_startidx = startidx;
    }

    *S->startidx = startidx;
    *S->endidx = endidx;

    FLINT_ASSERT(ends[0] >= startidx);

    Alen = 0;
    Di = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        while (Di < Dlen && mpoly_monomial_gt(Dexp + N*Di, exp, N, S->cmpmask))
        {
            _nmod_mpolyn_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, N, S->ctx);
            mpoly_monomial_set(Aexp + N*Alen, Dexp + N*Di, N);
            if (saveD)
                n_poly_set(Acoeff + Alen, Dcoeff + Di);
            else
                n_poly_swap(Acoeff + Alen, (n_poly_struct *)(Dcoeff + Di));
            Alen++;
            Di++;
        }

        _nmod_mpolyn_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, N, S->ctx);

        mpoly_monomial_set(Aexp + N*Alen, exp, N);

        if (Di < Dlen && mpoly_monomial_equal(Dexp + N*Di, exp, N))
        {
            n_poly_set(Acoeff + Alen, Dcoeff + Di);
            Di++;
        }
        else
        {
            n_poly_zero(Acoeff + Alen);
        }

        do
        {
            exp_list[--exp_next] = heap[1].exp;
            x = _mpoly_heap_pop(heap, &heap_len, N, S->cmpmask);

            do {
                hind[x->i] |= WORD(1);
                *store++ = x->i;
                *store++ = x->j;
                FLINT_ASSERT(startidx <= x->j);
                FLINT_ASSERT(x->j < ends[0]);
                n_poly_mod_mul(pp, Bcoeff + x->i, Ccoeff + x->j, S->ctx->mod);
                n_poly_mod_sub(Acoeff + Alen, Acoeff + Alen, pp, S->ctx->mod);
            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

        Alen += !n_poly_is_zero(Acoeff + Alen);

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

                mpoly_monomial_add(exp_list[exp_next], Bexp + N*x->i,
                                                          Cexp + N*x->j, N);

                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                          &next_loc, &heap_len, N, S->cmpmask);
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

                mpoly_monomial_add(exp_list[exp_next], Bexp + N*x->i,
                                                       Cexp + N*x->j, N);

                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                          &next_loc, &heap_len, N, S->cmpmask);
            }
        }
    }

    FLINT_ASSERT(Di <= Dlen);
    _nmod_mpolyn_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + Dlen - Di, N, S->ctx);

    if (saveD)
    {
        for (i = 0; i < Dlen - Di; i++)
            n_poly_set(Acoeff + Alen + i, Dcoeff + Di + i);
    }
    else
    {
        for (i = 0; i < Dlen - Di; i++)
            n_poly_swap(Acoeff + Alen + i, (n_poly_struct *)(Dcoeff + Di + i));
    }

    mpoly_copy_monomials(Aexp + N*Alen, Dexp + N*Di, Dlen - Di, N);
    Alen += Dlen - Di;

    *A_coeff = Acoeff;
    *A_exp = Aexp;
    *A_alloc = Aalloc;

    n_poly_clear(pp);

    return Alen;
}

/*
    Q = stripe of A/B (assume A != 0)
    return Qlen = 0 if exact division is impossible
*/
static slong _nmod_mpolyn_divides_stripe1(
            n_poly_struct ** Q_coeff,     ulong ** Q_exp, slong * Q_alloc,
        const n_poly_struct * Acoeff, const ulong * Aexp, slong Alen,
        const n_poly_struct * Bcoeff, const ulong * Bexp, slong Blen,
                                                  const nmod_mpolyn_stripe_t S)
{
    flint_bitcnt_t bits = S->bits;
    ulong emin = S->emin[0];
    ulong cmpmask = 0;
    ulong texp;
    int lt_divides;
    slong i, j, s;
    slong next_loc, heap_len;
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    slong Qlen;
    slong Qalloc = * Q_alloc;
    n_poly_struct * Qcoeff = * Q_coeff;
    ulong * Qexp = * Q_exp;
    ulong exp;
    ulong mask;
    slong * hind;
    n_poly_t acc_lg, pp;

    FLINT_ASSERT(S->cmpmask[0] == 0);

    FLINT_ASSERT(Alen > 0);
    FLINT_ASSERT(Blen > 0);
    FLINT_ASSERT(S->N == 1);
    FLINT_ASSERT(S->bits <= FLINT_BITS);

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

    n_poly_init(acc_lg);
    n_poly_init(pp);

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

        _nmod_mpolyn_fit_length(&Qcoeff, &Qexp, &Qalloc, Qlen + 1, 1, S->ctx);

        lt_divides = mpoly_monomial_divides1(Qexp + Qlen, exp, Bexp[0], mask);

        n_poly_zero(acc_lg);
        do
        {
            x = _mpoly_heap_pop1(heap, &heap_len, cmpmask);
            do
            {
                *store++ = x->i;
                *store++ = x->j;
                if (x->i != -WORD(1))
                    hind[x->i] |= WORD(1);

                if (x->i == -WORD(1))
                {
                    n_poly_mod_add(acc_lg, acc_lg, Acoeff + x->j, S->ctx->mod);
                }
                else
                {
                    n_poly_mod_mul(pp, Bcoeff + x->i, Qcoeff + x->j, S->ctx->mod);
                    n_poly_mod_sub(acc_lg, acc_lg, pp, S->ctx->mod);
                }
            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && heap[1].exp == exp);

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

                    FLINT_ASSERT(mpoly_monomial_cmp1(Aexp[x->j], emin, cmpmask)
                                                                         >= 0);

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

        if (n_poly_is_zero(acc_lg))
        {
            continue;
        }

        n_poly_mod_divrem(Qcoeff + Qlen, pp, acc_lg, Bcoeff + 0, S->ctx->mod);

        if (!n_poly_is_zero(pp))
        {
            goto not_exact_division;
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

cleanup:

    n_poly_clear(acc_lg);
    n_poly_clear(pp);

    *Q_alloc = Qalloc;
    *Q_coeff = Qcoeff;
    *Q_exp = Qexp;

    return Qlen;

not_exact_division:
    Qlen = 0;
    goto cleanup;
}

static slong _nmod_mpolyn_divides_stripe(
                n_poly_struct ** Q_coeff,     ulong ** Q_exp, slong * Q_alloc,
            const n_poly_struct * Acoeff, const ulong * Aexp, slong Alen,
            const n_poly_struct * Bcoeff, const ulong * Bexp, slong Blen,
                                                  const nmod_mpolyn_stripe_t S)
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
    slong Qalloc = * Q_alloc;
    n_poly_struct * Qcoeff = * Q_coeff;
    ulong * Qexp = * Q_exp;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    ulong mask;
    slong * hind;
    n_poly_t acc_lg, pp;

    FLINT_ASSERT(Alen > 0);
    FLINT_ASSERT(Blen > 0);
    FLINT_ASSERT(S->bits <= FLINT_BITS);

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
    FLINT_ASSERT(i <= S->big_mem_alloc);

    n_poly_init(acc_lg);
    n_poly_init(pp);

    exp_next = 0;
    for (i = 0; i < Blen; i++)
        exp_list[i] = exps + i*N;

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
    heap[1].exp = exp_list[exp_next++];

    FLINT_ASSERT(mpoly_monomial_cmp(Aexp + N*0, S->emin, N, S->cmpmask) >= 0);

    mpoly_monomial_set(heap[1].exp, Aexp + N*0, N);

    while (heap_len > 1)
    {
        exp = heap[1].exp;

        if (mpoly_monomial_overflows(exp, N, mask))
        {
            goto not_exact_division;
        }

        FLINT_ASSERT(mpoly_monomial_cmp(exp, S->emin, N, S->cmpmask) >= 0);

        _nmod_mpolyn_fit_length(&Qcoeff, &Qexp, &Qalloc, Qlen + 1, N, S->ctx);

        lt_divides = mpoly_monomial_divides(Qexp + N*Qlen, exp,
                                                         Bexp + N*0, N, mask);

        n_poly_zero(acc_lg);
        do
        {
            exp_list[--exp_next] = heap[1].exp;
            x = _mpoly_heap_pop(heap, &heap_len, N, S->cmpmask);
            do
            {
                *store++ = x->i;
                *store++ = x->j;
                if (x->i != -WORD(1))
                    hind[x->i] |= WORD(1);

                if (x->i == -WORD(1))
                {
                    n_poly_mod_add(acc_lg, acc_lg, Acoeff + x->j, S->ctx->mod);
                }
                else
                {
                    n_poly_mod_mul(pp, Bcoeff + x->i, Qcoeff + x->j, S->ctx->mod);
                    n_poly_mod_sub(acc_lg, acc_lg, pp, S->ctx->mod);
                }
            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

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

                    FLINT_ASSERT(mpoly_monomial_cmp(exp_list[exp_next],
                                                 S->emin, N, S->cmpmask) >= 0);

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

                    mpoly_monomial_add(exp_list[exp_next], Bexp + N*x->i,
                                                              Qexp + N*x->j, N);

                    if (mpoly_monomial_cmp(exp_list[exp_next], S->emin, N,
                                                              S->cmpmask) >= 0)
                    {
                        exp_next += _mpoly_heap_insert(heap, exp_list[exp_next],
                                       x, &next_loc, &heap_len, N, S->cmpmask);
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

                    mpoly_monomial_add(exp_list[exp_next], Bexp + N*x->i,
                                                              Qexp + N*x->j, N);

                    if (mpoly_monomial_cmp(exp_list[exp_next], S->emin, N,
                                                              S->cmpmask) >= 0)
                    {
                        exp_next += _mpoly_heap_insert(heap, exp_list[exp_next],
                                       x, &next_loc, &heap_len, N, S->cmpmask);
                    }
                    else
                    {
                        hind[x->i] |= 1;
                    }
                }
            }
        }

        if (n_poly_is_zero(acc_lg))
        {
            continue;
        }

        n_poly_mod_divrem(Qcoeff + Qlen, pp, acc_lg, Bcoeff + 0, S->ctx->mod);

        if (!n_poly_is_zero(pp))
        {
            goto not_exact_division;
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

            mpoly_monomial_add(exp_list[exp_next], Bexp + N*x->i,
                                                             Qexp + N*x->j, N);

            if (mpoly_monomial_cmp(exp_list[exp_next], S->emin, N, S->cmpmask)
                                                                          >= 0)
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


cleanup:

    n_poly_clear(acc_lg);
    n_poly_clear(pp);

    *Q_alloc = Qalloc;
    *Q_coeff = Qcoeff;
    *Q_exp = Qexp;

    return Qlen;

not_exact_division:
    Qlen = 0;
    goto cleanup;
}




static slong _chunk_find_exp(ulong * exp, slong a, const divides_heap_base_t H)
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

static void chunk_mulsub(worker_arg_t W, divides_heap_chunk_t L, slong q_prev_length)
{
    divides_heap_base_struct * H = W->H;
    slong N = H->N;
    nmod_mpolyn_struct * C = L->polyC;
    const nmod_mpolyn_struct * B = H->polyB;
    const nmod_mpolyn_struct * A = H->polyA;
    nmod_mpolyn_ts_struct * Q = H->polyQ;
    nmod_mpolyn_struct * T1 = W->polyT1;
    nmod_mpolyn_stripe_struct * S = W->S;

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
            T1->length = _nmod_mpolyn_mulsub_stripe1(
                    &T1->coeffs, &T1->exps, &T1->alloc,
                    C->coeffs, C->exps, C->length, 1,
                    Q->coeffs + L->mq, Q->exps + N*L->mq, q_prev_length - L->mq,
                    B->coeffs, B->exps, B->length, S);
        }
        else
        {
            T1->length = _nmod_mpolyn_mulsub_stripe(
                    &T1->coeffs, &T1->exps, &T1->alloc,
                    C->coeffs, C->exps, C->length, 1,
                    Q->coeffs + L->mq, Q->exps + N*L->mq, q_prev_length - L->mq,
                    B->coeffs, B->exps, B->length, S);
        }
        nmod_mpolyn_swap(C, T1);
    }
    else
    {
        slong startidx, stopidx;
        if (L->upperclosed)
        {
            startidx = 0;
            stopidx = _chunk_find_exp(L->emin, 1, H);
        }
        else
        {
            startidx = _chunk_find_exp(L->emax, 1, H);
            stopidx = _chunk_find_exp(L->emin, startidx, H);
        }

        L->Cinited = 1;
        nmod_mpolyn_init(C, H->bits, H->ctx);
        nmod_mpolyn_fit_length(C, 16 + stopidx - startidx, H->ctx);

        if (N == 1)
        {
            C->length = _nmod_mpolyn_mulsub_stripe1(
                    &C->coeffs, &C->exps, &C->alloc,
                    A->coeffs + startidx, A->exps + N*startidx, stopidx - startidx, 1,
                    Q->coeffs + L->mq, Q->exps + N*L->mq, q_prev_length - L->mq,
                    B->coeffs, B->exps, B->length, S);
        }
        else
        {
            C->length = _nmod_mpolyn_mulsub_stripe(
                    &C->coeffs, &C->exps, &C->alloc,
                    A->coeffs + startidx, A->exps + N*startidx, stopidx - startidx, 1,
                    Q->coeffs + L->mq, Q->exps + N*L->mq, q_prev_length - L->mq,
                    B->coeffs, B->exps, B->length, S);
        }
    }

    L->mq = q_prev_length;
}

static void trychunk(worker_arg_t W, divides_heap_chunk_t L)
{
    divides_heap_base_struct * H = W->H;
    slong N = H->N;
    nmod_mpolyn_struct * C = L->polyC;
    slong q_prev_length;
    const nmod_mpolyn_struct * B = H->polyB;
    const nmod_mpolyn_struct * A = H->polyA;
    nmod_mpolyn_ts_struct * Q = H->polyQ;
    nmod_mpolyn_struct * T2 = W->polyT2;

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

        chunk_mulsub(W, L, q_prev_length);
    }

    if (L->producer == 1)
    {
        divides_heap_chunk_struct * next;
        n_poly_struct * Rcoeff;
        ulong * Rexp;
        slong Rlen;

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
                stopidx = _chunk_find_exp(L->emin, 1, H);
            }
            else
            {
                startidx = _chunk_find_exp(L->emax, 1, H);
                stopidx = _chunk_find_exp(L->emin, startidx, H);
            }
            Rlen = stopidx - startidx;
            Rcoeff = A->coeffs + startidx;
            Rexp = A->exps + N*startidx;
        }

        /* if we have remaining terms, add to quotient  */
        if (Rlen > 0)
        {
            nmod_mpolyn_stripe_struct * S = W->S;
            S->startidx = &L->startidx;
            S->endidx = &L->endidx;
            S->emin = L->emin;
            S->emax = L->emax;
            S->upperclosed = L->upperclosed;
            if (N == 1)
            {
                T2->length = _nmod_mpolyn_divides_stripe1(
                                    &T2->coeffs, &T2->exps, &T2->alloc,
                                       Rcoeff, Rexp, Rlen,
                                       B->coeffs, B->exps, B->length,  S);
            }
            else
            {
                T2->length = _nmod_mpolyn_divides_stripe(
                                    &T2->coeffs, &T2->exps, &T2->alloc,
                                       Rcoeff, Rexp, Rlen,
                                       B->coeffs, B->exps, B->length,  S);
            }
            if (T2->length == 0)
            {
                H->failed = 1;
                return;
            }
            else
            {
                nmod_mpolyn_ts_append(H->polyQ, T2->coeffs, T2->exps,
                                                        T2->length, N, H->ctx);
            }
        }

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
    nmod_mpolyn_stripe_struct * S = W->S;
    const nmod_mpolyn_struct * B = H->polyB;
    nmod_mpolyn_struct * T1 = W->polyT1;
    nmod_mpolyn_struct * T2 = W->polyT2;
    slong N = H->N;
    slong Blen = B->length;

    /* initialize stripe working memory */
    S->ctx = H->ctx;
    S->N = N;
    S->bits = H->bits;
    S->cmpmask = H->cmpmask;
    S->big_mem_alloc = 0;
    S->big_mem = NULL;

    stripe_fit_length(S, Blen);

    nmod_mpolyn_init(T1, H->bits, H->ctx);
    nmod_mpolyn_fit_length(T1, 16, H->ctx);

    nmod_mpolyn_init(T2, H->bits, H->ctx);
    nmod_mpolyn_fit_length(T2, 16, H->ctx);

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

    nmod_mpolyn_clear(T1, H->ctx);
    nmod_mpolyn_clear(T2, H->ctx);
    flint_free(S->big_mem);

    return;
}

/* return 1 if quotient is exact */
int nmod_mpolyn_divides_threaded_pool(
    nmod_mpolyn_t Q,
    const nmod_mpolyn_t A,
    const nmod_mpolyn_t B,
    const nmod_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    ulong mask;
    int divides;
    fmpz_mpoly_ctx_t zctx;
    fmpz_mpoly_t S;
    slong i, k, N;
    flint_bitcnt_t bits = A->bits;
    ulong * cmpmask;
    worker_arg_struct * worker_args;
    n_poly_t qcoeff, r;
    ulong * texps, * qexps;
    divides_heap_base_t H;
    TMP_INIT;

    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(A->bits == bits);
    FLINT_ASSERT(B->bits == bits);
    FLINT_ASSERT(Q->bits == bits);

    if (B->length < 2 || A->length < 2)
    {
        return nmod_mpolyn_divides(Q, A, B, ctx);
    }

    TMP_START;

    n_poly_init(qcoeff);
    n_poly_init(r);

    N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, bits, ctx->minfo);

    fmpz_mpoly_ctx_init(zctx, ctx->minfo->nvars, ctx->minfo->ord);
    fmpz_mpoly_init(S, zctx);

    n_poly_mod_rem(r, A->coeffs + 0, B->coeffs + 0, ctx->mod);
    if (!n_poly_is_zero(r))
    {
        divides = 0;
        nmod_mpolyn_zero(Q, ctx);
        goto cleanup1;
    }

    if (mpoly_divides_select_exps(S, zctx, num_handles,
                                 A->exps, A->length, B->exps, B->length, bits))
    {
        divides = 0;
        nmod_mpolyn_zero(Q, ctx);
        goto cleanup1;
    }

    /*
        At this point A and B both have at least two terms
            and the leading coefficients and monomials divide
            and the exponent selection did not give an easy exit
    */

    divides_heap_base_init(H);

    *H->polyA = *A;
    *H->polyB = *B;
    H->ctx = ctx;
    H->bits = bits;
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
        divides_heap_base_add_chunk(H, L);
    }

    H->head->upperclosed = 1;
    H->head->producer = 1;
    H->cur = H->head;

    /* generate at least the first quotient term */

    texps = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    qexps = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    mpoly_monomial_sub(qexps + N*0, A->exps + N*0, B->exps + N*0, N);
    n_poly_mod_div(qcoeff, A->coeffs + 0, B->coeffs + 0, ctx->mod); /* already checked */

    nmod_mpolyn_ts_init(H->polyQ, qcoeff, qexps, 1, H->bits, H->N, ctx);

    mpoly_monomial_add(texps, qexps + N*0, B->exps + N*1, N);

    mask = mpoly_overflow_mask_sp(bits);

    k = 1;
    while (k < A->length && mpoly_monomial_gt(A->exps + N*k, texps, N, cmpmask))
    {
        if (!mpoly_monomial_divides(qexps, A->exps + N*k, B->exps + N*0, N, mask))
        {
            H->failed = 1;
            break;
        }

        n_poly_mod_divrem(qcoeff, r, A->coeffs + k, B->coeffs + 0, ctx->mod);
        if (!n_poly_is_zero(r))
        {
            H->failed = 1;
            break;
        }

        nmod_mpolyn_ts_append(H->polyQ, qcoeff, qexps, 1, H->N, ctx);
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
    {
        thread_pool_wait(global_thread_pool, handles[i]);
    }

    flint_free(worker_args);

#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&H->mutex);
#endif

    divides = divides_heap_base_clear(Q, H);

cleanup1:
    fmpz_mpoly_clear(S, zctx);
    fmpz_mpoly_ctx_clear(zctx);

    n_poly_clear(qcoeff);
    n_poly_clear(r);

    TMP_END;

    return divides;
}
#else
typedef int this_file_is_empty;
#endif
