/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"
#include "thread_pool.h"


typedef struct _fmpz_mpolyu_stripe_struct
{
    char * big_mem;
    slong big_mem_alloc;
    flint_bitcnt_t bits;
    const ulong * cmpmask;
    slong * startidx;
    slong * endidx;
    ulong emin;
    ulong emax;

    const fmpz_mpoly_ctx_struct * ctx;
    flint_bitcnt_t minor_bits;
    slong minor_N;
    ulong * minor_cmpmask;

    ulong main_overflow_mask;


    int upperclosed;

} fmpz_mpolyu_stripe_struct;

typedef fmpz_mpolyu_stripe_struct fmpz_mpolyu_stripe_t[1];


/*
    a thread safe mpolyu supports three mutating operations
    - init from an array of terms
    - append an array of terms
    - clear out contents to a normal mpolyu
*/
typedef struct _fmpz_mpolyu_ts_struct
{
    fmpz_mpoly_struct * volatile coeffs; /* this is coeff_array[idx] */
    ulong * volatile exps;       /* this is exp_array[idx] */
    volatile slong length;
    slong alloc;
    flint_bitcnt_t bits;
    flint_bitcnt_t idx;    
    ulong * exp_array[FLINT_BITS];
    fmpz_mpoly_struct * coeff_array[FLINT_BITS];
} fmpz_mpolyu_ts_struct;

typedef fmpz_mpolyu_ts_struct fmpz_mpolyu_ts_t[1];

/* Bcoeff is clobbered */
void fmpz_mpolyu_ts_init(
    fmpz_mpolyu_ts_t A,
    fmpz_mpoly_struct * Bcoeff,
    ulong * Bexp,
    slong Blen,
    flint_bitcnt_t Bbits,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    flint_bitcnt_t idx = FLINT_BIT_COUNT(Blen);
    idx = (idx <= 8) ? 0 : idx - 8;
    for (i = 0; i < FLINT_BITS; i++)
    {
        A->exp_array[i] = NULL;
        A->coeff_array[i] = NULL;
    }
    A->idx = idx;
    A->bits = Bbits,
    A->alloc = WORD(256) << idx;
    A->exps = A->exp_array[idx]
            = (ulong *) flint_malloc(A->alloc*sizeof(ulong));
    A->coeffs = A->coeff_array[idx]
              = (fmpz_mpoly_struct *) flint_malloc(A->alloc*sizeof(fmpz_mpoly_struct));
    for (i = 0; i < A->alloc; i++)
    {
        fmpz_mpoly_init3(A->coeffs + i, 0, Bbits, ctx);
    }

    A->length = Blen;
    for (i = 0; i < Blen; i++)
    {
        fmpz_mpoly_swap(A->coeffs + i, Bcoeff + i, ctx);
        A->exps[i] = Bexp[i];
    }
}

void fmpz_mpolyuu_ts_print(
    const fmpz_mpolyu_ts_t B,
    const char ** x,
    slong nmainvars,
    const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpolyu_t A;
    A->length = B->length;
    A->alloc = B->alloc;
    A->coeffs = B->coeffs;
    A->exps = B->exps;
    A->bits = B->bits;
    fmpz_mpolyuu_print_pretty(A, x, nmainvars, ctx);
    FLINT_ASSERT(fmpz_mpolyu_is_canonical(A, ctx));
}

void fmpz_mpolyu_ts_clear(fmpz_mpolyu_ts_t A, const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    for (i = 0; i < A->alloc; i++)
    {
        fmpz_mpoly_clear(A->coeffs + i, ctx);
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

void fmpz_mpolyu_ts_clear_poly(
    fmpz_mpolyu_t Q,
    fmpz_mpolyu_ts_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    if (Q->alloc != 0)
    {

        FLINT_ASSERT(Q->exps != NULL);
        FLINT_ASSERT(Q->coeffs != NULL);

        for (i = 0; i < Q->alloc; i++)
        {
            fmpz_mpoly_clear(Q->coeffs + i, ctx);
        }
        flint_free(Q->exps);
        flint_free(Q->coeffs);
    }

    Q->exps = A->exps;
    Q->coeffs = A->coeffs;
    Q->bits = A->bits;
    Q->alloc = A->alloc;
    Q->length = A->length;

    A->coeff_array[A->idx] = NULL;
    A->exp_array[A->idx] = NULL;

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


/* put B on the end of A - Bcoeff is changed*/
void fmpz_mpolyu_ts_append(
    fmpz_mpolyu_ts_t A,
    fmpz_mpoly_struct * Bcoeff,
    ulong * Bexps,
    slong Blen,
    const fmpz_mpoly_ctx_t ctx)
{
/* TODO: this needs barriers on non-x86 */
    flint_bitcnt_t bits = A->bits;
    slong i;
    ulong * oldexps = A->exps;
    fmpz_mpoly_struct * oldcoeffs = A->coeffs;
    slong oldlength = A->length;
    slong newlength = A->length + Blen;

    if (newlength <= A->alloc)
    {
        /* write new terms first */
        for (i = 0; i < Blen; i++)
        {
            FLINT_ASSERT((Bcoeff + i)->bits == bits);
            FLINT_ASSERT((oldcoeffs + oldlength + i)->bits == bits);
            fmpz_mpoly_swap(oldcoeffs + oldlength + i, Bcoeff + i, ctx);
            oldexps[oldlength + i] = Bexps[i];
        }
    }
    else
    {
        slong newalloc;
        ulong * newexps;
        fmpz_mpoly_struct * newcoeffs;
        flint_bitcnt_t newidx;
        newidx = FLINT_BIT_COUNT(newlength - 1);
        newidx = (newidx > 8) ? newidx - 8 : 0;
        FLINT_ASSERT(newidx > A->idx);

        newalloc = UWORD(256) << newidx;
        FLINT_ASSERT(newlength <= newalloc);
        newexps = A->exp_array[newidx]
                = (ulong *) flint_malloc(newalloc*sizeof(ulong));
        newcoeffs = A->coeff_array[newidx]
                  = (fmpz_mpoly_struct *) flint_malloc(newalloc*sizeof(fmpz_mpoly_struct));

        for (i = 0; i < oldlength; i++)
        {
            newcoeffs[i] = oldcoeffs[i]; /* just copy the bits */
            newexps[i] = oldexps[i];
        }
        for (i = oldlength; i < newalloc; i++)
        {
            fmpz_mpoly_init3(newcoeffs + i, 0, bits, ctx);
        }

        for (i = 0; i < Blen; i++)
        {
            FLINT_ASSERT((Bcoeff + i)->bits == bits);
            FLINT_ASSERT((newcoeffs + oldlength + i)->bits == bits);
            fmpz_mpoly_swap(newcoeffs + oldlength + i, Bcoeff + i, ctx);
            newexps[oldlength + i] = Bexps[i];
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
    fmpz_mpolyu_t polyC;
    struct _divides_heap_chunk_struct * next;
    ulong emin;
    ulong emax;
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
    volatile int failed;
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif
    divides_heap_chunk_struct * head;
    divides_heap_chunk_struct * tail;
    divides_heap_chunk_struct * volatile cur;
    const fmpz_mpolyu_t polyA;
    const fmpz_mpolyu_t polyB;
    fmpz_mpolyu_ts_t polyQ;
    slong length;

    const fmpz_mpoly_ctx_struct * ctx;
    slong minor_N;
    ulong * minor_cmpmask;
    flint_bitcnt_t minor_bits;

    slong main_nvars;
    flint_bitcnt_t main_bits;
    ulong main_overflow_mask;

} divides_heap_base_struct;

typedef divides_heap_base_struct divides_heap_base_t[1];

/*
    the worker stuct has a big chunk of memory in the stripe_t
    and two polys for work space
*/
typedef struct _worker_arg_struct
{
    divides_heap_base_struct * H;
    fmpz_mpolyu_stripe_t S;
    fmpz_mpolyu_t polyT1;
    fmpz_mpolyu_t polyT2;
} worker_arg_struct;

typedef worker_arg_struct worker_arg_t[1];


static void divides_heap_base_init(divides_heap_base_t H)
{
    H->head = NULL;
    H->tail = NULL;
    H->cur = NULL;
    H->ctx = NULL;
    H->length = 0;
    H->minor_N = 0;
    H->minor_bits = 0;
    H->minor_cmpmask = NULL;
}

static void divides_heap_chunk_clear(divides_heap_chunk_t L, divides_heap_base_t H)
{
    if (L->Cinited)
    {
        fmpz_mpolyu_clear(L->polyC, H->ctx);
    }
}
static int divides_heap_base_clear(fmpz_mpolyu_t Q, divides_heap_base_t H)
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
    H->minor_N = 0;
    H->main_bits = 0;
    H->minor_bits = 0;
    H->minor_cmpmask = NULL;

    if (H->failed)
    {
        fmpz_mpolyu_zero(Q, H->ctx);
        fmpz_mpolyu_ts_clear(H->polyQ, H->ctx);
        return 0;
    }
    else
    {
        fmpz_mpolyu_ts_clear_poly(Q, H->polyQ, H->ctx);
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


void _fmpz_mpolyu_fit_length(fmpz_mpoly_struct ** coeffs,
                              ulong ** exps, slong * alloc, slong length,
                                  flint_bitcnt_t bits, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    slong old_alloc = *alloc;
    slong new_alloc = FLINT_MAX(length, 2*old_alloc);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            *exps = (ulong *) flint_malloc(new_alloc*sizeof(ulong));
            *coeffs = (fmpz_mpoly_struct *) flint_malloc(
                                          new_alloc*sizeof(fmpz_mpoly_struct));
        }
        else
        {
            *exps = (ulong *) flint_realloc(*exps,
                                                      new_alloc*sizeof(ulong));
            *coeffs = (fmpz_mpoly_struct *) flint_realloc(*coeffs,
                                          new_alloc*sizeof(fmpz_mpoly_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fmpz_mpoly_init3(*coeffs + i, 0, bits, ctx);
        }
        *alloc = new_alloc;
    }
}


/*
void fmpz_mpolyu_fit_length(fmpz_mpolyu_t A, slong length,
                                                   const fmpz_mpoly_ctx_t uctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->exps = (ulong *) flint_malloc(new_alloc*sizeof(ulong));
            A->coeffs = (fmpz_mpoly_struct *) flint_malloc(
                                          new_alloc*sizeof(fmpz_mpoly_struct));
        } else
        {
            A->exps = (ulong *) flint_realloc(A->exps,
                                                      new_alloc*sizeof(ulong));
            A->coeffs = (fmpz_mpoly_struct *) flint_realloc(A->coeffs,
                                          new_alloc*sizeof(fmpz_mpoly_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fmpz_mpoly_init(A->coeffs + i, uctx);
            fmpz_mpoly_fit_bits(A->coeffs + i, A->bits, uctx);
            (A->coeffs + i)->bits = A->bits;
        }
        A->alloc = new_alloc;
    }
}
*/


slong _fmpz_mpoly_mulsub(
                fmpz ** A_coeff, ulong ** A_exp, slong * A_alloc,
                 const fmpz * Dcoeff, const ulong * Dexp, slong Dlen, int saveD,
                 const fmpz * Bcoeff, const ulong * Bexp, slong Blen,
                 const fmpz * Ccoeff, const ulong * Cexp, slong Clen,
                              flint_bitcnt_t bits, slong N, const ulong * cmpmask);


/*
    A = D - (a stripe of B * C)
    S->startidx and S->endidx are assumed to be correct
        that is, we expect and successive calls to keep
            B decreasing
            C the same

    saveD: 0 means we can modify coeffs of input D
           1 means we must not modify coeffs of input D
*/
slong _fmpz_mpolyuu_mulsub_stripe(fmpz_mpoly_struct ** A_coeff, ulong ** A_exp, slong * A_alloc,
                         const fmpz_mpoly_struct * Dcoeff, const ulong * Dexp, slong Dlen, int saveD,
                         const fmpz_mpoly_struct * Bcoeff, const ulong * Bexp, slong Blen,
                         const fmpz_mpoly_struct * Ccoeff, const ulong * Cexp, slong Clen,
                                                   const fmpz_mpolyu_stripe_t S)
{
    int upperclosed;
    slong startidx, endidx;
    ulong prev_startidx;
    ulong maskhi = 0;
    const ulong emax = S->emax;
    const ulong emin = S->emin;
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
    fmpz_mpoly_struct * Acoeff = *A_coeff;
    ulong * Aexp = *A_exp;
    ulong exp;
    slong * ends;
    ulong texp;
    slong * hind;
    fmpz_mpoly_t T;

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

    fmpz_mpoly_init3(T, 16, S->minor_bits, S->ctx);

    Alen = 0;
    Di = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        while (Di < Dlen && mpoly_monomial_gt1(Dexp[Di], exp, maskhi))
        {
            _fmpz_mpolyu_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, S->minor_bits, S->ctx);
            Aexp[Alen] = Dexp[Di];
            if (saveD)
                fmpz_mpoly_set(Acoeff + Alen, Dcoeff + Di, S->ctx);
            else
                fmpz_mpoly_swap(Acoeff + Alen, (fmpz_mpoly_struct *)(Dcoeff + Di), S->ctx);
            Alen++;
            Di++;
        }

        _fmpz_mpolyu_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, S->minor_bits, S->ctx);

        Aexp[Alen] = exp;

        {
            if (Di < Dlen && Dexp[Di] == exp)
            {
                fmpz_mpoly_set(Acoeff + Alen, Dcoeff + Di, S->ctx);
                Di++;
            }
            else
            {
                fmpz_mpoly_zero(Acoeff + Alen, S->ctx);
            }

            FLINT_ASSERT((Acoeff + Alen)->bits == S->minor_bits);

            do
            {
                x = _mpoly_heap_pop1(heap, &heap_len, maskhi);

                do
                {
                    hind[x->i] |= WORD(1);
                    *store++ = x->i;
                    *store++ = x->j;

                    T->length = _fmpz_mpoly_mulsub(
                                    &T->coeffs, &T->exps, &T->alloc,
                                    (Acoeff + Alen)->coeffs, (Acoeff + Alen)->exps, (Acoeff + Alen)->length, 1,
                                    (Bcoeff + x->i)->coeffs, (Bcoeff + x->i)->exps, (Bcoeff + x->i)->length,
                                    (Ccoeff + x->j)->coeffs, (Ccoeff + x->j)->exps, (Ccoeff + x->j)->length,
                                                         S->minor_bits, S->minor_N, S->minor_cmpmask);
                    fmpz_mpoly_swap(T, Acoeff + Alen, S->ctx);

                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && heap[1].exp == exp);
        }

        Alen += !fmpz_mpoly_is_zero(Acoeff + Alen, S->ctx);

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
    _fmpz_mpolyu_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + Dlen - Di, S->minor_bits, S->ctx);
    for (i = 0; i < Dlen - Di; i++)
    {
        if (saveD)
        {
            fmpz_mpoly_set(Acoeff + Alen + i, Dcoeff + Di + i, S->ctx);
        }
        else
        {
            fmpz_mpoly_swap(Acoeff + Alen + i, (fmpz_mpoly_struct *)(Dcoeff + Di + i), S->ctx);
        }
    }
    mpoly_copy_monomials(Aexp + 1*Alen, Dexp + 1*Di, Dlen - Di, 1);
    Alen += Dlen - Di;

    *A_coeff = Acoeff;
    *A_exp = Aexp;
    *A_alloc = Aalloc;

    fmpz_mpoly_clear(T, S->ctx);

    return Alen;
}


/*
    Q = stripe of A/B (assume A != 0)
    return Qlen = 0 if exact division is impossible
*/
slong _fmpz_mpolyuu_divides_stripe(
                        fmpz_mpoly_struct ** Q_coeff,     ulong ** Q_exp, slong * Q_alloc,
                    const fmpz_mpoly_struct * Acoeff, const ulong * Aexp, slong Alen,
                    const fmpz_mpoly_struct * Bcoeff, const ulong * Bexp, slong Blen,
                                                   const fmpz_mpolyu_stripe_t S)
{
    const ulong emin = S->emin;
    const ulong cmpmask = 0;
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
    fmpz_mpoly_struct * Qcoeff = * Q_coeff;
    ulong * Qexp = * Q_exp;
    ulong exp;
    slong * hind;
    fmpz_mpoly_t tT, tS;
    const fmpz_mpoly_struct * a, * b;
    fmpz_mpoly_struct * q;

    FLINT_ASSERT(Alen > 0);
    FLINT_ASSERT(Blen > 0);

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

    fmpz_mpoly_init3(tT, 16, S->minor_bits, S->ctx);
    fmpz_mpoly_init3(tS, 16, S->minor_bits, S->ctx);

    for (i = 0; i < Blen; i++)
        hind[i] = 1;

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

        if (mpoly_monomial_overflows1(exp, S->main_overflow_mask))
        {
            goto not_exact_division;
        }

        FLINT_ASSERT(mpoly_monomial_cmp1(exp, emin, cmpmask) >= 0);

        _fmpz_mpolyu_fit_length(&Qcoeff, &Qexp, &Qalloc, Qlen + 1, S->minor_bits, S->ctx);

        lt_divides = mpoly_monomial_divides1(Qexp + Qlen, exp, Bexp[0], S->main_overflow_mask);

        tT->length = 0;

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
                    a = Acoeff + x->j;

                    fmpz_mpoly_fit_length(tS, tT->length + a->length, S->ctx);
                    tS->length = _fmpz_mpoly_add(
                                    tS->coeffs, tS->exps,
                                    tT->coeffs, tT->exps, tT->length,
                                    a->coeffs, a->exps, a->length,
                                                 S->minor_N, S->minor_cmpmask);
                }
                else
                {
                    b = Bcoeff + x->i;
                    q = Qcoeff + x->j;
                    tS->length = _fmpz_mpoly_mulsub(
                                    &tS->coeffs, &tS->exps, &tS->alloc,
                                    tT->coeffs, tT->exps, tT->length, 0,
                                    b->coeffs, b->exps, b->length,
                                    q->coeffs, q->exps, q->length,
                                  S->minor_bits, S->minor_N, S->minor_cmpmask);
                }
                fmpz_mpoly_swap(tS, tT, S->ctx);

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
                /* should we go right */
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

        if (tT->length == 0)
        {
            continue;
        }

        if (mpoly_monomials_overflow_test(tT->exps, tT->length, S->minor_bits, S->ctx->minfo))
        {
            goto not_exact_division;            
        }

        q = Qcoeff + Qlen;
        FLINT_ASSERT(q->bits == S->minor_bits);
        b = Bcoeff + 0;

        q->length = _fmpz_mpoly_divides_monagan_pearce(
                            &q->coeffs, &q->exps, &q->alloc,
                            tT->coeffs, tT->exps, tT->length,
                            b->coeffs, b->exps, b->length,
                                  S->minor_bits, S->minor_N, S->minor_cmpmask);
        if (q->length == 0)
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

    fmpz_mpoly_clear(tT, S->ctx);
    fmpz_mpoly_clear(tS, S->ctx);

    *Q_alloc = Qalloc;
    *Q_coeff = Qcoeff;
    *Q_exp = Qexp;

    return Qlen;

not_exact_division:

    Qlen = 0;
    goto cleanup;
}


static slong chunk_find_exp(ulong exp, slong a, const divides_heap_base_t H)
{
    slong b = H->polyA->length;
    const ulong * Aexp = H->polyA->exps;
    ulong cmpmask = 0;

try_again:
    FLINT_ASSERT(b >= a);

    FLINT_ASSERT(a > 0);
    FLINT_ASSERT(mpoly_monomial_cmp1(Aexp[a - 1], exp, cmpmask) >= 0);
    FLINT_ASSERT(b >= H->polyA->length
                  ||  mpoly_monomial_cmp1(Aexp[b], exp, cmpmask) < 0);

    if (b - a < 5)
    {
        slong i = a;
        while (i < b
                && mpoly_monomial_cmp1(Aexp[i], exp, cmpmask) >= 0)
        {
            i++;
        }
        return i;
    }
    else
    {
        slong c = a + (b - a)/2;
        if (mpoly_monomial_cmp1(Aexp[c], exp, cmpmask) < 0)
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

static void stripe_fit_length(fmpz_mpolyu_stripe_struct * S, slong new_len)
{
    slong N = 1; /* main exponents in one word */
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
    slong N = 1; /* main exponents take one word */
    fmpz_mpolyu_struct * C = L->polyC;
    const fmpz_mpolyu_struct * B = H->polyB;
    const fmpz_mpolyu_struct * A = H->polyA;
    fmpz_mpolyu_ts_struct * Q = H->polyQ;
    fmpz_mpolyu_struct * T1 = W->polyT1;
    fmpz_mpolyu_stripe_struct * S = W->S;

    S->startidx = &L->startidx;
    S->endidx = &L->endidx;
    S->emin = L->emin;
    S->emax = L->emax;
    S->upperclosed = L->upperclosed;
    stripe_fit_length(S, q_prev_length - L->mq);

    if (L->Cinited)
    {
        T1->length = _fmpz_mpolyuu_mulsub_stripe(
                &T1->coeffs, &T1->exps, &T1->alloc,
                C->coeffs, C->exps, C->length, 1,
                Q->coeffs + L->mq, Q->exps + N*L->mq, q_prev_length - L->mq,
                B->coeffs, B->exps, B->length, S);
        fmpz_mpolyu_swap(C, T1, H->ctx);
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
        fmpz_mpolyu_init(C, H->minor_bits, H->ctx); /*any is OK*/
        fmpz_mpolyu_fit_length(C, 16 + stopidx - startidx, H->ctx);

        C->length = _fmpz_mpolyuu_mulsub_stripe(
                &C->coeffs, &C->exps, &C->alloc,
                A->coeffs + startidx, A->exps + N*startidx, stopidx - startidx, 1,
                Q->coeffs + L->mq, Q->exps + N*L->mq, q_prev_length - L->mq,
                B->coeffs, B->exps, B->length, S);
    }

    FLINT_ASSERT(fmpz_mpolyu_is_canonical(C, S->ctx));

    L->mq = q_prev_length;
}

static void trychunk(worker_arg_t W, divides_heap_chunk_t L)
{
    divides_heap_base_struct * H = W->H;
    slong N = 1; /* main exponents take one word */
    fmpz_mpolyu_struct * C = L->polyC;
    slong q_prev_length;
    const fmpz_mpolyu_struct * B = H->polyB;
    const fmpz_mpolyu_struct * A = H->polyA;
    fmpz_mpolyu_ts_struct * Q = H->polyQ;
    fmpz_mpolyu_struct * T2 = W->polyT2;

    /* return if this section has already finished processing */
    if (L->mq < 0)
    {
        return;
    }

    /* process more quotient terms if available */
    q_prev_length = Q->length;
    if (q_prev_length > L->mq)
    {
        if (L->producer == 0 && q_prev_length - L->mq < 2)
            return;

        chunk_mulsub(W, L, q_prev_length);
    }

    if (L->producer == 1)
    {
        divides_heap_chunk_struct * next;
        fmpz_mpoly_struct * Rcoeff;
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
            fmpz_mpolyu_stripe_struct * S = W->S;
            S->startidx = &L->startidx;
            S->endidx = &L->endidx;
            S->emin = L->emin;
            S->emax = L->emax;
            S->upperclosed = L->upperclosed;

            T2->length = _fmpz_mpolyuu_divides_stripe(
                                &T2->coeffs, &T2->exps, &T2->alloc,
                                   Rcoeff, Rexp, Rlen,
                                   B->coeffs, B->exps, B->length,  S);
            if (T2->length == 0)
            {
                H->failed = 1;
                return;
            }
            else
            {
                fmpz_mpolyu_ts_append(H->polyQ, T2->coeffs, T2->exps,
                                                           T2->length, H->ctx);
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
    fmpz_mpolyu_stripe_struct * S = W->S;
    const fmpz_mpolyu_struct * B = H->polyB;
    fmpz_mpolyu_struct * T1 = W->polyT1;
    fmpz_mpolyu_struct * T2 = W->polyT2;
    slong Blen = B->length;

    /* initialize stripe working memory */
    S->big_mem_alloc = 0;
    S->big_mem = NULL;

    S->ctx = H->ctx;
    S->minor_bits = H->minor_bits;
    S->minor_N = H->minor_N;
    S->minor_cmpmask = H->minor_cmpmask;

    S->main_overflow_mask = H->main_overflow_mask;

    stripe_fit_length(S, Blen);

    fmpz_mpolyu_init(T1, H->minor_bits, H->ctx);
    fmpz_mpolyu_init(T2, H->minor_bits, H->ctx);

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

    fmpz_mpolyu_clear(T1, H->ctx);
    fmpz_mpolyu_clear(T2, H->ctx);
    flint_free(S->big_mem);

    return;
}


/* return 1 if quotient is exact */
int fmpz_mpolyuu_divides_threaded_pool(
    fmpz_mpolyu_t Q,
    const fmpz_mpolyu_t A,
    const fmpz_mpolyu_t B,
    slong main_nvars,
    const fmpz_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    flint_bitcnt_t minor_bits = A->bits;
    int divides;
    fmpz_mpoly_ctx_t zctx;
    fmpz_mpoly_t S;
    slong i, k, minor_N;
    ulong * minor_cmpmask;
    flint_bitcnt_t exp_bits;
    worker_arg_struct * worker_args;
    fmpz_mpoly_t qcoeff;
    ulong texp, qexp;
    divides_heap_base_t H;
    TMP_INIT;

    FLINT_ASSERT(Q->bits == minor_bits);
    FLINT_ASSERT(A->bits == minor_bits);
    FLINT_ASSERT(B->bits == minor_bits);

#if !FLINT_KNOW_STRONG_ORDER
    return fmpz_mpolyuu_divides(Q, A, B, main_nvars, ctx);
#endif

    if (B->length < 2 || A->length < 2)
    {
        if (B->length == 0)
        {
            flint_throw(FLINT_DIVZERO,
                       "Divide by zero in fmpz_mpolyuu_divides_heap_threaded");
        }

        if (A->length == 0)
        {
            fmpz_mpolyu_zero(Q, ctx);
            return 1;
        }
        return fmpz_mpolyuu_divides(Q, A, B, main_nvars, ctx);
    }

    TMP_START;

    fmpz_mpoly_init3(qcoeff, 8, minor_bits, ctx);

    fmpz_mpoly_ctx_init(zctx, main_nvars, ORD_LEX);
    fmpz_mpoly_init(S, zctx);

    minor_N = mpoly_words_per_exp(minor_bits, ctx->minfo);
    minor_cmpmask = (ulong *) TMP_ALLOC(minor_N*sizeof(ulong));
    mpoly_get_cmpmask(minor_cmpmask, minor_N, minor_bits, ctx->minfo);

    qcoeff->length = _fmpz_mpoly_divides_monagan_pearce(
                        &qcoeff->coeffs, &qcoeff->exps, &qcoeff->alloc,
                        (A->coeffs + 0)->coeffs, (A->coeffs + 0)->exps, (A->coeffs + 0)->length,
                        (B->coeffs + 0)->coeffs, (B->coeffs + 0)->exps, (B->coeffs + 0)->length,
                                          minor_bits, minor_N, minor_cmpmask);
    if (qcoeff->length == 0)
    {
        divides = 0;
        fmpz_mpolyu_zero(Q, ctx);
        goto cleanup1;
    }

    exp_bits = FLINT_BITS/main_nvars;
    if (mpoly_divides_select_exps(S, zctx, num_handles,
                             A->exps, A->length, B->exps, B->length, exp_bits))
    {
        divides = 0;
        fmpz_mpolyu_zero(Q, ctx);
        goto cleanup1;
    }
    FLINT_ASSERT(S->bits == exp_bits);
    FLINT_ASSERT(1 == mpoly_words_per_exp(S->bits, zctx->minfo));

    /*
        At this point A and B both have at least two terms
            and the leading coefficients divide
            and the exponent selection did not give an easy exit
                in particular the leading exponents divide
    */

    divides_heap_base_init(H);

    *(fmpz_mpolyu_struct *) H->polyA = *A;
    *(fmpz_mpolyu_struct *) H->polyB = *B;

    H->failed = 0;

    H->ctx = ctx;
    H->minor_bits = minor_bits;
    H->minor_N = minor_N;
    H->minor_cmpmask = minor_cmpmask;

    H->main_bits = FLINT_BITS/main_nvars;
    H->main_nvars = main_nvars;
    H->main_overflow_mask = 0;
    for (i = 0; i < main_nvars; i++)
        H->main_overflow_mask = (H->main_overflow_mask << H->main_bits)
                                            + (UWORD(1) << (H->main_bits - 1));

    for (i = 0; i + 1 < S->length; i++)
    {
        divides_heap_chunk_struct * L;
        L = (divides_heap_chunk_struct *) flint_malloc(
                                            sizeof(divides_heap_chunk_struct));
        L->ma = 0;
        L->mq = 0;
        L->emax = S->exps[i];
        L->emin = S->exps[i + 1];
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

    qexp = A->exps[0] - B->exps[0];
    /* exponent selection was supposed to guarantee divisibility of lt */
    FLINT_ASSERT((qexp & H->main_overflow_mask) == 0);

    fmpz_mpolyu_ts_init(H->polyQ, qcoeff, &qexp, 1, H->minor_bits, H->ctx);

    /* texp might overflow into sign bit, but the comparison still works */
    texp = qexp + B->exps[1];

    k = 1;
    while (k < A->length && A->exps[k] > texp)
    {
        const fmpz_mpoly_struct * a, * b;

        qexp = A->exps[k] - B->exps[0];
        if ((qexp & H->main_overflow_mask) != 0)
        {
            H->failed = 1;
            break;
        }

        a = A->coeffs + k;
        b = B->coeffs + 0;
        qcoeff->length = _fmpz_mpoly_divides_monagan_pearce(
                            &qcoeff->coeffs, &qcoeff->exps, &qcoeff->alloc,
                                a->coeffs, a->exps, a->length,
                                b->coeffs, b->exps, b->length,
                                  H->minor_bits, H->minor_N, H->minor_cmpmask);
        if (qcoeff->length == 0)
        {
            H->failed = 1;
            break;
        }

        fmpz_mpolyu_ts_append(H->polyQ, qcoeff, &qexp, 1, H->ctx);
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

    fmpz_mpoly_clear(qcoeff, ctx);

    TMP_END;

    return divides;
}
