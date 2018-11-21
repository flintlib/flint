/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "thread_pool.h"
#include "nmod_mpoly.h"


/*
    set A = the part of B*C with exps in [start, end)
    this functions reallocates A and returns the length of A
    version for N == 1
*/
slong _nmod_mpoly_mul_heap_part1(mp_limb_t ** A_coeff, ulong ** A_exp, slong * A_alloc,
              const mp_limb_t * Bcoeff, const ulong * Bexp, slong Blen,
              const mp_limb_t * Ccoeff, const ulong * Cexp, slong Clen,
         slong * start, slong * end, slong * hind, const nmod_mpoly_stripe_t S)
{
    const ulong cmpmask = S->cmpmask[0];
    slong i, j;
    ulong exp;
    mpoly_heap_t * x;
    slong next_loc;
    slong heap_len;
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    slong Alen;
    ulong * Aexp = *A_exp;
    slong Aalloc = *A_alloc;
    mp_limb_t * Acoeff = *A_coeff;
    mp_limb_t acc0, acc1, acc2, pp0, pp1;

    FLINT_ASSERT(S->N == 1);

    /* tmp allocs from S->big_mem */
    i = 0;
    store = store_base = (slong *) (S->big_mem + i);
    i += 2*Blen*sizeof(slong);
    heap = (mpoly_heap1_s *)(S->big_mem + i);
    i += (Blen + 1)*sizeof(mpoly_heap1_s);
    chain = (mpoly_heap_t *)(S->big_mem + i);
    i += Blen*sizeof(mpoly_heap_t);
    FLINT_ASSERT(i <= S->big_mem_alloc);

    /* put all the starting nodes on the heap */
    heap_len = 1; /* heap zero index unused */
    next_loc = Blen + 4;   /* something bigger than heap can ever be */
    for (i = 0; i < Blen; i++)
    {
        hind[i] = 2*start[i] + 1;
    }
    for (i = 0; i < Blen; i++)
    {
        if (  (start[i] < end[i])
           && (  (i == 0)
              || (start[i] < start[i - 1])
              )
           )
        {
            x = chain + i;
            x->i = i;
            x->j = start[i];
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;
            _mpoly_heap_insert1(heap, Bexp[x->i] + Cexp[x->j], x,
                                                &next_loc, &heap_len, cmpmask);
        }
    }

    Alen = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;
      
        _nmod_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, 1);

        Aexp[Alen] = exp;

        acc0 = acc1 = acc2 = 0;
        do
        {
            x = _mpoly_heap_pop1(heap, &heap_len, cmpmask);

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
        Alen += (Acoeff[Alen] != UWORD(0));

        /* for each node temporarily stored */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            /* should we go right? */
            if (  (i + 1 < Blen)
               && (j + 0 < end[i + 1])
               && (hind[i + 1] == 2*j + 1)
               )
            {
                x = chain + i + 1;
                x->i = i + 1;
                x->j = j;
                x->next = NULL;
                hind[x->i] = 2*(x->j+1) + 0;
                _mpoly_heap_insert1(heap, Bexp[x->i] + Cexp[x->j], x,
                                                &next_loc, &heap_len, cmpmask);
            }

            /* should we go up? */
            if (  (j + 1 < end[i + 0])
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
                hind[x->i] = 2*(x->j+1) + 0;
                _mpoly_heap_insert1(heap, Bexp[x->i] + Cexp[x->j], x,
                                                &next_loc, &heap_len, cmpmask);
            }
        }
    }

    *A_coeff = Acoeff;
    *A_exp = Aexp;
    *A_alloc = Aalloc;

    return Alen;
}


slong _nmod_mpoly_mul_heap_part(mp_limb_t ** A_coeff, ulong ** A_exp, slong * A_alloc,
              const mp_limb_t * Bcoeff, const ulong * Bexp, slong Blen,
              const mp_limb_t * Ccoeff, const ulong * Cexp, slong Clen,
         slong * start, slong * end, slong * hind, const nmod_mpoly_stripe_t S)
{
    mp_bitcnt_t bits = S->bits;
    slong N = S->N;
    const ulong * cmpmask = S->cmpmask;
    slong i, j;
    slong next_loc;
    slong heap_len;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    mpoly_heap_t * x;
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    slong Alen;
    ulong * Aexp = *A_exp;
    slong Aalloc = *A_alloc;
    mp_limb_t * Acoeff = *A_coeff;
    ulong acc0, acc1, acc2, pp0, pp1;

    /* tmp allocs from S->big_mem */
    i = 0;
    store = store_base = (slong *) (S->big_mem + i);
    i += 2*Blen*sizeof(slong);
    exp_list = (ulong **) (S->big_mem + i);
    i += Blen*sizeof(ulong *);
    exps = (ulong *) (S->big_mem + i);
    i += Blen*N*sizeof(ulong);
    heap = (mpoly_heap_s *) (S->big_mem + i);
    i += (Blen + 1)*sizeof(mpoly_heap_s);
    chain = (mpoly_heap_t *) (S->big_mem + i);
    i += Blen*sizeof(mpoly_heap_t);
    FLINT_ASSERT(i <= S->big_mem_alloc);

    /* put all the starting nodes on the heap */
    heap_len = 1; /* heap zero index unused */
    next_loc = Blen + 4;   /* something bigger than heap can ever be */
    exp_next = 0;
    for (i = 0; i < Blen; i++)
        exp_list[i] = exps + N*i;
    for (i = 0; i < Blen; i++)
        hind[i] = 2*start[i] + 1;
    for (i = 0; i < Blen; i++)
    {
        if (  (start[i] < end[i])
           && (  (i == 0)
              || (start[i] < start[i - 1])
              )
           )
        {
            x = chain + i;
            x->i = i;
            x->j = start[i];
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;

            if (bits <= FLINT_BITS)
                mpoly_monomial_add(exp_list[exp_next], Bexp + N*x->i,
                                                       Cexp + N*x->j, N);
            else
                mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i,
                                                          Cexp + N*x->j, N);

            exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
        }
    }

    Alen = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        _nmod_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, N);

        mpoly_monomial_set(Aexp + N*Alen, exp, N);

        acc0 = acc1 = acc2 = 0;
        do
        {
            exp_list[--exp_next] = heap[1].exp;

            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);

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
        Alen += (Acoeff[Alen] != UWORD(0));

        /* for each node temporarily stored */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            /* should we go right? */
            if (  (i + 1 < Blen)
               && (j + 0 < end[i + 1])
               && (hind[i + 1] == 2*j + 1)
               )
            {
                x = chain + i + 1;
                x->i = i + 1;
                x->j = j;
                x->next = NULL;

                hind[x->i] = 2*(x->j + 1) + 0;

                if (bits <= FLINT_BITS)
                    mpoly_monomial_add(exp_list[exp_next], Bexp + N*x->i,
                                                           Cexp + N*x->j, N);
                else
                    mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i,
                                                              Cexp + N*x->j, N);

                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
            }

            /* should we go up? */
            if (  (j + 1 < end[i + 0])
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

                if (bits <= FLINT_BITS)
                    mpoly_monomial_add(exp_list[exp_next], Bexp + N*x->i,
                                                           Cexp + N*x->j, N);
                else
                    mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i,
                                                              Cexp + N*x->j, N);

                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
            }
        }
    }

    *A_coeff = Acoeff;
    *A_exp = Aexp;
    *A_alloc = Aalloc;

    return Alen;
}


/*
    The workers calculate product terms from 4*n divisions, where n is the
    number of threads.
*/

typedef struct
{
    pthread_mutex_t mutex;
    const nmod_mpoly_ctx_struct * ctx;
    slong nthreads;
    slong ndivs;
    const mp_limb_t * Bcoeff; const ulong * Bexp; slong Blen;
    const mp_limb_t * Ccoeff; const ulong * Cexp; slong Clen;
    slong N;
    mp_bitcnt_t bits;
    const ulong * cmpmask;
    volatile int idx;
}
mul_heap_threaded_base_t;

typedef struct
{
    slong lower;
    slong upper;
    slong Alen;
    slong Aalloc;
    ulong * Aexp;
    mp_limb_t * Acoeff;
}
mul_heap_threaded_div_t;

typedef struct
{
    nmod_mpoly_stripe_t S;
    slong idx;
    slong time;
    mul_heap_threaded_base_t * basep;
    mul_heap_threaded_div_t * divp;
    pthread_mutex_t mutex;
    pthread_cond_t cond;
    slong * t1, * t2, * t3, * t4;
    ulong * exp;
}
mul_heap_threaded_arg_t;


/*
    The workers simply take the next available division and calculate all
    product terms in this division.
*/

#define SWAP_PTRS(xx, yy) \
   do { \
      tt = xx; \
      xx = yy; \
      yy = tt; \
   } while (0)

void _nmod_mpoly_mul_heap_threaded_worker(void * arg_ptr)
{
    mul_heap_threaded_arg_t * arg = (mul_heap_threaded_arg_t *) arg_ptr;
    nmod_mpoly_stripe_struct * S = arg->S;
    mul_heap_threaded_div_t * divs = arg->divp;
    mul_heap_threaded_base_t * base = arg->basep;
    slong Blen = base->Blen;
    slong N = base->N;
    slong i, j;
    ulong * exp;
    slong score;
    slong * start, * end, * t1, * t2, * t3, * t4, * tt;

    exp = (ulong *) flint_malloc(N*sizeof(ulong));
    t1 = (slong *) flint_malloc(Blen*sizeof(slong));
    t2 = (slong *) flint_malloc(Blen*sizeof(slong));
    t3 = (slong *) flint_malloc(Blen*sizeof(slong));
    t4 = (slong *) flint_malloc(Blen*sizeof(slong));

    S->N = N;
    S->bits = base->bits;
    S->cmpmask = base->cmpmask;
    S->ctx = base->ctx;
    S->mod = base->ctx->ffinfo->mod;

    S->big_mem_alloc = 0;
    if (N == 1)
    {
        S->big_mem_alloc += 2*Blen*sizeof(slong);
        S->big_mem_alloc += (Blen + 1)*sizeof(mpoly_heap1_s);
        S->big_mem_alloc += Blen*sizeof(mpoly_heap_t);
    }
    else
    {
        S->big_mem_alloc += 2*Blen*sizeof(slong);
        S->big_mem_alloc += (Blen + 1)*sizeof(mpoly_heap_s);
        S->big_mem_alloc += Blen*sizeof(mpoly_heap_t);
        S->big_mem_alloc += Blen*S->N*sizeof(ulong);
        S->big_mem_alloc += Blen*sizeof(ulong *);
    }
    S->big_mem = (char *) flint_malloc(S->big_mem_alloc);

    /* get index to start working on */
    if (arg->idx + 1 < base->nthreads)
    {
        pthread_mutex_lock(&base->mutex);
        i = base->idx - 1;
        base->idx = i;
        pthread_mutex_unlock(&base->mutex);
    }
    else
    {
        i = base->ndivs - 1;
    }

    while (i >= 0)
    {
        /* calculate start */
        if (i + 1 < base-> ndivs)
        {
            mpoly_search_monomials(
                &start, exp, &score, t1, t2, t3,
                            divs[i].lower, divs[i].lower,
                            base->Bexp, base->Blen, base->Cexp, base->Clen,
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
            for (j = 0; j < base->Blen; j++)
                start[j] = 0;
        }

        /* calculate end */
        if (i > 0)
        {
            mpoly_search_monomials(
                &end, exp, &score, t2, t3, t4,
                            divs[i - 1].lower, divs[i - 1].lower,
                            base->Bexp, base->Blen, base->Cexp, base->Clen,
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
            for (j = 0; j < base->Blen; j++)
                end[j] = base->Clen;
        }
        /* t3 and t4 are free for workspace at this point */

        /* calculate products in [start, end) */
        _nmod_mpoly_fit_length(&divs[i].Acoeff, &divs[i].Aexp, &divs[i].Aalloc, 256, N);
        if (N == 1)
        {
            divs[i].Alen = _nmod_mpoly_mul_heap_part1(
                         &divs[i].Acoeff, &divs[i].Aexp, &divs[i].Aalloc,
                                      base->Bcoeff,  base->Bexp,  base->Blen,
                                      base->Ccoeff,  base->Cexp,  base->Clen,
                                                            start, end, t3, S);
        }
        else
        {
            divs[i].Alen = _nmod_mpoly_mul_heap_part(
                         &divs[i].Acoeff, &divs[i].Aexp, &divs[i].Aalloc,
                                      base->Bcoeff,  base->Bexp,  base->Blen,
                                      base->Ccoeff,  base->Cexp,  base->Clen,
                                                            start, end, t3, S);
        }

        /* get next index to work on */
        pthread_mutex_lock(&base->mutex);
        i = base->idx - 1;
        base->idx = i;
        pthread_mutex_unlock(&base->mutex);
    }

    /* clean up */
    flint_free(S->big_mem);
    flint_free(t4);
    flint_free(t3);
    flint_free(t2);
    flint_free(t1);
    flint_free(exp);
}



slong _nmod_mpoly_mul_heap_threaded(mp_limb_t ** A_coeff, ulong ** A_exp, slong * A_alloc,
                 const mp_limb_t * Bcoeff, const ulong * Bexp, slong Blen,
                 const mp_limb_t * Ccoeff, const ulong * Cexp, slong Clen,
                            mp_bitcnt_t bits, slong N, const ulong * cmpmask,
                                                   const nmod_mpoly_ctx_t ctx)
{
    slong i, j, ndivs2;
    slong max_num_workers, num_workers;
    thread_pool_handle * handles;
    mul_heap_threaded_arg_t * args;
    mul_heap_threaded_base_t * base;
    mul_heap_threaded_div_t * divs;
    slong Alen;
    mp_limb_t * Acoeff;
    ulong * Aexp;

    FLINT_ASSERT(global_thread_pool_initialized);

    max_num_workers = thread_pool_get_size(global_thread_pool);
    handles = (thread_pool_handle *) flint_malloc(max_num_workers
                                                  *sizeof(thread_pool_handle));

    num_workers = thread_pool_request(global_thread_pool,
                                                     handles, max_num_workers);

    base = flint_malloc(sizeof(mul_heap_threaded_base_t));
    base->nthreads = num_workers + 1;
    base->ndivs    = base->nthreads*4;  /* number of divisons */
    base->Bcoeff = Bcoeff;
    base->Bexp = Bexp;
    base->Blen = Blen;
    base->Ccoeff = Ccoeff;
    base->Cexp = Cexp;
    base->Clen = Clen;
    base->bits = bits;
    base->N = N;
    base->cmpmask = cmpmask;
    base->idx = base->ndivs - 1;    /* decremented by worker threads */
    base->ctx = ctx;

    ndivs2 = base->ndivs*base->ndivs;

    divs = flint_malloc(sizeof(mul_heap_threaded_div_t) * base->ndivs);
    args = flint_malloc(sizeof(mul_heap_threaded_arg_t) * base->nthreads);

    /* allocate space and set the boundary for each division */
    for (i = base->ndivs - 1; i >= 0; i--)
    {
        /* divisions decrease in size so that no worker finishes too early */
        divs[i].lower = (i + 1)*(i + 1)*Blen*Clen/ndivs2;
        divs[i].upper = divs[i].lower;

        divs[i].Alen = 0;
        if (i == base->ndivs - 1)
        {
            /* highest division writes to original poly */
            divs[i].Aalloc = *A_alloc;
            divs[i].Aexp = *A_exp;
            divs[i].Acoeff = *A_coeff;
        } else
        {
            /* lower divisions write to a new worker poly */
            divs[i].Aalloc = 0;
            divs[i].Aexp = NULL;
            divs[i].Acoeff = NULL;
        }
    }

    /* compute each chunk in parallel */
    pthread_mutex_init(&base->mutex, NULL);
    for (i = 0; i + 1 < base->nthreads; i++)
    {
        /* start ith worker */
        args[i].idx = i;
        args[i].basep = base;
        args[i].divp = divs;
        thread_pool_wake(global_thread_pool, handles[i],
                               _nmod_mpoly_mul_heap_threaded_worker, &args[i]);
    }
    /* main thread starts on highest index */
    i = base->nthreads - 1;
    args[i].idx = i;
    args[i].basep = base;
    args[i].divp = divs;
    _nmod_mpoly_mul_heap_threaded_worker(&args[i]);

    /* wait for workers to finish */
    for (i = 0; i + 1 < base->nthreads; i++)
    {
        thread_pool_wait(global_thread_pool, handles[i]);
        thread_pool_give_back(global_thread_pool, handles[i]);
    }
    pthread_mutex_destroy(&base->mutex);

    /* start concatenating the outputs */
    i = base->ndivs - 1;
    Alen = divs[i].Alen;
    Acoeff = divs[i].Acoeff;
    Aexp = divs[i].Aexp;
    *A_alloc = divs[i].Aalloc;

    /* make space for all coeffs in one call to fit length */
    j = Alen;
    for (i = base->ndivs - 2; i >= 0; i--)
        j += divs[i].Alen;
    _nmod_mpoly_fit_length(&Acoeff, &Aexp, A_alloc, j, N);

    for (i = base->ndivs - 2; i >= 0; i--)
    {
        /* transfer from worker poly to original poly */
        FLINT_ASSERT(divs[i].Acoeff != NULL);
        FLINT_ASSERT(divs[i].Aexp != NULL);
        flint_mpn_copyi(Acoeff + Alen, divs[i].Acoeff, divs[i].Alen);
        flint_mpn_copyi(Aexp + N*Alen, divs[i].Aexp, N*divs[i].Alen);
        Alen += divs[i].Alen;
        flint_free(divs[i].Acoeff);
        flint_free(divs[i].Aexp);
    }

    flint_free(handles);
    flint_free(args);
    flint_free(divs);
    flint_free(base);

    *A_coeff = Acoeff;
    *A_exp  = Aexp;
    return Alen;
}

void nmod_mpoly_mul_heap_threaded(nmod_mpoly_t A, const nmod_mpoly_t B,
                              const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)
{
    slong i, N;
    mp_bitcnt_t Abits;
    ulong * cmpmask;
    fmpz * maxBfields, * maxCfields;
    ulong * Bexp, * Cexp;
    int freeBexp, freeCexp;
    TMP_INIT;

    if (B->length == 0 || C->length == 0)
    {
        nmod_mpoly_zero(A, ctx);
        return;
    }

    TMP_START;

    /* work out bits required for A */
    maxBfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    maxCfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_init(maxBfields + i);
        fmpz_init(maxCfields + i);
    }
    mpoly_max_fields_fmpz(maxBfields, B->exps, B->length, B->bits, ctx->minfo);
    mpoly_max_fields_fmpz(maxCfields, C->exps, C->length, C->bits, ctx->minfo);
    _fmpz_vec_add(maxBfields, maxBfields, maxCfields, ctx->minfo->nfields);
    Abits = _fmpz_vec_max_bits(maxBfields, ctx->minfo->nfields);
    Abits = FLINT_MAX(MPOLY_MIN_BITS, Abits + 1);
    Abits = FLINT_MAX(Abits, B->bits);
    Abits = FLINT_MAX(Abits, C->bits);
    Abits = mpoly_fix_bits(Abits, ctx->minfo);
    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_clear(maxBfields + i);
        fmpz_clear(maxCfields + i);
    }

    N = mpoly_words_per_exp(Abits, ctx->minfo);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, Abits, ctx->minfo);

    /* ensure input exponents are packed into same sized fields as output */
    freeBexp = 0;
    Bexp = B->exps;
    if (Abits > B->bits)
    {
        freeBexp = 1;
        Bexp = (ulong *) flint_malloc(N*B->length*sizeof(ulong));
        mpoly_repack_monomials(Bexp, Abits, B->exps, B->bits, B->length, ctx->minfo);
    }

    freeCexp = 0;
    Cexp = C->exps;
    if (Abits > C->bits)
    {
        freeCexp = 1;
        Cexp = (ulong *) flint_malloc(N*C->length*sizeof(ulong));
        mpoly_repack_monomials(Cexp, Abits, C->exps, C->bits, C->length, ctx->minfo);
    }

    /* deal with aliasing and do multiplication */
    if (A == B || A == C)
    {
        nmod_mpoly_t T;
        nmod_mpoly_init2(T, B->length + C->length, ctx);
        nmod_mpoly_fit_bits(T, Abits, ctx);
        T->bits = Abits;

        /* algorithm more efficient if smaller poly first */
        if (B->length > C->length)
            T->length = _nmod_mpoly_mul_heap_threaded(
                                             &T->coeffs, &T->exps, &T->alloc,
                                                  C->coeffs, Cexp, C->length,
                                                  B->coeffs, Bexp, B->length,
                                                       Abits, N, cmpmask, ctx);
        else
            T->length = _nmod_mpoly_mul_heap_threaded(
                                             &T->coeffs, &T->exps, &T->alloc,
                                                  B->coeffs, Bexp, B->length,
                                                  C->coeffs, Cexp, C->length,
                                                       Abits, N, cmpmask, ctx);

        nmod_mpoly_swap(T, A, ctx);
        nmod_mpoly_clear(T, ctx);
    }
    else
    {
        nmod_mpoly_fit_length(A, B->length + C->length, ctx);
        nmod_mpoly_fit_bits(A, Abits, ctx);
        A->bits = Abits;

        /* algorithm more efficient if smaller poly first */
        if (B->length > C->length)
            A->length = _nmod_mpoly_mul_heap_threaded(
                                             &A->coeffs, &A->exps, &A->alloc,
                                                  C->coeffs, Cexp, C->length,
                                                  B->coeffs, Bexp, B->length,
                                                       Abits, N, cmpmask, ctx);
        else
            A->length = _nmod_mpoly_mul_heap_threaded(
                                             &A->coeffs, &A->exps, &A->alloc,
                                                  B->coeffs, Bexp, B->length,
                                                  C->coeffs, Cexp, C->length,
                                                       Abits, N, cmpmask, ctx);
    }

    if (freeBexp)
        flint_free(Bexp);

    if (freeCexp)
        flint_free(Cexp);

    TMP_END;
}
