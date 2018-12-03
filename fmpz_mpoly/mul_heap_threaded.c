/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "thread_pool.h"
#include "fmpz_mpoly.h"


typedef struct _fmpz_mpoly_stripe_struct
{
    char * big_mem;
    slong big_mem_alloc;
    slong N;
    mp_bitcnt_t bits;
    const ulong * cmpmask;
    int flint_small;
} fmpz_mpoly_stripe_struct;

typedef fmpz_mpoly_stripe_struct fmpz_mpoly_stripe_t[1];



slong _fmpz_mpoly_mul_heap_part1(fmpz ** A_coeff, ulong ** A_exp, slong * A_alloc,
              const fmpz * Bcoeff, const ulong * Bexp, slong Blen,
              const fmpz * Ccoeff, const ulong * Cexp, slong Clen,
         slong * start, slong * end, slong * hind, const fmpz_mpoly_stripe_t S)
{
    const int flint_small = S->flint_small;
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
    fmpz * Acoeff = *A_coeff;
    ulong * Aexp = *A_exp;
    slong Aalloc = *A_alloc;
    ulong acc[3], p[3];
    int first_prod;

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

        _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, 1);

        Aexp[Alen] = exp;

        acc[0] = acc[1] = acc[2] = 0;
        first_prod = 1;
        while (heap_len > 1 && heap[1].exp == exp)
        {
            x = _mpoly_heap_pop1(heap, &heap_len, cmpmask);

            hind[x->i] |= WORD(1);
            *store++ = x->i;
            *store++ = x->j;

            if (flint_small)
            {
                smul_ppmm(p[1], p[0], Bcoeff[x->i], Ccoeff[x->j]);
                p[2] = FLINT_SIGN_EXT(p[1]);
                add_sssaaaaaa(acc[2], acc[1], acc[0], acc[2], acc[1], acc[0],
                                                         p[2], p[1], p[0]);
                first_prod = 0;
                while ((x = x->next) != NULL)
                {          
                    smul_ppmm(p[1], p[0], Bcoeff[x->i], Ccoeff[x->j]);
                    p[2] = FLINT_SIGN_EXT(p[1]);
                    add_sssaaaaaa(acc[2], acc[1], acc[0], acc[2], acc[1], acc[0],
                                                             p[2], p[1], p[0]);
                    hind[x->i] |= WORD(1);
                    *store++ = x->i;
                    *store++ = x->j;
                }
            }
            else /* output coeffs require multiprecision */
            {
                if (first_prod)
                    fmpz_mul(Acoeff + Alen, Bcoeff + x->i, Ccoeff + x->j);
                else
                    fmpz_addmul(Acoeff + Alen, Bcoeff + x->i, Ccoeff + x->j);
                first_prod = 0; 
                while ((x = x->next) != NULL)
                {
                    fmpz_addmul(Acoeff + Alen, Bcoeff + x->i, Ccoeff + x->j);
                    hind[x->i] |= WORD(1);
                    *store++ = x->i;
                    *store++ = x->j;
                }
            }
        }

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

        /* set output poly coeff from temporary accumulation, if not multiprec */
        if (flint_small)
        {
            fmpz_set_signed_uiuiui(Acoeff + Alen, acc[2], acc[1], acc[0]);
        }

        Alen += !fmpz_is_zero(Acoeff + Alen);
    }

    *A_coeff = Acoeff;
    *A_exp = Aexp;
    *A_alloc = Aalloc;
    return Alen;
}


slong _fmpz_mpoly_mul_heap_part(fmpz ** A_coeff, ulong ** A_exp, slong * A_alloc,
                 const fmpz * Bcoeff, const ulong * Bexp, slong Blen,
                 const fmpz * Ccoeff, const ulong * Cexp, slong Clen,
         slong * start, slong * end, slong * hind, const fmpz_mpoly_stripe_t S)
{
    const int flint_small = S->flint_small;
    mp_bitcnt_t bits = S->bits;
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
    slong * store, * store_base;
    slong Alen;
    ulong * Aexp = *A_exp;
    slong Aalloc = *A_alloc;
    fmpz * Acoeff = *A_coeff;
    ulong acc[3], p[3];
    int first_prod;

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
        exp_list[i] = exps + i*N;
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
                mpoly_monomial_add(exp_list[exp_next], Bexp + x->i*N,
                                                       Cexp + x->j*N, N);
            else
                mpoly_monomial_add_mp(exp_list[exp_next], Bexp + x->i*N,
                                                          Cexp + x->j*N, N);

            exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
        }
    }

    Alen = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, N);

        mpoly_monomial_set(Aexp + N*Alen, exp, N);

        acc[0] = acc[1] = acc[2] = 0;
        first_prod = 1;
        while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N))
        {
            exp_list[--exp_next] = heap[1].exp;

            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);

            hind[x->i] |= WORD(1);
            *store++ = x->i;
            *store++ = x->j;

            /* if output coeffs will fit in three words */
            if (flint_small)
            {
                smul_ppmm(p[1], p[0], Bcoeff[x->i], Ccoeff[x->j]);
                p[2] = FLINT_SIGN_EXT(p[1]);
                add_sssaaaaaa(acc[2], acc[1], acc[0], acc[2], acc[1], acc[0],
                                                         p[2], p[1], p[0]);
                first_prod = 0;
                while ((x = x->next) != NULL)
                {
                    smul_ppmm(p[1], p[0], Bcoeff[x->i], Ccoeff[x->j]);
                    p[2] = FLINT_SIGN_EXT(p[1]);
                    add_sssaaaaaa(acc[2], acc[1], acc[0], acc[2], acc[1], acc[0],
                                                             p[2], p[1], p[0]);
                    hind[x->i] |= WORD(1);
                    *store++ = x->i;
                    *store++ = x->j;
                }
            }
            else /* output coeffs require multiprecision */
            {
                if (first_prod)
                    fmpz_mul(Acoeff + Alen, Bcoeff + x->i, Ccoeff + x->j);
                else
                    fmpz_addmul(Acoeff + Alen, Bcoeff + x->i, Ccoeff + x->j);
                first_prod = 0;
                while ((x = x->next) != NULL)
                {
                    fmpz_addmul(Acoeff + Alen, Bcoeff + x->i, Ccoeff + x->j);
                    hind[x->i] |= WORD(1);
                    *store++ = x->i;
                    *store++ = x->j;
                }
            }
        }
      
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
                    mpoly_monomial_add(exp_list[exp_next], Bexp + x->i*N,
                                                           Cexp + x->j*N, N);
                else
                    mpoly_monomial_add_mp(exp_list[exp_next], Bexp + x->i*N,
                                                              Cexp + x->j*N, N);

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
                    mpoly_monomial_add(exp_list[exp_next], Bexp + x->i*N,
                                                           Cexp + x->j*N, N);
                else
                    mpoly_monomial_add_mp(exp_list[exp_next], Bexp + x->i*N,
                                                              Cexp + x->j*N, N);

                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
            }
        }     

        /* set output poly coeff from temporary accumulation, if not multiprec */
        if (flint_small)
        {
            fmpz_set_signed_uiuiui(Acoeff + Alen, acc[2], acc[1], acc[0]);
        }

        Alen += !fmpz_is_zero(Acoeff + Alen);
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
    volatile int idx;
    pthread_mutex_t mutex;
    slong nthreads;
    slong ndivs;
    const fmpz * Bcoeff;
    const ulong * Bexp;
    slong Blen;
    const fmpz * Ccoeff;
    const ulong * Cexp;
    slong Clen;
    slong N;
    mp_bitcnt_t bits;
    const ulong * cmpmask;
    int flint_small;
}
_base_struct;

typedef _base_struct _base_t[1];

typedef struct
{
    slong lower;
    slong upper;
    slong Alen;
    slong Aalloc;
    ulong * Aexp;
    fmpz * Acoeff;
}
_div_struct;

typedef struct
{
    fmpz_mpoly_stripe_t S;
    slong idx;
    slong time;
    _base_struct * base;
    _div_struct * divs;
    pthread_mutex_t mutex;
    pthread_cond_t cond;
    slong * t1, * t2, * t3, * t4;
    ulong * exp;
}
_worker_arg_struct;


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

void _fmpz_mpoly_mul_heap_threaded_worker(void * varg)
{
    _worker_arg_struct * arg = (_worker_arg_struct *) varg;
    fmpz_mpoly_stripe_struct * S = arg->S;
    _div_struct * divs = arg->divs;
    _base_struct * base = arg->base;
    slong Blen = base->Blen;
    slong N = base->N;
    slong i, j;
    ulong *exp;
    slong score;
    slong *start, *end, *t1, *t2, *t3, *t4, *tt;

    exp = (ulong *) flint_malloc(N*sizeof(ulong));
    t1 = (slong *) flint_malloc(Blen*sizeof(slong));
    t2 = (slong *) flint_malloc(Blen*sizeof(slong));
    t3 = (slong *) flint_malloc(Blen*sizeof(slong));
    t4 = (slong *) flint_malloc(Blen*sizeof(slong));

    S->N = N;
    S->bits = base->bits;
    S->cmpmask = base->cmpmask;
    S->flint_small = base->flint_small;

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
        if (N == 1)
        {
            divs[i].Alen = _fmpz_mpoly_mul_heap_part1(
                         &divs[i].Acoeff, &divs[i].Aexp, &divs[i].Aalloc,
                              base->Bcoeff, base->Bexp, base->Blen,
                              base->Ccoeff, base->Cexp, base->Clen,
                                                          start, end, t3, S);
        }
        else
        {
            divs[i].Alen = _fmpz_mpoly_mul_heap_part(
                         &divs[i].Acoeff, &divs[i].Aexp, &divs[i].Aalloc,
                              base->Bcoeff, base->Bexp, base->Blen,
                              base->Ccoeff, base->Cexp, base->Clen,
                                                          start, end, t3, S);
        }

        /* get next index to work on */
        pthread_mutex_lock(&base->mutex);
        i = base->idx - 1;
        base->idx = i;
        pthread_mutex_unlock(&base->mutex);
    }

    flint_free(S->big_mem);
    flint_free(t4);
    flint_free(t3);
    flint_free(t2);
    flint_free(t1);
    flint_free(exp);
}



slong _fmpz_mpoly_mul_heap_threaded(fmpz ** A_coeff, ulong ** A_exp, slong * A_alloc,
                 const fmpz * Bcoeff, const ulong * Bexp, slong Blen,
                 const fmpz * Ccoeff, const ulong * Cexp, slong Clen,
                              mp_bitcnt_t bits, slong N, const ulong * cmpmask)
{
    slong i, j, ndivs2;
    _base_t base;
    _div_struct * divs;
    _worker_arg_struct * args;
    slong Alen;
    fmpz * Acoeff;
    ulong * Aexp;
    slong max_num_workers, num_workers;
    thread_pool_handle * handles;

    /* bail here if no workers */
    FLINT_ASSERT(global_thread_pool_initialized);
    max_num_workers = thread_pool_get_size(global_thread_pool);
    max_num_workers = FLINT_MIN(max_num_workers, Clen/32);
    if (max_num_workers == 0)
    {
        return _fmpz_mpoly_mul_johnson(A_coeff, A_exp, A_alloc,
                                         Bcoeff, Bexp, Blen,
                                         Ccoeff, Cexp, Clen, bits, N, cmpmask);
    }
    handles = (thread_pool_handle *) flint_malloc(max_num_workers
                                                  *sizeof(thread_pool_handle));
    num_workers = thread_pool_request(global_thread_pool,
                                                     handles, max_num_workers);
    if (num_workers == 0)
    {
        flint_free(handles);
        return _fmpz_mpoly_mul_johnson(A_coeff, A_exp, A_alloc,
                                         Bcoeff, Bexp, Blen,
                                         Ccoeff, Cexp, Clen, bits, N, cmpmask);
    }

    base->nthreads = num_workers + 1;
    base->ndivs = base->nthreads*4;  /* number of divisons */
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
    base->flint_small =   _fmpz_mpoly_fits_small(Bcoeff, Blen)
                       && _fmpz_mpoly_fits_small(Ccoeff, Clen);

    ndivs2 = base->ndivs*base->ndivs;

    divs = (_div_struct *) flint_malloc(base->ndivs*sizeof(_div_struct));
    args = (_worker_arg_struct *) flint_malloc(base->nthreads
                                                  *sizeof(_worker_arg_struct));

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
        }
        else
        {
            /* lower divisions write to a new worker poly */
            divs[i].Aalloc = Blen + Clen/base->ndivs;
            divs[i].Aexp = (ulong *) flint_malloc(divs[i].Aalloc*N*sizeof(ulong)); 
            divs[i].Acoeff = (fmpz *) flint_calloc(divs[i].Aalloc, sizeof(fmpz));
        }
    }

    /* compute each chunk in parallel */
    pthread_mutex_init(&base->mutex, NULL);
    for (i = 0; i < num_workers; i++)
    {
        args[i].idx = i;
        args[i].base = base;
        args[i].divs = divs;
        thread_pool_wake(global_thread_pool, handles[i],
                               _fmpz_mpoly_mul_heap_threaded_worker, &args[i]);
    }
    i = num_workers;
    args[i].idx = i;
    args[i].base = base;
    args[i].divs = divs;
    _fmpz_mpoly_mul_heap_threaded_worker(&args[i]);
    for (i = 0; i < num_workers; i++)
    {
        thread_pool_wait(global_thread_pool, handles[i]);
        thread_pool_give_back(global_thread_pool, handles[i]);
    }
    pthread_mutex_destroy(&base->mutex);

    flint_free(handles);

    /* concatenate the outputs */

    i = base->ndivs - 1;
    Alen = divs[i].Alen;
    Acoeff = divs[i].Acoeff;
    Aexp = divs[i].Aexp;
    *A_alloc = divs[i].Aalloc;

    /* make space for all coeffs in one call to fit length */
    j = Alen;
    for (i = base->ndivs - 2; i >= 0; i--)
        j += divs[i].Alen;
    _fmpz_mpoly_fit_length(&Acoeff, &Aexp, A_alloc, j, N);

    for (i = base->ndivs - 2; i >= 0; i--)
    {
        /* transfer from worker poly to original poly */
        /* swap coeffs */
        for (j = 0; j < divs[i].Alen; j++)
        {
            fmpz_swap(Acoeff + Alen + j, divs[i].Acoeff + j);
            fmpz_clear(divs[i].Acoeff + j);
        }
        /* clear remaining coeffs */
        for ( ; j < divs[i].Aalloc; j++)
        {
            fmpz_clear(divs[i].Acoeff + j);
        }
        /* copy exps */
        flint_mpn_copyi(Aexp + N*Alen, divs[i].Aexp, N*divs[i].Alen);

        Alen += divs[i].Alen;
        flint_free(divs[i].Acoeff);
        flint_free(divs[i].Aexp);
    }

    flint_free(args);
    flint_free(divs);

    *A_coeff = Acoeff;
    *A_exp = Aexp;
    return Alen;
}



/* maxBfields gets clobbered */
void _fmpz_mpoly_mul_heap_threaded_maxfields(fmpz_mpoly_t A,
                                 const fmpz_mpoly_t B, fmpz * maxBfields,
                                 const fmpz_mpoly_t C, fmpz * maxCfields,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong N, Alen;
    mp_bitcnt_t exp_bits;
    ulong * cmpmask;
    ulong * Bexp, * Cexp;
    int freeBexp, freeCexp;
    TMP_INIT;

    if (!global_thread_pool_initialized)
    {
        _fmpz_mpoly_mul_johnson_maxfields(A, B, maxBfields, C, maxCfields, ctx);
        return;
    }

    TMP_START;

    _fmpz_vec_add(maxBfields, maxBfields, maxCfields, ctx->minfo->nfields);

    exp_bits = _fmpz_vec_max_bits(maxBfields, ctx->minfo->nfields);
    exp_bits = FLINT_MAX(MPOLY_MIN_BITS, exp_bits + 1);
    exp_bits = FLINT_MAX(exp_bits, B->bits);
    exp_bits = FLINT_MAX(exp_bits, C->bits);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);

    N = mpoly_words_per_exp(exp_bits, ctx->minfo);
    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, exp_bits, ctx->minfo);

    /* ensure input exponents are packed into same sized fields as output */
    freeBexp = 0;
    Bexp = B->exps;
    if (exp_bits > B->bits)
    {
        freeBexp = 1;
        Bexp = (ulong *) flint_malloc(N*B->length*sizeof(ulong));
        mpoly_repack_monomials(Bexp, exp_bits, B->exps, B->bits,
                                                        B->length, ctx->minfo);
    }

    freeCexp = 0;
    Cexp = C->exps;
    if (exp_bits > C->bits)
    {
        freeCexp = 1;
        Cexp = (ulong *) flint_malloc(N*C->length*sizeof(ulong));
        mpoly_repack_monomials(Cexp, exp_bits, C->exps, C->bits,
                                                        C->length, ctx->minfo);
    }

    /* deal with aliasing and do multiplication */
    if (A == B || A == C)
    {
        fmpz_mpoly_t T;
        fmpz_mpoly_init2(T, B->length + C->length - 1, ctx);
        fmpz_mpoly_fit_bits(T, exp_bits, ctx);
        T->bits = exp_bits;

        /* algorithm more efficient if smaller poly first */
        if (B->length >= C->length)
        {
            Alen = _fmpz_mpoly_mul_heap_threaded(&T->coeffs, &T->exps, &T->alloc,
                                                  C->coeffs, Cexp, C->length,
                                                  B->coeffs, Bexp, B->length,
                                                         exp_bits, N, cmpmask);
        }
        else
        {
            Alen = _fmpz_mpoly_mul_heap_threaded(&T->coeffs, &T->exps, &T->alloc,
                                                  B->coeffs, Bexp, B->length,
                                                  C->coeffs, Cexp, C->length,
                                                         exp_bits, N, cmpmask);
        }

        fmpz_mpoly_swap(T, A, ctx);
        fmpz_mpoly_clear(T, ctx);
    }
    else
    {
        fmpz_mpoly_fit_length(A, B->length + C->length - 1, ctx);
        fmpz_mpoly_fit_bits(A, exp_bits, ctx);
        A->bits = exp_bits;

        /* algorithm more efficient if smaller poly first */
        if (B->length > C->length)
        {
            Alen = _fmpz_mpoly_mul_heap_threaded(&A->coeffs, &A->exps, &A->alloc,
                                                  C->coeffs, Cexp, C->length,
                                                  B->coeffs, Bexp, B->length,
                                                         exp_bits, N, cmpmask);
        }
        else
        {
            Alen = _fmpz_mpoly_mul_heap_threaded(&A->coeffs, &A->exps, &A->alloc,
                                                  B->coeffs, Bexp, B->length,
                                                  C->coeffs, Cexp, C->length,
                                                         exp_bits, N, cmpmask);
        }
    }

    if (freeBexp)
        flint_free(Bexp);

    if (freeCexp)
        flint_free(Cexp);

    _fmpz_mpoly_set_length(A, Alen, ctx);

    TMP_END;
}


void fmpz_mpoly_mul_heap_threaded(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                              const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    fmpz * maxBfields, * maxCfields;
    TMP_INIT;

    if (B->length == 0 || C->length == 0)
    {
        fmpz_mpoly_zero(A, ctx);
        return;
    }

    TMP_START;

    maxBfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    maxCfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_init(maxBfields + i);
        fmpz_init(maxCfields + i);
    }
    mpoly_max_fields_fmpz(maxBfields, B->exps, B->length, B->bits, ctx->minfo);
    mpoly_max_fields_fmpz(maxCfields, C->exps, C->length, C->bits, ctx->minfo);

    _fmpz_mpoly_mul_heap_threaded_maxfields(A, B, maxBfields, C, maxCfields, ctx);

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_clear(maxBfields + i);
        fmpz_clear(maxCfields + i);
    }

    TMP_END;
}
