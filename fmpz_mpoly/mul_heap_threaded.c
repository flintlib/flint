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
#include <pthread.h>

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
    slong next_loc = Blen + 4;   /* something bigger than heap can ever be */
    slong Q_len = 0, heap_len = 1; /* heap zero index unused */
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * Q;
    mpoly_heap_t * x;
    slong Alen;
    fmpz * Acoeff = *A_coeff;
    ulong * Aexp = *A_exp;
    slong Aalloc = *A_alloc;
    ulong exp, cy;
    ulong c[3], p[2]; /* for accumulating coefficients */
    int first;

    i = 0;
    Q = (slong *) (S->big_mem + i);
    i += 2*Blen*sizeof(slong);
    heap = (mpoly_heap1_s *)(S->big_mem + i);
    i += (Blen + 1)*sizeof(mpoly_heap1_s);
    chain = (mpoly_heap_t *)(S->big_mem + i);
    i += Blen*sizeof(mpoly_heap_t);
    FLINT_ASSERT(i <= S->big_mem_alloc);


    /* put all the starting nodes on the heap */
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
            hind[x->i] = 2*(x->j+1) + 0;
            _mpoly_heap_insert1(heap, Bexp[x->i] + Cexp[x->j], x,
                                                &next_loc, &heap_len, cmpmask);
        }
    }

    Alen = -WORD(1);
    while (heap_len > 1)
    {
        exp = heap[1].exp;
      
        Alen++;

        _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, 1);

        /* whether we are on first coeff product for this output exponent */
        first = 1;

        c[0] = c[1] = c[2] = 0;

        /* while heap nonempty and contains chain with current output exponent */
        while (heap_len > 1 && heap[1].exp == exp)
        {
            /* pop chain from heap */
            x = _mpoly_heap_pop1(heap, &heap_len, cmpmask);

            /* temporarily store indicies of this node */
            hind[x->i] |= WORD(1);
            Q[Q_len++] = x->i;
            Q[Q_len++] = x->j;

            /* if output coeffs will fit in three words */
            if (flint_small)
            {
                /* compute product of input poly coeffs */
                if (first)
                {
                    smul_ppmm(c[1], c[0], Bcoeff[x->i], Ccoeff[x->j]);
                    c[2] = -(c[1] >> (FLINT_BITS - 1));

                    /* set output monomial */
                    Aexp[Alen] = exp;
                    first = 0; 
                } else /* addmul product of input poly coeffs */
                {
                    smul_ppmm(p[1], p[0], Bcoeff[x->i], Ccoeff[x->j]);
                    add_sssaaaaaa(cy, c[1], c[0], 0, c[1], c[0], 0, p[1], p[0]);
                    c[2] += (0 <= (slong) p[1]) ? cy : cy - 1;
                }

                /* for every node in this chain */
                while ((x = x->next) != NULL)
                {          
                    /* addmul product of input poly coeffs */
                    smul_ppmm(p[1], p[0], Bcoeff[x->i], Ccoeff[x->j]);
                    add_sssaaaaaa(cy, c[1], c[0], 0, c[1], c[0], 0, p[1], p[0]);
                    c[2] += (0 <= (slong) p[1]) ? cy : cy - 1;

                    /* take node out of heap and put into store */
                    hind[x->i] |= WORD(1);
                    Q[Q_len++] = x->i;
                    Q[Q_len++] = x->j;
                }
           }
            else /* output coeffs require multiprecision */
           {
                if (first) /* compute product of input poly coeffs */
                {
                    fmpz_mul(Acoeff + Alen, Bcoeff + x->i, Ccoeff + x->j);
               
                    Aexp[Alen] = exp;
                    first = 0; 
                } else
                {
                    fmpz_addmul(Acoeff + Alen, Bcoeff + x->i, Ccoeff + x->j);
                }

                /* for each node in this chain */
                while ((x = x->next) != NULL)
                {
                    fmpz_addmul(Acoeff + Alen, Bcoeff + x->i, Ccoeff + x->j);

                    hind[x->i] |= WORD(1);
                    Q[Q_len++] = x->i;
                    Q[Q_len++] = x->j;
                }
            }
        }

        /* for each node temporarily stored */
        while (Q_len > 0)
        {
            /* take node from store */
            j = Q[--Q_len];
            i = Q[--Q_len];

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
                  || (hind[i - 1] >  2*(j + 2) + 1)
                  || (hind[i - 1] == 2*(j + 2) + 1) /* gcc should fuse */
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
            fmpz_set_signed_uiuiui(Acoeff + Alen, c[2], c[1], c[0]);

        if (fmpz_is_zero(Acoeff + Alen))
            Alen--;
    }

    Alen++;

    *A_coeff = Acoeff;
    *A_exp = Aexp;
    *A_alloc = Aalloc;

    return Alen;
}



slong _fmpz_mpoly_mul_heap_part(fmpz ** poly1, ulong ** exp1, slong * alloc,
                 const fmpz * poly2, const ulong * exp2, slong len2,
                 const fmpz * poly3, const ulong * exp3, slong len3,
         slong * start, slong * end, slong * hind, const fmpz_mpoly_stripe_t S)
{
    const int flint_small = S->flint_small;
    mp_bitcnt_t bits = S->bits;
    slong N = S->N;
    const ulong * cmpmask = S->cmpmask;
    slong i, j, k;
    slong next_loc = len2 + 4;   /* something bigger than heap can ever be */
    slong Q_len = 0, heap_len = 1; /* heap zero index unused */
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * Q;
    mpoly_heap_t * x;
    fmpz * p1 = *poly1;
    ulong * e1 = *exp1;
    ulong cy;
    ulong c[3], p[2]; /* for accumulating coefficients */
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    int first;
    slong Blen = len2;


    i = 0;

/*    store = store_base = (slong *) (S->big_mem + i);*/
    Q = (slong *) (S->big_mem + i);

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

    for (i = 0; i < len2; i++)
        exp_list[i] = exps + i*N;

    /* heap indices */
    for (i = 0; i < len2; i++)
        hind[i] = 2*start[i] + 1;

    /* start with no heap nodes and no exponent vectors in use */
    exp_next = 0;

    /* put all the starting nodes on the heap */
    for (i = 0; i < len2; i++)
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
                mpoly_monomial_add(exp_list[exp_next], exp2 + x->i*N, exp3 + x->j*N, N);
            else
                mpoly_monomial_add_mp(exp_list[exp_next], exp2 + x->i*N, exp3 + x->j*N, N);

            if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                      &next_loc, &heap_len, N, cmpmask))
               exp_next--;
        }
    }

    /* output poly index starts at -1, will be immediately updated to 0 */
    k = -WORD(1);

    /* while heap is nonempty */
    while (heap_len > 1)
    {
        /* get pointer to exponent field of heap top */
        exp = heap[1].exp;

        /* realloc output poly ready for next product term */
        k++;
        _fmpz_mpoly_fit_length(&p1, &e1, alloc, k + 1, N);

        /* whether we are on first coeff product for this output exponent */
        first = 1;

        /* set temporary coeff to zero */
        c[0] = c[1] = c[2] = 0;

        /* while heap nonempty and contains chain with current output exponent */
        while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N))
        {
            /* pop chain from heap and set exponent field to be reused */
            exp_list[--exp_next] = heap[1].exp;

            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);

            /* take node out of heap and put into store */
            hind[x->i] |= WORD(1);
            Q[Q_len++] = x->i;
            Q[Q_len++] = x->j;

            /* if output coeffs will fit in three words */
            if (flint_small)
            {
                /* compute product of input poly coeffs */
                if (first)
                {
                    smul_ppmm(c[1], c[0], poly2[x->i], poly3[x->j]);
                    c[2] = -(c[1] >> (FLINT_BITS - 1));

                    /* set output monomial */
                    mpoly_monomial_set(e1 + k*N, exp, N);

                    first = 0; 
                } else /* addmul product of input poly coeffs */
                {
                    smul_ppmm(p[1], p[0], poly2[x->i], poly3[x->j]);
                    add_sssaaaaaa(cy, c[1], c[0], 0, c[1], c[0], 0, p[1], p[0]);
                    c[2] += (0 <= (slong) p[1]) ? cy : cy - 1;
                }

                /* for every node in this chain */
                while ((x = x->next) != NULL)
                {
                    /* addmul product of input poly coeffs */
                    smul_ppmm(p[1], p[0], poly2[x->i], poly3[x->j]);
                    add_sssaaaaaa(cy, c[1], c[0], 0, c[1], c[0], 0, p[1], p[0]);
                    c[2] += (0 <= (slong) p[1]) ? cy : cy - 1;

                    /* take node out of heap and put into store */
                    hind[x->i] |= WORD(1);
                    Q[Q_len++] = x->i;
                    Q[Q_len++] = x->j;
                }
            } else /* output coeffs require multiprecision */
            {
                if (first) /* compute product of input poly coeffs */
                {
                    fmpz_mul(p1 + k, poly2 + x->i, poly3 + x->j);

                    /* set output monomial */
                    mpoly_monomial_set(e1 + k*N, exp, N);

                    first = 0; 
                } else
                {   /* addmul product of input poly coeffs */
                    fmpz_addmul(p1 + k, poly2 + x->i, poly3 + x->j);
                }

                /* for each node in this chain */
                while ((x = x->next) != NULL)
                {
                    /* addmul product of input poly coeffs */
                    fmpz_addmul(p1 + k, poly2 + x->i, poly3 + x->j);

                    /* take node out of heap and put into store */
                    hind[x->i] |= WORD(1);
                    Q[Q_len++] = x->i;
                    Q[Q_len++] = x->j;
                }
            }
        }
      
        /* for each node temporarily stored */
        while (Q_len > 0)
        {
            /* take node from store */
            j = Q[--Q_len];
            i = Q[--Q_len];

            /* should we go right? */
            if (  (i + 1 < len2)
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
                    mpoly_monomial_add(exp_list[exp_next], exp2 + x->i*N, exp3 + x->j*N, N);
                else
                    mpoly_monomial_add_mp(exp_list[exp_next], exp2 + x->i*N, exp3 + x->j*N, N);

                if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                      &next_loc, &heap_len, N, cmpmask))
                    exp_next--;
            }

            /* should we go up? */
            if (  (j + 1 < end[i + 0])
               && ((hind[i] & 1) == 1)
               && (  (i == 0)
                  || (hind[i - 1] >  2*(j + 2) + 1)
                  || (hind[i - 1] == 2*(j + 2) + 1) /* gcc should fuse */
                  )
               )
            {
                x = chain + i;
                x->i = i;
                x->j = j + 1;
                x->next = NULL;

                hind[x->i] = 2*(x->j + 1) + 0;

                if (bits <= FLINT_BITS)
                    mpoly_monomial_add(exp_list[exp_next], exp2 + x->i*N, exp3 + x->j*N, N);
                else
                    mpoly_monomial_add_mp(exp_list[exp_next], exp2 + x->i*N, exp3 + x->j*N, N);

                if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                      &next_loc, &heap_len, N, cmpmask))
                    exp_next--;
            }
        }     

        /* set output poly coeff from temporary accumulation, if not multiprec */
        if (flint_small)
            fmpz_set_signed_uiuiui(p1 + k, c[2], c[1], c[0]);

        if (fmpz_is_zero(p1 + k))
            k--;

    }

    k++;

    (*poly1) = p1;
    (*exp1) = e1;

    return k;
}


/*
    The workers calculate product terms from 4*n divisions, where n is the
    number of threads. 
    This contains the address of a mul_heap_threaded_div_t structure
*/

typedef struct
{
    volatile int idx;
    pthread_mutex_t mutex;
    slong nthreads;
    slong ndivs;
    const fmpz * coeff2; const ulong * exp2; slong len2;
    const fmpz * coeff3; const ulong * exp3; slong len3;
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
    slong len1;
    slong alloc1;
    ulong * exp1;
    fmpz * coeff1;
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
    slong Blen = base->len2;
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
    } else
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
                            base->exp2, base->len2, base->exp3, base->len3,
                                          base->N, base->cmpmask);
            if (start == t2)
            {
                SWAP_PTRS(t1, t2);
            } else if (start == t3)
            {
                SWAP_PTRS(t1, t3);
            }

        } else {
            start = t1;
            for (j = 0; j < base->len2; j++)
                start[j] = 0;
        }

        /* calculate end */
        if (i > 0) {
            mpoly_search_monomials(
                &end, exp, &score, t2, t3, t4,
                            divs[i - 1].lower, divs[i - 1].lower,
                            base->exp2, base->len2, base->exp3, base->len3,
                                          base->N, base->cmpmask);
            if (end == t3)
            {
                SWAP_PTRS(t2, t3);
            } else if (end == t4)
            {
                SWAP_PTRS(t2, t4);
            }

        } else {
            end = t2;
            for (j = 0; j < base->len2; j++)
                end[j] = base->len3;
        }
        /* t3 and t4 are free for workspace at this point */

        /* calculate products in [start, end) */
        if (N == 1)
        {
            divs[i].len1 = _fmpz_mpoly_mul_heap_part1(
                         &divs[i].coeff1, &divs[i].exp1, &divs[i].alloc1,
                          base->coeff2,  base->exp2,  base->len2,
                          base->coeff3,  base->exp3,  base->len3,
                                                          start, end, t3, S);
        }
        else
        {
            divs[i].len1 = _fmpz_mpoly_mul_heap_part(
                         &divs[i].coeff1, &divs[i].exp1, &divs[i].alloc1,
                          base->coeff2,  base->exp2,  base->len2,
                          base->coeff3,  base->exp3,  base->len3,
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



slong _fmpz_mpoly_mul_heap_threaded(fmpz ** poly1, ulong ** exp1, slong * alloc,
                 const fmpz * coeff2, const ulong * exp2, slong len2,
                 const fmpz * coeff3, const ulong * exp3, slong len3,
                              mp_bitcnt_t bits, slong N, const ulong * cmpmask)
{
    slong i, j, k, ndivs2;
    _worker_arg_struct * args;
    _base_t base;
    _div_struct * divs;
    fmpz * p1;
    ulong * e1;
    slong max_num_workers, num_workers;
    thread_pool_handle * handles;

    /* bail here if no workers */
    FLINT_ASSERT(global_thread_pool_initialized);
    max_num_workers = thread_pool_get_size(global_thread_pool);
    max_num_workers = FLINT_MIN(max_num_workers, len3/32);
    if (max_num_workers == 0)
    {
        return _fmpz_mpoly_mul_johnson(poly1, exp1, alloc, coeff2, exp2, len2,
                                         coeff3, exp3, len3, bits, N, cmpmask);
    }
    handles = (thread_pool_handle *) flint_malloc(max_num_workers
                                                  *sizeof(thread_pool_handle));
    num_workers = thread_pool_request(global_thread_pool,
                                                     handles, max_num_workers);
    if (num_workers == 0)
    {
        flint_free(handles);
        return _fmpz_mpoly_mul_johnson(poly1, exp1, alloc, coeff2, exp2, len2,
                                         coeff3, exp3, len3, bits, N, cmpmask);
    }

    base->nthreads = num_workers + 1;
    base->ndivs    = base->nthreads*4;  /* number of divisons */
    base->coeff2 = coeff2;
    base->exp2 = exp2;
    base->len2 = len2;
    base->coeff3 = coeff3;
    base->exp3 = exp3;
    base->len3 = len3;
    base->bits = bits;
    base->N = N;
    base->cmpmask = cmpmask;
    base->idx = base->ndivs - 1;    /* decremented by worker threads */
    base->flint_small =    _fmpz_mpoly_fits_small(coeff2, len2)
                        && _fmpz_mpoly_fits_small(coeff3, len3);

    ndivs2 = base->ndivs*base->ndivs;

    divs = (_div_struct *) flint_malloc(base->ndivs*sizeof(_div_struct));
    args = (_worker_arg_struct *) flint_malloc(base->nthreads
                                                  *sizeof(_worker_arg_struct));

    /* allocate space and set the boundary for each division */
    for (i = base->ndivs - 1; i >= 0; i--)
    {
        /* divisions decrease in size so that no worker finishes too early */
        divs[i].lower = (i + 1)*(i + 1)*len2*len3/ndivs2;
        divs[i].upper = divs[i].lower;

        divs[i].len1 = 0;
        if (i == base->ndivs - 1)
        {
            /* highest division writes to original poly */
            divs[i].alloc1 = *alloc;
            divs[i].exp1 = *exp1;
            divs[i].coeff1 = *poly1;
        } else
        {
            /* lower divisions write to a new worker poly */
            divs[i].alloc1 = len2 + len3/base->ndivs;
            divs[i].exp1 = (ulong *) flint_malloc(divs[i].alloc1*N*sizeof(ulong)); 
            divs[i].coeff1 = (fmpz *) flint_calloc(divs[i].alloc1, sizeof(fmpz));
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
    k = divs[i].len1;
    p1 = divs[i].coeff1;
    e1 = divs[i].exp1;
    *alloc = divs[i].alloc1;

    /* make space for all coeffs in one call to fit length */
    j = k;
    for (i = base->ndivs - 2; i >= 0; i--)
        j += divs[i].len1;
    _fmpz_mpoly_fit_length(&p1, &e1, alloc, j, N);

    for (i = base->ndivs - 2; i >= 0; i--)
    {
        /* transfer from worker poly to original poly */
        /* swap coeffs */
        for (j = 0; j < divs[i].len1; j++)
        {
            fmpz_swap(p1 + k + j, divs[i].coeff1 + j);
            fmpz_clear(divs[i].coeff1 + j);
        }
        /* clear remaining coeffs */
        for ( ; j < divs[i].alloc1; j++)
        {
            fmpz_clear(divs[i].coeff1 + j);
        }
        /* copy exps */
        flint_mpn_copyi(e1 + N*k, divs[i].exp1, N*divs[i].len1);

        k += divs[i].len1;
        flint_free(divs[i].coeff1);
        flint_free(divs[i].exp1);
    }

    flint_free(args);
    flint_free(divs);

    *poly1 = p1;
    *exp1  = e1;
    return k;
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
