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
#include "nmod_mpoly.h"


/*
    Set poly1 to poly2*poly3 using Johnson's heap method. The function
    realocates its output and returns the length of the product. This
    version of the function assumes the exponent vectors all fit in a
    single word. Assumes input polys are nonzero.
    Only terms t with start >= t > end are written;
    "start" and "end" are not monomials but arrays of indicies into exp3
*/
slong _nmod_mpoly_mul_heap_part1(mp_limb_t ** coeff1, ulong ** exp1, slong * alloc,
              const mp_limb_t * coeff2, const ulong * exp2, slong len2,
              const mp_limb_t * coeff3, const ulong * exp3, slong len3,
               slong * start, slong * end, slong * hind, ulong maskhi,
                                                        const nmodf_ctx_t fctx)
{
    slong i, j, len1;
    slong next_loc = len2 + 4;   /* something bigger than heap can ever be */
    slong Q_len = 0, heap_len = 1; /* heap zero index unused */
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * Q;
    mpoly_heap_t * x;
    mp_limb_t * p1 = * coeff1;
    ulong * e1 = *exp1;
/*  slong * hind;   */
    ulong exp;
    mp_limb_t acc0, acc1, acc2, pp0, pp1;
    TMP_INIT;

    TMP_START;

    next_loc = len2 + 4;   /* something bigger than heap can ever be */
    heap = (mpoly_heap1_s *) TMP_ALLOC((len2 + 1)*sizeof(mpoly_heap1_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(len2*sizeof(mpoly_heap_t));
    Q = (slong *) TMP_ALLOC(2*len2*sizeof(slong));
   
    /* space for heap indices */
/*  hind = (slong *) TMP_ALLOC(len2*sizeof(slong)); */
    for (i = 0; i < len2; i++)
        hind[i] = 2*start[i] + 1;

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
            hind[x->i] = 2*(x->j+1) + 0;
            _mpoly_heap_insert1(heap, exp2[x->i] + exp3[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
        }
    }

    len1 = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;
      
        _nmod_mpoly_fit_length(&p1, &e1, alloc, len1 + 1, 1);

        e1[len1] = exp;

        acc0 = acc1 = acc2 = 0;
        do
        {
            x = _mpoly_heap_pop1(heap, &heap_len, maskhi);

            hind[x->i] |= WORD(1);
            Q[Q_len++] = x->i;
            Q[Q_len++] = x->j;
            umul_ppmm(pp1, pp0, coeff2[x->i], coeff3[x->j]);
            add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);

            while ((x = x->next) != NULL)
            {
                hind[x->i] |= WORD(1);
                Q[Q_len++] = x->i;
                Q[Q_len++] = x->j;
                umul_ppmm(pp1, pp0, coeff2[x->i], coeff3[x->j]);
                add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);
            }
        } while (heap_len > 1 && heap[1].exp == exp);

        NMOD_RED3(p1[len1], acc2, acc1, acc0, fctx->mod);
        len1 += (p1[len1] != 0);

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

                hind[x->i] = 2*(x->j+1) + 0;
                _mpoly_heap_insert1(heap, exp2[x->i] + exp3[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
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
                _mpoly_heap_insert1(heap, exp2[x->i] + exp3[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
            }
        }
    }

    (*coeff1) = p1;
    (*exp1) = e1;
   
    TMP_END;

    return len1;
}


/*
    Set poly1 to poly2*poly3 using Johnson's heap method. The function
    realocates its output and returns the length of the product. This
    version of the function assumes the exponent vectors take N words.
    Only terms t with start >= t > end are written;
    "start" and "end" are not monomials but arrays of indicies into exp3
*/
slong _nmod_mpoly_mul_heap_part(mp_limb_t ** coeff1, ulong ** exp1, slong * alloc,
                 const mp_limb_t * coeff2, const ulong * exp2, slong len2,
                 const mp_limb_t * coeff3, const ulong * exp3, slong len3,
                  slong * start, slong * end, slong * hind,
      mp_bitcnt_t bits, slong N, const ulong * cmpmask, const nmodf_ctx_t fctx)
{
    slong i, j, len1;
    slong next_loc = len2 + 4;   /* something bigger than heap can ever be */
    slong Q_len = 0, heap_len = 1; /* heap zero index unused */
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * Q;
    mpoly_heap_t * x;
    mp_limb_t * p1 = * coeff1;
    ulong * e1 = *exp1;
    ulong acc0, acc1, acc2, pp0, pp1;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    TMP_INIT;

    /* if exponent vectors fit in single word, call special version */
    if (N == 1)
        return _nmod_mpoly_mul_heap_part1(coeff1, exp1, alloc,
                                          coeff2, exp2, len2,
                                          coeff3, exp3, len3,
                                           start, end, hind, cmpmask[0], fctx);

    TMP_START;

    heap = (mpoly_heap_s *) TMP_ALLOC((len2 + 1)*sizeof(mpoly_heap_s));
    chain = (mpoly_heap_t *) TMP_ALLOC(len2*sizeof(mpoly_heap_t));
    Q = (slong *) TMP_ALLOC(2*len2*sizeof(slong));
    exps = (ulong *) TMP_ALLOC(len2*N*sizeof(ulong));
    exp_list = (ulong **) TMP_ALLOC(len2*sizeof(ulong *));
    exp_next = 0;
    for (i = 0; i < len2; i++)
        exp_list[i] = exps + i*N;

    /* heap indices */
    for (i = 0; i < len2; i++)
        hind[i] = 2*start[i] + 1;

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

    len1 = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        _nmod_mpoly_fit_length(&p1, &e1, alloc, len1 + 1, N);

        mpoly_monomial_set(e1 + len1*N, exp, N);

        acc0 = acc1 = acc2 = 0;
        do
        {
            exp_list[--exp_next] = heap[1].exp;

            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);

            hind[x->i] |= WORD(1);
            Q[Q_len++] = x->i;
            Q[Q_len++] = x->j;
            umul_ppmm(pp1, pp0, coeff2[x->i], coeff3[x->j]);
            add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);

            while ((x = x->next) != NULL)
            {
                hind[x->i] |= WORD(1);
                Q[Q_len++] = x->i;
                Q[Q_len++] = x->j;
                umul_ppmm(pp1, pp0, coeff2[x->i], coeff3[x->j]);
                add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);
            }
        } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

        NMOD_RED3(p1[len1], acc2, acc1, acc0, fctx->mod);
        len1 += (p1[len1] != 0);

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
                    mpoly_monomial_add(exp_list[exp_next], exp2 + x->i*N, exp3 + x->j*N, N);
                else
                    mpoly_monomial_add_mp(exp_list[exp_next], exp2 + x->i*N, exp3 + x->j*N, N);

                if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                      &next_loc, &heap_len, N, cmpmask))
                    exp_next--;
            }
        }
    }

    (*coeff1) = p1;
    (*exp1) = e1;

    TMP_END;

    return len1;
}


/*
    The workers calculate product terms from 4*n divisions, where n is the
    number of threads. 
    This contains the address of a mul_heap_threaded_div_t structure
*/

typedef struct
{
    pthread_mutex_t mutex;
    const nmodf_ctx_struct * fctx;
    slong nthreads;
    slong ndivs;
    const mp_limb_t * coeff2; const ulong * exp2; slong len2;
    const mp_limb_t * coeff3; const ulong * exp3; slong len3;
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
    slong len1;
    slong alloc1;
    ulong * exp1;
    mp_limb_t * coeff1;
}
mul_heap_threaded_div_t;

typedef struct
{
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

void * _nmod_mpoly_mul_heap_threaded_worker(void * arg_ptr)
{
    mul_heap_threaded_arg_t * arg = (mul_heap_threaded_arg_t *) arg_ptr;

    mul_heap_threaded_div_t * divs;
    mul_heap_threaded_base_t * base;
    slong i, j;

    ulong *exp;
    slong score;
    slong *start, *end, *t1, *t2, *t3, *t4, *tt;

    divs = arg->divp;
    base = arg->basep;
    exp = (ulong *) flint_malloc(base->N*sizeof(ulong));
    t1 = (slong *) flint_malloc(base->len2*sizeof(slong));
    t2 = (slong *) flint_malloc(base->len2*sizeof(slong));
    t3 = (slong *) flint_malloc(base->len2*sizeof(slong));
    t4 = (slong *) flint_malloc(base->len2*sizeof(slong));

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


    while (i >= 0) {

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
        divs[i].len1 = _nmod_mpoly_mul_heap_part(
                     &divs[i].coeff1, &divs[i].exp1, &divs[i].alloc1,
                      base->coeff2,  base->exp2,  base->len2,
                      base->coeff3,  base->exp3,  base->len3,
                       start, end, t3, base->bits, base->N,
                                       base->cmpmask, base->fctx);

        /* get next index to work on */
        pthread_mutex_lock(&base->mutex);
        i = base->idx - 1;
        base->idx = i;
        pthread_mutex_unlock(&base->mutex);
    }

    /* clean up */
    flint_free(t4);
    flint_free(t3);
    flint_free(t2);
    flint_free(t1);
    flint_free(exp);
    if (arg->idx + 1 < base->nthreads)
    {
        flint_cleanup();
    }
    return NULL;
}



slong _nmod_mpoly_mul_heap_threaded(mp_limb_t ** coeff1, ulong ** exp1, slong * alloc,
                 const mp_limb_t * coeff2, const ulong * exp2, slong len2,
                 const mp_limb_t * coeff3, const ulong * exp3, slong len3,
      mp_bitcnt_t bits, slong N, const ulong * cmpmask, const nmodf_ctx_t fctx)
{
    slong i, j, k, ndivs2;
    pthread_t * threads;
    mul_heap_threaded_arg_t * args;
    mul_heap_threaded_base_t * base;
    mul_heap_threaded_div_t * divs;
    mp_limb_t * p1;
    ulong * e1;

    base = flint_malloc(sizeof(mul_heap_threaded_base_t));
    base->nthreads = flint_get_num_threads();
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
    base->idx = base->ndivs-1;    /* decremented by worker threads */
    base->fctx = fctx;

    ndivs2 = base->ndivs*base->ndivs;

    divs    = flint_malloc(sizeof(mul_heap_threaded_div_t) * base->ndivs);
    threads = flint_malloc(sizeof(pthread_t) * base->nthreads);
    args    = flint_malloc(sizeof(mul_heap_threaded_arg_t) * base->nthreads);

    /* allocate space and set the boundary for each division */
    for (i = base->ndivs - 1; i >= 0; i--)
    {
        /* divisions decrease in size so that no worker finishes too early */
        divs[i].lower = (i + 1)*(i + 1)*len2*len3/ndivs2;
        divs[i].upper = divs[i].lower;

        divs[i].len1 = 0;
        k = 0; /* avoid bogus warning */
        if (i == base->ndivs - 1)
        {
            /* highest division writes to original poly */
            divs[i].alloc1 = *alloc;
            divs[i].exp1 = *exp1;
            divs[i].coeff1 = *coeff1;
        } else
        {
            /* lower divisions write to a new worker poly */
            divs[i].alloc1 = len2 + len3/base->ndivs;
            divs[i].exp1 = (ulong *) flint_malloc(divs[i].alloc1*N*sizeof(ulong)); 
            divs[i].coeff1 = (ulong *) flint_malloc(divs[i].alloc1*sizeof(ulong));
        }
    }

    /* compute each chunk in parallel */
    pthread_mutex_init(&base->mutex, NULL);
    for (i = 0; i < base->nthreads; i++)
    {
        args[i].idx = i;
        args[i].basep = base;
        args[i].divp = divs;
        if (i + 1 < base->nthreads)
        {
            pthread_create(&threads[i], NULL, _nmod_mpoly_mul_heap_threaded_worker, &args[i]);
        } else
        {
            _nmod_mpoly_mul_heap_threaded_worker(&args[i]);
        }
    }
    for (i = base->nthreads-1; i >= 0; i--)
    {
        if (i + 1 < base->nthreads)
        {
            pthread_join(threads[i], NULL);
        }
    }
    pthread_mutex_destroy(&base->mutex);

    /* concatenate the outputs */ 
    k = 0; /* avoid bogus warning */
    for (i = base->ndivs - 1; i >= 0; i--)
    {
        if (i + 1 < base->ndivs)
        {
            /* transfer from worker poly to original poly */
            for (j = 0; j < divs[i].len1; j++)
            {
                _nmod_mpoly_fit_length(&p1, &e1, alloc, k + 1, N);
                p1[k] = divs[i].coeff1[j];
                mpoly_monomial_set(e1 + N*k, divs[i].exp1 + N*j, N);
                k++;
            }
            flint_free(divs[i].coeff1);
            flint_free(divs[i].exp1);
        } else
        {
            /* highest thread used original poly */
            k = divs[i].len1;
            p1 = divs[i].coeff1;
            e1 = divs[i].exp1;
            *alloc = divs[i].alloc1;
        }
    }

    flint_free(args);
    flint_free(threads);
    flint_free(divs);
    flint_free(base);

    *coeff1 = p1;
    *exp1  = e1;
    return k;
}

void nmod_mpoly_mul_heap_threaded(nmod_mpoly_t poly1, const nmod_mpoly_t poly2,
                          const nmod_mpoly_t poly3, const nmod_mpoly_ctx_t ctx)
{
    slong i, N, len1;
    mp_bitcnt_t exp_bits;
    fmpz * max_fields2, * max_fields3;
    ulong * cmpmask;
    ulong * exp2 = poly2->exps, * exp3 = poly3->exps;
    int free2 = 0, free3 = 0;
    TMP_INIT;

    if (poly2->length == 0 || poly3->length == 0)
    {
        nmod_mpoly_zero(poly1, ctx);
        return;
    }

    TMP_START;

    max_fields2 = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    max_fields3 = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_init(max_fields2 + i);
        fmpz_init(max_fields3 + i);
    }
    mpoly_max_fields_fmpz(max_fields2, poly2->exps, poly2->length,
                                                      poly2->bits, ctx->minfo);
    mpoly_max_fields_fmpz(max_fields3, poly3->exps, poly3->length,
                                                      poly3->bits, ctx->minfo);
    _fmpz_vec_add(max_fields2, max_fields2, max_fields3, ctx->minfo->nfields);

    exp_bits = _fmpz_vec_max_bits(max_fields2, ctx->minfo->nfields);
    exp_bits = FLINT_MAX(MPOLY_MIN_BITS, exp_bits + 1);
    exp_bits = FLINT_MAX(exp_bits, poly2->bits);
    exp_bits = FLINT_MAX(exp_bits, poly3->bits);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_clear(max_fields2 + i);
        fmpz_clear(max_fields3 + i);
    }

    N = mpoly_words_per_exp(exp_bits, ctx->minfo);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, exp_bits, ctx->minfo);

    /* ensure input exponents are packed into same sized fields as output */
    if (exp_bits > poly2->bits)
    {
        free2 = 1;
        exp2 = (ulong *) flint_malloc(N*poly2->length*sizeof(ulong));
        mpoly_repack_monomials(exp2, exp_bits, poly2->exps, poly2->bits,
                                                    poly2->length, ctx->minfo);
    }

    if (exp_bits > poly3->bits)
    {
        free3 = 1;
        exp3 = (ulong *) flint_malloc(N*poly3->length*sizeof(ulong));
        mpoly_repack_monomials(exp3, exp_bits, poly3->exps, poly3->bits,
                                                    poly3->length, ctx->minfo);
    }

    /* deal with aliasing and do multiplication */
    if (poly1 == poly2 || poly1 == poly3)
    {
        nmod_mpoly_t temp;

        nmod_mpoly_init2(temp, poly2->length + poly3->length - 1, ctx);
        nmod_mpoly_fit_bits(temp, exp_bits, ctx);
        temp->bits = exp_bits;

        /* algorithm more efficient if smaller poly first */
        if (poly2->length >= poly3->length)
            len1 = _nmod_mpoly_mul_heap_threaded(
                                    &temp->coeffs, &temp->exps, &temp->alloc,
                                      poly3->coeffs, exp3, poly3->length,
                                      poly2->coeffs, exp2, poly2->length,
                                            exp_bits, N, cmpmask, ctx->ffinfo);
        else
            len1 = _nmod_mpoly_mul_heap_threaded(
                                   &temp->coeffs, &temp->exps, &temp->alloc,
                                      poly2->coeffs, exp2, poly2->length,
                                      poly3->coeffs, exp3, poly3->length,
                                            exp_bits, N, cmpmask, ctx->ffinfo);

        nmod_mpoly_swap(temp, poly1, ctx);

        nmod_mpoly_clear(temp, ctx);
    } else
    {
        nmod_mpoly_fit_length(poly1, poly2->length + poly3->length - 1, ctx);
        nmod_mpoly_fit_bits(poly1, exp_bits, ctx);
        poly1->bits = exp_bits;

        /* algorithm more efficient if smaller poly first */
        if (poly2->length > poly3->length)
            len1 = _nmod_mpoly_mul_heap_threaded(
                                &poly1->coeffs, &poly1->exps, &poly1->alloc,
                                      poly3->coeffs, exp3, poly3->length,
                                      poly2->coeffs, exp2, poly2->length,
                                            exp_bits, N, cmpmask, ctx->ffinfo);
        else
            len1 = _nmod_mpoly_mul_heap_threaded(
                                &poly1->coeffs, &poly1->exps, &poly1->alloc,
                                      poly2->coeffs, exp2, poly2->length,
                                      poly3->coeffs, exp3, poly3->length,
                                            exp_bits, N, cmpmask, ctx->ffinfo);
    }

    if (free2)
        flint_free(exp2);

    if (free3)
        flint_free(exp3);

    _nmod_mpoly_set_length(poly1, len1, ctx);

    TMP_END;
}
