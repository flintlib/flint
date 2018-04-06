/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/


#include "fmpz_mpoly.h"

/* improve locality */
#define BLOCK 128
#define MAX_ARRAY_SIZE (WORD(300000))
#define MAX_LEX_SIZE (WORD(300))


typedef struct
{
    slong idx;
    slong work;
    slong len;
    fmpz_mpoly_t poly;
}
mul_array_threaded_chunk_t;

typedef struct
{
    pthread_mutex_t mutex;
    slong nthreads;
    slong Al, Bl, Pl;
    fmpz * Acoeffs, * Bcoeffs;
    slong * Amax, * Bmax, * Asum, * Bsum;
    slong * Amain, * Bmain;
    ulong * Apexp, * Bpexp;
    slong * perm;
    slong nvars;
    const ulong * mults;
    slong array_size;
    slong degb;
    mul_array_threaded_chunk_t * Pchunks;
    int rev;
    volatile int idx;
}
mul_array_threaded_base_t;

typedef struct
{
    slong idx;
    slong time;
    mul_array_threaded_base_t * basep;
    ulong * exp;
}
mul_array_threaded_arg_t;







/******************
    LEX
******************/

void * _fmpz_mpoly_mul_array_threaded_worker_LEX(void * arg_ptr)
{
    slong i, j, Pi;
    mul_array_threaded_arg_t * arg = (mul_array_threaded_arg_t *) arg_ptr;
    mul_array_threaded_base_t * base;
    ulong * coeff_array;
    TMP_INIT;

    base = arg->basep;

    TMP_START;
    coeff_array = (ulong *) TMP_ALLOC(3*base->array_size*sizeof(ulong));
    for (j = 0; j < 3*base->array_size; j++)
        coeff_array[j] = 0;;

    pthread_mutex_lock(&base->mutex);
    Pi = base->idx;
    base->idx = Pi + 1;
    pthread_mutex_unlock(&base->mutex);

    while (Pi < base->Pl)
    {
        /* work out bit counts for this chunk */
        slong Abits = 0;
        slong Bbits = 0;
        slong Pbits = 0;
        slong number = 0;
        for (i = 0, j = base->perm[Pi]; i < base->Al && j >= 0; i++, j--)
        {
            if (j < base->Bl)
            {
                number++;
                Abits = FLINT_MAX(Abits, base->Amax[i]);
                Bbits = FLINT_MAX(Bbits, base->Bmax[j]);
                Pbits = FLINT_MAX(Pbits,
                            FLINT_MIN(base->Asum[i] + base->Bmax[j],
                                      base->Amax[i] + base->Bsum[j]));
            }
        }
        Pbits += FLINT_BIT_COUNT(number) + 1; /* includes one bit for sign */

        if (Abits <= FLINT_BITS - 2 && Bbits <= FLINT_BITS - 2)
        {
            if (Pbits <= FLINT_BITS)
            {
                for (i = 0, j = base->perm[Pi]; i < base->Al && j >= 0; i++, j--)
                {
                    if (j < base->Bl)
                    {
                        _fmpz_mpoly_addmul_array1_slong1((ulong *)coeff_array, 
                            base->Acoeffs + base->Amain[i],
                                base->Apexp + base->Amain[i],
                                base->Amain[i + 1] - base->Amain[i],
                            base->Bcoeffs + base->Bmain[j],
                                base->Bpexp + base->Bmain[j],
                                base->Bmain[j + 1] - base->Bmain[j]);
                    }
                }
                (base->Pchunks + base->perm[Pi])->len = 
                    fmpz_mpoly_append_array_sm1_LEX(
                        (base->Pchunks + base->perm[Pi])->poly, 0,
                        (ulong *)coeff_array, base->mults, base->nvars - 1,
                        base->array_size, base->Pl - base->perm[Pi] - 1);

            } else if (Pbits <= 2*FLINT_BITS)
            {
                for (i = 0, j = base->perm[Pi]; i < base->Al && j >= 0; i++, j--)
                {
                    if (j < base->Bl)
                    {
                        _fmpz_mpoly_addmul_array1_slong2((ulong *)coeff_array, 
                            base->Acoeffs + base->Amain[i],
                                base->Apexp + base->Amain[i],
                                base->Amain[i + 1] - base->Amain[i],
                            base->Bcoeffs + base->Bmain[j],
                                base->Bpexp + base->Bmain[j],
                                base->Bmain[j + 1] - base->Bmain[j]);
                    }
                }
                (base->Pchunks + base->perm[Pi])->len = 
                    fmpz_mpoly_append_array_sm2_LEX(
                        (base->Pchunks + base->perm[Pi])->poly, 0,
                        (ulong *)coeff_array, base->mults, base->nvars - 1,
                        base->array_size, base->Pl - base->perm[Pi] - 1);
            } else
            {
                for (i = 0, j = base->perm[Pi]; i < base->Al && j >= 0; i++, j--)
                {
                    if (j < base->Bl)
                    {
                        _fmpz_mpoly_addmul_array1_slong((ulong *)coeff_array, 
                            base->Acoeffs + base->Amain[i],
                                base->Apexp + base->Amain[i],
                                base->Amain[i + 1] - base->Amain[i],
                            base->Bcoeffs + base->Bmain[j],
                                base->Bpexp + base->Bmain[j],
                                base->Bmain[j + 1] - base->Bmain[j]);
                    }
                }
                (base->Pchunks + base->perm[Pi])->len = 
                    fmpz_mpoly_append_array_sm3_LEX(
                        (base->Pchunks + base->perm[Pi])->poly, 0,
                        (ulong *)coeff_array, base->mults, base->nvars - 1,
                        base->array_size, base->Pl - base->perm[Pi] - 1);
            }
        } else
        {
            for (i = 0, j = base->perm[Pi]; i < base->Al && j >= 0; i++, j--)
            {
                if (j < base->Bl)
                {
                    _fmpz_mpoly_addmul_array1_fmpz((fmpz *)coeff_array, 
                        base->Acoeffs + base->Amain[i],
                            base->Apexp + base->Amain[i],
                            base->Amain[i + 1] - base->Amain[i],
                        base->Bcoeffs + base->Bmain[j],
                            base->Bpexp + base->Bmain[j],
                            base->Bmain[j + 1] - base->Bmain[j]);
                }
            }
            (base->Pchunks + base->perm[Pi])->len = 
                fmpz_mpoly_append_array_fmpz_LEX(
                    (base->Pchunks + base->perm[Pi])->poly, 0,
                    (fmpz *)coeff_array, base->mults, base->nvars - 1,
                    base->array_size, base->Pl - base->perm[Pi] - 1);
        }

        pthread_mutex_lock(&base->mutex);
        Pi = base->idx;
        base->idx = Pi + 1;
        pthread_mutex_unlock(&base->mutex);
    }

    if (arg->idx > 0)
    {
        flint_cleanup();
    }

    TMP_END;
    return NULL;
}


void _fmpz_mpoly_mul_array_chunked_threaded_LEX(fmpz_mpoly_t P,
                             const fmpz_mpoly_t A, const fmpz_mpoly_t B, 
                               const ulong * mults, const fmpz_mpoly_ctx_t ctx)
{
    slong nvars = ctx->minfo->nvars;
    slong Pi, i, j, Plen, Pl, Al, Bl, array_size;
    slong * Asum, * Amax, * Bsum, * Bmax;
    slong * Amain, * Bmain;
    ulong * Apexp, * Bpexp;
    pthread_t * threads;
    mul_array_threaded_arg_t * args;
    mul_array_threaded_base_t * base;
    mul_array_threaded_chunk_t * Pchunks;
    slong * perm;
    TMP_INIT;

    array_size = 1;
    for (i = 0; i < nvars - 1; i++) {
        array_size *= mults[i];
    }

    /* compute lengths of poly2 and poly3 in chunks */
    Al = 1 + (slong) (A->exps[0] >> (A->bits*(nvars - 1)));
    Bl = 1 + (slong) (B->exps[0] >> (B->bits*(nvars - 1)));

    TMP_START;

    /* compute indices and lengths of coefficients of polys in main variable */
    Amain = (slong *) TMP_ALLOC((Al + 1)*sizeof(slong));
    Bmain = (slong *) TMP_ALLOC((Bl + 1)*sizeof(slong));
    Asum  = (slong *) TMP_ALLOC(Al*sizeof(slong));
    Amax  = (slong *) TMP_ALLOC(Al*sizeof(slong));
    Bsum  = (slong *) TMP_ALLOC(Bl*sizeof(slong));
    Bmax  = (slong *) TMP_ALLOC(Bl*sizeof(slong));
    Apexp = (ulong *) TMP_ALLOC(A->length*sizeof(ulong));
    Bpexp = (ulong *) TMP_ALLOC(B->length*sizeof(ulong));
    mpoly_main_variable_split_LEX(Amain, Apexp, A->exps, Al, A->length,
                                                    mults, nvars - 1, A->bits);
    mpoly_main_variable_split_LEX(Bmain, Bpexp, B->exps, Bl, B->length,
                                                    mults, nvars - 1, B->bits);

    /* work out bit counts for each chunk */
    for (i = 0; i < Al; i++)
    {
        _fmpz_vec_sum_max_bits(&Asum[i], &Amax[i],
                                A->coeffs + Amain[i], Amain[i + 1] - Amain[i]);
    }
    for (j = 0; j < Bl; j++)
    {
        _fmpz_vec_sum_max_bits(&Bsum[j], &Bmax[j],
                                B->coeffs + Bmain[j], Bmain[j + 1] - Bmain[j]);
    }

    Pl = Al + Bl - 1;

    /* work out data for each chunk of the output */
    Pchunks = (mul_array_threaded_chunk_t *) TMP_ALLOC(Pl
                                          *sizeof(mul_array_threaded_chunk_t));
    perm = (slong *) TMP_ALLOC(Pl*sizeof(slong));
    for (Pi = 0; Pi < Pl; Pi++)
    {
        fmpz_mpoly_init2((Pchunks + Pi)->poly, 8, ctx);
        fmpz_mpoly_fit_bits((Pchunks + Pi)->poly, P->bits, ctx);
        (Pchunks + Pi)->work = 0;
        perm[Pi] = Pi;
        for (i = 0, j = Pi; i < Al && j >= 0; i++, j--)
        {
            if (j < Bl)
            {
                (Pchunks + Pi)->work += (Amain[i + 1] - Amain[i])
                                       *(Bmain[j + 1] - Bmain[j]);
            }
        }
    }

    for (i = 0; i < Pl; i++)
    {
        for (j = i; j > 0 && (Pchunks + perm[j - 1])->work
                                   < (Pchunks + perm[j])->work; j--)
        {
            slong t = perm[j - 1];
            perm[j - 1] = perm[j];
            perm[j] = t;
        }
    }

    base = (mul_array_threaded_base_t *) TMP_ALLOC(sizeof(
                                                   mul_array_threaded_base_t));
    base->nthreads = flint_get_num_threads();
    base->Al = Al;
    base->Bl = Bl;
    base->Pl = Pl;
    base->Amax = Amax;
    base->Bmax = Bmax;
    base->Asum = Asum;
    base->Bsum = Bsum;
    base->Acoeffs = A->coeffs;
    base->Amain = Amain;
    base->Apexp = Apexp;
    base->Bcoeffs = B->coeffs;
    base->Bmain = Bmain;
    base->Bpexp = Bpexp;
    base->idx = 0;
    base->perm = perm;
    base->nvars = nvars;
    base->Pchunks = Pchunks;
    base->array_size = array_size;
    base->mults = mults;

    args    = (mul_array_threaded_arg_t *) TMP_ALLOC(
                            sizeof(mul_array_threaded_arg_t) * base->nthreads);
    threads = (pthread_t *) TMP_ALLOC(sizeof(pthread_t) * base->nthreads);

    pthread_mutex_init(&base->mutex, NULL);
    for (i = base->nthreads - 1; i >= 0; i--)
    {
        args[i].idx = i;
        args[i].basep = base;
        if (i > 0)
        {
            pthread_create(&threads[i], NULL,
                          _fmpz_mpoly_mul_array_threaded_worker_LEX, &args[i]);
        } else
        {
            _fmpz_mpoly_mul_array_threaded_worker_LEX(&args[i]);
        }
    }
    for (i = base->nthreads - 1; i > 0; i--)
    {
        pthread_join(threads[i], NULL);
    }
    pthread_mutex_destroy(&base->mutex);

    Plen = 0;
    for (Pi = 0; Pi < Pl; Pi++)
    {
        _fmpz_mpoly_fit_length(&P->coeffs, &P->exps, &P->alloc,
                                                Plen + (Pchunks + Pi)->len, 1);
        for (i = 0; i < (Pchunks + Pi)->len; i++)
        {
            P->exps[Plen] = (Pchunks + Pi)->poly->exps[i];
            fmpz_swap(P->coeffs + Plen, (Pchunks + Pi)->poly->coeffs + i);
            fmpz_clear((Pchunks + Pi)->poly->coeffs + i);
            Plen++;
        }

        flint_free((Pchunks + Pi)->poly->coeffs);
        flint_free((Pchunks + Pi)->poly->exps);
    }

    TMP_END;
    _fmpz_mpoly_set_length(P, Plen, ctx);
}


int fmpz_mpoly_mul_array_threaded_LEX(fmpz_mpoly_t poly1,
                          const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i, exp_bits, array_size;
    ulong max, * max_fields2, * max_fields3;
    int success = 1;
    TMP_INIT;

    /* input poly is zero */
    if (poly2->length == 0 || poly3->length == 0)
    {
        fmpz_mpoly_zero(poly1, ctx);
        return 1;
    }
    /* lets only work with exponents packed into 1 word */
    if (    1 != mpoly_words_per_exp(poly2->bits, ctx->minfo)
         || 1 != mpoly_words_per_exp(poly3->bits, ctx->minfo))
    {
        return 0;
    }

    TMP_START;

    /* compute maximum exponents for each variable */
    max_fields2 = (ulong *) TMP_ALLOC(ctx->minfo->nfields*sizeof(ulong));
    max_fields3 = (ulong *) TMP_ALLOC(ctx->minfo->nfields*sizeof(ulong));
    mpoly_max_fields_ui(max_fields2, poly2->exps, poly2->length,
                                                      poly2->bits, ctx->minfo);
    mpoly_max_fields_ui(max_fields3, poly3->exps, poly3->length,
                                                      poly3->bits, ctx->minfo);

    /* the field of index n-1 is the one that wil be pulled out */
    i = ctx->minfo->nfields - 1;
    max_fields2[i] += max_fields3[i] + 1;
    max = max_fields2[i];
    if (((slong) max_fields2[i]) <= 0 || max_fields2[i] > MAX_LEX_SIZE)
    {
        success = 0;
        goto cleanup;
    }

    /* the fields of index n-2...0, contribute to the array size */
    array_size = WORD(1);
    for (i--; i >= 0; i--)
    {
        ulong hi;
        max_fields2[i] += max_fields3[i] + 1;
        max |= max_fields2[i];
        umul_ppmm(hi, array_size, array_size, max_fields2[i]);
        if (hi != WORD(0) || (array_size | (slong) max_fields2[i]) <= 0
                          || array_size > MAX_ARRAY_SIZE)
        {
            success = 0;
            goto cleanup;
        }
    }

    exp_bits = FLINT_MAX(WORD(8), FLINT_BIT_COUNT(max) + 1);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);

    /* array multiplication assumes result fit into 1 word */
    if (ctx->minfo->ord != ORD_LEX ||
            1 != mpoly_words_per_exp(exp_bits, ctx->minfo))
    {
        success = 0;
        goto cleanup;
    }

    /* handle aliasing and do array multiplication */
    success = 1;
    if (poly1 == poly2 || poly1 == poly3)
    {
        fmpz_mpoly_t temp;
        fmpz_mpoly_init2(temp, poly2->length + poly3->length - 1, ctx);
        fmpz_mpoly_fit_bits(temp, exp_bits, ctx);
        temp->bits = exp_bits;
        _fmpz_mpoly_mul_array_chunked_threaded_LEX(temp, poly3, poly2,
                                                             max_fields2, ctx);
        fmpz_mpoly_swap(temp, poly1, ctx);
        fmpz_mpoly_clear(temp, ctx);
    } else
    {
        fmpz_mpoly_fit_length(poly1, poly2->length + poly3->length - 1, ctx);
        fmpz_mpoly_fit_bits(poly1, exp_bits, ctx);
        poly1->bits = exp_bits;
        _fmpz_mpoly_mul_array_chunked_threaded_LEX(poly1, poly3, poly2,
                                                             max_fields2, ctx);
    }

cleanup:

    TMP_END;

    return success;
}






/*****************************
    DEGLEX and DEGREVLEX
*****************************/


void * _fmpz_mpoly_mul_array_threaded_worker_DEG(void * arg_ptr)
{
    slong i, j, Pi;
    mul_array_threaded_arg_t * arg = (mul_array_threaded_arg_t *) arg_ptr;
    mul_array_threaded_base_t * base;
    ulong * coeff_array;
    slong (* upack_sm1)(fmpz_mpoly_t, slong, ulong *, slong, slong, slong); 
    slong (* upack_sm2)(fmpz_mpoly_t, slong, ulong *, slong, slong, slong); 
    slong (* upack_sm3)(fmpz_mpoly_t, slong, ulong *, slong, slong, slong); 
    slong (* upack_fmpz)(fmpz_mpoly_t, slong, fmpz *, slong, slong, slong); 
    TMP_INIT;

    base = arg->basep;

    upack_sm1  = &fmpz_mpoly_append_array_sm1_DEGLEX;
    upack_sm2  = &fmpz_mpoly_append_array_sm2_DEGLEX;
    upack_sm3  = &fmpz_mpoly_append_array_sm3_DEGLEX;
    upack_fmpz = &fmpz_mpoly_append_array_fmpz_DEGLEX;
    if (base->rev) {
        upack_sm1  = &fmpz_mpoly_append_array_sm1_DEGREVLEX;
        upack_sm2  = &fmpz_mpoly_append_array_sm2_DEGREVLEX;
        upack_sm3  = &fmpz_mpoly_append_array_sm3_DEGREVLEX;
        upack_fmpz = &fmpz_mpoly_append_array_fmpz_DEGREVLEX;
    }

    TMP_START;
    coeff_array = (ulong *) TMP_ALLOC(3*base->array_size*sizeof(ulong));
    for (j = 0; j < 3*base->array_size; j++)
        coeff_array[j] = 0;;

    pthread_mutex_lock(&base->mutex);
    Pi = base->idx;
    base->idx = Pi + 1;
    pthread_mutex_unlock(&base->mutex);

    while (Pi < base->Pl)
    {
        /* work out bit counts for this chunk */
        slong Abits = 0;
        slong Bbits = 0;
        slong Pbits = 0;
        slong number = 0;
        for (i = 0, j = base->perm[Pi]; i < base->Al && j >= 0; i++, j--)
        {
            if (j < base->Bl)
            {
                number++;
                Abits = FLINT_MAX(Abits, base->Amax[i]);
                Bbits = FLINT_MAX(Bbits, base->Bmax[j]);
                Pbits = FLINT_MAX(Pbits,
                            FLINT_MIN(base->Asum[i] + base->Bmax[j],
                                      base->Amax[i] + base->Bsum[j]));
            }
        }
        Pbits += FLINT_BIT_COUNT(number) + 1; /* includes one bit for sign */

        if (Abits <= FLINT_BITS - 2 && Bbits <= FLINT_BITS - 2)
        {
            if (Pbits <= FLINT_BITS)
            {
                for (i = 0, j = base->perm[Pi]; i < base->Al && j >= 0; i++, j--)
                {
                    if (j < base->Bl)
                    {
                        _fmpz_mpoly_addmul_array1_slong1((ulong *)coeff_array, 
                            base->Acoeffs + base->Amain[i],
                                base->Apexp + base->Amain[i],
                                base->Amain[i + 1] - base->Amain[i],
                            base->Bcoeffs + base->Bmain[j],
                                base->Bpexp + base->Bmain[j],
                                base->Bmain[j + 1] - base->Bmain[j]);
                    }
                }
                (base->Pchunks + base->perm[Pi])->len = 
                    upack_sm1((base->Pchunks + base->perm[Pi])->poly, 0,
                        (ulong *)coeff_array, base->Pl - base->perm[Pi] - 1,
                                                      base->nvars, base->degb);

            } else if (Pbits <= 2*FLINT_BITS)
            {
                for (i = 0, j = base->perm[Pi]; i < base->Al && j >= 0; i++, j--)
                {
                    if (j < base->Bl)
                    {
                        _fmpz_mpoly_addmul_array1_slong2((ulong *)coeff_array, 
                            base->Acoeffs + base->Amain[i],
                                base->Apexp + base->Amain[i],
                                base->Amain[i + 1] - base->Amain[i],
                            base->Bcoeffs + base->Bmain[j],
                                base->Bpexp + base->Bmain[j],
                                base->Bmain[j + 1] - base->Bmain[j]);
                    }
                }
                (base->Pchunks + base->perm[Pi])->len = 
                    upack_sm2((base->Pchunks + base->perm[Pi])->poly, 0,
                        (ulong *)coeff_array, base->Pl - base->perm[Pi] - 1,
                                                      base->nvars, base->degb);
            } else
            {
                for (i = 0, j = base->perm[Pi]; i < base->Al && j >= 0; i++, j--)
                {
                    if (j < base->Bl)
                    {
                        _fmpz_mpoly_addmul_array1_slong((ulong *)coeff_array, 
                            base->Acoeffs + base->Amain[i],
                                base->Apexp + base->Amain[i],
                                base->Amain[i + 1] - base->Amain[i],
                            base->Bcoeffs + base->Bmain[j],
                                base->Bpexp + base->Bmain[j],
                                base->Bmain[j + 1] - base->Bmain[j]);
                    }
                }
                (base->Pchunks + base->perm[Pi])->len = 
                    upack_sm3((base->Pchunks + base->perm[Pi])->poly, 0,
                        (ulong *)coeff_array, base->Pl - base->perm[Pi] - 1,
                                                      base->nvars, base->degb);
            }
        } else
        {
            for (i = 0, j = base->perm[Pi]; i < base->Al && j >= 0; i++, j--)
            {
                if (j < base->Bl)
                {
                    _fmpz_mpoly_addmul_array1_fmpz((fmpz *)coeff_array, 
                        base->Acoeffs + base->Amain[i],
                            base->Apexp + base->Amain[i],
                            base->Amain[i + 1] - base->Amain[i],
                        base->Bcoeffs + base->Bmain[j],
                            base->Bpexp + base->Bmain[j],
                            base->Bmain[j + 1] - base->Bmain[j]);
                }
            }
            (base->Pchunks + base->perm[Pi])->len = 
                upack_fmpz((base->Pchunks + base->perm[Pi])->poly, 0,
                    (fmpz *)coeff_array, base->Pl - base->perm[Pi] - 1,
                                                      base->nvars, base->degb);
        }

        pthread_mutex_lock(&base->mutex);
        Pi = base->idx;
        base->idx = Pi + 1;
        pthread_mutex_unlock(&base->mutex);
    }

    if (arg->idx > 0)
    {
        flint_cleanup();
    }

    TMP_END;
    return NULL;
}



void _fmpz_mpoly_mul_array_chunked_threaded_DEG(fmpz_mpoly_t P,
                             const fmpz_mpoly_t A, const fmpz_mpoly_t B, 
                                        ulong degb, const fmpz_mpoly_ctx_t ctx)
{
    slong nvars = ctx->minfo->nvars;
    slong Pi, i, j, Plen, Pl, Al, Bl, array_size;
    slong * Asum, * Amax, * Bsum, * Bmax;
    slong * Amain, * Bmain;
    ulong * Apexp, * Bpexp;

    pthread_t * threads;
    mul_array_threaded_arg_t * args;
    mul_array_threaded_base_t * base;
    mul_array_threaded_chunk_t * Pchunks;

    slong * perm;
    TMP_INIT;

    TMP_START;

    /* compute lengths of poly2 and poly3 in chunks */
    Al = 1 + (slong) (A->exps[0] >> (A->bits*nvars));
    Bl = 1 + (slong) (B->exps[0] >> (B->bits*nvars));

    array_size = 1;
    for (i = 0; i < nvars-1; i++) {
        array_size *= degb;
    }

    /* compute indices and lengths of coefficients of polys in main variable */
    Amain = (slong *) TMP_ALLOC((Al + 1)*sizeof(slong));
    Bmain = (slong *) TMP_ALLOC((Bl + 1)*sizeof(slong));
    Asum  = (slong *) TMP_ALLOC(Al*sizeof(slong));
    Amax  = (slong *) TMP_ALLOC(Al*sizeof(slong));
    Bsum  = (slong *) TMP_ALLOC(Bl*sizeof(slong));
    Bmax  = (slong *) TMP_ALLOC(Bl*sizeof(slong));
    Apexp = (ulong *) TMP_ALLOC(A->length*sizeof(ulong));
    Bpexp = (ulong *) TMP_ALLOC(B->length*sizeof(ulong));
    mpoly_main_variable_split_DEG(Amain, Apexp, A->exps, Al, A->length,
                                                         degb, nvars, A->bits);
    mpoly_main_variable_split_DEG(Bmain, Bpexp, B->exps, Bl, B->length,
                                                         degb, nvars, B->bits);

    /* work out bit counts for each chunk */
    for (i = 0; i < Al; i++)
    {
        _fmpz_vec_sum_max_bits(&Asum[i], &Amax[i],
                                A->coeffs + Amain[i], Amain[i + 1] - Amain[i]);
    }
    for (j = 0; j < Bl; j++)
    {
        _fmpz_vec_sum_max_bits(&Bsum[j], &Bmax[j],
                                B->coeffs + Bmain[j], Bmain[j + 1] - Bmain[j]);
    }

    Pl = Al + Bl - 1;
    FLINT_ASSERT(Pl == degb);

    /* work out data for each chunk of the output */
    Pchunks = (mul_array_threaded_chunk_t *) TMP_ALLOC(Pl
                                          *sizeof(mul_array_threaded_chunk_t));
    perm = (slong *) TMP_ALLOC(Pl*sizeof(slong));
    for (Pi = 0; Pi < Pl; Pi++)
    {
        fmpz_mpoly_init2((Pchunks + Pi)->poly, 8, ctx);
        fmpz_mpoly_fit_bits((Pchunks + Pi)->poly, P->bits, ctx);
        (Pchunks + Pi)->work = 0;
        perm[Pi] = Pi;
        for (i = 0, j = Pi; i < Al && j >= 0; i++, j--)
        {
            if (j < Bl)
            {
                (Pchunks + Pi)->work += (Amain[i + 1] - Amain[i])
                                       *(Bmain[j + 1] - Bmain[j]);
            }
        }
    }

    for (i = 0; i < Pl; i++)
    {
        for (j = i; j > 0 && (Pchunks + perm[j-1])->work
                                         < (Pchunks + perm[j])->work; j--)
        {
            slong t = perm[j - 1];
            perm[j - 1] = perm[j];
            perm[j] = t;
        }
    }

    base = (mul_array_threaded_base_t *) TMP_ALLOC(sizeof(
                                                   mul_array_threaded_base_t));
    base->nthreads = flint_get_num_threads();
    base->Al = Al;
    base->Bl = Bl;
    base->Pl = Pl;
    base->Amax = Amax;
    base->Bmax = Bmax;
    base->Asum = Asum;
    base->Bsum = Bsum;
    base->Acoeffs = A->coeffs;
    base->Amain = Amain;
    base->Apexp = Apexp;
    base->Bcoeffs = B->coeffs;
    base->Bmain = Bmain;
    base->Bpexp = Bpexp;
    base->idx = 0;
    base->perm = perm;
    base->nvars = nvars;
    base->Pchunks = Pchunks;
    base->array_size = array_size;
    base->degb = degb;
    base->rev = (ctx->minfo->ord == ORD_DEGREVLEX);

    args    = (mul_array_threaded_arg_t *) TMP_ALLOC(sizeof(
                                   mul_array_threaded_arg_t) * base->nthreads);
    threads = (pthread_t *) TMP_ALLOC(sizeof(pthread_t) * base->nthreads);

    pthread_mutex_init(&base->mutex, NULL);
    for (i = base->nthreads - 1; i >= 0; i--)
    {
        args[i].idx = i;
        args[i].basep = base;
        if (i > 0)
        {
            pthread_create(&threads[i], NULL,
                          _fmpz_mpoly_mul_array_threaded_worker_DEG, &args[i]);
        } else
        {
            _fmpz_mpoly_mul_array_threaded_worker_DEG(&args[i]);
        }
    }
    for (i = base->nthreads - 1; i > 0; i--)
    {
        pthread_join(threads[i], NULL);
    }
    pthread_mutex_destroy(&base->mutex);

    Plen = 0;
    for (Pi = 0; Pi < Pl; Pi++)
    {
        _fmpz_mpoly_fit_length(&P->coeffs, &P->exps, &P->alloc,
                                                Plen + (Pchunks + Pi)->len, 1);
        for (i = 0; i < (Pchunks + Pi)->len; i++)
        {
            P->exps[Plen] = (Pchunks + Pi)->poly->exps[i];
            fmpz_swap(P->coeffs + Plen, (Pchunks + Pi)->poly->coeffs + i);
            fmpz_clear((Pchunks + Pi)->poly->coeffs + i);
            Plen++;
        }

        flint_free((Pchunks + Pi)->poly->coeffs);
        flint_free((Pchunks + Pi)->poly->exps);
    }

    TMP_END;
    _fmpz_mpoly_set_length(P, Plen, ctx);
}


int fmpz_mpoly_mul_array_threaded_DEG(fmpz_mpoly_t poly1,
                         const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i, exp_bits, array_size;
    ulong deg, * max_fields2, * max_fields3;
    int success = 1;
    TMP_INIT;

    /* input poly is zero */
    if (poly2->length == 0 || poly3->length == 0)
    {
        fmpz_mpoly_zero(poly1, ctx);
        return 1;
    }

    /* lets only work with exponents packed into 1 word */
    if ((     ctx->minfo->ord != ORD_DEGREVLEX 
           && ctx->minfo->ord != ORD_DEGLEX)
        || 1 != mpoly_words_per_exp(poly2->bits, ctx->minfo)
        || 1 != mpoly_words_per_exp(poly3->bits, ctx->minfo)
       )
    {
        return 0;
    }

    TMP_START;

    /* compute maximum exponents for each variable */
    max_fields2 = (ulong *) TMP_ALLOC(ctx->minfo->nfields*sizeof(ulong));
    max_fields3 = (ulong *) TMP_ALLOC(ctx->minfo->nfields*sizeof(ulong));
    mpoly_max_fields_ui(max_fields2, poly2->exps, poly2->length,
                                                      poly2->bits, ctx->minfo);
    mpoly_max_fields_ui(max_fields3, poly3->exps, poly3->length,
                                                      poly3->bits, ctx->minfo);

    /* the field of index n-1 is the one that wil be pulled out */
    i = ctx->minfo->nfields - 1;
    deg = max_fields2[i] + max_fields3[i] + 1;
    if (((slong) deg) <= 0 || deg > MAX_ARRAY_SIZE)
    {
        success = 0;
        goto cleanup;
    }

    /* the fields of index n-2...1, contribute to the array size */
    array_size = WORD(1);
    for (i--; i >= 1; i--)
    {
        ulong hi;
        umul_ppmm(hi, array_size, array_size, deg);
        if (hi != WORD(0) || array_size <= 0
                          || array_size > MAX_ARRAY_SIZE)
        {
            success = 0;
            goto cleanup;
        }
    }

    exp_bits = FLINT_MAX(WORD(8), FLINT_BIT_COUNT(deg) + 1);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);

    /* array multiplication assumes result fit into 1 word */
    if (1 != mpoly_words_per_exp(exp_bits, ctx->minfo))
    {
        success = 0;
        goto cleanup;
    }

    /* handle aliasing and do array multiplication */
    success = 1;
    if (poly1 == poly2 || poly1 == poly3)
    {
        fmpz_mpoly_t temp;
        fmpz_mpoly_init2(temp, poly2->length + poly3->length - 1, ctx);
        fmpz_mpoly_fit_bits(temp, exp_bits, ctx);
        temp->bits = exp_bits;
        _fmpz_mpoly_mul_array_chunked_threaded_DEG(temp, poly3, poly2, deg, ctx);
        fmpz_mpoly_swap(temp, poly1, ctx);
        fmpz_mpoly_clear(temp, ctx);
    } else
    {
        fmpz_mpoly_fit_length(poly1, poly2->length + poly3->length - 1, ctx);
        fmpz_mpoly_fit_bits(poly1, exp_bits, ctx);
        poly1->bits = exp_bits;
        _fmpz_mpoly_mul_array_chunked_threaded_DEG(poly1, poly3, poly2, deg, ctx);
    }

cleanup:

    TMP_END;

    return success;
}



int fmpz_mpoly_mul_array_threaded(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                          const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)
{
    switch (ctx->minfo->ord)
    {
        case ORD_LEX:
            return fmpz_mpoly_mul_array_threaded_LEX(poly1, poly2, poly3, ctx);
        case ORD_DEGREVLEX:
        case ORD_DEGLEX:
            return fmpz_mpoly_mul_array_threaded_DEG(poly1, poly2, poly3, ctx);
        default:
            return 0;
    }
}
