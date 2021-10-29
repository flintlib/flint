/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "string.h"
#include "thread_pool.h"
#include "nmod_mpoly.h"

/* improve locality */
#define BLOCK 128
#define MAX_ARRAY_SIZE (WORD(300000))
#define MAX_LEX_SIZE (WORD(300))


typedef struct
{
    slong idx;
    slong work;
    slong len;
    nmod_mpoly_t poly;
}
_chunk_struct;


typedef struct
{
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif
    volatile int idx;
    slong nthreads;
    slong Al, Bl, Pl;
    mp_limb_t * Acoeffs, * Bcoeffs;
    slong * Amain, * Bmain;
    ulong * Apexp, * Bpexp;
    slong * perm;
    slong nvars;
    const ulong * mults;
    slong array_size;
    slong degb;
    const nmod_mpoly_ctx_struct * ctx;
    _chunk_struct * Pchunks;
    int rev;
}
_base_struct;

typedef _base_struct _base_t[1];


typedef struct
{
    slong idx;
    slong time;
    _base_struct * base;
    ulong * exp;
}
_worker_arg_struct;


/******************
    LEX
******************/

static void _nmod_mpoly_mul_array_threaded_worker_LEX(void * varg)
{
    slong i, j, Pi;
    _worker_arg_struct * arg = (_worker_arg_struct *) varg;
    _base_struct * base = arg->base;
    slong  Al = base->Al;
    slong Bl = base->Bl;
    slong Pl = base->Pl;
    slong * Amain = base->Amain;
    slong * Bmain = base->Bmain;
    ulong * coeff_array;
    
    TMP_INIT;

    TMP_START;
    coeff_array = (ulong *) TMP_ALLOC(3*base->array_size*sizeof(ulong));
    for (j = 0; j < 3*base->array_size; j++)
        coeff_array[j] = 0;

#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&base->mutex);
#endif
    Pi = base->idx;
    base->idx = Pi + 1;
#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&base->mutex);
#endif

    while (Pi < Pl)
    {
        slong len;
        mp_limb_t t2, t1, t0, u1, u0;

        Pi = base->perm[Pi];

        /* work out bit counts for this chunk */
        len = 0;
        for (i = 0, j = Pi; i < Al && j >= 0; i++, j--)
        {
            if (j < Bl)
            {
                len += FLINT_MIN(Amain[i + 1] - Amain[i],
                                 Bmain[j + 1] - Bmain[j]);
            }
        }

        umul_ppmm(t1, t0, base->ctx->mod.n - 1, base->ctx->mod.n - 1);
        umul_ppmm(t2, t1, t1, len);
        umul_ppmm(u1, u0, t0, len);
        add_sssaaaaaa(t2, t1, t0,  t2, t1, UWORD(0),  UWORD(0), u1, u0);

        (base->Pchunks + Pi)->len = 0;

        if (t2 != 0)
        {
            /* need three words */
            for (i = 0, j = Pi; i < Al && j >= 0; i++, j--)
            {
                if (j >= Bl)
                    continue;

                _nmod_mpoly_addmul_array1_ulong3(coeff_array, 
                        base->Acoeffs + base->Amain[i],
                            base->Apexp + base->Amain[i],
                            base->Amain[i + 1] - base->Amain[i],
                        base->Bcoeffs + base->Bmain[j],
                            base->Bpexp + base->Bmain[j],
                            base->Bmain[j + 1] - base->Bmain[j]);
            }

            (base->Pchunks + Pi)->len = 
                nmod_mpoly_append_array_sm3_LEX(
                    (base->Pchunks + Pi)->poly, 0,
                    coeff_array, base->mults, base->nvars - 1,
                    base->array_size, Pl - Pi - 1, base->ctx);
        }
        else if (t1 != 0)
        {
            /* fits into two words */
            for (i = 0, j = Pi; i < Al && j >= 0; i++, j--)
            {
                if (j >= Bl)
                    continue;

                _nmod_mpoly_addmul_array1_ulong2(coeff_array, 
                        base->Acoeffs + base->Amain[i],
                            base->Apexp + base->Amain[i],
                            base->Amain[i + 1] - base->Amain[i],
                        base->Bcoeffs + base->Bmain[j],
                            base->Bpexp + base->Bmain[j],
                            base->Bmain[j + 1] - base->Bmain[j]);
            }

            (base->Pchunks + Pi)->len = 
                nmod_mpoly_append_array_sm2_LEX(
                    (base->Pchunks + Pi)->poly, 0,
                    coeff_array, base->mults, base->nvars - 1,
                    base->array_size, Pl - Pi - 1, base->ctx);
            
        }
        else if (t0 != 0)
        {
            /* fits into one word */
            for (i = 0, j = Pi; i < Al && j >= 0; i++, j--)
            {
                if (j >= Bl)
                    continue;

                _nmod_mpoly_addmul_array1_ulong1(coeff_array, 
                        base->Acoeffs + base->Amain[i],
                            base->Apexp + base->Amain[i],
                            base->Amain[i + 1] - base->Amain[i],
                        base->Bcoeffs + base->Bmain[j],
                            base->Bpexp + base->Bmain[j],
                            base->Bmain[j + 1] - base->Bmain[j]);
            }

            (base->Pchunks + Pi)->len = 
                nmod_mpoly_append_array_sm1_LEX(
                    (base->Pchunks + Pi)->poly, 0,
                    coeff_array, base->mults, base->nvars - 1,
                    base->array_size, Pl - Pi - 1, base->ctx);
        }

#if FLINT_USES_PTHREAD
        pthread_mutex_lock(&base->mutex);
#endif
	Pi = base->idx;
        base->idx = Pi + 1;
#if FLINT_USES_PTHREAD
        pthread_mutex_unlock(&base->mutex);
#endif
    }

    TMP_END;
}

void _nmod_mpoly_mul_array_chunked_threaded_LEX(
    nmod_mpoly_t P,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const ulong * mults,
    const nmod_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    slong nvars = ctx->minfo->nvars;
    slong Pi, i, j, Plen, Pl, Al, Bl, array_size;
    slong * Amain, * Bmain;
    ulong * Apexp, * Bpexp;
    _base_t base;
    _worker_arg_struct * args;
    _chunk_struct * Pchunks;
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
    Apexp = (ulong *) flint_malloc(A->length*sizeof(ulong));
    Bpexp = (ulong *) flint_malloc(B->length*sizeof(ulong));
    mpoly_main_variable_split_LEX(Amain, Apexp, A->exps, Al, A->length,
                                                    mults, nvars - 1, A->bits);
    mpoly_main_variable_split_LEX(Bmain, Bpexp, B->exps, Bl, B->length,
                                                    mults, nvars - 1, B->bits);

    Pl = Al + Bl - 1;

    /* work out data for each chunk of the output */
    Pchunks = (_chunk_struct *) TMP_ALLOC(Pl*sizeof(_chunk_struct));
    perm = (slong *) TMP_ALLOC(Pl*sizeof(slong));
    for (Pi = 0; Pi < Pl; Pi++)
    {
        nmod_mpoly_init3((Pchunks + Pi)->poly, 8, P->bits, ctx);
        (Pchunks + Pi)->work = 0;
        perm[Pi] = Pi;
        for (i = 0, j = Pi; i < Al && j >= 0; i++, j--)
        {
            if (j >= Bl)
                continue;
            
            (Pchunks + Pi)->work += (Amain[i + 1] - Amain[i])
                                   *(Bmain[j + 1] - Bmain[j]);
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

    base->nthreads = num_handles + 1;
    base->Al = Al;
    base->Bl = Bl;
    base->Pl = Pl;
    base->Acoeffs = A->coeffs;
    base->Amain = Amain;
    base->Apexp = Apexp;
    base->Bcoeffs = B->coeffs;
    base->Bmain = Bmain;
    base->Bpexp = Bpexp;
    base->idx = 0;
    base->perm = perm;
    base->nvars = nvars;
    base->ctx = ctx;
    base->Pchunks = Pchunks;
    base->array_size = array_size;
    base->mults = mults;

    args = (_worker_arg_struct *) TMP_ALLOC(base->nthreads
                                                  *sizeof(_worker_arg_struct));

#if FLINT_USES_PTHREAD
    pthread_mutex_init(&base->mutex, NULL);
#endif
    for (i = 0; i < num_handles; i++)
    {
        args[i].idx = i;
        args[i].base = base;
        thread_pool_wake(global_thread_pool, handles[i], 0,
                          _nmod_mpoly_mul_array_threaded_worker_LEX, &args[i]);
    }
    i = num_handles;
    args[i].idx = i;
    args[i].base = base;
    _nmod_mpoly_mul_array_threaded_worker_LEX(&args[i]);
    for (i = 0; i < num_handles; i++)
    {
        thread_pool_wait(global_thread_pool, handles[i]);
    }
#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&base->mutex);
#endif

    /* join answers */
    Plen = 0;
    for (Pi = 0; Pi < Pl; Pi++)
    {
        _nmod_mpoly_fit_length(&P->coeffs, &P->coeffs_alloc,
                      &P->exps, &P->exps_alloc, 1, Plen + (Pchunks + Pi)->len);

        FLINT_ASSERT((Pchunks + Pi)->poly->coeffs != NULL);
        FLINT_ASSERT((Pchunks + Pi)->poly->exps != NULL);

        memcpy(P->exps + Plen, (Pchunks + Pi)->poly->exps, (Pchunks + Pi)->len*sizeof(ulong));
        memcpy(P->coeffs + Plen, (Pchunks + Pi)->poly->coeffs, (Pchunks + Pi)->len*sizeof(mp_limb_t));

        Plen += (Pchunks + Pi)->len;

        flint_free((Pchunks + Pi)->poly->coeffs);
        flint_free((Pchunks + Pi)->poly->exps);
    }

    P->length = Plen;

    flint_free(Apexp);
    flint_free(Bpexp);
    TMP_END;
}


int _nmod_mpoly_mul_array_threaded_pool_LEX(
    nmod_mpoly_t A,
    const nmod_mpoly_t B, fmpz * maxBfields,
    const nmod_mpoly_t C, fmpz * maxCfields,
    const nmod_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    slong i, exp_bits, array_size;
    ulong max, * mults;
    int success;
    TMP_INIT;

    FLINT_ASSERT(B->length != 0);
    FLINT_ASSERT(C->length != 0);

    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    FLINT_ASSERT(1 == mpoly_words_per_exp(B->bits, ctx->minfo));
    FLINT_ASSERT(1 == mpoly_words_per_exp(C->bits, ctx->minfo));

    TMP_START;

    /* compute maximum exponents for each variable */
    mults = (ulong *) TMP_ALLOC(ctx->minfo->nfields*sizeof(ulong));

    /* the field of index n-1 is the one that wil be pulled out */
    i = ctx->minfo->nfields - 1;
    FLINT_ASSERT(fmpz_fits_si(maxBfields + i));
    FLINT_ASSERT(fmpz_fits_si(maxCfields + i));
    mults[i] = 1 + fmpz_get_ui(maxBfields + i) + fmpz_get_ui(maxCfields + i);
    max = mults[i];
    if (((slong) mults[i]) <= 0 || mults[i] > MAX_LEX_SIZE)
    {
        success = 0;
        goto cleanup;
    }

    /* the fields of index n-2...0, contribute to the array size */
    array_size = WORD(1);
    for (i--; i >= 0; i--)
    {
        ulong hi;
        FLINT_ASSERT(fmpz_fits_si(maxBfields + i));
        FLINT_ASSERT(fmpz_fits_si(maxCfields + i));
        mults[i] = 1 + fmpz_get_ui(maxBfields + i) + fmpz_get_ui(maxCfields + i);
        max |= mults[i];
        umul_ppmm(hi, array_size, array_size, mults[i]);
        if (hi != 0 || (slong) mults[i] <= 0
                    || array_size <= 0
                    || array_size > MAX_ARRAY_SIZE)
        {
            success = 0;
            goto cleanup;
        }
    }

    exp_bits = FLINT_MAX(MPOLY_MIN_BITS, FLINT_BIT_COUNT(max) + 1);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);

    /* array multiplication assumes result fits into 1 word */
    if (1 != mpoly_words_per_exp(exp_bits, ctx->minfo))
    {
        success = 0;
        goto cleanup;
    }

    /* handle aliasing and do array multiplication */
    if (A == B || A == C)
    {
        nmod_mpoly_t T;
        nmod_mpoly_init3(T, B->length + C->length - 1, exp_bits, ctx);
        _nmod_mpoly_mul_array_chunked_threaded_LEX(T, C, B, mults, ctx,
                                                         handles, num_handles);
        nmod_mpoly_swap(T, A, ctx);
        nmod_mpoly_clear(T, ctx);
    }
    else
    {
        nmod_mpoly_fit_length_reset_bits(A, B->length + C->length - 1, exp_bits, ctx);
        _nmod_mpoly_mul_array_chunked_threaded_LEX(A, C, B, mults, ctx,
                                                         handles, num_handles);
    }
    success = 1;

cleanup:

    TMP_END;

    return success;
}




/*****************************
    DEGLEX and DEGREVLEX
*****************************/

static void _nmod_mpoly_mul_array_threaded_worker_DEG(void * varg)
{
    slong i, j, Pi;
    _worker_arg_struct * arg = (_worker_arg_struct *) varg;
    _base_struct * base = arg->base;
    slong  Al = base->Al;
    slong Bl = base->Bl;
    slong Pl = base->Pl;
    slong * Amain = base->Amain;
    slong * Bmain = base->Bmain;
    ulong * coeff_array;
    slong (* upack_sm1)(nmod_mpoly_t, slong, ulong *, slong, slong, slong, const nmod_mpoly_ctx_t);
    slong (* upack_sm2)(nmod_mpoly_t, slong, ulong *, slong, slong, slong, const nmod_mpoly_ctx_t);
    slong (* upack_sm3)(nmod_mpoly_t, slong, ulong *, slong, slong, slong, const nmod_mpoly_ctx_t);
    TMP_INIT;

    upack_sm1 = &nmod_mpoly_append_array_sm1_DEGLEX;
    upack_sm2 = &nmod_mpoly_append_array_sm2_DEGLEX;
    upack_sm3 = &nmod_mpoly_append_array_sm3_DEGLEX;
    if (base->rev)
    {
        upack_sm1 = &nmod_mpoly_append_array_sm1_DEGREVLEX;
        upack_sm2 = &nmod_mpoly_append_array_sm2_DEGREVLEX;
        upack_sm3 = &nmod_mpoly_append_array_sm3_DEGREVLEX;
    }

    TMP_START;
    coeff_array = (ulong *) TMP_ALLOC(3*base->array_size*sizeof(ulong));
    for (j = 0; j < 3*base->array_size; j++)
        coeff_array[j] = 0;

#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&base->mutex);
#endif
    Pi = base->idx;
    base->idx = Pi + 1;
#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&base->mutex);
#endif

    while (Pi < Pl)
    {
        slong len;
        mp_limb_t t2, t1, t0, u1, u0;

        Pi = base->perm[Pi];

        /* work out bit counts for this chunk */
        len = 0;
        for (i = 0, j = Pi; i < Al && j >= 0; i++, j--)
        {
            if (j < Bl)
            {
                len += FLINT_MIN(Amain[i + 1] - Amain[i],
                                 Bmain[j + 1] - Bmain[j]);
            }
        }

        umul_ppmm(t1, t0, base->ctx->mod.n - 1, base->ctx->mod.n - 1);
        umul_ppmm(t2, t1, t1, len);
        umul_ppmm(u1, u0, t0, len);
        add_sssaaaaaa(t2, t1, t0,  t2, t1, UWORD(0),  UWORD(0), u1, u0);

        (base->Pchunks + Pi)->len = 0;

        if (t2 != 0)
        {
            /* need three words */
            for (i = 0, j = Pi; i < Al && j >= 0; i++, j--)
            {
                if (j >= Bl)
                    continue;

                _nmod_mpoly_addmul_array1_ulong3(coeff_array, 
                        base->Acoeffs + base->Amain[i],
                            base->Apexp + base->Amain[i],
                            base->Amain[i + 1] - base->Amain[i],
                        base->Bcoeffs + base->Bmain[j],
                            base->Bpexp + base->Bmain[j],
                            base->Bmain[j + 1] - base->Bmain[j]);
            }

            (base->Pchunks + Pi)->len = upack_sm3((base->Pchunks + Pi)->poly, 0,
                       coeff_array, Pl - Pi - 1, base->nvars, base->degb, base->ctx);
        }
        else if (t1 != 0)
        {
            /* fits into two words */
            for (i = 0, j = Pi; i < Al && j >= 0; i++, j--)
            {
                if (j >= Bl)
                    continue;

                _nmod_mpoly_addmul_array1_ulong2(coeff_array, 
                        base->Acoeffs + base->Amain[i],
                            base->Apexp + base->Amain[i],
                            base->Amain[i + 1] - base->Amain[i],
                        base->Bcoeffs + base->Bmain[j],
                            base->Bpexp + base->Bmain[j],
                            base->Bmain[j + 1] - base->Bmain[j]);
            }

            (base->Pchunks + Pi)->len = upack_sm2((base->Pchunks + Pi)->poly, 0,
                       coeff_array, Pl - Pi - 1, base->nvars, base->degb, base->ctx);
        }
        else if (t0 != 0)
        {
            /* fits into one word */
            for (i = 0, j = Pi; i < Al && j >= 0; i++, j--)
            {
                if (j >= Bl)
                    continue;

                _nmod_mpoly_addmul_array1_ulong1(coeff_array, 
                        base->Acoeffs + base->Amain[i],
                            base->Apexp + base->Amain[i],
                            base->Amain[i + 1] - base->Amain[i],
                        base->Bcoeffs + base->Bmain[j],
                            base->Bpexp + base->Bmain[j],
                            base->Bmain[j + 1] - base->Bmain[j]);
            }

            (base->Pchunks + Pi)->len = upack_sm1((base->Pchunks + Pi)->poly, 0,
                       coeff_array, Pl - Pi - 1, base->nvars, base->degb, base->ctx);
        }

#if FLINT_USES_PTHREAD
	pthread_mutex_lock(&base->mutex);
#endif
        Pi = base->idx;
        base->idx = Pi + 1;
#if FLINT_USES_PTHREAD
        pthread_mutex_unlock(&base->mutex);
#endif
    }

    TMP_END;
}



void _nmod_mpoly_mul_array_chunked_threaded_DEG(
    nmod_mpoly_t P,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    ulong degb,
    const nmod_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    slong nvars = ctx->minfo->nvars;
    slong Pi, i, j, Plen, Pl, Al, Bl, array_size;
    slong * Amain, * Bmain;
    ulong * Apexp, * Bpexp;
    _base_t base;
    _worker_arg_struct * args;
    _chunk_struct * Pchunks;
    slong * perm;
    TMP_INIT;

    /* compute lengths of poly2 and poly3 in chunks */
    Al = 1 + (slong) (A->exps[0] >> (A->bits*nvars));
    Bl = 1 + (slong) (B->exps[0] >> (B->bits*nvars));

    array_size = 1;
    for (i = 0; i < nvars-1; i++) {
        array_size *= degb;
    }

    TMP_START;

    /* compute indices and lengths of coefficients of polys in main variable */
    Amain = (slong *) TMP_ALLOC((Al + 1)*sizeof(slong));
    Bmain = (slong *) TMP_ALLOC((Bl + 1)*sizeof(slong));
    Apexp = (ulong *) flint_malloc(A->length*sizeof(ulong));
    Bpexp = (ulong *) flint_malloc(B->length*sizeof(ulong));
    mpoly_main_variable_split_DEG(Amain, Apexp, A->exps, Al, A->length,
                                                         degb, nvars, A->bits);
    mpoly_main_variable_split_DEG(Bmain, Bpexp, B->exps, Bl, B->length,
                                                         degb, nvars, B->bits);

    Pl = Al + Bl - 1;
    FLINT_ASSERT(Pl == degb);

    /* work out data for each chunk of the output */
    Pchunks = (_chunk_struct *) TMP_ALLOC(Pl*sizeof(_chunk_struct));
    perm = (slong *) TMP_ALLOC(Pl*sizeof(slong));
    for (Pi = 0; Pi < Pl; Pi++)
    {
        nmod_mpoly_init3((Pchunks + Pi)->poly, 8, P->bits, ctx);
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

    base->nthreads = num_handles + 1;
    base->Al = Al;
    base->Bl = Bl;
    base->Pl = Pl;
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
    base->ctx = ctx;
    base->array_size = array_size;
    base->degb = degb;
    base->rev = (ctx->minfo->ord == ORD_DEGREVLEX);

    args = (_worker_arg_struct *) TMP_ALLOC(base->nthreads
                                                  *sizeof(_worker_arg_struct));

#if FLINT_USES_PTHREAD
    pthread_mutex_init(&base->mutex, NULL);
#endif
    for (i = 0; i < num_handles; i++)
    {
        args[i].idx = i;
        args[i].base = base;

        thread_pool_wake(global_thread_pool, handles[i], 0,
                          _nmod_mpoly_mul_array_threaded_worker_DEG, &args[i]);
    }
    i = num_handles;
    args[i].idx = i;
    args[i].base = base;
    _nmod_mpoly_mul_array_threaded_worker_DEG(&args[i]);
    for (i = 0; i < num_handles; i++)
    {
        thread_pool_wait(global_thread_pool, handles[i]);
    }
#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&base->mutex);
#endif

    /* join answers */
    Plen = 0;
    for (Pi = 0; Pi < Pl; Pi++)
    {
        _nmod_mpoly_fit_length(&P->coeffs, &P->coeffs_alloc,
                      &P->exps, &P->exps_alloc, 1, Plen + (Pchunks + Pi)->len);

        FLINT_ASSERT((Pchunks + Pi)->poly->coeffs != NULL);
        FLINT_ASSERT((Pchunks + Pi)->poly->exps != NULL);

        memcpy(P->exps + Plen, (Pchunks + Pi)->poly->exps, (Pchunks + Pi)->len*sizeof(ulong));
        memcpy(P->coeffs + Plen, (Pchunks + Pi)->poly->coeffs, (Pchunks + Pi)->len*sizeof(mp_limb_t));

        Plen += (Pchunks + Pi)->len;

        flint_free((Pchunks + Pi)->poly->coeffs);
        flint_free((Pchunks + Pi)->poly->exps);
    }

    P->length = Plen;

    flint_free(Apexp);
    flint_free(Bpexp);
    TMP_END;
}

int _nmod_mpoly_mul_array_threaded_pool_DEG(
    nmod_mpoly_t A,
    const nmod_mpoly_t B, fmpz * maxBfields,
    const nmod_mpoly_t C, fmpz * maxCfields,
    const nmod_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    slong i, exp_bits, array_size;
    ulong deg;
    int success;

    FLINT_ASSERT(B->length != 0);
    FLINT_ASSERT(C->length != 0);

    FLINT_ASSERT(  ctx->minfo->ord == ORD_DEGREVLEX
                || ctx->minfo->ord == ORD_DEGLEX);

    FLINT_ASSERT(1 == mpoly_words_per_exp(B->bits, ctx->minfo));
    FLINT_ASSERT(1 == mpoly_words_per_exp(C->bits, ctx->minfo));

    /* the field of index n-1 is the one that wil be pulled out */
    i = ctx->minfo->nfields - 1;
    FLINT_ASSERT(fmpz_fits_si(maxBfields + i));
    FLINT_ASSERT(fmpz_fits_si(maxCfields + i));
    deg = 1 + fmpz_get_ui(maxBfields + i) + fmpz_get_ui(maxCfields + i);
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

    exp_bits = FLINT_MAX(MPOLY_MIN_BITS, FLINT_BIT_COUNT(deg) + 1);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);

    /* array multiplication assumes result fit into 1 word */
    if (1 != mpoly_words_per_exp(exp_bits, ctx->minfo))
    {
        success = 0;
        goto cleanup;
    }

    /* handle aliasing and do array multiplication */
    if (A == B || A == C)
    {
        nmod_mpoly_t T;
        nmod_mpoly_init3(T, B->length + C->length - 1, exp_bits, ctx);
        _nmod_mpoly_mul_array_chunked_threaded_DEG(T, C, B, deg, ctx,
                                                         handles, num_handles);
        nmod_mpoly_swap(T, A, ctx);
        nmod_mpoly_clear(T, ctx);
    }
    else
    {
        nmod_mpoly_fit_length_reset_bits(A, B->length + C->length - 1, exp_bits, ctx);
        _nmod_mpoly_mul_array_chunked_threaded_DEG(A, C, B, deg, ctx,
                                                         handles, num_handles);
    }
    success = 1;

cleanup:

    return success;
}


int nmod_mpoly_mul_array_threaded(
    nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_t C,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    int success;
    fmpz * maxBfields, * maxCfields;
    thread_pool_handle * handles;
    slong num_handles;
    slong thread_limit = FLINT_MIN(B->length, C->length)/16;
    TMP_INIT;

    if (B->length == 0 || C->length == 0)
    {
        nmod_mpoly_zero(A, ctx);
        return 1;
    }

    if (ctx->minfo->nvars < 1 ||
        1 != mpoly_words_per_exp(B->bits, ctx->minfo) ||
        1 != mpoly_words_per_exp(C->bits, ctx->minfo))
    {
        return 0;
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

    num_handles = flint_request_threads(&handles, thread_limit);

    switch (ctx->minfo->ord)
    {
        case ORD_LEX:
        {
            success = _nmod_mpoly_mul_array_threaded_pool_LEX(A,
                      B, maxBfields, C, maxCfields, ctx, handles, num_handles);
            break;
        }
        case ORD_DEGREVLEX:
        case ORD_DEGLEX:
        {
            success = _nmod_mpoly_mul_array_threaded_pool_DEG(A,
                      B, maxBfields, C, maxCfields, ctx, handles, num_handles);
            break;
        }
        default:
        {
            success = 0;
            break;
        }
    }

    flint_give_back_threads(handles, num_handles);

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_clear(maxBfields + i);
        fmpz_clear(maxCfields + i);
    }

    TMP_END;
    return success;
}
