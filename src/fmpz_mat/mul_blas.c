/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdint.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "ulong_extras.h"
#include "n_poly.h"
#include "mpn_extras.h"

#if FLINT_USES_BLAS && FLINT_BITS == 64

#include "cblas.h"


typedef struct {
    slong m;
    slong k;
    slong n;
    slong Astartrow;
    slong Astoprow;
    slong Bstartrow;
    slong Bstoprow;
    fmpz ** Arows;
    fmpz ** Brows;
    double * dA;
    double * dB;
} _red_worker_arg;

static void _red_worker(void * varg)
{
    _red_worker_arg * arg = (_red_worker_arg *) varg;
    slong i, j;
    slong k = arg->k;
    slong n = arg->n;
    slong Astartrow = arg->Astartrow;
    slong Astoprow = arg->Astoprow;
    slong Bstartrow = arg->Bstartrow;
    slong Bstoprow = arg->Bstoprow;
    fmpz ** Arows = arg->Arows;
    fmpz ** Brows = arg->Brows;
    double * dA = arg->dA;
    double * dB = arg->dB;

    for (i = Astartrow; i < Astoprow; i++)
        for (j = 0; j < k; j++)
            dA[k*i + j] = (double)(Arows[i][j]);

    for (i = Bstartrow; i < Bstoprow; i++)
        for (j = 0; j < n; j++)
            dB[n*i + j] = (double)(Brows[i][j]);
}

static int _fmpz_mat_mul_blas_direct(
    fmpz_mat_t C,
    const fmpz_mat_t A,
    const fmpz_mat_t B)
{
    slong i, j, start, stop;
    slong m = A->r;
    slong k = A->c;
    slong n = B->c;
    double * dC, * dA, * dB;
    slong limit;
    _red_worker_arg mainarg;
    _red_worker_arg * args;
    slong num_workers;
    thread_pool_handle * handles;

    dA = FLINT_ARRAY_ALLOC(m*k, double);
    dB = FLINT_ARRAY_ALLOC(k*n, double);
    dC = (double *) flint_calloc(m*n, sizeof(double));

    mainarg.m = m = A->r;
    mainarg.k = k = A->c;
    mainarg.n = n = B->c;
    mainarg.Arows = A->rows;
    mainarg.Brows = B->rows;
    mainarg.dA = dA;
    mainarg.dB = dB;

    limit = m + k + n;
    limit = (limit < 300) ? 0 : (limit - 300)/128;

    if (limit < 1)
    {
red_single:
        mainarg.Astartrow = 0;
        mainarg.Astoprow = m;
        mainarg.Bstartrow = 0;
        mainarg.Bstoprow = k;
        _red_worker(&mainarg);
    }
    else
    {
        num_workers = flint_request_threads(&handles, limit);
        if (num_workers < 1)
        {
            flint_give_back_threads(handles, num_workers);
            goto red_single;
        }

        args = FLINT_ARRAY_ALLOC(num_workers, _red_worker_arg);
        for (start = 0, i = 0; i < num_workers; start = stop, i++)
        {
            args[i] = mainarg;
            stop = _thread_pool_find_work_2(m, k, k, n, i + 1, num_workers + 1);
            _thread_pool_distribute_work_2(start, stop,
                                     &args[i].Astartrow, &args[i].Astoprow, m,
                                     &args[i].Bstartrow, &args[i].Bstoprow, k);
        }

        _thread_pool_distribute_work_2(start, m + k,
                                     &mainarg.Astartrow, &mainarg.Astoprow, m,
                                     &mainarg.Bstartrow, &mainarg.Bstoprow, k);

        for (i = 0; i < num_workers; i++)
            thread_pool_wake(global_thread_pool, handles[i], 0, _red_worker, &args[i]);
        _red_worker(&mainarg);
        for (i = 0; i < num_workers; i++)
            thread_pool_wait(global_thread_pool, handles[i]);

        flint_give_back_threads(handles, num_workers);
        flint_free(args);
    }

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                                       m, n, k, 1.0, dA, k, dB, n, 0.0, dC, n);

    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            fmpz_set_si(&C->rows[i][j], (slong)(dC[n*i + j]));

    flint_free(dA);
    flint_free(dB);
    flint_free(dC);

    return 1;
}


/* Stuff that should be shared with nmod_mat_mul_blas */

#define MAX_BLAS_DP_INT  (UWORD(1) << 53)

/* mod/lift helpers */

static void _lift_vec(double * a, const uint32_t * b, slong len, uint32_t n)
{
    slong i;
    for (i = 0; i < len; i++)
        a[i] = (int32_t)(b[i] - (n & (-(uint32_t)((int32_t)(n/2 - b[i]) < 0))));
}

static uint32_t _reduce_uint32(mp_limb_t a, nmod_t mod)
{
    mp_limb_t r;
    NMOD_RED(r, a, mod);
    return (uint32_t)r;
}

static void fmpz_multi_mod_uint32_stride(
    uint32_t * out, slong stride,
    const fmpz_t input,
    const fmpz_comb_t C,
    fmpz_comb_temp_t CT)
{
    slong i, k, l;
    fmpz * A = CT->A;
    mod_lut_entry * lu;
    slong * offsets;
    slong klen = C->mod_klen;
    fmpz_t ttt;

    /* high level split */
    if (klen == 1)
    {
        *ttt = A[0];
        A[0] = *input;
    }
    else
    {
        _fmpz_multi_mod_precomp(A, C->mod_P, input, -1, CT->T);
    }

    offsets = C->mod_offsets;
    lu = C->mod_lu;

    for (k = 0, i = 0, l = 0; k < klen; k++)
    {
        slong j = offsets[k];

        for ( ; i < j; i++)
        {
            /* mid level split: depends on FMPZ_MOD_UI_CUTOFF */
            mp_limb_t t = fmpz_get_nmod(A + k, lu[i].mod);

            /* low level split: 1, 2, or 3 small primes */
            if (lu[i].mod2.n != 0)
            {
                FLINT_ASSERT(l + 3 <= C->num_primes);
                out[l*stride] = _reduce_uint32(t, lu[i].mod0); l++;
                out[l*stride] = _reduce_uint32(t, lu[i].mod1); l++;
                out[l*stride] = _reduce_uint32(t, lu[i].mod2); l++;
            }
            else if (lu[i].mod1.n != 0)
            {
                FLINT_ASSERT(l + 2 <= C->num_primes);
                out[l*stride] = _reduce_uint32(t, lu[i].mod0); l++;
                out[l*stride] = _reduce_uint32(t, lu[i].mod1); l++;
            }
            else
            {
                FLINT_ASSERT(l + 1 <= C->num_primes);
                out[l*stride] = (uint32_t)(t); l++;
            }
        }
    }

    FLINT_ASSERT(l == C->num_primes);

    if (klen == 1)
        A[0] = *ttt;
}

/* workers */

typedef struct {
    mp_limb_t prime;
    slong l;
    slong num_primes;
    slong m;
    slong k;
    slong n;
    slong Astartrow;
    slong Astoprow;
    slong Bstartrow;
    slong Bstoprow;
    slong Cstartrow;
    slong Cstoprow;
    uint32_t * bigA;
    uint32_t * bigB;
    uint32_t * bigC;
    double * dA;
    double * dB;
    double * dC;
    fmpz ** Arows;
    fmpz ** Brows;
    fmpz ** Crows;
    const fmpz_comb_struct * comb;
    int sign;
} _worker_arg;

static void _mod_worker(void * arg_ptr)
{
    _worker_arg * arg = (_worker_arg *) arg_ptr;
    slong i, j;
    slong num_primes = arg->num_primes;
    slong k = arg->k;
    slong n = arg->n;
    slong Astartrow = arg->Astartrow;
    slong Astoprow = arg->Astoprow;
    slong Bstartrow = arg->Bstartrow;
    slong Bstoprow = arg->Bstoprow;
    uint32_t * bigA = arg->bigA;
    uint32_t * bigB = arg->bigB;
    fmpz ** Arows = arg->Arows;
    fmpz ** Brows = arg->Brows;
    const fmpz_comb_struct * comb = arg->comb;
    fmpz_comb_temp_t comb_temp;

    fmpz_comb_temp_init(comb_temp, comb);

    for (i = Astartrow; i < Astoprow; i++)
        for (j = 0; j < k; j++)
            fmpz_multi_mod_uint32_stride(bigA + i*k*num_primes + j, k,
                                                &Arows[i][j], comb, comb_temp);

    for (i = Bstartrow; i < Bstoprow; i++)
        for (j = 0; j < n; j++)
            fmpz_multi_mod_uint32_stride(bigB + i*n*num_primes + j, n,
                                                &Brows[i][j], comb, comb_temp);

    fmpz_comb_temp_clear(comb_temp);
}

void _tod_worker(void * arg_ptr)
{
    _worker_arg * arg = (_worker_arg *) arg_ptr;
    slong i;
    slong l = arg->l;
    slong num_primes = arg->num_primes;
    slong k = arg->k;
    slong n = arg->n;
    slong Astartrow = arg->Astartrow;
    slong Astoprow = arg->Astoprow;
    slong Bstartrow = arg->Bstartrow;
    slong Bstoprow = arg->Bstoprow;
    const uint32_t * bigA = arg->bigA;
    const uint32_t * bigB = arg->bigB;
    double * dA = arg->dA;
    double * dB = arg->dB;
    uint32_t prime = arg->prime;

    for (i = Astartrow; i < Astoprow; i++)
        _lift_vec(dA + i*k, bigA + l*k + i*k*num_primes, k, prime);

    for (i = Bstartrow; i < Bstoprow; i++)
        _lift_vec(dB + i*n, bigB + l*n + i*n*num_primes, n, prime);
}

void _fromd_worker(void * arg_ptr)
{
    _worker_arg * arg = (_worker_arg *) arg_ptr;
    slong i, j;
    slong l = arg->l;
    slong num_primes = arg->num_primes;
    slong n = arg->n;
    slong Cstartrow = arg->Cstartrow;
    slong Cstoprow = arg->Cstoprow;
    uint32_t * bigC = arg->bigC;
    double * dC = arg->dC;
    ulong shift;
    nmod_t mod;

    nmod_init(&mod, arg->prime);

    shift = ((2*MAX_BLAS_DP_INT)/mod.n)*mod.n;

    for (i = Cstartrow; i < Cstoprow; i++)
    {
        for (j = 0; j < n; j++)
        {
            mp_limb_t r;
            slong a = (slong) dC[i*n + j];
            mp_limb_t b = (a < 0) ? a + shift : a;
            NMOD_RED(r, b, mod);
            bigC[n*(num_primes*i + l) + j] = r;
        }
    }
}

void _crt_worker(void * arg_ptr)
{
    _worker_arg * arg = (_worker_arg *) arg_ptr;
    slong i, j, k;
    slong num_primes = arg->num_primes;
    slong n = arg->n;
    slong Cstartrow = arg->Cstartrow;
    slong Cstoprow = arg->Cstoprow;
    uint32_t * bigC = arg->bigC;
    fmpz ** Crows = arg->Crows;
    const fmpz_comb_struct * comb = arg->comb;
    fmpz_comb_temp_t comb_temp;
    mp_limb_t * r;
    int sign = arg->sign;

    fmpz_comb_temp_init(comb_temp, comb);
    r = FLINT_ARRAY_ALLOC(num_primes, mp_limb_t);

    for (i = Cstartrow; i < Cstoprow; i++)
    {
        for (j = 0; j < n; j++)
        {
            for (k = 0; k < num_primes; k++)
                r[k] = bigC[(i*num_primes + k)*n + j];

            fmpz_multi_CRT_ui(&Crows[i][j], r, comb, comb_temp, sign);
        }
    }

    flint_free(r);
    fmpz_comb_temp_clear(comb_temp);
}

static mp_limb_t * _calculate_primes(
    slong * num_primes_,
    flint_bitcnt_t bits,
    slong k)
{
    slong num_primes, primes_alloc;
    mp_limb_t * primes;
    mp_limb_t p;
    fmpz_t prod;

    p = 2 + 2*n_sqrt((MAX_BLAS_DP_INT - 1)/(ulong)k);
    if (bits > 200)
    {
        /* if mod is the bottleneck, ensure p1*p2*p3 < 2^62 */
        p = FLINT_MIN(p, UWORD(1664544));
    }

    primes_alloc = 1 + bits/FLINT_BIT_COUNT(p);
    primes = FLINT_ARRAY_ALLOC(primes_alloc, mp_limb_t);
    num_primes = 0;

    fmpz_init_set_ui(prod, 1);

    do {
        do {
            if (p < 1000)
            {
                fmpz_clear(prod);
                flint_free(primes);
                *num_primes_ = 0;
                return NULL;
            }
            p--;
        } while (!n_is_prime(p));

        if (num_primes + 1 > primes_alloc)
        {
            primes_alloc = FLINT_MAX(num_primes + 1, primes_alloc*5/4);
            primes = FLINT_ARRAY_REALLOC(primes, primes_alloc, mp_limb_t);
        }

        primes[num_primes] = p;
        num_primes++;

        fmpz_mul_ui(prod, prod, p);

    } while (fmpz_bits(prod) <= bits);

    fmpz_clear(prod);

    *num_primes_ = num_primes;
    return primes;
}

/*
    max|A| < 2^Abits
    max|B| < 2^Bbits
    max|C| < 2^Cbits

    sign = 1:   either A or B could have negative entries.
    sign = 0:   all entries are >= 0.
*/
int _fmpz_mat_mul_blas(
    fmpz_mat_t C,
    const fmpz_mat_t A, flint_bitcnt_t Abits,
    const fmpz_mat_t B, flint_bitcnt_t Bbits,
    int sign,
    flint_bitcnt_t Cbits)
{
    slong i, l, start, stop;
    slong m = A->r;
    slong k = A->c;
    slong n = B->c;
    uint32_t * bigC, * bigA, * bigB;
    double * dC, * dA, * dB;
    mp_limb_t * primes;
    slong num_primes;
    fmpz_comb_t comb;
    thread_pool_handle * handles;
    slong num_workers;
    _worker_arg * args;

    FLINT_ASSERT(sign == 0 || sign == 1);
    FLINT_ASSERT(m == A->r && m == C->r);
    FLINT_ASSERT(k == A->c && k == B->r);
    FLINT_ASSERT(n == B->c && n == C->c);

    if (m < 1 || k < 1 || n < 1 || m > INT_MAX || k > INT_MAX || n > INT_MAX)
        return 0;

    if (Abits + Bbits + FLINT_BIT_COUNT(k) <= 53)
        return _fmpz_mat_mul_blas_direct(C, A, B);

    primes = _calculate_primes(&num_primes, Cbits + sign, k);
    if (primes == NULL)
        return 0;

    fmpz_comb_init(comb, primes, num_primes);

    /*
        To allow the 3D transpose to and from blas to be at least slightly
        cache friendly, matrices M mod p[l] are stored with the prime index
        in the middle:
            M[i,j] mod p[l] is at bigM[(i*num_primes + l)*M->c + j]
    */
    bigA = FLINT_ARRAY_ALLOC(m*k*num_primes, uint32_t);
    bigB = FLINT_ARRAY_ALLOC(k*n*num_primes, uint32_t);
    bigC = FLINT_ARRAY_ALLOC(m*n*num_primes, uint32_t);
    dA = FLINT_ARRAY_ALLOC(m*k, double);
    dB = FLINT_ARRAY_ALLOC(k*n, double);
    dC = (double *) flint_calloc(m*n, sizeof(double));

    num_workers = flint_request_threads(&handles, INT_MAX);

    args = FLINT_ARRAY_ALLOC(num_workers + 1, _worker_arg);
    for (start = 0, i = 0; i <= num_workers; start = stop, i++)
    {
        args[i].l = -1;
        args[i].prime = 0;
        args[i].num_primes = num_primes;
        args[i].m = m;
        args[i].k = k;
        args[i].n = n;
        args[i].bigA = bigA;
        args[i].bigB = bigB;
        args[i].bigC = bigC;
        args[i].Arows = A->rows;
        args[i].Brows = B->rows;
        args[i].Crows = C->rows;
        args[i].dA = dA;
        args[i].dB = dB;
        args[i].dC = dC;
        args[i].comb = comb;
        args[i].sign = sign;

        /* distribute rows of C evenly */
        args[i].Cstartrow = ((i + 0)*m)/(num_workers + 1);
        args[i].Cstoprow  = ((i + 1)*m)/(num_workers + 1);

        /* distribute rows of A and B evenly as weighted by their columns */
        stop = _thread_pool_find_work_2(m, k, k, n, i + 1, num_workers + 1);
        _thread_pool_distribute_work_2(start, stop,
                          &args[i].Astartrow, &args[i].Astoprow, m,
                          &args[i].Bstartrow, &args[i].Bstoprow, k);
    }

    /* mod inputs */
    for (i = 0; i < num_workers; i++)
        thread_pool_wake(global_thread_pool, handles[i], 0, _mod_worker, &args[i]);
    _mod_worker(&args[num_workers]);
    for (i = 0; i < num_workers; i++)
        thread_pool_wait(global_thread_pool, handles[i]);

    /* multiply and reduce answers mod the primes */
    for (l = 0; l < num_primes; l++)
    {
        for (i = 0; i <= num_workers; i++)
        {
            args[i].l = l;
            args[i].prime = primes[l];
        }

        for (i = 0; i < num_workers; i++)
            thread_pool_wake(global_thread_pool, handles[i], 0, _tod_worker, &args[i]);
        _tod_worker(&args[num_workers]);
        for (i = 0; i < num_workers; i++)
            thread_pool_wait(global_thread_pool, handles[i]);

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                                       m, n, k, 1.0, dA, k, dB, n, 0.0, dC, n);

        for (i = 0; i < num_workers; i++)
            thread_pool_wake(global_thread_pool, handles[i], 0, _fromd_worker, &args[i]);
        _fromd_worker(&args[num_workers]);
        for (i = 0; i < num_workers; i++)
            thread_pool_wait(global_thread_pool, handles[i]);
    }

    /* crt output */
    for (i = 0; i < num_workers; i++)
        thread_pool_wake(global_thread_pool, handles[i], 0, _crt_worker, &args[i]);
    _crt_worker(&args[num_workers]);
    for (i = 0; i < num_workers; i++)
        thread_pool_wait(global_thread_pool, handles[i]);

    flint_give_back_threads(handles, num_workers);

    /* cleanup */
    fmpz_comb_clear(comb);
    flint_free(primes);

    flint_free(args);
    flint_free(bigA);
    flint_free(bigB);
    flint_free(bigC);
    flint_free(dA);
    flint_free(dB);
    flint_free(dC);

    return 1;
}

#else

int _fmpz_mat_mul_blas(
    fmpz_mat_t C,
    const fmpz_mat_t A, flint_bitcnt_t Abits,
    const fmpz_mat_t B, flint_bitcnt_t Bbits,
    int sign,
    flint_bitcnt_t Cbits)
{
    return 0;
}

#endif


int fmpz_mat_mul_blas(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
{
    slong Abits = fmpz_mat_max_bits(A);
    slong Bbits = fmpz_mat_max_bits(B);
    flint_bitcnt_t Cbits;
    int sign = 0;

    if (Abits < 0)
    {
        sign = 1;
        Abits = -Abits;
    }

    if (Bbits < 0)
    {
        sign = 1;
        Bbits = -Bbits;
    }

    Cbits = Abits + Bbits + FLINT_BIT_COUNT(A->c);

    return _fmpz_mat_mul_blas(C, A, Abits, B, Bbits, sign, Cbits);
}

