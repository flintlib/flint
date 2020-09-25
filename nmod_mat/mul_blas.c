/*
    Copyright (C) 2020 Daniel Schultz
    This file is part of FLINT.
    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mat.h"
#include "thread_support.h"

#if FLINT_USES_BLAS && FLINT_BITS == 64

#include "cblas.h"

/*
    This code is on the edge of disaster. Blas is used for dot products

        sum_{i=1}^k a_i*b_i, |a_i|, |b_i| <= n

    and we **assume** that if k*n^2 < MAX_BLAS_{DP|SP}_INT, then blas
    will calculate this dot product exactly.

    Test code fails easily with MAX_BLAS_DP_INT = 2^54 and also fails
    with MAX_BLAS_SP_INT = 2^25.
*/

/* for dgemm */
#define MAX_BLAS_DP_INT  (UWORD(1) << 53)

/* for sgemm */
#define MAX_BLAS_SP_INT  (UWORD(1) << 24)


/* helper for distributing the input conversion work */
static void _distribute_rows(
    slong * Astartrow, slong * Astoprow,    /* limits for A */
    slong * Bstartrow, slong * Bstoprow,    /* limits for B */
    slong m,                                /* number of rows of A */
    slong start, slong stop)
{
    FLINT_ASSERT(stop <= stop);

    if (start >= m)
    {
        *Astartrow = 0;
        *Astoprow  = 0;
        *Bstartrow = start - m;
        *Bstoprow  = stop - m;
    }
    else if (stop <= m)
    {
        *Astartrow = start;
        *Astoprow  = stop;
        *Bstartrow = 0;
        *Bstoprow  = 0;
    }
    else
    {
        *Astartrow = start;
        *Astoprow  = m;
        *Bstartrow = 0;
        *Bstoprow  = stop - m;
    }
}


/************ small enough that a single sgemm suffices **********************/

static void _lift_vec_sp(float * a, ulong * b, slong len, ulong n)
{
    slong i;
    for (i = 0; i < len; i++)
        a[i] = (int)(b[i] - (n & FLINT_SIGN_EXT(n/2 - b[i])));
}

typedef struct {
    slong m;
    slong n;
    slong k;
    slong Astartrow;
    slong Astoprow;
    slong Bstartrow;
    slong Bstoprow;
    mp_limb_t ctxn;
    float * dA;
    float * dB;
    mp_limb_t ** Arows;
    mp_limb_t ** Brows;
} _lift_sp_worker_arg_struct;

void _lift_sp_worker(void * arg_ptr)
{
    _lift_sp_worker_arg_struct * arg = (_lift_sp_worker_arg_struct *) arg_ptr;
    slong n = arg->n;
    slong k = arg->k;
    slong Astartrow = arg->Astartrow;
    slong Astoprow = arg->Astoprow;
    slong Bstartrow = arg->Bstartrow;
    slong Bstoprow = arg->Bstoprow;
    mp_limb_t ctxn = arg->ctxn;
    float * dA = arg->dA;
    float * dB = arg->dB;
    mp_limb_t ** Arows = arg->Arows;
    mp_limb_t ** Brows = arg->Brows;
    slong i;

    for (i = Astartrow; i < Astoprow; i++)
        _lift_vec_sp(dA + i*k, Arows[i], k, ctxn);

    for (i = Bstartrow; i < Bstoprow; i++)
        _lift_vec_sp(dB + i*n, Brows[i], n, ctxn);
}

typedef struct {
    slong n;
    slong Cstartrow;
    slong Cstoprow;
    nmod_t * ctx;
    mp_limb_t shift;
    float * dC;
    mp_limb_t ** Crows;
} _reduce_sp_worker_arg_struct;

void _reduce_sp_worker(void * arg_ptr)
{
    _reduce_sp_worker_arg_struct * arg = (_reduce_sp_worker_arg_struct *) arg_ptr;
    slong n = arg->n;
    slong Cstartrow = arg->Cstartrow;
    slong Cstoprow = arg->Cstoprow;
    nmod_t ctx = *arg->ctx;
    mp_limb_t shift = arg->shift;
    float * dC = arg->dC;
    mp_limb_t ** Crows = arg->Crows;
    slong i, j;

    for (i = Cstartrow; i < Cstoprow; i++)
    {
        for (j = 0; j < n; j++)
        {
            slong a = (slong) dC[i*n + j];
            mp_limb_t b = (a < 0) ? a + shift : a;
            NMOD_RED(Crows[i][j], b, ctx);
        }
    }
}

static int _nmod_mat_mul_blas_sp(nmod_mat_t C,
                                        const nmod_mat_t A, const nmod_mat_t B)
{
    slong i;
    slong m = A->r;
    slong k = A->c;
    slong n = B->c;
    float * dC, * dA, * dB;
    ulong shift;
    nmod_t ctx = C->mod;
    slong num_workers;
    thread_pool_handle * handles;
    void * tmp;

    dA = flint_malloc(m*k*sizeof(float));
    dB = flint_malloc(k*n*sizeof(float));
    dC = flint_calloc(m*n, sizeof(float));

    num_workers = flint_request_threads(&handles, INT_MAX);

    tmp = flint_malloc((num_workers + 1)*FLINT_MAX(
                                        sizeof(_lift_sp_worker_arg_struct),
                                        sizeof(_reduce_sp_worker_arg_struct)));
    /* convert inputs */
    {
        _lift_sp_worker_arg_struct * args = (_lift_sp_worker_arg_struct *) tmp;

        for (i = 0; i <= num_workers; i++)
        {
            args[i].m = m;
            args[i].n = n;
            args[i].k = k;
            args[i].ctxn = ctx.n;
            args[i].dA = dA;
            args[i].dB = dB;
            args[i].Arows = A->rows;
            args[i].Brows = B->rows;
            _distribute_rows(&args[i].Astartrow, &args[i].Astoprow,
                             &args[i].Bstartrow, &args[i].Bstoprow, m,
                                          ((m + k)*(i + 0))/(num_workers + 1),
                                          ((m + k)*(i + 1))/(num_workers + 1));
        }

        for (i = 0; i < num_workers; i++)
            thread_pool_wake(global_thread_pool, handles[i], 0, _lift_sp_worker, &args[i]);
        _lift_sp_worker(&args[num_workers]);
        for (i = 0; i < num_workers; i++)
            thread_pool_wait(global_thread_pool, handles[i]);
    }

    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
                                                1.0, dA, k, dB, n, 0.0, dC, n);

    /* convert output */

    /*
        the shift for negative outputs must satisfy
            ctx.n divides shift, and
            2^24 <= shift <= 2^64
    */
    shift = ((2*MAX_BLAS_SP_INT)/ctx.n)*ctx.n;
    FLINT_ASSERT(MAX_BLAS_SP_INT <= shift);

    {
        _reduce_sp_worker_arg_struct * args = (_reduce_sp_worker_arg_struct *) tmp;

        for (i = 0; i <= num_workers; i++)
        {
            args[i].n = n;
            args[i].Cstartrow = ((i + 0)*m)/(num_workers + 1);
            args[i].Cstoprow  = ((i + 1)*m)/(num_workers + 1);
            args[i].ctx = &ctx;
            args[i].shift = shift;
            args[i].dC = dC;
            args[i].Crows = C->rows;
        }

        for (i = 0; i < num_workers; i++)
            thread_pool_wake(global_thread_pool, handles[i], 0, _reduce_sp_worker, &args[i]);
        _reduce_sp_worker(&args[num_workers]);
        for (i = 0; i < num_workers; i++)
            thread_pool_wait(global_thread_pool, handles[i]);
    }

    flint_free(tmp);

    flint_give_back_threads(handles, num_workers);

    flint_free(dA);
    flint_free(dB);
    flint_free(dC);

    return 1;
}


/******** handle larger larger moduli via several dgemm's and crt ************/

#define MAX_CRT_NUM 12

static void _lift_vec_crt(double * a, ulong * b, slong len, nmod_t ctx)
{
    slong i;
    for (i = 0; i < len; i++)
    {
        mp_limb_t bn;
        NMOD_RED(bn, b[i], ctx);
        a[i] = (int)(bn - (ctx.n & FLINT_SIGN_EXT(ctx.n/2 - bn)));
    }
}

typedef struct {
    slong m;
    slong n;
    slong k;
    slong Astartrow;
    slong Astoprow;
    slong Bstartrow;
    slong Bstoprow;
    nmod_t crtmod;
    double * dA;
    double * dB;
    mp_limb_t ** Arows;
    mp_limb_t ** Brows;
} _lift_crt_worker_arg_struct;

void _lift_crt_worker(void * arg_ptr)
{
    _lift_crt_worker_arg_struct * arg = (_lift_crt_worker_arg_struct *) arg_ptr;
    slong n = arg->n;
    slong k = arg->k;
    slong Astartrow = arg->Astartrow;
    slong Astoprow = arg->Astoprow;
    slong Bstartrow = arg->Bstartrow;
    slong Bstoprow = arg->Bstoprow;
    nmod_t crtmod = arg->crtmod;
    double * dA = arg->dA;
    double * dB = arg->dB;
    mp_limb_t ** Arows = arg->Arows;
    mp_limb_t ** Brows = arg->Brows;
    slong i;

    for (i = Astartrow; i < Astoprow; i++)
        _lift_vec_crt(dA + i*k, Arows[i], k, crtmod);

    for (i = Bstartrow; i < Bstoprow; i++)
        _lift_vec_crt(dB + i*n, Brows[i], n, crtmod);
}

typedef struct {
    slong m;
    slong n;
    slong Cstartrow;
    slong Cstoprow;
    slong crtnum;
    nmod_t * crtmod;
    nmod_t * ctx;
    double * dC;
    mp_limb_t ** Crows;
} _reduce_crt_worker_arg_struct;

void _reduce_crt_worker(void * arg_ptr)
{
    _reduce_crt_worker_arg_struct * arg = (_reduce_crt_worker_arg_struct *) arg_ptr;
    slong Cstartrow = arg->Cstartrow;
    slong Cstoprow = arg->Cstoprow;
    slong i, j, pi, pj;
    slong m = arg->m;
    slong n = arg->n;
    double * dC = arg->dC;
    nmod_t ctx = *arg->ctx;
    ulong s, t, hi, lo, reshi, reslo;
    slong crtnum = arg->crtnum;
    mp_limb_t ** Crows = arg->Crows;
    nmod_t crtmod[MAX_CRT_NUM];
    mp_limb_t q[MAX_CRT_NUM], v[MAX_CRT_NUM], u[MAX_CRT_NUM];
    mp_limb_t shifts[MAX_CRT_NUM], pmodinv[MAX_CRT_NUM*MAX_CRT_NUM];

    for (i = 0; i < crtnum; i++)
        crtmod[i] = arg->crtmod[i];

    /*
        set p_i = crtmod[i].n
        for finding u given its image u[i] mod p_i, first solve for the v[i]:
            u = v[0] + v[1]*p_0 + v[2]*p_0*p_1 + ... + v[crtnum-1]*p_0*...*p_{crtnum-1}
        then evaluate this dot product modulo ctx.n to find u mod ctx.n.

        The lower triangular matrix pmodinv is used for calculting the v[i].
    */

    for (pi = 0; pi < crtnum; pi++)
    {
        t = 1;
        s = 1;
        for (pj = pi - 1; pj >= 0; pj--)
        {
            t = nmod_mul(t, crtmod[pj].n, crtmod[pi]);
            FLINT_ASSERT(crtmod[pj].n < ctx.n);
            s = nmod_mul(s, crtmod[pj].n, ctx);
            pmodinv[MAX_CRT_NUM*pi + pj] = nmod_neg(
                                          nmod_inv(t, crtmod[pi]), crtmod[pi]);
        }
        q[pi] = s; /* q[i] = p_0 * ... * p[i-1] mod ctx.n */

        /* for double -> nmod conversion */
        shifts[pi] = ((2*MAX_BLAS_DP_INT)/crtmod[pi].n)*crtmod[pi].n;
    }

    for (i = Cstartrow; i < Cstoprow; i++)
    {
        for (j = 0; j < n; j++)
        {
            for (pi = 0; pi < crtnum; pi++)
            {
                slong a = (slong) dC[i*n + j + pi*m*n];
                mp_limb_t b = (a < 0) ? a + shifts[pi] : a;
                NMOD_RED(u[pi], b, crtmod[pi]);
            }

            reslo = u[0];
            reshi = 0;
            for (pi = 1; pi < crtnum; pi++)
            {
                FLINT_ASSERT(u[pi] < crtmod[pi].n);
                FLINT_ASSERT(u[0] < crtmod[pi].n);
                t = pmodinv[MAX_CRT_NUM*pi + 0]*nmod_sub(u[0], u[pi], crtmod[pi]);
                for (pj = 1; pj < pi; pj++)
                    t += pmodinv[MAX_CRT_NUM*pi + pj]*v[pj];
                NMOD_RED(v[pi], t, crtmod[pi]);
                umul_ppmm(hi, lo, v[pi], q[pi]);
                add_ssaaaa(reshi, reslo, reshi, reslo, hi, lo);
            }

            if (reshi < ctx.n)
                NMOD_RED2(Crows[i][j], reshi, reslo, ctx);
            else
                NMOD2_RED2(Crows[i][j], reshi, reslo, ctx);
        }
    }
}

static int _nmod_mat_mul_blas_crt(nmod_mat_t C,
                                        const nmod_mat_t A, const nmod_mat_t B)
{
    slong i, pi, crtnum;
    slong m = A->r;
    slong k = A->c;
    slong n = B->c;
    double * dC, * dA, * dB;
    nmod_t ctx = C->mod;
    ulong t;
    nmod_t crtmod[MAX_CRT_NUM]; /* not nec prime */
    fmpz_t prodmod, maxentry;
    void * tmp;
    slong num_workers;
    thread_pool_handle * handles;

    fmpz_init_set_ui(maxentry, k);
    fmpz_mul_ui(maxentry, maxentry, ctx.n - 1);
    fmpz_mul_ui(maxentry, maxentry, ctx.n - 1);

    t = n_sqrt((4*MAX_BLAS_DP_INT)/k - 1);
    fmpz_init_set_ui(prodmod, t);
    nmod_init(crtmod + 0, t);
    crtnum = 1;

    do {
        t = crtmod[crtnum - 1].n;
        do {
            if (crtnum >= MAX_CRT_NUM || t < 100)
            {
                fmpz_clear(maxentry);
                fmpz_clear(prodmod);
                return 0;
            }
            t--;
        } while (n_gcd(fmpz_fdiv_ui(prodmod, t), t) != 1);
        fmpz_mul_ui(prodmod, prodmod, t);
        nmod_init(crtmod + crtnum, t);
        crtnum++;
    } while (fmpz_cmp(prodmod, maxentry) <= 0);

    /* note that if k is sufficiently big, i.e. k >= 4, then crtnum >= 3 */

    /* arange the crt moduli in increasing order */
    for (pi = 0; pi < crtnum/2; pi++)
    {
        nmod_t tmp = crtmod[pi];
        crtmod[pi] = crtmod[crtnum - 1 - pi];
        crtmod[crtnum - 1 - pi] = tmp;
    }

    dA = flint_malloc(m*k*sizeof(double));
    dB = flint_malloc(k*n*sizeof(double));
    dC = flint_calloc(crtnum*m*n, sizeof(double));

    num_workers = flint_request_threads(&handles, INT_MAX);

    tmp = flint_malloc((num_workers + 1)*FLINT_MAX(
                                       sizeof(_lift_crt_worker_arg_struct),
                                       sizeof(_reduce_crt_worker_arg_struct)));

    for (pi = 0; pi < crtnum; pi++)
    {
        _lift_crt_worker_arg_struct * args = (_lift_crt_worker_arg_struct *) tmp;

        for (i = 0; i <= num_workers; i++)
        {
            args[i].m = m;
            args[i].n = n;
            args[i].k = k;
            args[i].crtmod = crtmod[pi];
            args[i].dA = dA;
            args[i].dB = dB;
            args[i].Arows = A->rows;
            args[i].Brows = B->rows;
            _distribute_rows(&args[i].Astartrow, &args[i].Astoprow,
                             &args[i].Bstartrow, &args[i].Bstoprow, m,
                                          ((m + k)*(i + 0))/(num_workers + 1),
                                          ((m + k)*(i + 1))/(num_workers + 1));
        }

        for (i = 0; i < num_workers; i++)
            thread_pool_wake(global_thread_pool, handles[i], 0, _lift_crt_worker, &args[i]);
        _lift_crt_worker(&args[num_workers]);
        for (i = 0; i < num_workers; i++)
            thread_pool_wait(global_thread_pool, handles[i]);

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
                                       1.0, dA, k, dB, n, 0.0, dC + pi*m*n, n);
    }

    {
        _reduce_crt_worker_arg_struct * args = (_reduce_crt_worker_arg_struct *) tmp;

        for (i = 0; i <= num_workers; i++)
        {
            args[i].m = m;
            args[i].n = n;
            args[i].Cstartrow = ((i + 0)*m)/(num_workers + 1);
            args[i].Cstoprow  = ((i + 1)*m)/(num_workers + 1);
            args[i].crtnum = crtnum;
            args[i].crtmod = crtmod;
            args[i].ctx = &ctx;
            args[i].dC = dC;
            args[i].Crows = C->rows;
        }

        for (i = 0; i < num_workers; i++)
            thread_pool_wake(global_thread_pool, handles[i], 0, _reduce_crt_worker, &args[i]);
        _reduce_crt_worker(&args[num_workers]);
        for (i = 0; i < num_workers; i++)
            thread_pool_wait(global_thread_pool, handles[i]);
    }

    flint_free(tmp);

    flint_give_back_threads(handles, num_workers);

    flint_free(dA);
    flint_free(dB);
    flint_free(dC);

    fmpz_clear(maxentry);
    fmpz_clear(prodmod);

    return 1;
}

/********** try with a single dgemm if possible ******************************/

static void _lift_vec_dp(double * a, ulong * b, slong len, ulong n)
{
    slong i;
    for (i = 0; i < len; i++)
        a[i] = (int)(b[i] - (n & FLINT_SIGN_EXT(n/2 - b[i])));
}

typedef struct {
    slong m;
    slong n;
    slong k;
    slong Astartrow;
    slong Astoprow;
    slong Bstartrow;
    slong Bstoprow;
    mp_limb_t ctxn;
    double * dA;
    double * dB;
    mp_limb_t ** Arows;
    mp_limb_t ** Brows;
} _lift_dp_worker_arg_struct;

void _lift_dp_worker(void * arg_ptr)
{
    _lift_dp_worker_arg_struct * arg = (_lift_dp_worker_arg_struct *) arg_ptr;
    slong n = arg->n;
    slong k = arg->k;
    slong Astartrow = arg->Astartrow;
    slong Astoprow = arg->Astoprow;
    slong Bstartrow = arg->Bstartrow;
    slong Bstoprow = arg->Bstoprow;
    mp_limb_t ctxn = arg->ctxn;
    double * dA = arg->dA;
    double * dB = arg->dB;
    mp_limb_t ** Arows = arg->Arows;
    mp_limb_t ** Brows = arg->Brows;
    slong i;

    for (i = Astartrow; i < Astoprow; i++)
        _lift_vec_dp(dA + i*k, Arows[i], k, ctxn);

    for (i = Bstartrow; i < Bstoprow; i++)
        _lift_vec_dp(dB + i*n, Brows[i], n, ctxn);
}

typedef struct {
    slong n;
    slong Cstartrow;
    slong Cstoprow;
    nmod_t * ctx;
    mp_limb_t shift;
    double * dC;
    mp_limb_t ** Crows;
} _reduce_dp_worker_arg_struct;

void _reduce_dp_worker(void * arg_ptr)
{
    _reduce_dp_worker_arg_struct * arg = (_reduce_dp_worker_arg_struct *) arg_ptr;
    slong n = arg->n;
    slong Cstartrow = arg->Cstartrow;
    slong Cstoprow = arg->Cstoprow;
    nmod_t ctx = *arg->ctx;
    mp_limb_t shift = arg->shift;
    double * dC = arg->dC;
    mp_limb_t ** Crows = arg->Crows;
    slong i, j;

    for (i = Cstartrow; i < Cstoprow; i++)
    {
        for (j = 0; j < n; j++)
        {
            slong a = (slong) dC[i*n + j];
            mp_limb_t b = (a < 0) ? a + shift : a;
            NMOD_RED(Crows[i][j], b, ctx);
        }
    }
}

int nmod_mat_mul_blas(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    slong i;
    slong m = A->r;
    slong k = A->c;
    slong n = B->c;
    double * dC, * dA, * dB;
    ulong hi, lo, shift;
    nmod_t ctx = C->mod;
    slong num_workers;
    thread_pool_handle * handles;
    void * tmp;

    if (m < 1 || k < 1 || n < 1 || m > INT_MAX || k > INT_MAX || n > INT_MAX)
        return 0;

    /*
        Each nmod has |_lift()| <= floor(mod.n/2) in the signed representation.
        Want k*floor(mod.n/2)^2 < MAX_BLAS_DP_INT.
    */

    umul_ppmm(hi, lo, ctx.n/2, ctx.n/2);
    if (hi != 0)
        return _nmod_mat_mul_blas_crt(C, A, B);

    umul_ppmm(hi, lo, lo, k);
    if (hi != 0 || lo >= MAX_BLAS_DP_INT)
        return _nmod_mat_mul_blas_crt(C, A, B);

    if (lo < MAX_BLAS_SP_INT)
        return _nmod_mat_mul_blas_sp(C, A, B);

    dA = flint_malloc(m*k*sizeof(double));
    dB = flint_malloc(k*n*sizeof(double));
    dC = flint_calloc(m*n, sizeof(double));

    num_workers = flint_request_threads(&handles, INT_MAX);

    tmp = flint_malloc((num_workers + 1)*FLINT_MAX(
                                        sizeof(_lift_dp_worker_arg_struct),
                                        sizeof(_reduce_dp_worker_arg_struct)));
    /* convert inputs */
    {
        _lift_dp_worker_arg_struct * args = (_lift_dp_worker_arg_struct *) tmp;

        for (i = 0; i <= num_workers; i++)
        {
            args[i].m = m;
            args[i].n = n;
            args[i].k = k;
            args[i].ctxn = ctx.n;
            args[i].dA = dA;
            args[i].dB = dB;
            args[i].Arows = A->rows;
            args[i].Brows = B->rows;
            _distribute_rows(&args[i].Astartrow, &args[i].Astoprow,
                             &args[i].Bstartrow, &args[i].Bstoprow, m,
                                          ((m + k)*(i + 0))/(num_workers + 1),
                                          ((m + k)*(i + 1))/(num_workers + 1));
        }

        for (i = 0; i < num_workers; i++)
            thread_pool_wake(global_thread_pool, handles[i], 0, _lift_dp_worker, &args[i]);
        _lift_dp_worker(&args[num_workers]);
        for (i = 0; i < num_workers; i++)
            thread_pool_wait(global_thread_pool, handles[i]);
    }

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
                                                1.0, dA, k, dB, n, 0.0, dC, n);

    /* convert output */

    /*
        the shift for negative outputs must satisfy
            ctx.n divides shift, and
            2^53 <= shift <= 2^64
    */
    shift = ((2*MAX_BLAS_DP_INT)/ctx.n)*ctx.n;
    FLINT_ASSERT(MAX_BLAS_DP_INT <= shift);

    {
        _reduce_dp_worker_arg_struct * args = (_reduce_dp_worker_arg_struct *) tmp;

        for (i = 0; i <= num_workers; i++)
        {
            args[i].n = n;
            args[i].Cstartrow = ((i + 0)*m)/(num_workers + 1);
            args[i].Cstoprow  = ((i + 1)*m)/(num_workers + 1);
            args[i].ctx = &ctx;
            args[i].shift = shift;
            args[i].dC = dC;
            args[i].Crows = C->rows;
        }

        for (i = 0; i < num_workers; i++)
            thread_pool_wake(global_thread_pool, handles[i], 0, _reduce_dp_worker, &args[i]);
        _reduce_dp_worker(&args[num_workers]);
        for (i = 0; i < num_workers; i++)
            thread_pool_wait(global_thread_pool, handles[i]);
    }

    flint_free(tmp);

    flint_give_back_threads(handles, num_workers);

    flint_free(dA);
    flint_free(dB);
    flint_free(dC);

    return 1;
}

#else

int nmod_mat_mul_blas(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    return 0;
}

#endif
