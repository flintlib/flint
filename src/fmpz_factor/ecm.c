/*
    Copyright (C) 2015 Kushagra Singh

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Outer wrapper for ECM
   makes calls to stage I and stage II (one) */

#include "ulong_extras.h"
#include "mpn_extras.h"
#include "fmpz.h"
#include "fmpz_factor.h"
#include "thread_pool.h"
#include "thread_support.h"
#if FLINT_USES_PTHREAD
#include <pthread.h>
#endif

static
ulong n_ecm_primorial[] =
{
#ifdef FLINT64

    UWORD(2), UWORD(6), UWORD(30), UWORD(210), UWORD(2310), UWORD(30030),
    UWORD(510510), UWORD(9699690), UWORD(223092870), UWORD(6469693230),
    UWORD(200560490130), UWORD(7420738134810), UWORD(304250263527210),
    UWORD(13082761331670030), UWORD(614889782588491410)
    /* 15 values */

#else

    UWORD(2), UWORD(6), UWORD(30), UWORD(210), UWORD(2310), UWORD(30030),
    UWORD(510510), UWORD(9699690)
    /* 9 values */

#endif
};

#ifdef FLINT64
#define num_n_ecm_primorials 15
#else
#define num_n_ecm_primorials 9
#endif

/*
   One worker tries a strided subset of the curves: curve indices
   start, start + stride, start + 2*stride, ....  All curves are
   independent, so this parallelises with near-linear speedup until a
   factor is found.  Shared state (n, prime_array, GCD_table, prime_table,
   the precomputed sigmas) is read-only; each worker has its own scratch
   (its own ecm_s buffers and gcd output buffer `fac`).  The first worker
   to find a factor sets `found`; the others stop after their current curve.
*/
typedef struct
{
    /* shared, read-only */
    nn_srcptr n;
    const ulong * prime_array;
    nn_srcptr mpsig_all;        /* curves * n_size precomputed sigmas */
    ulong num, B1, B2, P;
    ulong curves, n_size;
    ulong start, stride;

    /* private scratch */
    ecm_s * ecm;                /* own buffers; tables/ninv/one aliased/copied */
    nn_ptr fac;                 /* own gcd output buffer (n_size limbs) */

    /* shared result + synchronisation */
#if FLINT_USES_PTHREAD
    pthread_mutex_t * mutex;
#endif
    volatile int * found;
    volatile ulong * found_curve;
    nn_ptr result;              /* winning gcd (n_size limbs, normalised) */
    slong * result_size;
    int * result_ret;
}
_ecm_worker_arg;

static void
_ecm_worker(void * varg)
{
    _ecm_worker_arg * a = (_ecm_worker_arg *) varg;
    ulong j;

    for (j = a->start; j < a->curves; j += a->stride)
    {
        int gcdlimbs, stage_code = 0;
        nn_srcptr sig;

        if (*a->found)              /* another curve already succeeded */
            break;

        sig = a->mpsig_all + j * a->n_size;

        gcdlimbs = fmpz_factor_ecm_select_curve(a->fac, (nn_ptr) sig,
                                                (nn_ptr) a->n, a->ecm);

        if (gcdlimbs == -1)         /* degenerate curve: skip (cf. serial guard) */
            continue;
        else if (gcdlimbs > 0)
            stage_code = -1;        /* factor found while selecting curve */
        else
        {
            gcdlimbs = fmpz_factor_ecm_stage_I(a->fac, a->prime_array,
                                       a->num, a->B1, (nn_ptr) a->n, a->ecm);
            if (gcdlimbs > 0)
                stage_code = 1;     /* factor found in stage I */
            else
            {
                gcdlimbs = fmpz_factor_ecm_stage_II(a->fac, a->B1, a->B2,
                                       a->P, (nn_ptr) a->n, a->ecm);
                if (gcdlimbs > 0)
                    stage_code = 2; /* factor found in stage II */
            }
        }

        if (stage_code != 0)
        {
#if FLINT_USES_PTHREAD
            if (a->mutex != NULL)
                pthread_mutex_lock(a->mutex);
#endif
            /* keep the lowest-index winning curve for reproducibility */
            if (!*a->found || j < *a->found_curve)
            {
                flint_mpn_copyi(a->result, a->fac, gcdlimbs);
                *a->result_size = gcdlimbs;
                *a->result_ret = stage_code;
                *a->found_curve = j;
                *a->found = 1;
            }
#if FLINT_USES_PTHREAD
            if (a->mutex != NULL)
                pthread_mutex_unlock(a->mutex);
#endif
            break;
        }
    }
}

int
fmpz_factor_ecm(fmpz_t f, ulong curves, ulong B1, ulong B2,
                flint_rand_t state, const fmpz_t n_in)
{
    fmpz_t sig, nm8;
    ulong P, num, maxP, mmin, mmax, mdiff, prod, maxj, n_size, cy;
    ulong i, j;
    int ret;
    ecm_t ecm_inf;
    mpz_ptr fac, mptr;
    nn_ptr n;

    TMP_INIT;

    const ulong *prime_array;
    n_size = fmpz_size(n_in);

    if (n_size == 1)
    {
        ret = n_factor_ecm(&P, curves, B1, B2, state, fmpz_get_ui(n_in));
        fmpz_set_ui(f, P);
        return ret;
    }

    fmpz_factor_ecm_init(ecm_inf, n_size);

    TMP_START;

    n      = TMP_ALLOC(n_size * sizeof(ulong));

    if ((!COEFF_IS_MPZ(* n_in)))
    {
        ecm_inf->normbits = flint_clz(fmpz_get_ui(n_in));
        n[0] = fmpz_get_ui(n_in);
        n[0] <<= ecm_inf->normbits;
    }
    else
    {
        mptr = COEFF_TO_PTR(* n_in);
        ecm_inf->normbits = flint_clz(mptr->_mp_d[n_size - 1]);
        if (ecm_inf->normbits)
           mpn_lshift(n, mptr->_mp_d, n_size, ecm_inf->normbits);
        else
           flint_mpn_copyi(n, mptr->_mp_d, n_size);
    }

    flint_mpn_preinvn(ecm_inf->ninv, n, n_size);
    ecm_inf->one[0] = UWORD(1) << ecm_inf->normbits;

    fmpz_init(sig);
    fmpz_init(nm8);
    fmpz_sub_ui(nm8, n_in, 8);

    ret = 0;
    /* FIXME: Wait to promote f until after stage 2 precomputations? */
    fac = _fmpz_promote(f);
    {
        int alloc = fmpz_size(n_in);
        FLINT_MPZ_REALLOC(fac, alloc);
    }

    /************************ STAGE I PRECOMPUTATIONS ************************/

    num = n_prime_pi(B1);   /* number of primes under B1 */

    /* compute list of primes under B1 for stage I */
    prime_array = n_primes_arr_readonly(num);

    /************************ STAGE II PRECOMPUTATIONS ***********************/

    maxP = n_sqrt(B2);

    /* Selecting primorial */

    j = 1;
    while ((j < num_n_ecm_primorials) && (n_ecm_primorial[j] < maxP))
        j += 1;

    P = n_ecm_primorial[j - 1];

    mmin = (B1 + (P/2)) / P;
    mmax = ((B2 - P/2) + P - 1)/P;      /* ceil */
    if (mmax < mmin)
    {
       flint_throw(FLINT_ERROR, "Exception (ecm). B1 > B2 encountered.\n");
    }
    maxj = (P + 1)/2;
    mdiff = mmax - mmin + 1;

    /* compute GCD_table */

    ecm_inf->GCD_table = flint_malloc(maxj + 1);

    for (j = 1; j <= maxj; j += 2)
    {
        if ((j%2) && n_gcd(j, P) == 1)
            ecm_inf->GCD_table[j] = 1;
        else
            ecm_inf->GCD_table[j] = 0;
    }

    /* compute prime table */

    ecm_inf->prime_table = flint_malloc(mdiff * sizeof(unsigned char*));

    for (i = 0; i < mdiff; i++)
        ecm_inf->prime_table[i] = flint_malloc((maxj + 1) * sizeof(unsigned char));

    for (i = 0; i < mdiff; i++)
    {
        for (j = 1; j <= maxj; j += 2)
        {
            ecm_inf->prime_table[i][j] = 0;

            /* if (i + mmin)*P + j
               is prime, mark 1. Can be possibly prime
               only if gcd(j, P) = 1 */

            if (ecm_inf->GCD_table[j] == 1)
            {
                prod = (i + mmin)*P + j;
                if (n_is_prime(prod))
                    ecm_inf->prime_table[i][j] = 1;

                prod = (i + mmin)*P - j;
                if (n_is_prime(prod))
                    ecm_inf->prime_table[i][j] = 1;
            }
        }
    }


    /****************************** TRY "CURVES" *****************************/

    /* All curves are independent, so distribute them across the thread pool. */
    {
        thread_pool_handle * handles;
        slong num_handles, nworkers, w;
        ulong jj;
        nn_ptr mpsig_all, result;
        ecm_s * worker_ecm;
        nn_ptr * worker_fac;
        _ecm_worker_arg * args;
        volatile int found = 0;
        volatile ulong found_curve = 0;
        slong result_size = 0;
        int result_ret = 0;
#if FLINT_USES_PTHREAD
        pthread_mutex_t mutex;
#endif

        mpsig_all = flint_calloc(FLINT_MAX(curves, 1) * (ulong) ecm_inf->n_size,
                                 sizeof(ulong));

        for (jj = 0; jj < curves; jj++)
        {
            nn_ptr mpsig = mpsig_all + jj * (ulong) ecm_inf->n_size;

            fmpz_randm(sig, state, nm8);
            fmpz_add_ui(sig, sig, 7);

            if ((!COEFF_IS_MPZ(*sig)))
            {
                mpsig[0] = fmpz_get_ui(sig);
                if (ecm_inf->normbits)
                {
                    cy = mpn_lshift(mpsig, mpsig, 1, ecm_inf->normbits);
                    if (cy)
                        mpsig[1] = cy;
                }
            }
            else
            {
                mptr = COEFF_TO_PTR(*sig);

                if (ecm_inf->normbits)
                {
                    cy = mpn_lshift(mpsig, mptr->_mp_d, mptr->_mp_size, ecm_inf->normbits);
                    if (cy)
                        mpsig[mptr->_mp_size] = cy;
                }
                else
                {
                    flint_mpn_copyi(mpsig, mptr->_mp_d, mptr->_mp_size);
                }
            }
        }

        /* request threads and decide how many workers to use */
        num_handles = flint_request_threads(&handles, flint_get_num_threads());
        nworkers = num_handles + 1;
        if ((ulong) nworkers > curves)
            nworkers = (slong) curves;
        if (nworkers < 1)
            nworkers = 1;

        result      = flint_calloc(ecm_inf->n_size, sizeof(ulong));
        worker_ecm  = flint_malloc(nworkers * sizeof(ecm_s));
        worker_fac  = flint_malloc(nworkers * sizeof(nn_ptr));
        args        = flint_malloc(nworkers * sizeof(_ecm_worker_arg));

#if FLINT_USES_PTHREAD
        pthread_mutex_init(&mutex, NULL);
#endif

        for (w = 0; w < nworkers; w++)
        {
            /* private scratch; ninv/one copied, big tables aliased read-only */
            fmpz_factor_ecm_init(&worker_ecm[w], ecm_inf->n_size);
            flint_mpn_copyi(worker_ecm[w].ninv, ecm_inf->ninv, ecm_inf->n_size);
            worker_ecm[w].one[0]      = ecm_inf->one[0];
            worker_ecm[w].normbits    = ecm_inf->normbits;
            worker_ecm[w].n_size      = ecm_inf->n_size;
            worker_ecm[w].GCD_table   = ecm_inf->GCD_table;
            worker_ecm[w].prime_table = ecm_inf->prime_table;

            worker_fac[w] = flint_calloc(ecm_inf->n_size, sizeof(ulong));

            args[w].n           = n;
            args[w].prime_array = prime_array;
            args[w].mpsig_all   = mpsig_all;
            args[w].num         = num;
            args[w].B1          = B1;
            args[w].B2          = B2;
            args[w].P           = P;
            args[w].curves      = curves;
            args[w].n_size      = ecm_inf->n_size;
            args[w].start       = (ulong) w;
            args[w].stride      = (ulong) nworkers;
            args[w].ecm         = &worker_ecm[w];
            args[w].fac         = worker_fac[w];
#if FLINT_USES_PTHREAD
            args[w].mutex       = (nworkers > 1) ? &mutex : NULL;
#endif
            args[w].found       = &found;
            args[w].found_curve = &found_curve;
            args[w].result      = result;
            args[w].result_size = &result_size;
            args[w].result_ret  = &result_ret;
        }

        for (w = 0; w < nworkers - 1; w++)
            thread_pool_wake(global_thread_pool, handles[w], 0,
                             _ecm_worker, &args[w]);

        _ecm_worker(&args[nworkers - 1]);

        for (w = 0; w < nworkers - 1; w++)
            thread_pool_wait(global_thread_pool, handles[w]);

        /* combine: write the winning factor into f, denormalising */
        if (found)
        {
            flint_mpn_copyi(fac->_mp_d, result, result_size);
            if (ecm_inf->normbits)
                mpn_rshift(fac->_mp_d, fac->_mp_d, result_size, ecm_inf->normbits);
            MPN_NORM(fac->_mp_d, result_size);
            fac->_mp_size = result_size;
            _fmpz_demote_val(f);
            ret = result_ret;
        }
        else
            ret = 0;

        for (w = 0; w < nworkers; w++)
        {
            fmpz_factor_ecm_clear(&worker_ecm[w]);
            flint_free(worker_fac[w]);
        }
        flint_free(worker_ecm);
        flint_free(worker_fac);
        flint_free(args);
        flint_free(result);
        flint_free(mpsig_all);
#if FLINT_USES_PTHREAD
        pthread_mutex_destroy(&mutex);
#endif
        flint_give_back_threads(handles, num_handles);
    }

    flint_free(ecm_inf->GCD_table);
    for (i = 0; i < mdiff; i++)
        flint_free(ecm_inf->prime_table[i]);
    flint_free(ecm_inf->prime_table);

    fmpz_factor_ecm_clear(ecm_inf);

    fmpz_clear(nm8);
    fmpz_clear(sig);

    TMP_END;

    return ret;
}
