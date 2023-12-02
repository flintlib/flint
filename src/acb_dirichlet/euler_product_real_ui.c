/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_support.h"
#include "acb_dirichlet.h"

#ifdef __GNUC__
# define log __builtin_log
# define sqrt __builtin_sqrt
#else
# include <math.h>
#endif

#define ONE_OVER_LOG2 1.4426950408889634

typedef struct
{
    ulong s;
    int mod;
    const signed char * chi;
    mp_ptr primes;
    double * powmags;
    slong num_primes;
    slong wp;
    slong index;
    slong num_threads;
    arb_struct res;
}
euler_work_chunk_t;

static void
euler_worker(void * _work)
{
    euler_work_chunk_t * work = ((euler_work_chunk_t *) _work);
    slong i;

    slong powprec;
    double powmag;
    arb_t t, u;
    ulong p;

    arb_init(t);
    arb_init(u);

    for (i = work->index; i < work->num_primes; i += work->num_threads)
    {
        p = work->primes[i];
        powmag = work->powmags[i];

        powprec = FLINT_MAX(work->wp - powmag, 8);

        arb_ui_pow_ui(t, p, work->s, powprec);
        arb_set_round(u, &work->res, powprec);
        arb_div(t, u, t, powprec);

        if (work->mod == 1 || (work->chi[p % work->mod] == 1))
            arb_sub(&work->res, &work->res, t, work->wp);
        else
            arb_add(&work->res, &work->res, t, work->wp);
    }

    arb_clear(t);
    arb_clear(u);
}

void _acb_dirichlet_euler_product_real_ui(arb_t res, ulong s,
    const signed char * chi, int mod, int reciprocal, slong prec)
{
    slong wp, powprec;
    double logp, powmag, errmag, limit;
    arb_t t, u;
    ulong p;
    mag_t err;
    slong num_threads;

    if (s <= 1)
    {
        arb_indeterminate(res);
        return;
    }

    if (prec < 2)
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    /* L(s), 1/L(s) = 1 + ...  For s >= 3, zeta(s,2) < 2^(1-s). */
    if (s > (ulong) prec)
    {
        arf_one(arb_midref(res));
        mag_set_ui_2exp_si(arb_radref(res), 1, -prec);
        return;
    }

    /* L(s), 1/L(s) = 1 +/- chi(2) 2^(-s) + ...
       For s >= 2, zeta(s,3) < 2^(2-floor(3s/2)). */
    if (s > 0.7 * prec)
    {
        arb_one(res);

        if (chi[2 % mod] != 0)
        {
            arf_t t;
            arf_init(t);
            arf_set_si_2exp_si(t, chi[2 % mod], -s);
            if (reciprocal)
                arf_neg(t, t);
            arb_add_arf(res, res, t, prec);
            arf_clear(t);
        }

        arb_add_error_2exp_si(res, 2 - (3 * s) / 2);
        return;
    }

    /* Heuristic. */
    wp = prec + FLINT_BIT_COUNT(prec) + (prec / s) + 4;

    arb_init(t);
    arb_init(u);

    /* 1 - chi(2) 2^(-s) */
    arb_one(res);
    arf_set_ui_2exp_si(arb_midref(t), 1, -s);   /* -s cannot overflow */
    if (chi[2 % mod] == 1)
        arb_sub(res, res, t, wp);
    else if (chi[2 % mod] == -1)
        arb_add(res, res, t, wp);

    /* Cut at some point if this algorithm just isn't a good fit... */
    /* The C * prec / log(prec) cutoff in arb_zeta_ui implies that the limit
       should be at least prec ^ (log(2) / C) for the Riemann zeta function,
       which gives prec ^ 1.2956 here. */
    limit = 100 + prec * sqrt(prec);

    num_threads = flint_get_num_available_threads();

    if (num_threads > 1 && prec > 5000 && s > 5000)
    {
        n_primes_t iter;
        slong i;
        mp_ptr primes;
        double * powmags;
        slong num_primes = 0;
        slong alloc = 16;
        slong thread_limit, num_threads, num_workers;
        thread_pool_handle * handles;
        euler_work_chunk_t * work;

        ulong first_omitted_p = 3;

        n_primes_init(iter);
        n_primes_jump_after(iter, 3);

        primes = flint_malloc(alloc * sizeof(mp_limb_t));
        powmags = flint_malloc(alloc * sizeof(double));

        for (p = 3; p < limit; p = n_primes_next(iter))
        {
            first_omitted_p = p;

            if (mod == 1 || chi[p % mod] != 0)
            {
                /* p^s */
                logp = log(p);
                powmag = s * logp * ONE_OVER_LOG2;

                /* zeta(s,p) ~= 1/p^s + 1/((s-1) p^(s-1)) */
                errmag = (log(s - 1.0) + (s - 1.0) * logp) * ONE_OVER_LOG2;
                errmag = FLINT_MIN(powmag, errmag);

                if (errmag > prec + 2)
                    break;

                if (num_primes >= alloc)
                {
                    alloc *= 2;
                    primes = flint_realloc(primes, alloc * sizeof(mp_limb_t));
                    powmags = flint_realloc(powmags, alloc * sizeof(double));
                }

                primes[num_primes] = p;
                powmags[num_primes] = powmag;

                num_primes++;
            }
        }

        n_primes_clear(iter);

        thread_limit = flint_get_num_threads();
        thread_limit = FLINT_MIN(thread_limit, num_primes / 4);
        thread_limit = FLINT_MAX(thread_limit, 1);

        num_workers = flint_request_threads(&handles, thread_limit);
        num_threads = num_workers + 1;

        work = flint_malloc(num_threads * sizeof(euler_work_chunk_t));

        for (i = 0; i < num_threads; i++)
        {
            work[i].s = s;
            work[i].mod = mod;
            work[i].chi = chi;
            work[i].primes = primes;
            work[i].powmags = powmags;
            work[i].num_primes = num_primes;
            work[i].wp = wp;
            work[i].index = i;
            work[i].num_threads = num_threads;
            arb_init(&work[i].res);
            arb_one(&work[i].res);
        }

        for (i = 0; i < num_workers; i++)
            thread_pool_wake(global_thread_pool, handles[i], 0, euler_worker, &work[i]);

        euler_worker(&work[num_workers]);

        for (i = 0; i < num_workers; i++)
            thread_pool_wait(global_thread_pool, handles[i]);

        flint_give_back_threads(handles, num_workers);

        for (i = 0; i < num_threads; i++)
        {
            arb_mul(res, res, &work[i].res, wp);
            arb_clear(&work[i].res);
        }

        flint_free(work);

        flint_free(primes);
        flint_free(powmags);

        mag_init(err);
        mag_hurwitz_zeta_uiui(err, s, first_omitted_p);
        arb_add_error_mag(res, err);
        mag_clear(err);
    }
    else
    {
        /* todo: prime iterator here too? */

        for (p = 3; p < limit; p = n_nextprime(p, 1))
        {
            if (mod == 1 || chi[p % mod] != 0)
            {
                /* p^s */
                logp = log(p);
                powmag = s * logp * ONE_OVER_LOG2;

                /* zeta(s,p) ~= 1/p^s + 1/((s-1) p^(s-1)) */
                errmag = (log(s - 1.0) + (s - 1.0) * logp) * ONE_OVER_LOG2;
                errmag = FLINT_MIN(powmag, errmag);

                if (errmag > prec + 2)
                    break;

                powprec = FLINT_MAX(wp - powmag, 8);

                arb_ui_pow_ui(t, p, s, powprec);
                arb_set_round(u, res, powprec);
                arb_div(t, u, t, powprec);

                if (mod == 1 || (chi[p % mod] == 1))
                    arb_sub(res, res, t, wp);
                else
                    arb_add(res, res, t, wp);
            }
        }

        mag_init(err);
        mag_hurwitz_zeta_uiui(err, s, p);
        arb_add_error_mag(res, err);
        mag_clear(err);
    }

    if (reciprocal)
        arb_set_round(res, res, prec);
    else
        arb_inv(res, res, prec);

    arb_clear(t);
    arb_clear(u);
}

