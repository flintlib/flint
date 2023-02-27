/*
    Copyright (C) 2012-2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <pthread.h>
#include "acb_poly.h"

typedef struct
{
    acb_ptr z;
    acb_srcptr s;
    acb_srcptr a;
    acb_srcptr q;
    slong n0;
    slong n1;
    slong d0;
    slong len;
    slong prec;
}
powsum_arg_t;

void *
_acb_zeta_powsum_evaluator(void * arg_ptr)
{
    powsum_arg_t arg = *((powsum_arg_t *) arg_ptr);
    slong i, k;
    int q_one, s_int;

    acb_t t, u, v, ak, qpow, negs;
    arb_t f;

    acb_init(t);
    acb_init(u);
    acb_init(v);
    acb_init(ak);
    acb_init(qpow);
    acb_init(negs);
    arb_init(f);

    _acb_vec_zero(arg.z, arg.len);

    q_one = acb_is_one(arg.q);
    s_int = arb_is_int(acb_realref(arg.s)) && arb_is_zero(acb_imagref(arg.s));

    if (!q_one)
        acb_pow_ui(qpow, arg.q, arg.n0, arg.prec);

    acb_neg(negs, arg.s);
    arb_fac_ui(f, arg.d0, arg.prec);

    for (k = arg.n0; k < arg.n1; k++)
    {
        acb_add_ui(ak, arg.a, k, arg.prec);

        if (arg.d0 == 0 && arg.len == 1)
        {
            /* u = (a+k)^(-s) */
            acb_pow(u, ak, negs, arg.prec);
        }
        else
        {
            /* t = log(a+k) */
            acb_log(t, ak, arg.prec);

            /* u = (a+k)^(-s) */
            if (s_int)
            {
                acb_pow(u, ak, negs, arg.prec);
            }
            else
            {
                acb_mul(u, t, negs, arg.prec);
                acb_exp(u, u, arg.prec);
            }
        }

        /* u = u * q^k */
        if (!q_one)
        {
            acb_mul(u, u, qpow, arg.prec);
            if (k < arg.n1 - 1)
                acb_mul(qpow, qpow, arg.q, arg.prec);
        }

        /* forward: u *= (-1)^d * log(a+k)^d / d! */
        if (arg.d0 != 0)
        {
            acb_pow_ui(v, t, arg.d0, arg.prec);
            acb_mul(u, u, v, arg.prec);
            arb_div(acb_realref(u), acb_realref(u), f, arg.prec);
            arb_div(acb_imagref(u), acb_imagref(u), f, arg.prec);
            if (arg.d0 % 2)
                acb_neg(u, u);
        }

        acb_add(arg.z, arg.z, u, arg.prec);

        for (i = 1; i < arg.len; i++)
        {
            acb_mul(u, u, t, arg.prec);
            acb_div_si(u, u, -(arg.d0 + i), arg.prec);
            acb_add(arg.z + i, arg.z + i, u, arg.prec);
        }
    }

    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
    acb_clear(ak);
    acb_clear(qpow);
    acb_clear(negs);
    arb_clear(f);

    flint_cleanup();

    return NULL;
}

void
_acb_poly_powsum_series_naive_threaded(acb_ptr z,
    const acb_t s, const acb_t a, const acb_t q, slong n, slong len, slong prec)
{
    pthread_t * threads;
    powsum_arg_t * args;
    slong i, num_threads;
    int split_each_term;

    num_threads = flint_get_num_threads();

    threads = flint_malloc(sizeof(pthread_t) * num_threads);
    args = flint_malloc(sizeof(powsum_arg_t) * num_threads);

    split_each_term = (len > 1000);

    for (i = 0; i < num_threads; i++)
    {
        args[i].s = s;
        args[i].a = a;
        args[i].q = q;

        if (split_each_term)
        {
            slong n0, n1;
            n0 = (len * i) / num_threads;
            n1 = (len * (i + 1)) / num_threads;
            args[i].z = z + n0;
            args[i].n0 = 0;
            args[i].n1 = n;
            args[i].d0 = n0;
            args[i].len = n1 - n0;
        }
        else
        {
            args[i].z = _acb_vec_init(len);
            args[i].n0 = (n * i) / num_threads;
            args[i].n1 = (n * (i + 1)) / num_threads;
            args[i].d0 = 0;
            args[i].len = len;
        }

        args[i].prec = prec;
        pthread_create(&threads[i], NULL, _acb_zeta_powsum_evaluator, &args[i]);
    }

    for (i = 0; i < num_threads; i++)
    {
        pthread_join(threads[i], NULL);
    }

    if (!split_each_term)
    {
        _acb_vec_zero(z, len);
        for (i = 0; i < num_threads; i++)
        {
            _acb_vec_add(z, z, args[i].z, len, prec);
            _acb_vec_clear(args[i].z, len);
        }
    }

    flint_free(threads);
    flint_free(args);
}

