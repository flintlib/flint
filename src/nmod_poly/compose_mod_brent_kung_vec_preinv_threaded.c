/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Martin Lee
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_support.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod_mat.h"

typedef struct
{
    nmod_poly_struct * res;
    nmod_mat_struct * C;
    mp_srcptr h;
    mp_srcptr poly;
    mp_srcptr polyinv;
    nmod_t p;
    mp_ptr t;
    volatile slong * j;
    slong k;
    slong m;
    slong len;
    slong leninv;
    slong len2;
#if FLINT_USES_PTHREAD
    pthread_mutex_t * mutex;
#endif
} compose_vec_arg_t;

void
_nmod_poly_compose_mod_brent_kung_vec_preinv_worker(void * arg_ptr)
{
    compose_vec_arg_t arg = *((compose_vec_arg_t *) arg_ptr);
    slong i, j, k = arg.k, n = arg.len - 1;
    slong len = arg.len, leninv = arg.leninv;
    mp_ptr t = arg.t;
    mp_srcptr h = arg.h;
    mp_srcptr poly = arg.poly;
    mp_srcptr polyinv = arg.polyinv;
    nmod_poly_struct * res = arg.res;
    nmod_mat_struct * C = arg.C;
    nmod_t p = arg.p;

    while (1)
    {
#if FLINT_USES_PTHREAD
        pthread_mutex_lock(arg.mutex);
#endif
        j = *arg.j;
        *arg.j = j + 1;
#if FLINT_USES_PTHREAD
        pthread_mutex_unlock(arg.mutex);
#endif

        if (j >= arg.len2)
            return;

        _nmod_vec_set(res[j].coeffs, C->rows[(j + 1)*k - 1], n);

        if (n == 1) /* special case, constant polynomials */
        {
            for (i = 2; i <= k; i++)
            {
                t[0] = n_mulmod2_preinv(res[j].coeffs[0], h[0], p.n, p.ninv);
                res[j].coeffs[0] = n_addmod(t[0],
                                               C->rows[(j + 1)*k - i][0], p.n);
            }
        } else
        {
            for (i = 2; i <= k; i++)
            {
                _nmod_poly_mulmod_preinv(t, res[j].coeffs, n, h, n, poly,
                                                      len, polyinv, leninv, p);
                _nmod_poly_add(res[j].coeffs, t, n,
                                                 C->rows[(j + 1)*k - i], n, p);
            }
        }
    }
}

void
_nmod_poly_compose_mod_brent_kung_vec_preinv_threaded_pool(nmod_poly_struct * res,
                                             const nmod_poly_struct * polys,
                                             slong lenpolys, slong l,
                                             mp_srcptr g, slong glen,
                                             mp_srcptr poly, slong len,
                                             mp_srcptr polyinv, slong leninv,
                                             nmod_t mod,
                                             thread_pool_handle * threads,
                                             slong num_threads)
{
    nmod_mat_t A, B, C;
    slong i, j, n, m, k, len2 = l, len1, shared_j = 0;
    mp_ptr h;
    compose_vec_arg_t * args;
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif

    n = len - 1;

    m = n_sqrt(n*len2) + 1;

    h = _nmod_vec_init(n);

    k = len/m + 1;

    nmod_mat_init(A, m, n, mod.n);
    nmod_mat_init(B, k*len2, m, mod.n);
    nmod_mat_init(C, k*len2, n, mod.n);

    /* Set rows of B to the segments of polys */
    for (j = 0; j < len2; j++)
    {
        len1 = polys[j].length;

        for (i = 0; i < len1 / m; i++)
            _nmod_vec_set(B->rows[i + j * k], polys[j].coeffs + i * m, m);

        _nmod_vec_set(B->rows[i + j * k], polys[j].coeffs + i * m,
                      len1 % m);
    }

    /* Set rows of A to powers of g */
    _nmod_poly_powers_mod_preinv_threaded_pool(A->rows, g, glen,
                         m, poly, len, polyinv, leninv, mod, threads, num_threads);

    _nmod_mat_mul_classical_threaded_pool_op(C, NULL, B, A, 0,
                                                         threads, num_threads);

    /* Evaluate block composition using the Horner scheme */
    if (n == 1)
    {
        h[0] = n_mulmod2_preinv(A->rows[m - 1][0],
                                               A->rows[1][0], mod.n, mod.ninv);
    } else
    {

        _nmod_poly_mulmod_preinv(h, A->rows[m - 1], n, A->rows[1], n, poly,
                             len, polyinv, leninv, mod);
    }

    args = (compose_vec_arg_t *)
                   flint_malloc(sizeof(compose_vec_arg_t) * (num_threads + 1));

    for (i = 0; i < num_threads + 1; i++)
    {
        args[i].res     = res;
        args[i].C       = C;
        args[i].h       = h;
        args[i].k       = k;
        args[i].m       = m;
        args[i].j       = &shared_j;
        args[i].poly    = poly;
        args[i].t       = _nmod_vec_init(len);
        args[i].len     = len;
        args[i].polyinv = polyinv;
        args[i].leninv  = leninv;
        args[i].p       = mod;
        args[i].len2    = len2;
#if FLINT_USES_PTHREAD
        args[i].mutex   = &mutex;
#endif
    }

#if FLINT_USES_PTHREAD
    pthread_mutex_init(&mutex, NULL);
#endif

    for (i = 0; i < num_threads; i++)
    {
        thread_pool_wake(global_thread_pool, threads[i], 0,
                _nmod_poly_compose_mod_brent_kung_vec_preinv_worker, &args[i]);
    }

    _nmod_poly_compose_mod_brent_kung_vec_preinv_worker(&args[num_threads]);

    for (i = 0; i < num_threads; i++)
    {
        thread_pool_wait(global_thread_pool, threads[i]);
    }

#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&mutex);
#endif

    for (i = 0; i < num_threads + 1; i++)
       _nmod_vec_clear(args[i].t);

    flint_free(args);

    _nmod_vec_clear(h);

    nmod_mat_clear(A);
    nmod_mat_clear(B);
    nmod_mat_clear(C);
}

void
nmod_poly_compose_mod_brent_kung_vec_preinv_threaded_pool(nmod_poly_struct * res,
                                            const nmod_poly_struct * polys,
                                            slong len1, slong n,
                                            const nmod_poly_t g,
                                            const nmod_poly_t poly,
                                            const nmod_poly_t polyinv,
                                            thread_pool_handle * threads,
                                            slong num_threads)
{
    slong len2 = poly->length;
    slong i;

    if (n == 0)
        return;

    if (len2 == 1)
    {
        for (i = 0; i < n; i++)
            nmod_poly_zero(res + i);

        return;
    }

    if (len2 == 2)
    {
        for (i = 0; i < n; i++)
            nmod_poly_set(res + i, polys + i);

        return;
    }

    for (i = 0; i < n; i++)
    {
        nmod_poly_fit_length(res + i, len2 - 1);
        _nmod_poly_set_length(res + i, len2 - 1);
    }

    _nmod_poly_compose_mod_brent_kung_vec_preinv_threaded_pool(res, polys, len1, n,
                                                          g->coeffs, g->length,
                                                          poly->coeffs, len2,
                                                          polyinv->coeffs,
                                                          polyinv->length,
                                                          poly->mod,
                                                          threads,
                                                          num_threads);

    for (i = 0; i < n; i++)
        _nmod_poly_normalise(res + i);
}

void
nmod_poly_compose_mod_brent_kung_vec_preinv_threaded(nmod_poly_struct * res,
                                            const nmod_poly_struct * polys,
                                            slong len1, slong n,
                                            const nmod_poly_t g,
                                            const nmod_poly_t poly,
                                            const nmod_poly_t polyinv)
{
    slong i, len2 = poly->length;
    thread_pool_handle * threads;
    slong num_threads;

    for (i = 0; i < len1; i++)
    {
        if (polys[i].length >= len2)
        {
            flint_throw(FLINT_ERROR, "(nmod_poly_compose_mod_brent_kung_vec_preinv_threaded): "
                 "The degree of the first polynomial must be smaller than that of the modulus\n");
        }
    }

    if (n > len1)
    {
        flint_throw(FLINT_ERROR, "(nmod_poly_compose_mod_brent_kung_vec_preinv_threaded): "
                "n is larger than the length of polys\n");
    }

    if (n == 0)
        return;

    if (len2 == 1)
    {
        for (i = 0; i < n; i++)
            nmod_poly_zero(res + i);

        return;
    }

    if (len2 == 2)
    {
        for (i = 0; i < n; i++)
            nmod_poly_set(res + i, polys + i);

        return;
    }

    for (i = 0; i < n; i++)
    {
        nmod_poly_fit_length(res + i, len2 - 1);
        _nmod_poly_set_length(res + i, len2 - 1);
    }

    num_threads = flint_request_threads(&threads, flint_get_num_threads());

    _nmod_poly_compose_mod_brent_kung_vec_preinv_threaded_pool(res, polys, len1, n,
                                                          g->coeffs, g->length,
                                                          poly->coeffs, len2,
                                                          polyinv->coeffs,
                                                          polyinv->length,
                                                          poly->mod,
                                                          threads,
                                                          num_threads);

    flint_give_back_threads(threads, num_threads);

    for (i = 0; i < n; i++)
        _nmod_poly_normalise(res + i);
}
