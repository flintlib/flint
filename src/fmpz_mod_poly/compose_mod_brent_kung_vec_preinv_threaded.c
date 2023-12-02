/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
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
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod.h"
#include "fmpz_mod_mat.h"
#include "fmpz_mod_poly.h"

typedef struct
{
    fmpz_mod_poly_struct * res;
    fmpz_mod_mat_struct * C;
    const fmpz * h;
    const fmpz * poly;
    const fmpz * polyinv;
    const fmpz_mod_ctx_struct * ctx;
    fmpz * t;
    volatile slong * j;
    slong k;
    slong m;
    slong len;
    slong leninv;
    slong len2;
#if FLINT_USES_PTHREAD
    pthread_mutex_t * mutex;
#endif
}
compose_vec_arg_t;

void
_fmpz_mod_poly_compose_mod_brent_kung_vec_preinv_worker(void * arg_ptr)
{
    compose_vec_arg_t arg = *((compose_vec_arg_t *) arg_ptr);
    slong i, j, k = arg.k, n = arg.len - 1;
    slong len = arg.len, leninv = arg.leninv;
    fmpz * t = arg.t;
    const fmpz * h = arg.h;
    const fmpz * poly = arg.poly;
    const fmpz * polyinv = arg.polyinv;
    fmpz_mod_poly_struct * res = arg.res;
    fmpz_mat_struct * C = arg.C->mat;
    const fmpz_mod_ctx_struct * ctx = arg.ctx;

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

        _fmpz_vec_set(res[j].coeffs, C->rows[(j + 1)*k - 1], n);

        if (n == 1) /* special case, constant polynomials */
        {
            for (i = 2; i <= k; i++)
            {
                fmpz_mod_mul(t + 0, res[j].coeffs + 0, h + 0, ctx);
                fmpz_mod_add(res[j].coeffs + 0, t + 0,
                                                 C->rows[(j + 1)*k - i] + 0, ctx);
            }
        }
        else
        {
            for (i = 2; i <= k; i++)
            {
                _fmpz_mod_poly_mulmod_preinv(t, res[j].coeffs, n, h, n, poly,
                                                      len, polyinv, leninv, ctx);
                _fmpz_mod_poly_add(res[j].coeffs, t, n,
                                                 C->rows[(j + 1)*k - i], n, ctx);
            }
        }
    }
}

void
_fmpz_mod_poly_compose_mod_brent_kung_vec_preinv_threaded_pool(fmpz_mod_poly_struct * res,
                                                 const fmpz_mod_poly_struct * polys,
                                                 slong lenpolys, slong l,
                                                 const fmpz * g, slong glen,
                                                 const fmpz * poly, slong len,
                                                 const fmpz * polyinv, slong leninv,
                                                 const fmpz_mod_ctx_t ctx,
                                                 thread_pool_handle * threads,
                                                 slong num_threads)
{
    fmpz_mod_mat_t A, B, C;
    slong i, j, n, m, k, len2 = l, len1, shared_j = 0;
    fmpz * h;
    compose_vec_arg_t * args;
    const fmpz * p = fmpz_mod_ctx_modulus(ctx);
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif

    n = len - 1;

    m = n_sqrt(n*len2) + 1;

    h = _fmpz_vec_init(n);

    k = len/m + 1;

    fmpz_mod_mat_init(A, m, n, p);
    fmpz_mod_mat_init(B, k*len2, m, p);
    fmpz_mod_mat_init(C, k*len2, n, p);

    /* Set rows of B to the segments of polys */
    for (j = 0; j < len2; j++)
    {
        len1 = polys[j].length;

        for (i = 0; i < len1 / m; i++)
            _fmpz_vec_set(B->mat->rows[i + j*k], polys[j].coeffs + i*m, m);

        _fmpz_vec_set(B->mat->rows[i + j*k], polys[j].coeffs + i*m,
                      len1 % m);
    }

    /* Set rows of A to powers of last element of polys */
    _fmpz_mod_poly_powers_mod_preinv_threaded_pool(A->mat->rows, g, glen,
	               m, poly, len, polyinv, leninv, ctx, threads, num_threads);

    _fmpz_mod_mat_mul_classical_threaded_pool_op(C, NULL, B, A, 0,
                                                         threads, num_threads);

    /* Evaluate block composition using the Horner scheme */
    if (n == 1)
    {
        fmpz_mod_mul(h + 0, A->mat->rows[m - 1] + 0, A->mat->rows[1] + 0, ctx);
    }
    else
    {
        _fmpz_mod_poly_mulmod_preinv(h, A->mat->rows[m - 1], n, A->mat->rows[1],
                        n, poly, len, polyinv, leninv, ctx);
    }

    args = (compose_vec_arg_t *)
                     flint_malloc(sizeof(compose_vec_arg_t)*(num_threads + 1));

    for (i = 0; i < num_threads + 1; i++)
    {
        args[i].res     = res;
        args[i].C       = C;
        args[i].h       = h;
        args[i].k       = k;
        args[i].m       = m;
        args[i].j       = &shared_j;
        args[i].poly    = poly;
        args[i].t       = _fmpz_vec_init(len);
        args[i].len     = len;
        args[i].polyinv = polyinv;
        args[i].leninv  = leninv;
        args[i].ctx     = ctx;
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
           _fmpz_mod_poly_compose_mod_brent_kung_vec_preinv_worker, &args[i]);
    }

    _fmpz_mod_poly_compose_mod_brent_kung_vec_preinv_worker(&args[num_threads]);

    for (i = 0; i < num_threads; i++)
    {
        thread_pool_wait(global_thread_pool, threads[i]);
    }

#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&mutex);
#endif

    for (i = 0; i < num_threads + 1; i++)
       _fmpz_vec_clear(args[i].t, len);

    flint_free(args);

    _fmpz_vec_clear(h, n);

    fmpz_mod_mat_clear(A);
    fmpz_mod_mat_clear(B);
    fmpz_mod_mat_clear(C);
}

void
fmpz_mod_poly_compose_mod_brent_kung_vec_preinv_threaded_pool(fmpz_mod_poly_struct * res,
                                                const fmpz_mod_poly_struct * polys,
                                                slong len1, slong n,
                                                const fmpz_mod_poly_t g,
                                                const fmpz_mod_poly_t poly,
                                                const fmpz_mod_poly_t polyinv,
                                               const fmpz_mod_ctx_t ctx,
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
            fmpz_mod_poly_zero(res + i, ctx);
    }

    if (len2 == 2)
    {
        for (i = 0; i < n; i++)
            fmpz_mod_poly_set(res + i, polys + i, ctx);

        return;
    }

    for (i = 0; i < n; i++)
    {
        fmpz_mod_poly_fit_length(res + i, len2 - 1, ctx);
        _fmpz_mod_poly_set_length(res + i, len2 - 1);
    }

    _fmpz_mod_poly_compose_mod_brent_kung_vec_preinv_threaded_pool(res, polys,
                                                     len1, n,
                                                     g->coeffs, g->length,
                                                     poly->coeffs, len2,
                                                     polyinv->coeffs,
                                                     polyinv->length,
                                                     ctx,
                                                     threads,
                                                     num_threads);

    for (i = 0; i < n; i++)
        _fmpz_mod_poly_normalise(res + i);
}

void
fmpz_mod_poly_compose_mod_brent_kung_vec_preinv_threaded(fmpz_mod_poly_struct * res,
                                            const fmpz_mod_poly_struct * polys,
                                            slong len1, slong n,
                                            const fmpz_mod_poly_t g,
                                            const fmpz_mod_poly_t poly,
                                            const fmpz_mod_poly_t polyinv,
                                            const fmpz_mod_ctx_t ctx)
{
    slong i, len2 = poly->length;
    thread_pool_handle * threads;
    slong num_threads;

    for (i = 0; i < len1; i++)
    {
        if (polys[i].length >= len2)
        {
            flint_throw(FLINT_ERROR, "(fmpz_mod_poly_compose_mod_brent_kung_vec_preinv_threaded): "
                 "The degree of the first polynomial must be smaller than that of the modulus\n");
        }
    }

    if (n > len1)
    {
        flint_throw(FLINT_ERROR, "(fmpz_mod_poly_compose_mod_brent_kung_vec_preinv_threaded): "
             "n is larger than the length of polys\n");
    }

    if (n == 0)
        return;

    if (len2 == 1)
    {
        for (i = 0; i < n; i++)
            fmpz_mod_poly_zero(res + i, ctx);

        return;
    }

    if (len2 == 2)
    {
        for (i = 0; i < n; i++)
            fmpz_mod_poly_set(res + i, polys + i, ctx);

        return;
    }

    for (i = 0; i < n; i++)
    {
        fmpz_mod_poly_fit_length(res + i, len2 - 1, ctx);
        _fmpz_mod_poly_set_length(res + i, len2 - 1);
    }

    num_threads = flint_request_threads(&threads, flint_get_num_threads());

    _fmpz_mod_poly_compose_mod_brent_kung_vec_preinv_threaded_pool(res, polys,
                                                          len1, n,
                                                          g->coeffs, g->length,
                                                          poly->coeffs, len2,
                                                          polyinv->coeffs,
                                                          polyinv->length,
                                                          ctx,
                                                          threads,
                                                          num_threads);

    flint_give_back_threads(threads, num_threads);

    for (i = 0; i < n; i++)
        _fmpz_mod_poly_normalise(res + i);
}
