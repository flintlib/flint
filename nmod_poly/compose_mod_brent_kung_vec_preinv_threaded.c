/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Martin Lee

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <pthread.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod_mat.h"
#include "ulong_extras.h"

typedef struct
{
    nmod_poly_struct res;
    nmod_poly_struct g;
    nmod_mat_struct C;
    mp_srcptr h;
    mp_srcptr poly;
    mp_srcptr polyinv;
    nmod_t p;
    slong j;
    slong k;
    slong m;
    slong len;
    slong leninv;
}
compose_vec_arg_t;

void *
_nmod_poly_compose_mod_brent_kung_vec_preinv_worker(void * arg_ptr)
{
    compose_vec_arg_t arg= *((compose_vec_arg_t *) arg_ptr);
    slong i, n;
    mp_ptr t;

    n = arg.len - 1;
    t = _nmod_vec_init(n);

    _nmod_vec_set(arg.res.coeffs, arg.C.rows[(arg.j + 1) * arg.k - 1], n);
    for (i = 2; i <= arg.k; i++)
    {
        _nmod_poly_mulmod_preinv(t, arg.res.coeffs, n, arg.h, n, arg.poly,
                                     arg.len, arg.polyinv, arg.leninv, arg.p);
        _nmod_poly_add(arg.res.coeffs, t, n,
                           arg.C.rows[(arg.j + 1) * arg.k - i], n, arg.p);
    }

    _nmod_vec_clear(t);

    flint_cleanup();
    return NULL;
}

void
_nmod_poly_compose_mod_brent_kung_vec_preinv_threaded(nmod_poly_struct * res,
                                             const nmod_poly_struct * polys,
                                             slong lenpolys, slong l,
                                             mp_srcptr poly, slong len,
                                             mp_srcptr polyinv, slong leninv,
                                             nmod_t mod)
{
    nmod_mat_t A, B, C;
    slong i, j, n, m, k, len2 = l, len1, num_threads, c;
    mp_ptr h;
    pthread_t *threads;
    compose_vec_arg_t * args;

    n = len - 1;

    m = n_sqrt(n * len2) + 1;

    h = _nmod_vec_init(n);

    k = len / m + 1;

    nmod_mat_init(A, m, n, mod.n);
    nmod_mat_init(B, k * len2, m, mod.n);
    nmod_mat_init(C, k * len2, n, mod.n);

    /* Set rows of B to the segments of polys */
    for (j = 0; j < len2; j++)
    {
        len1 = (polys + j)->length;
        for (i = 0; i < len1 / m; i++)
            _nmod_vec_set(B->rows[i + j * k], (polys + j)->coeffs + i * m, m);
        _nmod_vec_set(B->rows[i + j * k], (polys + j)->coeffs + i * m,
                      len1 % m);
    }

    /* Set rows of A to powers of last element of polys */
    A->rows[0][0] = UWORD(1);
    _nmod_vec_set(A->rows[1], (polys + lenpolys - 1)->coeffs,
                  (polys + lenpolys - 1)->length);
    flint_mpn_zero(A->rows[1] + (polys + lenpolys - 1)->length,
                   n - (polys + lenpolys - 1)->length);
    for (i = 2; i < m; i++)
        _nmod_poly_mulmod_preinv(A->rows[i], A->rows[i - 1], n, A->rows[1],
                                     n, poly, len, polyinv, leninv, mod);

    nmod_mat_mul(C, B, A);

    /* Evaluate block composition using the Horner scheme */
    _nmod_poly_mulmod_preinv(h, A->rows[m - 1], n, A->rows[1], n, poly,
                             len, polyinv, leninv, mod);

    num_threads = flint_get_num_threads();

    threads = flint_malloc(sizeof(pthread_t) * num_threads);
    args = flint_malloc(sizeof(compose_vec_arg_t) * num_threads);

    for (j = 0; j < len2 / num_threads + 1; j++)
    {
        c = 0;
        for (i = 0; i < num_threads; i++)
        {
            if (i + j * num_threads < len2)
            {
                c++;
                args[i].res     = res[i + j * num_threads];
                args[i].C       = *C;
                args[i].g       = polys[i + j * num_threads];
                args[i].h       = h;
                args[i].k       = k;
                args[i].m       = m;
                args[i].j       = i + j * num_threads;
                args[i].poly    = poly;
                args[i].len     = len;
                args[i].polyinv = polyinv;
                args[i].leninv  = leninv;
                args[i].p       = mod;

                pthread_create(&threads[i], NULL,
                        _nmod_poly_compose_mod_brent_kung_vec_preinv_worker,
                        &args[i]);
            }
        }
        for (i = 0; i < c; i++)
            pthread_join(threads[i], NULL);
    }

    flint_free(threads);
    flint_free(args);

    _nmod_vec_clear(h);

    nmod_mat_clear(A);
    nmod_mat_clear(B);
    nmod_mat_clear(C);
}

void
nmod_poly_compose_mod_brent_kung_vec_preinv_threaded(nmod_poly_struct * res,
                                            const nmod_poly_struct * polys,
                                            slong len1, slong n,
                                            const nmod_poly_t poly,
                                            const nmod_poly_t polyinv)
{
    slong len2 = poly->length;
    slong len3, i;

    for (i = 0; i < len1; i++)
    {
        len3 = (polys + i)->length;
        if (len3 >= len2)
        {
            flint_printf
                ("Exception (nmod_poly_compose_mod_brent_kung_vec_preinv_threaded)."
                 "The degree of the first polynomial must be smaller than that of the "
                 " modulus\n");
            flint_abort();
        }
    }

    if (n > len1)
    {
        flint_printf
            ("Exception (nmod_poly_compose_mod_brent_kung_vec_preinv_threaded)."
             "n is larger than the length of polys\n");
        flint_abort();
    }

    if (n == 0)
        return;

    if (len2 == 1)
    {
        for (i = 0; i < n; i++)
        {
            nmod_poly_init_preinv(res + i, poly->mod.n, poly->mod.ninv);
            nmod_poly_zero(res + i);
        }
        return;
    }

    if (len2 == 2)
    {
        for (i = 0; i < n; i++)
        {
            nmod_poly_init_preinv(res + i, poly->mod.n, poly->mod.ninv);
            nmod_poly_set(res + i, polys + i);
        }
        return;
    }

    for (i = 0; i < n; i++)
    {
        nmod_poly_init2_preinv(res + i, poly->mod.n, poly->mod.ninv, len2 - 1);
        _nmod_poly_set_length(res + i, len2 - 1);
    }

    _nmod_poly_compose_mod_brent_kung_vec_preinv_threaded(res, polys, len1, n,
                                                          poly->coeffs, len2,
                                                          polyinv->coeffs,
                                                          polyinv->length,
                                                          poly->mod);

    for (i = 0; i < n; i++)
        _nmod_poly_normalise(res + i);
}

