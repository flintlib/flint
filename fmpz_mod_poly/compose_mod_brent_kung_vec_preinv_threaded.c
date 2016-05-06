/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
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
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"
#include "fmpz_mat.h"
#include "ulong_extras.h"

typedef struct
{
    fmpz_mod_poly_struct res;
    fmpz_mod_poly_struct g;
    fmpz_mat_struct C;
    fmpz * h;
    fmpz * poly;
    fmpz * polyinv;
    fmpz p;
    slong j;
    slong k;
    slong m;
    slong len;
    slong leninv;
}
compose_vec_arg_t;

void *
_fmpz_mod_poly_compose_mod_brent_kung_vec_preinv_worker(void * arg_ptr)
{
    compose_vec_arg_t arg= *((compose_vec_arg_t *) arg_ptr);
    slong i, n;
    fmpz * t;

    n = arg.len - 1;
    t = _fmpz_vec_init(n);

    _fmpz_vec_set(arg.res.coeffs, arg.C.rows[(arg.j + 1) * arg.k - 1], n);
    for (i = 2; i <= arg.k; i++)
    {
        _fmpz_mod_poly_mulmod_preinv(t, arg.res.coeffs, n, arg.h, n, arg.poly,
                                     arg.len, arg.polyinv, arg.leninv, &arg.p);
        _fmpz_mod_poly_add(arg.res.coeffs, t, n,
                           arg.C.rows[(arg.j + 1) * arg.k - i], n, &arg.p);
    }

    _fmpz_vec_clear(t, n);

    flint_cleanup();
    return NULL;
}

void
_fmpz_mod_poly_compose_mod_brent_kung_vec_preinv_threaded(fmpz_mod_poly_struct * res,
                                                 const fmpz_mod_poly_struct *
                                                 polys, slong lenpolys,
                                                 slong l, const fmpz * poly,
                                                 slong len,
                                                 const fmpz * polyinv,
                                                 slong leninv, const fmpz_t p)
{
    fmpz_mat_t A, B, C;
    slong i, j, n, m, k, len2 = l, len1, num_threads, c;
    fmpz *h;
    pthread_t *threads;
    compose_vec_arg_t * args;

    n = len - 1;

    m = n_sqrt(n * len2) + 1;

    h = _fmpz_vec_init(n);

    k = len / m + 1;

    fmpz_mat_init(A, m, n);
    fmpz_mat_init(B, k * len2, m);
    fmpz_mat_init(C, k * len2, n);

    /* Set rows of B to the segments of polys */
    for (j = 0; j < len2; j++)
    {
        len1 = (polys + j)->length;
        for (i = 0; i < len1 / m; i++)
            _fmpz_vec_set(B->rows[i + j * k], (polys + j)->coeffs + i * m, m);
        _fmpz_vec_set(B->rows[i + j * k], (polys + j)->coeffs + i * m,
                      len1 % m);
    }

    /* Set rows of A to powers of last element of polys */
    fmpz_one(A->rows[0]);
    _fmpz_vec_set(A->rows[1], (polys + lenpolys - 1)->coeffs,
                  (polys + lenpolys - 1)->length);
    _fmpz_vec_zero(A->rows[1] + (polys + lenpolys - 1)->length,
                   n - (polys + lenpolys - 1)->length);
    for (i = 2; i < m; i++)
        _fmpz_mod_poly_mulmod_preinv(A->rows[i], A->rows[i - 1], n, A->rows[1],
                                     n, poly, len, polyinv, leninv, p);

    fmpz_mat_mul(C, B, A);
    for (i = 0; i < k * len2; i++)
        for (j = 0; j < n; j++)
            fmpz_mod(C->rows[i] + j, C->rows[i] + j, p);

    /* Evaluate block composition using the Horner scheme */
    _fmpz_mod_poly_mulmod_preinv(h, A->rows[m - 1], n, A->rows[1], n, poly,
                                 len, polyinv, leninv, p);

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
                args[i].poly    = (fmpz *) poly;
                args[i].len     = len;
                args[i].polyinv = (fmpz *) polyinv;
                args[i].leninv  = leninv;
                args[i].p       = *p;

                pthread_create(&threads[i], NULL,
                        _fmpz_mod_poly_compose_mod_brent_kung_vec_preinv_worker,
                        &args[i]);
            }
        }
        for (i = 0; i < c; i++)
            pthread_join(threads[i], NULL);
    }

    flint_free(threads);
    flint_free(args);

    _fmpz_vec_clear(h, n);

    fmpz_mat_clear(A);
    fmpz_mat_clear(B);
    fmpz_mat_clear(C);
}

void
fmpz_mod_poly_compose_mod_brent_kung_vec_preinv_threaded(fmpz_mod_poly_struct * res,
                                                const fmpz_mod_poly_struct *
                                                polys, slong len1, slong n,
                                                const fmpz_mod_poly_t poly,
                                                const fmpz_mod_poly_t polyinv)
{
    slong len2 = poly->length;
    slong len3, i;

    for (i = 0; i < len1; i++)
    {
        len3 = (polys + i)->length;
        if (len3 >= len2)
        {
            flint_printf
                ("Exception (fmpz_mod_poly_compose_mod_brent_kung_vec_preinv)."
                 "The degree of the first polynomial must be smaller than that of the "
                 " modulus\n");
            flint_abort();
        }
    }

    if (n > len1)
    {
        flint_printf
            ("Exception (fmpz_mod_poly_compose_mod_brent_kung_vec_preinv)."
             "n is larger than the length of polys\n");
        flint_abort();
    }

    if (n == 0)
        return;

    if (len2 == 1)
    {
        for (i = 0; i < n; i++)
        {
            fmpz_mod_poly_init(res + i, &poly->p);
            fmpz_mod_poly_zero(res + i);
        }
        return;
    }

    if (len2 == 2)
    {
        for (i = 0; i < n; i++)
        {
            fmpz_mod_poly_init(res + i, &poly->p);
            fmpz_mod_poly_set(res + i, polys + i);
        }
        return;
    }

    for (i = 0; i < n; i++)
    {
        fmpz_mod_poly_init2(res + i, &poly->p, len2 - 1);
        _fmpz_mod_poly_set_length(res + i, len2 - 1);
    }

    _fmpz_mod_poly_compose_mod_brent_kung_vec_preinv_threaded(res, polys, len1,
                                                     n, poly->coeffs, len2,
                                                     polyinv->coeffs,
                                                     polyinv->length,
                                                     &poly->p);

    for (i = 0; i < n; i++)
        _fmpz_mod_poly_normalise(res + i);
}
