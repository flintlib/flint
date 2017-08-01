/*
    Copyright (C) 2012, 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <pthread.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "nmod_poly.h"

typedef struct
{
    fmpz * poly;
    const fmpz * c;
    slong len;
    slong num_threads;
    slong num_total_threads;
}
worker_t;

void _fmpz_poly_taylor_shift_dc(fmpz * poly,
    const fmpz_t c, slong len, slong num_total_threads);

static void *
_fmpz_poly_taylor_shift_dc_worker(void * arg_ptr)
{
    worker_t * data = (worker_t *) arg_ptr;
    flint_set_num_threads(data->num_threads);
    _fmpz_poly_taylor_shift_dc(data->poly, data->c, data->len,
                               data->num_total_threads);
    flint_cleanup();
    return NULL;
}

void
_fmpz_poly_taylor_shift_dc(fmpz * poly, const fmpz_t c, slong len,
                                        slong num_total_threads)
{
    fmpz *tmp, *tmp2;
    slong k, len1, len2;
    slong bits, max_horner;

    if (len < 64 || fmpz_is_zero(c))
    {
        _fmpz_poly_taylor_shift_horner(poly, c, len);
        return;
    }

    bits = _fmpz_vec_max_bits(poly, len);
    bits = FLINT_ABS(bits);

    /* A big problem for parallel tuning is that we have no
       multithreading in shift_horner. Moreover, shift_horner is
       currently implemented without good cache locality, which makes
       it perform much worse when run on several instances in parallel.
       Therefore, only use it for smaller lengths. */
    if (num_total_threads == 1)
        max_horner = 3000;
    else
        max_horner = 200;

    /* Horner is faster with huge coefficients. */
    if (len < max_horner && bits > pow(2.0, 7.0 + len * 0.005))
    {
        _fmpz_poly_taylor_shift_horner(poly, c, len);
        return;
    }

    len1 = len / 2;
    len2 = len - len1;

    if (len < 200 || len + bits < 2000 || flint_get_num_threads() == 1)
    {
        _fmpz_poly_taylor_shift_dc(poly, c, len1, num_total_threads);
        _fmpz_poly_taylor_shift_dc(poly + len1, c, len2, num_total_threads);
    }
    else
    {
        pthread_t threads[2];
        worker_t args[2];

        args[0].poly = poly;
        args[0].c = c;
        args[0].len = len1;
        args[0].num_threads = flint_get_num_threads() / 2;

        if (num_total_threads == 1)
            args[0].num_total_threads = flint_get_num_threads();
        else
            args[0].num_total_threads = num_total_threads;

        args[1].poly = poly + len1;
        args[1].c = c;
        args[1].len = len2;
        args[1].num_threads = args[0].num_threads;
        args[1].num_total_threads = args[0].num_total_threads;

        pthread_create(&threads[0], NULL,
            _fmpz_poly_taylor_shift_dc_worker, &args[0]);
        pthread_create(&threads[1], NULL,
            _fmpz_poly_taylor_shift_dc_worker, &args[1]);

        pthread_join(threads[0], NULL);
        pthread_join(threads[1], NULL);
    }

    tmp = _fmpz_vec_init(len1 + 1);
    tmp2 = _fmpz_vec_init(len);

    /* Now we generate (x+c)^len1 using binomial expansion. It's redundant
       to do this in all branches of the tree, but since it's just O(d),
       it's going to be cheap compared to the actual multiplications
       anyway. */
    fmpz_one(tmp);
    for (k = 1; k <= len1; k++)
    {
        if (k > len1 - k)
        {
            fmpz_set(tmp + k, tmp + len1 - k);
        }
        else
        {
            fmpz_mul_ui(tmp + k, tmp + k - 1, len1 + 1 - k);
            fmpz_divexact_ui(tmp + k, tmp + k, k);
        }
    }

    if (!fmpz_is_one(c))
    {
        if (fmpz_cmp_si(c, -1) == 0)
        {
            for (k = len1 - 1; k >= 0; k -= 2)
                fmpz_neg(tmp + k, tmp + k);
        }
        else
        {
            fmpz_set(tmp2, c);

            for (k = len1 - 1; k >= 0; k--)
            {
                fmpz_mul(tmp + k, tmp + k, tmp2);
                fmpz_mul(tmp2, tmp2, c);
            }
        }
    }

    _fmpz_poly_mul(tmp2, tmp, len1 + 1, poly + len1, len2);

    _fmpz_vec_add(poly, poly, tmp2, len1);
    _fmpz_vec_set(poly + len1, tmp2 + len1, len2);

    _fmpz_vec_clear(tmp, len1 + 1);
    _fmpz_vec_clear(tmp2, len);
}

void
_fmpz_poly_taylor_shift_divconquer(fmpz * poly, const fmpz_t c, slong len)
{
    _fmpz_poly_taylor_shift_dc(poly, c, len, 1);
}

void
fmpz_poly_taylor_shift_divconquer(fmpz_poly_t g, const fmpz_poly_t f,
    const fmpz_t c)
{
    if (f != g)
        fmpz_poly_set(g, f);

    _fmpz_poly_taylor_shift_divconquer(g->coeffs, c, g->length);
}

