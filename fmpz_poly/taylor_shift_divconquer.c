/*
    Copyright (C) 2012, 2016 Fredrik Johansson
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
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
}
worker_t;

static void
_fmpz_poly_taylor_shift_divconquer_worker(void * arg_ptr)
{
    worker_t * data = (worker_t *) arg_ptr;
    _fmpz_poly_taylor_shift_divconquer(data->poly, data->c, data->len);
}

void
_fmpz_poly_taylor_shift_divconquer(fmpz * poly, const fmpz_t c, slong len)
{
    fmpz *tmp, *tmp2;
    slong k, len1, len2;
    slong bits, cutoff;
    slong nt, nw;
    thread_pool_handle * threads;
    worker_t args[2];

    if (len < 50 || fmpz_is_zero(c))
    {
        _fmpz_poly_taylor_shift_horner(poly, c, len);
        return;
    }

    bits = _fmpz_vec_max_bits(poly, len);
    bits = FLINT_ABS(bits);

    nt = flint_get_num_threads();

    cutoff = 100 + 10 * n_sqrt(FLINT_MAX(bits - FLINT_BITS, 0));

    /* Parallel cutoff is set lower since shift_horner is serial. */
    if (nt == 1)
        cutoff = FLINT_MIN(cutoff, 1000);
    else
        cutoff = FLINT_MIN(cutoff, 300);

    if (len < cutoff)
    {
        _fmpz_poly_taylor_shift_horner(poly, c, len);
        return;
    }

    len1 = len / 2;
    len2 = len - len1;

    nw = flint_request_threads(&threads, FLINT_MIN(nt, 2)); /* request one extra worker */

    if (len < 200 || len + bits < 2000 || nw == 0)
    {
        _fmpz_poly_taylor_shift_divconquer(poly, c, len1);
        _fmpz_poly_taylor_shift_divconquer(poly + len1, c, len2);
    }
    else
    {
        slong nw_save;
        
        args[0].poly = poly;
        args[0].c = c;
        args[0].len = len1;

        args[1].poly = poly + len1;
        args[1].c = c;
        args[1].len = len2;

        FLINT_ASSERT(nt >= 2);

        nw_save = flint_set_num_workers(nt - nt/2 - 1);
        
        thread_pool_wake(global_thread_pool, threads[0], nt/2 - 1,
                          _fmpz_poly_taylor_shift_divconquer_worker, &args[1]);

        _fmpz_poly_taylor_shift_divconquer_worker(&args[0]);
        flint_reset_num_workers(nw_save);

        thread_pool_wait(global_thread_pool, threads[0]);
    }

    flint_give_back_threads(threads, nw);

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
fmpz_poly_taylor_shift_divconquer(fmpz_poly_t g, const fmpz_poly_t f,
    const fmpz_t c)
{
    if (f != g)
        fmpz_poly_set(g, f);

    _fmpz_poly_taylor_shift_divconquer(g->coeffs, c, g->length);
}

