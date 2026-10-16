/*
    Copyright (C) 2022, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "flint.h"
#include "ulong_extras.h"
#include "arb.h"
#include "mp_real.h"
#include "impl.h"

/* Relation tables for the diophantine argument reductions.

   A table for k primes holds integer relations

       sum_j d_ij alpha_j = epsilon_i,   |epsilon_i| decreasing,

   among alpha_j = log p_j for the first k primes (the exponential,
   exp_diophantine.c), or among alpha_0 = pi/2 and
   alpha_j = 2 arg(pi_j) for the first k nonreal Gaussian primes
   (the trigonometric functions, sin_cos_diophantine.c): the descent
   of _mp_real_log_reduce subtracts nearest multiples of the epsilon_i
   from its argument row by row.  The weights are log2 p_j (the
   exponential's proxy for the size of the prime powers) and
   log N(pi_j) (the norms, as in arb) respectively, with weight 0 on
   the free first element.

   Tables for k = 2, 4, 6, 8, 10, 12, 13, 16, 20, 24, 32, 40, 48 are
   precomputed (rel_tab_data.c); _mp_real_rel_table generates any other
   k on first use through _arb_log_precompute_reductions (fast up to
   ~20 primes, seconds around 48) and caches it per thread, and
   _mp_real_rel_table_is_cached tells whether a call would find the
   table ready.  The returned structure carries the reciprocals
   1/epsilon_i and the weights, formed on first use, and for the
   rational case the primes themselves. */

static FLINT_TLS_PREFIX mp_real_rel_struct *
    _mp_real_rel_cache[2][MP_REAL_REL_MAX + 1] = { { NULL }, { NULL } };
static FLINT_TLS_PREFIX int _mp_real_rel_cleanup_registered = 0;

static void
_mp_real_rel_free(mp_real_rel_struct * t)
{
    if (!t->is_static)
    {
        flint_free((void *) t->d);
        flint_free((void *) t->epsilon);
    }
    flint_free((void *) t->primes);
    flint_free((void *) t->weights);
    flint_free((void *) t->epsilon_inv);
    flint_free(t);
}

static void
_mp_real_rel_cleanup(void)
{
    slong g, i;
    for (g = 0; g < 2; g++)
        for (i = 0; i <= MP_REAL_REL_MAX; i++)
        {
            if (_mp_real_rel_cache[g][i] != NULL)
                _mp_real_rel_free(_mp_real_rel_cache[g][i]);
            _mp_real_rel_cache[g][i] = NULL;
        }
    _mp_real_rel_cleanup_registered = 0;
}

static const mp_real_rel_static_struct *
_mp_real_rel_find_static(int gaussian, slong num)
{
    slong i;
    for (i = 0; i < _mp_real_rel_static_num; i++)
        if (_mp_real_rel_static[i].gaussian == gaussian
            && _mp_real_rel_static[i].num == num)
            return _mp_real_rel_static + i;
    return NULL;
}

int
_mp_real_rel_table_is_cached(int gaussian, slong num)
{
    gaussian = (gaussian != 0);
    if (num < 2 || num > MP_REAL_REL_MAX)
        return 0;
    return _mp_real_rel_cache[gaussian][num] != NULL
        || _mp_real_rel_find_static(gaussian, num) != NULL;
}

/* the primes and weights of a table */
static void
_mp_real_rel_fill_primes(mp_real_rel_struct * t)
{
    slong i;
    float * weights = flint_malloc(t->num * sizeof(float));

    if (t->gaussian)
    {
        t->primes = NULL;
        for (i = 0; i < t->num; i++)
        {
            slong a = _mp_real_gaussian_primes[2 * i];
            slong b = _mp_real_gaussian_primes[2 * i + 1];
            weights[i] = (i == 0) ? 0.0f : (float) log((double) (a * a + b * b));
        }
    }
    else
    {
        ulong * primes = flint_malloc(t->num * sizeof(ulong));
        n_primes_t iter;
        n_primes_init(iter);
        for (i = 0; i < t->num; i++)
        {
            primes[i] = n_primes_next(iter);
            weights[i] = (i == 0) ? 0.0f : (float) log2((double) primes[i]);
        }
        n_primes_clear(iter);
        t->primes = primes;
    }
    t->weights = weights;
}

#define REL_MAX_ROWS 512

const mp_real_rel_struct *
_mp_real_rel_table(int gaussian, slong num)
{
    mp_real_rel_struct * t;
    const mp_real_rel_static_struct * st;
    double * eps_inv;
    slong i;

    gaussian = (gaussian != 0);
    FLINT_ASSERT(num >= 2 && num <= MP_REAL_REL_MAX);

    if (_mp_real_rel_cache[gaussian][num] != NULL)
        return _mp_real_rel_cache[gaussian][num];

    t = flint_malloc(sizeof(mp_real_rel_struct));
    t->num = num;
    t->gaussian = gaussian;

    st = _mp_real_rel_find_static(gaussian, num);
    if (st != NULL)
    {
        t->is_static = 1;
        t->rows = st->rows;
        t->d = st->d;
        t->epsilon = st->epsilon;
    }
    else
    {
        arb_ptr alpha;
        short * d;
        double * eps;
        slong rows;

        /* the generator rounds the values to 3 i + 100 bits at its
           i-th step; 3000 bits cover every row it can produce */
        alpha = _arb_vec_init(num);
        if (gaussian)
            _mp_real_atan_gauss_vec_arb(alpha, num, 3000);
        else
            arb_log_primes_vec_bsplit(alpha, num, 3000);

        d = flint_malloc((REL_MAX_ROWS + 1) * num * sizeof(short));
        eps = flint_malloc((REL_MAX_ROWS + 1) * sizeof(double));
        _arb_log_precompute_reductions(d, eps, alpha, num, REL_MAX_ROWS, 8.0);
        _arb_vec_clear(alpha, num);

        for (rows = 0; d[rows * num] != MP_REAL_REL_TERMINATOR; rows++)
            ;

        t->is_static = 0;
        t->rows = rows;
        t->d = d;
        t->epsilon = eps;
    }

    eps_inv = flint_malloc(t->rows * sizeof(double));
    for (i = 0; i < t->rows; i++)
        eps_inv[i] = 1.0 / t->epsilon[i];
    t->epsilon_inv = eps_inv;
    t->epsilon_min = (t->rows > 0) ? fabs(t->epsilon[t->rows - 1]) : 1.0;

    _mp_real_rel_fill_primes(t);

    _mp_real_rel_cache[gaussian][num] = t;

    if (!_mp_real_rel_cleanup_registered)
    {
        flint_register_cleanup_function(_mp_real_rel_cleanup);
        _mp_real_rel_cleanup_registered = 1;
    }

    return t;
}
