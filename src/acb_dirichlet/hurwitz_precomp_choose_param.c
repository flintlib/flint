/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

#ifdef __GNUC__
# define fabs __builtin_fabs
# define log __builtin_log
# define sqrt __builtin_sqrt
#else
# include <math.h>
#endif

void
acb_dirichlet_hurwitz_precomp_choose_param(ulong * _A, ulong *_K, ulong *_N,
        const acb_t s, double num_eval, slong prec)
{
    double MUL_COST, POW_COST, ZETA_COST, height, abss, cost, best_cost;
    ulong A, K, N, best_A, best_K, best_N;
    mag_t err;

    /* Default algorithm: separate evaluations with no precomputation. */
    best_A = best_K = best_N = 0;
    *_A = *_K = *_N = 0;

    /* This is pointless. */
    if (num_eval < 10)
        return;

    /* Precomputation fails at nonpositive integers. */
    if (acb_contains_int(s) && !arb_is_positive(acb_realref(s)))
        return;

    prec = FLINT_MAX(prec, 8);
    height = arf_get_d(arb_midref(acb_imagref(s)), ARF_RND_DOWN);
    height = fabs(height);
    abss = arf_get_d(arb_midref(acb_realref(s)), ARF_RND_DOWN);
    abss = sqrt(abss*abss + height*height);

    /* Relative evaluation time, estimated empirically. */
    MUL_COST = 1.0;
    POW_COST = 10.0 + FLINT_MIN(0.005 * prec, 200.0);
    ZETA_COST = 200.0 + 2.0*height + (3.0*prec + 0.0002*height*prec)*prec;

    /* Cost for the default algorithm. */
    best_cost = num_eval * ZETA_COST;
    /* Give the default algorithm some more leeway. */
    best_cost *= 0.5;
    if (acb_is_int(s))
        best_cost *= 0.5;

    mag_init(err);

    for (N = 1; N <= FLINT_MIN(num_eval, 2048); N = FLINT_MAX(N+1, N*1.2))
    {
        /* AN should be at least as large as |s| */
        A = FLINT_MAX(abss / N + 1.0, 1);

        /* Estimate K based on asymptotics for the error term. We will
           have to adjust up by actually computing the error bound. */
        K = FLINT_MAX(prec * log(2) / (log(2*A*N) + 1.0) + 1.0, 1);

        for ( ; K < num_eval; K = FLINT_MAX(K+1, K*1.2))
        {
            /* Too much space. */
            if (_acb_vec_estimate_allocated_bytes(N * K, prec) > 1e9)
                break;

            acb_dirichlet_hurwitz_precomp_bound(err, s, A, K, N);

            cost = (K * num_eval) * MUL_COST +     /* polynomial evaluation */
                   (A * num_eval) * POW_COST +     /* power sum */
                   (N * K) * ZETA_COST;            /* precomputation cost */

            /* The accuracy is good enough. */
            if (mag_cmp_2exp_si(err, -prec) <= 0)
            {
                if (cost < best_cost)
                {
                    best_cost = cost;
                    best_A = A;
                    best_K = K;
                    best_N = N;
                }
                break;
            }

            /* It will only get worse from here. */
            if (cost > best_cost)
                break;
        }
    }

    *_A = best_A;
    *_K = best_K;
    *_N = best_N;

    mag_clear(err);
}

