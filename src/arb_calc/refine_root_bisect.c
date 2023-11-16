/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_calc.h"

/* TODO: refactor/combine some of this code with isolate_roots.c */

int arb_calc_partition(arf_interval_t L, arf_interval_t R,
    arb_calc_func_t func, void * param, const arf_interval_t block, slong prec)
{
    arb_t t, m;
    arf_t u;
    int msign;

    arb_init(t);
    arb_init(m);
    arf_init(u);

    /* Compute the midpoint (TODO: try other points) */
    arf_add(u, &block->a, &block->b, ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_mul_2exp_si(u, u, -1);

    /* Evaluate and get sign at midpoint */
    arb_set_arf(m, u);
    func(t, m, param, 1, prec);
    msign = arb_sgn_nonzero(t);

    /* L, R = block, split at midpoint */
    arf_set(&L->a, &block->a);
    arf_set(&R->b, &block->b);
    arf_set(&L->b, u);
    arf_set(&R->a, u);

    arb_clear(t);
    arb_clear(m);
    arf_clear(u);

    return msign;
}

int arb_calc_refine_root_bisect(arf_interval_t r, arb_calc_func_t func,
    void * param, const arf_interval_t start, slong iter, slong prec)
{
    int asign, bsign, msign, result;
    slong i;
    arf_interval_t t, u;
    arb_t m, v;

    arf_interval_init(t);
    arf_interval_init(u);
    arb_init(m);
    arb_init(v);

    arb_set_arf(m, &start->a);
    func(v, m, param, 1, prec);
    asign = arb_sgn_nonzero(v);

    arb_set_arf(m, &start->b);
    func(v, m, param, 1, prec);
    bsign = arb_sgn_nonzero(v);

    /* must have proper sign changes */
    if (asign == 0 || bsign == 0 || asign == bsign)
    {
        result = ARB_CALC_IMPRECISE_INPUT;
    }
    else
    {
        arf_interval_set(r, start);

        result = ARB_CALC_SUCCESS;

        for (i = 0; i < iter; i++)
        {
            msign = arb_calc_partition(t, u, func, param, r, prec);

            /* the algorithm fails if the value at the midpoint cannot
               be distinguished from zero */
            if (msign == 0)
            {
                result = ARB_CALC_NO_CONVERGENCE;
                break;
            }

            if (msign == asign)
                arf_interval_swap(r, u);
            else
                arf_interval_swap(r, t);
        }
    }

    arf_interval_clear(t);
    arf_interval_clear(u);
    arb_clear(m);
    arb_clear(v);

    return result;
}

