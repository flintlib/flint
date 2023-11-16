/*
    Copyright (C) 2016 Pascal Molin
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

void
acb_dirichlet_roots_init(acb_dirichlet_roots_t t, ulong order, slong num, slong prec)
{
    slong k, size, depth, best_depth, wp;
    ulong reduced_order;
    double cost, best_cost;

    /* exploit 90 deg symmetries */
    if (order % 4 == 0)
        reduced_order = order / 8 + 1;
    else if (order % 2 == 0)
        reduced_order = order / 4 + 1;
    else
        reduced_order = order / 2 + 1;

    wp = prec + 6 + 2 * FLINT_BIT_COUNT(reduced_order);

    t->order = order;
    t->reduced_order = reduced_order;
    t->use_pow = 0;

    if (reduced_order <= 2 || num <= 2)
    {
        depth = 0;
        size = 0;
    }
    else
    {
        /* At depth = d and reduced_order = n, for k evaluations we need about
              log_2(n) * k               muls if d == 0
              d * n^(1/d) + (d-1) * k    muls if d > 0     */
        best_cost = FLINT_BIT_COUNT(reduced_order) * (double) num;
        best_depth = 0;

        for (depth = 1; depth <= 4; depth++)
        {
            size = n_root(reduced_order, depth) + 1;

            /* limit memory usage */
            if (depth * _acb_vec_estimate_allocated_bytes(size, wp) > 1e9)
                continue;

            cost = depth * (double) size + (depth - 1) * (double) num;

            if (cost < best_cost)
            {
                best_depth = depth;
                best_cost = cost;
            }
        }

        depth = best_depth;
        size = n_root(reduced_order, depth) + 1;
    }

    t->size = size;
    t->depth = depth;
    acb_init(t->z);

    if (depth != 0)
    {
        acb_struct * z;
        acb_unit_root(t->z, order, wp);
        z = t->z;
        t->Z = flint_malloc(depth * sizeof(acb_ptr));

        /* todo: at the last level, we could avoid computing
           entries that will never be reached */
        for (k = 0; k < depth; k++)
        {
            t->Z[k] = _acb_vec_init(size + 1);
            _acb_vec_set_powers(t->Z[k], z, size + 1, wp);
            z = t->Z[k] + size;
        }
    }
    else
    {
        /* this tuning could be improved */
        if (reduced_order < 30)
            t->use_pow = 1;
        else if (reduced_order < 100)
            t->use_pow = (prec >= 512);
        else if (reduced_order < 10000)
            t->use_pow = (prec >= 4096);
        else
            t->use_pow = (prec >= 16384);

        if (t->use_pow)
            acb_unit_root(t->z, order, wp);

        t->Z = NULL;
    }
}

