/*
    Copyright (C) 2025, 2026 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"
#include "acb_poly.h"
#include "acb_ode.h"

void
acb_ode_sol_init(acb_ode_sol_t sol, slong nshifts,
                 slong nlogs, slong npts, slong nder)
{
    sol->series = _acb_poly_vec_init(nlogs);
    sol->sums = _acb_vec_init(npts * nlogs * nder);
    acb_mat_init(sol->extini, nshifts, nlogs);
    sol->tb = _mag_vec_init(nder);

    acb_ode_cvest_init(sol->cvest);

    sol->alloc_logs = nlogs;
    sol->npts = npts;
    sol->nder = nder;

    sol->future_logs = nlogs;
    sol->nlogs = 0;

    acb_ode_sol_zero(sol);
}

void
acb_ode_sol_set_ini(acb_ode_sol_t sol, acb_srcptr ini,
                    const acb_ode_group_t group)
{
    acb_mat_zero(sol->extini);

    for (slong j = 0, i = 0; j < group->nshifts; j++)
        for (slong k = group->shifts[j].mult - 1; k >= 0; k--, i++)
            acb_set(acb_mat_entry(sol->extini, j, k), ini + i);
}

void
acb_ode_sol_clear(acb_ode_sol_t sol)
{
    acb_mat_clear(sol->extini);
    _acb_vec_clear(sol->sums, sol->npts * sol->alloc_logs * sol->nder);
    _acb_poly_vec_clear(sol->series, sol->alloc_logs);
    _mag_vec_clear(sol->tb, sol->nder);
    acb_ode_cvest_clear(sol->cvest);
}

void
acb_ode_sol_fit_length(acb_ode_sol_t sol, slong len)
{
    for (slong k = 0; k < sol->alloc_logs; k++)
    {
        acb_poly_struct * f = sol->series + k;
        acb_poly_fit_length(f, len);
        _acb_poly_set_length(f, f->alloc);  /* for printing */
    }
}

void
acb_ode_sol_zero(acb_ode_sol_t sol)
{
    for (slong k = 0; k < sol->alloc_logs; k++)
        _acb_poly_set_length(sol->series + k, 0);

    sol->future_logs = sol->nlogs + sol->future_logs;
    sol->nlogs = 0;

    _acb_vec_zero(sol->sums, sol->npts * sol->alloc_logs * sol->nder);

    for (slong i = 0; i < sol->nder; i++)
        mag_inf(sol->tb + i);
    sol->done = 0;

    /* some estimates could be kept, but not prec_wanted */
    acb_ode_cvest_clear(sol->cvest);
    acb_ode_cvest_init(sol->cvest);
}
