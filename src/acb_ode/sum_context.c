/*
    Copyright (C) 2025, 2026 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_types.h"
#include "acb.h"
#include "acb_mat.h"
#include "acb_poly.h"
#include "acb_ode.h"
#include "arb_poly.h"
#include "fmpz_vec.h"

void
acb_ode_sum_zero(acb_ode_sum_t sum)
{
    sum->n = 0;
    sum->n0 = 0;

    for (slong m = 0; m < sum->nsols; m++)
        acb_ode_sol_zero(sum->sol + m);

    if (sum->nder > 0)
        for (slong i = 0; i < sum->npts; i++)
            acb_one(sum->pows + i);

    mag_one(sum->magpow);
}

void
acb_ode_sum_init(acb_ode_sum_t sum, slong dop_len, slong npts,
                 slong nsols, slong nder)
{
    slong dop_order = dop_len - 1;

    sum->dop = _acb_poly_vec_init(dop_len);
    sum->dop_len = dop_len;
    sum->dop_clen = 0;
    mag_init(sum->cvrad);

    acb_ode_group_init(sum->group, dop_order);
    for (slong s = 0; s < dop_order; s++)
        sum->group->shifts[s].n = -1;
    acb_poly_init(sum->ind);

    sum->sol = flint_malloc(nsols * sizeof(acb_ode_sol_struct));
    for (slong m = 0; m < nsols; m++)
        acb_ode_sol_init(sum->sol + m, 0, 0, 0, 0);  /* TBI */
    sum->nsols = nsols;

    sum->pts = _acb_vec_init(npts);
    sum->pows = _acb_vec_init(npts * nder);
    sum->shifted_sums = flint_malloc(npts);
    sum->npts = npts;
    sum->nder = nder;

    sum->wp = 0;
    sum->sum_wp = 0;

    mag_init(sum->mag);
    mag_init(sum->magpow);
    arb_poly_init(sum->itg_pol);
    arb_poly_init(sum->itg_num);

    sum->flags = 0;
    sum->data = NULL;

    sum->binom_n = _fmpz_vec_init(nder);

    acb_ode_sum_zero(sum);

    sum->have_precomputed = 0;
}

void
acb_ode_sum_clear(acb_ode_sum_t sum)
{
    _acb_poly_vec_clear(sum->dop, sum->dop_len);

    acb_ode_group_clear(sum->group);
    acb_poly_clear(sum->ind);

    for (slong m = 0; m < sum->nsols; m++)
        acb_ode_sol_clear(sum->sol + m);
    flint_free(sum->sol);

    _acb_vec_clear(sum->pts, sum->npts);
    _acb_vec_clear(sum->pows, sum->npts * sum->nder);
    flint_free(sum->shifted_sums);

    mag_clear(sum->mag);
    mag_clear(sum->magpow);
    mag_clear(sum->cvrad);
    arb_poly_clear(sum->itg_pol);
    arb_poly_clear(sum->itg_num);

    _fmpz_vec_clear(sum->binom_n, sum->nder);
}

void
acb_ode_sum_fit_length(acb_ode_sum_t sum, slong len)
{
    for (slong m = 0; m < sum->nsols; m++)
        acb_ode_sol_fit_length(sum->sol + m, len);
}

void
_acb_ode_sum_discard_head(acb_ode_sum_struct * sum, slong n0)
{
    slong len = n0 - sum->n0;

    if (len <= 0)
        return;

    for (slong m = 0; m < sum->nsols; m++)
    {
        for (slong k = 0; k < sum->sol[m].nlogs; k++)
        {
            acb_poly_struct * f = sum->sol[m].series + k;
            /* the sage version uses slightly different sizes near max_terms */
            _acb_poly_shift_right(f->coeffs, f->coeffs, f->length, len);
            _acb_vec_zero(f->coeffs + f->length - len, len);
        }
    }

    sum->n0 = n0;
}

/* Operator */

void
acb_ode_sum_set_diffop(acb_ode_sum_t sum, const acb_poly_struct * dop, slong dop_len,
                       const mag_t cvrad)
{
    FLINT_ASSERT(dop_len > 0);
    FLINT_ASSERT(dop_len == sum->dop_len);

    _acb_poly_vec_set(sum->dop, dop, dop_len);
    mag_set(sum->cvrad, cvrad);
    sum->have_precomputed = 0;
}

/* Solution group */

void
acb_ode_sum_set_group(acb_ode_sum_t sum, const acb_ode_group_t group)
{
    FLINT_ASSERT(sum->dop_len < 2 || sum->group->shifts->n == -1);
    FLINT_ASSERT(!sum->have_precomputed);

    /* Upper bound, but some tests currently depend on it.
       Will need updating to support inhomogeneous equations. */
    slong nlogs = sum->dop_len - 1;

    for (slong m = 0; m < sum->nsols; m++)
    {
        acb_ode_sol_clear(sum->sol + m);
        acb_ode_sol_init(sum->sol + m, group->nshifts, nlogs, sum->npts, sum->nder);
    }

    acb_ode_group_set(sum->group, group);
}

/* Initial values */

static void
acb_ode_sol_unit_ini(acb_ode_sol_t sol, slong i0,
                     const acb_ode_shift_struct * shifts, slong nshifts)
{
    for (slong i = 0, s = 0; s < nshifts; s++)
    {
        slong mult = shifts ? shifts[s].mult : 1;
        for (slong k = mult - 1; k >= 0; k--, i++)
            acb_set_si(acb_mat_entry(sol->extini, s, k), i == i0 ? 1 : 0);
    }
}

void
acb_ode_sum_set_ini_echelon(acb_ode_sum_t sum)
{
    for (slong m = 0; m < sum->nsols; m++)
        acb_ode_sol_unit_ini(sum->sol + m, m, sum->group->shifts,
                             sum->group->nshifts);
}

void
acb_ode_sum_set_ini_highest(acb_ode_sum_t sum)
{
    for (slong m = 0; m < sum->nsols; m++)
    {
        acb_mat_zero(sum->sol[m].extini);
        acb_one(acb_mat_entry(sum->sol[m].extini, m,
                              sum->group->shifts[m].mult - 1));
    }
}

/* Evaluation points */

void
acb_ode_sum_set_points(acb_ode_sum_t sum, acb_srcptr pts, slong npts)
{
    FLINT_ASSERT(sum->npts == npts);
    _acb_vec_set(sum->pts, pts, npts);
    sum->have_precomputed = 0;
}

/* Main precomputation function */

void
acb_ode_sum_precompute(acb_ode_sum_t sum)
{
    if (sum->have_precomputed)
        return;

    slong dop_len = sum->dop_len;

    sum->dop_clen = _acb_poly_vec_length(sum->dop, dop_len);

    /* XXX _gr_ore_poly_indicial_polynomial? */
    acb_poly_fit_length(sum->ind, dop_len);
    _acb_poly_set_length(sum->ind, dop_len);
    for (slong i = 0; i < dop_len; i++)
        acb_poly_get_coeff_acb((sum->ind)->coeffs + i, sum->dop + i, 0);
    _acb_poly_normalise(sum->ind);

    for (slong j = 0; j < sum->npts; j++)  /* XXX move? */
    {
        /* We can trade some full-width muls for muls by small integers by
           multiplying everyone by ξ^n, but this optimization is only
           valid when ξ does not contain 0. */
        sum->shifted_sums[j] = !acb_contains_zero(sum->pts + j);
    }

    _acb_vec_get_mag(sum->mag, sum->pts, sum->npts);
    FLINT_ASSERT(mag_cmp(sum->mag, sum->cvrad) <= 0);

    sum->have_precomputed = 1;
}

/* Misc */

slong
acb_ode_sum_max_nlogs(acb_ode_sum_struct * sum)
{
    slong nlogs = 0;
    for (slong m = 0; m < sum->nsols; m++)
        nlogs = FLINT_MAX(nlogs, sum->sol[m].nlogs);
    return nlogs;
}

slong
acb_ode_sum_max_nlogs_xn(acb_ode_sum_struct * sum, slong n)
{
    slong off = n - sum->n0;

    for (slong nlogs = acb_ode_sum_max_nlogs(sum); nlogs > 0; nlogs--)
    {
        for (slong m = 0; m < sum->nsols; m++)
        {
            if (nlogs <= sum->sol[m].nlogs)
            {
                acb_poly_struct * p = sum->sol[m].series + nlogs - 1;
                if (off < p->length && !acb_is_zero(p->coeffs + off))
                    return nlogs;
            }
        }
    }
    return 0;
}
