/*
    Copyright (C) 2025, 2026 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_types.h"
#include "acb_poly.h"
#include "acb_mat.h"
#include "acb_ode.h"
#include "gr.h"
#include "gr_ore_poly.h"
#include "gr_poly.h"

static void
acb_ode_sum_set_ordinary(acb_ode_sum_t sum)
{
    /* probably later:
       sum->flags |= ACB_ODE_ORDINARY; */

    for (slong m = 0; m < sum->nsols; m++)
    {
        acb_ode_sol_struct * sol = sum->sol + m;
        sol->future_logs = 1 - sol->nlogs;
    }
}

static slong
col(const acb_ode_group_struct * grp, slong s, slong k)
{
    slong j = 0;
    for (slong s1 = 0; s1 <= s; s1++)
        j += grp->shifts[s1].mult;
    return j - 1 - k;
}

static void
fill_column(acb_mat_t mat,
            const acb_ode_group_struct * grp, slong s, slong k,
            const acb_poly_struct * val)
{
    slong j = col(grp, s, k);
    slong len = FLINT_MIN(acb_mat_nrows(mat), val->length);

    slong i = 0;
    for (; i < len; i++)
    {
        acb_ptr a = acb_mat_entry(mat, i, j);
        acb_swap(a, val->coeffs + i);
    }
    for (; i < acb_mat_nrows(mat); i++)
        acb_zero(acb_mat_entry(mat, i, j));
}

static void
fix_column_echelon(acb_mat_struct * mat,
                   const acb_ode_group_struct * grp, slong s, slong k,
                   const acb_mat_t extini, slong prec)
{
    slong mult = grp->shifts[s].mult;
    slong delta = mult - 1 - k;
    slong j = col(grp, s, k);

    for (slong s1 = s + 1; s1 < grp->nshifts; s1++)
    {
        slong mult1 =  grp->shifts[s1].mult;
        /* At most the topmost (highest k) delta elements of the echelonized
           basis associated to the shift s1 have a nonzero coefficient in the
           decomposition on that basis of the Frobenius shift by delta of the
           solution we computed in full. */
        slong k1_0 = FLINT_MAX(0, mult1 - delta);
        slong j1 = col(grp, s1, k1_0);
        for (slong k1 = k1_0; k1 < mult1; k1++)
        {
            FLINT_ASSERT(j1 == col(grp, s1, k1));
            /* Get the coefficient by 'undoing' the Frobenius shift */
            acb_ptr cc = acb_mat_entry(extini, s1, k1 + delta);
            for (slong i = 0; i < acb_mat_nrows(mat); i++)
            {
                acb_ptr a  = acb_mat_entry(mat, i, j);
                acb_ptr a1 = acb_mat_entry(mat, i, j1);
                acb_submul(a, cc, a1, prec);
            }
            j1--;
        }
    }

}

static void
fill_group(acb_mat_t mat, const acb_ode_sum_t sum,
           const acb_ode_group_struct * grp, slong p,
           acb_ode_basis_t basis, slong prec)
{
    for (slong s = grp->nshifts - 1; s >= 0; s--)
    {
        acb_ode_sol_struct * sol = sum->sol + s;
        slong mult = grp->shifts[s].mult;

        acb_poly_struct * val = _acb_poly_vec_init(mult);

        acb_ode_sol_jet(val, sum->group->leader, sol, p, sum->pts + p,
                        acb_mat_nrows(mat), mult, prec);

        for (slong k = 0; k < mult; k++)
        {
            slong delta = mult - 1 - k;
            fill_column(mat, grp, s, k, val + delta);

            switch (basis)
            {
                case ACB_ODE_BASIS_CASCADE:
                    break;
                case ACB_ODE_BASIS_ECHELON:
                    fix_column_echelon(mat, grp, s, k, sum->sol[s].extini, prec);
                    break;
                default:
                    FLINT_ASSERT(0);
            }
        }

        _acb_poly_vec_clear(val, mult);
    }
}

void
_acb_ode_fundamental_matrix_vec(
        acb_mat_struct * mat,
        const acb_poly_struct * dop, slong dop_len,
        const acb_ode_exponents_t expos,
        acb_srcptr lcrt,
        acb_ode_bound_t bound,  /* mutable */
        acb_srcptr pts, slong npts,
        acb_ode_basis_t basis,
        acb_ode_sum_worker_t sum_worker,
        void * worker_data,
        ulong flags,
        slong prec)
{
    mag_t cvrad, mag;

    slong nder = 0;
    for (slong p = 0; p < npts; p++)
        nder = FLINT_MAX(nder, acb_mat_nrows(mat + p));

    if (nder == 0)
        return;

    mag_init(cvrad);
    mag_init(mag);

    _acb_vec_get_mag_lower(cvrad, lcrt, dop[dop_len - 1].length - 1);

    /* XXX maybe not the best place to test this;
       redundant with sum_set_points */
    _acb_vec_get_mag(mag, pts, npts);
    if (mag_cmp(mag, cvrad) >= 0)
    {
        for (slong p = 0; p < npts; p++)
            for (slong i = 0; i < mat[p].r; i++)
                for (slong j = 0; j < dop_len - 1; j++)
                    acb_indeterminate(acb_mat_entry_ptr(mat + p, i, j));
        goto cleanup;
    }

    for (slong g = 0, j = 0; g < expos->ngroups; g++)
    {
        acb_ode_sum_t sum;
        acb_ode_group_bound_t gbound;

        acb_ode_group_struct * group = expos->groups + g;
        slong glen = acb_ode_group_length(expos->groups + g);

        acb_ode_group_bound_init(gbound);
        acb_ode_group_bound_precompute(gbound, dop, dop_len, expos, g, bound, bound->prec);

        acb_ode_sum_init(sum, dop_len, npts, group->nshifts, nder);
        acb_ode_sum_set_diffop(sum, dop, dop_len, cvrad);
        acb_ode_sum_set_group(sum, group);
        /* todo: only avoids part of the redundant computations between
           solutions of the same group in the presence of logs */
        acb_ode_sum_set_ini_highest(sum);
        acb_ode_sum_set_points(sum, pts, npts);
        sum->data = worker_data;

        if (flags & ACB_ODE_ORDINARY)
        {
            FLINT_ASSERT(g == 0);
            acb_ode_sum_set_ordinary(sum);
        }

        sum_worker(sum, -1, bound, gbound, prec);

        for (slong p = 0; p < npts; p++)
        {
            acb_mat_t win;
            acb_mat_window_init(win, mat + p, 0, j, acb_mat_nrows(mat + p), j + glen);

            fill_group(win, sum, expos->groups + g, p, basis, prec);

            acb_mat_window_clear(win);
        }

        acb_ode_sum_clear(sum);
        acb_ode_group_bound_clear(gbound);

        j += glen;
    }

cleanup:
    mag_clear(cvrad);
    mag_clear(mag);
}

/* XXX Redundant with work done in _acb_ode_fundamental_matrix_vec and
   acb_ode_sum_precompute. */
static slong
choose_wp(gr_ore_poly_t dop, acb_srcptr lcrt, acb_srcptr pts, slong npts,
          slong prec, gr_ctx_t Dop)
{
    slong wp;
    mag_t cvrad, mag;

    if (dop->length == 0)
        return 2;

    acb_poly_struct * lc = gr_ore_poly_coeff_ptr(dop, dop->length - 1, Dop);

    if (acb_is_zero(lc->coeffs))
        return 2;

    mag_init(cvrad);
    mag_init(mag);

    _acb_vec_get_mag_lower(cvrad, lcrt, lc->length - 1);
    _acb_vec_get_mag(mag, pts, npts);
    acb_ode_choose_prec(&wp, dop->coeffs, dop->length, mag, cvrad, prec);

    mag_clear(cvrad);
    mag_clear(mag);

    return wp;
}

int
acb_ode_fundamental_matrix_vec(
        acb_mat_struct * mat,
        const gr_ore_poly_t dop, gr_ore_poly_ctx_t Dop,
        const acb_ode_exponents_t expos,
        acb_srcptr lcrt,
        acb_srcptr pts, slong npts,
        acb_ode_basis_t basis,
        ulong flags,
        slong prec)
{
    gr_ctx_struct * Cst, * Pol;
    gr_ctx_t Pol1, Euler, wCC, wCCx, wCCxT;
    acb_ode_exponents_t _expos;
    gr_ore_poly_struct * ldop, * acb_ldop;
    acb_struct * _lcrt = NULL;
    acb_ode_bound_t bound;
    slong wp = prec;

    int status = GR_SUCCESS;

    status |= gr_ore_poly_ctx_over_gr_poly_base_ptrs(&Cst, &Pol, Dop);
    if (status != GR_SUCCESS)
        return status;

    if (GR_ORE_POLY_CTX(Dop)->which_algebra != ORE_ALGEBRA_DERIVATIVE
        && GR_ORE_POLY_CTX(Dop)->which_algebra != ORE_ALGEBRA_EULER_DERIVATIVE)
        return GR_UNABLE;

    if (dop->length < 1)
        return GR_DOMAIN;

    slong order = dop->length - 1;
    gr_poly_struct * lc = (gr_poly_struct *) (dop->coeffs) + order;

    if (lcrt == NULL)  /* XXX make it mandatory? */
    {
        gr_ctx_t bCC;
        gr_ctx_init_complex_acb(bCC, MAG_BITS);

        _lcrt = _acb_vec_init(lc->length - 1);
        status |= acb_ode_lc_roots(_lcrt, dop, Dop, 0, bCC);
        lcrt = _lcrt;

        gr_ctx_clear(bCC);
    }

    acb_ode_local_diffop_ctx_init(Pol1, Euler, Dop, Cst);
    GR_TMP_INIT(ldop, Euler);

    flags |= acb_ode_analyze(dop, flags, Dop);
    status |= acb_ode_local_diffop(ldop, dop, Dop, Euler);
    flags |= acb_ode_analyze(ldop, flags, Euler);

    if (expos == NULL)
    {
        acb_ode_exponents_init(_expos);

        if (flags & ACB_ODE_ORDINARY)
            acb_ode_exponents_ordinary(_expos, order);
        else
        {
            gr_ctx_t CC;
            gr_ctx_init_complex_acb(CC, prec);

            status |= acb_ode_exponents(_expos, ldop, Euler, CC);

            gr_ctx_clear(CC);
        }

        expos = _expos;
    }

    /* When starting with an operator over an exact domain, try to ensure that
       the precision of the operator with ball coefficients we pass downstream
       is not the limiting factor for the computation of the coefficients.
       Eventually we may want to allow generic operators in sum_divconquer and
       related functions so that the approximation can be computed after the
       working precision is chosen for real. */
    if (gr_ctx_is_exact(Cst) == T_TRUE)
    {
        gr_ctx_t bCC, bCCx, bCCxT;
        gr_ctx_init_complex_acb(bCC, MAG_BITS);
        acb_ode_local_diffop_ctx_init(bCCx, bCCxT, Dop, bCC);
        GR_TMP_INIT(acb_ldop, bCCxT);

        status |= gr_poly_set_gr_poly_other((gr_ptr) acb_ldop, (gr_ptr) ldop, Pol1, bCCx);
        if (status == GR_SUCCESS)
            wp = choose_wp(acb_ldop, lcrt, pts, npts, prec, bCCxT);

        GR_TMP_CLEAR(acb_ldop, bCCxT);
        gr_ctx_clear(bCCxT);
        gr_ctx_clear(bCCx);
        gr_ctx_clear(bCC);
    }

    gr_ctx_init_complex_acb(wCC, wp);
    acb_ode_local_diffop_ctx_init(wCCx, wCCxT, Dop, wCC);
    GR_TMP_INIT(acb_ldop, wCCxT);

    if (Cst->which_ring == GR_CTX_CC_ACB)
        /* Avoid rounding to the destination precision.
           todo: something similar for RR_ARB. */
        status |= gr_poly_set((gr_ptr) acb_ldop, (gr_ptr) ldop, wCCx);
    else
        status |= gr_poly_set_gr_poly_other((gr_ptr) acb_ldop, (gr_ptr) ldop, Pol1, wCCx);

    gr_poly_struct * llc = (gr_poly_struct *) (acb_ldop->coeffs) + order;
    if (!(flags & ACB_ODE_NONZERO)
        || (status == GR_SUCCESS && gr_is_zero(llc->coeffs, wCC) == T_TRUE))
        status = GR_UNABLE;

    if (status != GR_SUCCESS)
        goto cleanup;

    if (dop->length == 1)
        goto cleanup;

    acb_ode_bound_init(bound);
    slong dop_clen = _acb_poly_vec_length(acb_ldop->coeffs, acb_ldop->length);
    slong pol_part_len = dop_clen / 2 + 3;
    acb_ode_bound_precompute(bound, acb_ldop->coeffs, acb_ldop->length, lcrt,
                             pol_part_len, MAG_BITS);

    _acb_ode_fundamental_matrix_vec(mat, acb_ldop->coeffs, acb_ldop->length,
                                    expos, lcrt, bound,
                                    pts, npts,
                                    basis,
                                    acb_ode_sum_divconquer, NULL,
                                    flags,
                                    prec);

    acb_ode_bound_clear(bound);
cleanup:
    if (_lcrt == lcrt)
        _acb_vec_clear(_lcrt, lc->length - 1);
    if (expos == _expos)
        acb_ode_exponents_clear(_expos);
    GR_TMP_CLEAR(ldop, Euler);
    GR_TMP_CLEAR(acb_ldop, wCCxT);
    gr_ctx_clear(wCCxT);
    gr_ctx_clear(wCCx);
    gr_ctx_clear(wCC);
    gr_ctx_clear(Pol1);
    gr_ctx_clear(Euler);

    return status;
}

int
acb_ode_fundamental_matrix(
        acb_mat_t mat,
        const gr_ore_poly_t dop, gr_ore_poly_ctx_t Dop,
        const acb_ode_exponents_t expos,
        acb_srcptr lcrt,
        const acb_t pt,
        acb_ode_basis_t basis,
        ulong flags,
        slong prec)
{
    return acb_ode_fundamental_matrix_vec(mat, dop, Dop, expos, lcrt, pt, 1,
                                          basis, flags, prec);
}
