/*
    Copyright (C) 2026 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_types.h"
#include "acb.h"
#include "acb_poly.h"
#include "acb_ode.h"
#include "arb_poly.h"

static int
group_has_shift_past(const acb_ode_group_t group, slong n)
{
    for (slong s = 0; s < group->nshifts; s++)
        if (group->shifts[s].n >= n)
            return 1;

    return 0;
}

static void
_mag_vec_inf(mag_ptr mag, slong len)
{
    for (slong i = 0; i < len; i++)
        mag_inf(mag + i);
}

/* todo: support the case when only part of the residual is known */
static void
bound_normalized_residuals(arb_poly_struct * nres_maj,
                           acb_ode_sum_struct * sum,
                           slong prec)
{
    acb_t cst;
    arf_t tmparf;
    acb_poly_t ind_nd;
    TMP_INIT;
    TMP_START;

    acb_init(cst);
    arf_init(tmparf);

    acb_poly_struct * dop_lc = sum->dop + sum->dop_len - 1;

    acb_poly_init(ind_nd);
    slong nlogs0 = acb_ode_sum_max_nlogs(sum);
    acb_ptr rhs_term = TMP_ALLOC(sizeof(acb_struct) * nlogs0);
    acb_ptr nres_term = _acb_vec_init(nlogs0);

    for (slong m = 0; m < sum->nsols; m++)
    {
        arb_poly_fit_length(nres_maj + m, sum->dop_clen);
        _arb_poly_set_length(nres_maj + m, sum->dop_clen);
    }

    for (slong d = 0, s = 0; d < sum->dop_clen - 1; d++)
    {
        slong n = sum->n;

        /* Multiplicity of λ+n+d as an indicial root */
        while (s < sum->group->nshifts && sum->group->shifts[s].n < n + d)
            s++;
        slong mult = 0;
        if (s < sum->group->nshifts && sum->group->shifts[s].n == n + d)
            mult =  sum->group->shifts[s].mult;

        /* ind_nd/cst = ± jet of monic indicial polynomial at λ+n+d */
        /* XXX change acb_ode_sum_max_nlogs to take k (=n+d here) or k - n instead of k - n0 */
        slong nlogs_d = acb_ode_sum_max_nlogs_xn(sum, n + d);

        acb_ode_poly_taylor_shift_aps_trunc(ind_nd, sum->ind, sum->group->leader,
                                            n + d, nlogs_d + mult, prec);
        acb_set_round(cst, dop_lc->coeffs, prec);
        acb_neg(cst, cst);

        acb_poly_fit_length(ind_nd, nlogs_d + mult);
        for (slong m = 0; m < sum->nsols; m++)  /* XXX skip converged? */
        {
            slong nlogs_dm = FLINT_MIN(sum->sol[m].nlogs, nlogs_d);  /* XXX tighten? */

            /* nres_term = coeff of x^{λ+n+d} in normalized residual of sol m */
            for (slong k = 0; k < nlogs_dm; k++)
                rhs_term[k] = (sum->sol[m].series + k)->coeffs[n + d - sum->n0];
            _acb_ode_poly_negdivrevhigh(nres_term, ind_nd->coeffs + mult,
                                        cst, rhs_term, nlogs_dm, prec);

            /* nres_maj[m] = common majorant of the coefficients of log(x)^k/k!
               in the normalized residual of solution m */
            arf_ptr c = arb_midref((nres_maj + m)->coeffs + d);
            for (slong k = 0; k < nlogs_dm; k++)
            {
                acb_get_abs_ubound_arf(tmparf, nres_term + k, prec);
                if (arf_is_nan(tmparf))
                {
                    arf_nan(c);
                    break;
                }
                arf_max(c, c, tmparf);
            }
        }
    }

    for (slong m = 0; m < sum->nsols; m++)
        _arb_poly_normalise(nres_maj + m);

    _acb_vec_clear(nres_term, nlogs0);
    acb_poly_clear(ind_nd);
    acb_clear(cst);
    arf_clear(tmparf);
    TMP_END;
}

double
_acb_ode_sum_update_tail_bounds(acb_ode_sum_struct * sum,
                                acb_ode_bound_t bound,
                                acb_ode_group_bound_t gbound,
                                slong prec)
{
    arb_t rad;
    arb_poly_t jet_tb;
    mag_t tb1;
    TMP_INIT;

    /* tail_bound_jet wants n past the largest shift */
    if (group_has_shift_past(sum->group, sum->n))
    {
        for (slong m = 0; m < sum->nsols; m++)
        {
            mag_ptr tb = (sum->sol + m)->tb;
            _mag_vec_inf(tb, sum->nder);
        }
        return 0.;
    }

    TMP_START;
    arb_init(rad);
    arb_poly_init(jet_tb);
    arb_poly_struct * nres_maj = TMP_ALLOC(sizeof(arb_poly_t) * sum->nsols);
    for (slong m = 0; m < sum->nsols; m++)
        arb_poly_init(nres_maj + m);
    mag_init(tb1);

    /* ore_algebra uses a bound on nlogs:
    slong nlogs = acb_ode_group_nlogs(sum->group, sum->n);
    */

    /* todo: Only call precompute_integrand if the current bounds are not tight
       enough or nlogs has changed. (We can also use the larger nlogs from the
       start, but this leads to worse bounds before nlogs increases.) In fact,
       for large n, the result will barely depend on n at all... */

    acb_ode_bound_precompute_integrand(sum->itg_pol, sum->itg_num, sum->group,
                                       bound, gbound, sum->n,
                                       acb_ode_sum_max_nlogs(sum), MAG_BITS);

    bound_normalized_residuals(nres_maj, sum, prec);

    /* ore_algebra uses a single nres_maj for all solutions:

    for (slong m = 1; m < sum->nsols; m++)
    {
        arb_poly_fit_length(nres_maj, nres_maj[m].length);
        for (slong i = 0; i < nres_maj[m].length; i++)
            arb_max(nres_maj->coeffs + i, nres_maj->coeffs + i,
                    nres_maj[m].coeffs + i, prec);
        nres_maj->length = FLINT_MAX(nres_maj->length, nres_maj[m].length);
    }
    _arb_poly_normalise(nres_maj);
    */

    double progress = 0.;

    for (slong m = 0; m < sum->nsols; m++)
    {
        slong i;
        mag_ptr tb = (sum->sol + m)->tb;

        if (sum->sol[m].done)
            continue;

        arf_set_mag(arb_midref(rad), sum->mag);
        acb_ode_tail_bound_jet_precomp(jet_tb, bound, sum->n, sum->itg_pol,
                                       sum->itg_num, nres_maj + m, rad,
                                       sum->nder, MAG_BITS);
        for (i = 0; i < jet_tb->length; i++)
        {
            arb_get_mag(tb1, jet_tb->coeffs + i);
            if (i == 0)
            {
                double om0 = mag_get_d_log2_approx(tb + i);
                double om1 = mag_get_d_log2_approx(tb1);
                progress += om0 - om1;
            }
            mag_min(tb + i, tb + i, tb1);
        }
        for (; i < sum->nder; i++)
            mag_zero(tb + i);
    }

    mag_clear(tb1);
    arb_poly_clear(jet_tb);
    arb_clear(rad);
    for (slong m = 0; m < sum->nsols; m++)
        arb_poly_clear(nres_maj + m);
    TMP_END;

    return progress;
}
