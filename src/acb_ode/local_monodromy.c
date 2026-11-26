/*
    Copyright (C) 2026 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* generated using Claude Code (Claude Fable 5, claude-opus-5-5) */

#include "acb_types.h"
#include "acb.h"
#include "acb_mat.h"
#include "acb_ode.h"
#include "gr.h"
#include "gr_ore_poly.h"
#include "gr_poly.h"
#include "ulong_extras.h"

/* On a group with leader λ, a counterclockwise loop around the origin acts by
   x^λ ↦ e^{2πiλ} x^λ and log x ↦ log x + 2πi, so that the corresponding
   diagonal block of the monodromy matrix is e^{2πiλ} exp(2πi N), where N is
   the (nilpotent) matrix in the echelon basis of the log-lowering map
   f ↦ Σ_k f_{k+1} log(x)^k/k!. Column c of N consists of the coefficients of
   x^{λ+n} log(x)^{k+1}/(k+1)! of the c-th basis element, for the pairs (n, k)
   associated with the rows of the block, which the summation code stores in
   extini. */
void
_acb_ode_local_monodromy(acb_mat_t mono,
                         const acb_poly_struct * dop, slong dop_len,
                         const acb_ode_exponents_t expos, slong prec)
{
    acb_t twopii, f;
    mag_t cvrad;

    acb_init(twopii);
    acb_init(f);
    mag_init(cvrad);  /* unused when there are no evaluation points */

    acb_const_pi(twopii, prec);
    acb_mul_onei(twopii, twopii);
    acb_mul_2exp_si(twopii, twopii, 1);

    acb_mat_zero(mono);

    for (slong g = 0, j0 = 0; g < expos->ngroups; g++)
    {
        const acb_ode_group_struct * grp = expos->groups + g;
        slong glen = acb_ode_group_length(grp);
        acb_ode_sum_t sum;
        acb_mat_t blk, tmp;

        /* compute the echelon solutions until all shifts have been reached */
        acb_ode_sum_init(sum, dop_len, 0, glen, 0);
        acb_ode_sum_set_diffop(sum, dop, dop_len, cvrad);
        acb_ode_sum_set_group(sum, grp);
        acb_ode_sum_set_ini_echelon(sum);
        acb_ode_sum_divconquer(sum, grp->shifts[grp->nshifts - 1].n + 1,
                               NULL, NULL, prec);

        acb_mat_init(tmp, glen, glen);

        /* 2πi N; the columns of the block and the solutions sum->sol[c] are
           both ordered by increasing shift and decreasing log degree */
        for (slong c = 0; c < glen; c++)
        {
            const acb_mat_struct * extini = sum->sol[c].extini;
            for (slong s = 0, r = 0; s < grp->nshifts; s++)
                for (slong k = grp->shifts[s].mult - 1; k >= 0; k--, r++)
                    if (k + 1 < acb_mat_ncols(extini))
                        acb_mul(acb_mat_entry(tmp, r, c),
                                acb_mat_entry(extini, s, k + 1), twopii, prec);
        }

        /* exact truncation since N is (structurally) nilpotent */
        acb_mat_exp(tmp, tmp, prec);

        acb_mul_2exp_si(f, grp->leader, 1);
        acb_exp_pi_i(f, f, prec);
        acb_mat_window_init(blk, mono, j0, j0, j0 + glen, j0 + glen);
        acb_mat_scalar_mul_acb(blk, tmp, f, prec);
        acb_mat_window_clear(blk);

        acb_mat_clear(tmp);
        acb_ode_sum_clear(sum);

        j0 += glen;
    }

    mag_clear(cvrad);
    acb_clear(f);
    acb_clear(twopii);
}

int
acb_ode_local_monodromy(acb_mat_t mono,
                        const gr_ore_poly_t dop, gr_ore_poly_ctx_t Dop,
                        const acb_ode_exponents_t expos,
                        ulong flags, slong prec)
{
    gr_ctx_struct * Cst, * Pol;
    gr_ctx_t Pol1, Euler, wCC, wCCx, wCCxT;
    gr_ore_poly_struct * ldop, * acb_ldop;
    acb_ode_exponents_t _expos;

    int status = GR_SUCCESS;

    status |= gr_ore_poly_ctx_over_gr_poly_base_ptrs(&Cst, &Pol, Dop);
    if (status != GR_SUCCESS)
        return status;

    if (dop->length < 1)
        return GR_DOMAIN;

    slong order = dop->length - 1;

    if (acb_mat_nrows(mono) != order || acb_mat_ncols(mono) != order)
        return GR_DOMAIN;

    if (order == 0)
        return GR_SUCCESS;

    acb_ode_exponents_init(_expos);
    acb_ode_local_diffop_ctx_init(Pol1, Euler, Dop, Cst);
    GR_TMP_INIT(ldop, Euler);

    flags |= acb_ode_analyze(dop, flags, Dop);
    status |= acb_ode_local_diffop(ldop, dop, Dop, Euler);
    flags |= acb_ode_analyze(ldop, flags, Euler);

    if (status != GR_SUCCESS)
        goto cleanup;

    if (flags & ACB_ODE_ORDINARY)
    {
        acb_mat_one(mono);
        goto cleanup;
    }

    if (!(flags & ACB_ODE_REGULAR))
    {
        status = GR_UNABLE;
        goto cleanup;
    }

    if (expos == NULL)
    {
        gr_ctx_t CC;
        gr_ctx_init_complex_acb(CC, prec);
        status |= acb_ode_exponents(_expos, ldop, Euler, CC);
        gr_ctx_clear(CC);
        expos = _expos;
    }

    if (status != GR_SUCCESS)
        goto cleanup;

    if (acb_ode_exponents_length(expos) != order)
    {
        status = GR_DOMAIN;
        goto cleanup;
    }

    /* Precision used for the recurrence by acb_ode_sum_divconquer, so that
       rounding the operator when converting it from an exact ring is not the
       limiting factor. */
    slong dop_clen = 0;
    for (slong i = 0; i < ldop->length; i++)
    {
        const gr_poly_struct * c = gr_ore_poly_coeff_srcptr(ldop, i, Euler);
        dop_clen = FLINT_MAX(dop_clen, c->length);
    }
    slong wp = 8 + (prec + 2 * ldop->length)
                   * (1 + n_clog(FLINT_MAX(1, dop_clen - 3), 2));

    gr_ctx_init_complex_acb(wCC, wp);
    acb_ode_local_diffop_ctx_init(wCCx, wCCxT, Dop, wCC);
    GR_TMP_INIT(acb_ldop, wCCxT);

    if (Cst->which_ring == GR_CTX_CC_ACB)
        status |= gr_poly_set((gr_ptr) acb_ldop, (gr_ptr) ldop, wCCx);
    else
        status |= gr_poly_set_gr_poly_other((gr_ptr) acb_ldop, (gr_ptr) ldop,
                                            Pol1, wCCx);

    if (status == GR_SUCCESS && acb_ldop->length != ldop->length)
        status = GR_UNABLE;

    if (status == GR_SUCCESS)
        _acb_ode_local_monodromy(mono, acb_ldop->coeffs, acb_ldop->length,
                                 expos, prec);

    GR_TMP_CLEAR(acb_ldop, wCCxT);
    gr_ctx_clear(wCCxT);
    gr_ctx_clear(wCCx);
    gr_ctx_clear(wCC);

cleanup:
    GR_TMP_CLEAR(ldop, Euler);
    gr_ctx_clear(Euler);
    gr_ctx_clear(Pol1);
    acb_ode_exponents_clear(_expos);

    return status;
}
