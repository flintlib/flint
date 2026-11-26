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
#include "acb_hypgeom.h"
#include "acb_mat.h"
#include "acb_ode.h"
#include "acb_poly.h"
#include "arb.h"
#include "fmpz.h"
#include "gr.h"
#include "gr_poly.h"
#include "gr_ore_poly.h"
#include "gr_vec.h"

/* Sets dop (in Dop = R[x][d/dx]) to an operator annihilating
   x ↦ pFq(a; b; 1 + x), namely the hypergeometric operator
       θ (θ + b_1 - 1) ··· (θ + b_q - 1) - z (θ + a_1) ··· (θ + a_p),
   θ = z d/dz, written in terms of d/dz and recentered at z = 1. */
static int
hypgeom_operator_at_1(gr_ore_poly_t dop, const gr_vec_t a, const gr_vec_t b,
                      gr_ctx_t R, gr_ore_poly_ctx_t Dop)
{
    gr_ctx_struct * Pol = GR_ORE_POLY_CTX(Dop)->base_ring;
    slong p = a->length, q = b->length;
    slong order = q + 1;
    gr_ore_poly_ctx_t Theta;
    gr_ore_poly_t L;
    gr_poly_t P1, P2;
    gr_vec_t rts;
    gr_ptr t;
    slong sh;

    int status = GR_SUCCESS;

    FLINT_ASSERT(p <= order);

    gr_ore_poly_ctx_init(Theta, Pol, 0, ORE_ALGEBRA_EULER_DERIVATIVE);
    gr_ore_poly_init(L, Theta);
    gr_poly_init(P1, R);
    gr_poly_init(P2, R);
    gr_vec_init(rts, order, R);
    GR_TMP_INIT(t, R);

    /* P1 = θ ∏ (θ + b_j - 1), P2 = ∏ (θ + a_i) as polynomials in θ */

    status |= gr_zero(gr_vec_entry_ptr(rts, 0, R), R);
    for (slong j = 0; j < q; j++)
    {
        gr_ptr r = gr_vec_entry_ptr(rts, 1 + j, R);
        status |= gr_neg(r, gr_vec_entry_srcptr(b, j, R), R);
        status |= gr_add_si(r, r, 1, R);
    }
    status |= gr_poly_product_roots(P1, rts, R);

    gr_vec_set_length(rts, p, R);
    for (slong i = 0; i < p; i++)
        status |= gr_neg(gr_vec_entry_ptr(rts, i, R),
                         gr_vec_entry_srcptr(a, i, R), R);
    status |= gr_poly_product_roots(P2, rts, R);

    /* L = P1(θ) - z P2(θ) */

    gr_ore_poly_fit_length(L, order + 1, Theta);
    for (slong k = 0; k <= order; k++)
    {
        gr_poly_struct * c = gr_ore_poly_coeff_ptr(L, k, Theta);
        status |= gr_poly_zero(c, R);
        if (k < P1->length)
            status |= gr_poly_set_coeff_scalar(c, 0,
                    gr_poly_coeff_srcptr(P1, k, R), R);
        if (k < P2->length)
        {
            status |= gr_neg(t, gr_poly_coeff_srcptr(P2, k, R), R);
            status |= gr_poly_set_coeff_scalar(c, 1, t, R);
        }
    }
    _gr_ore_poly_set_length(L, order + 1, Theta);
    _gr_ore_poly_normalise(L, Theta);

    /* rewrite in terms of d/dz and substitute z = 1 + x */

    status |= gr_ore_poly_convert(dop, &sh, L, Dop, Theta);

    status |= gr_one(t, R);
    for (slong i = 0; i < dop->length; i++)
    {
        gr_poly_struct * c = gr_ore_poly_coeff_ptr(dop, i, Dop);
        status |= gr_poly_taylor_shift(c, c, t, R);
    }

    GR_TMP_CLEAR(t, R);
    gr_vec_clear(rts, R);
    gr_poly_clear(P2, R);
    gr_poly_clear(P1, R);
    gr_ore_poly_clear(L, Theta);
    gr_ore_poly_ctx_clear(Theta);

    return status;
}

/* Exponent structure at 1 of the hypergeometric equation. When p = q + 1, the
   exponents are {0, ..., q-1} ∪ {s} with s = Σ b_j - Σ a_i; otherwise, 1 is
   an ordinary point. The shift arrays of grps[0] and grps[1] must have room
   for q + 1 and 1 elements respectively. */
static int
hypgeom_exponents_at_1(acb_ode_exponents_t expos, acb_ode_group_struct * grps,
                       const gr_vec_t a, const gr_vec_t b, gr_ctx_t R,
                       gr_ctx_t CC)
{
    slong p = a->length, q = b->length;
    acb_ode_group_struct * grp = grps;
    int status = GR_SUCCESS;

    expos->groups = grps;
    expos->ngroups = 1;

    if (p < q + 1)
    {
        acb_zero(grp->leader);
        for (slong i = 0; i <= q; i++)
        {
            grp->shifts[i].n = i;
            grp->shifts[i].mult = 1;
        }
        grp->nshifts = q + 1;
        return GR_SUCCESS;
    }

    gr_ptr t;
    acb_t s;
    fmpz_t m;
    GR_TMP_INIT(t, R);
    acb_init(s);
    fmpz_init(m);

    status |= gr_zero(t, R);
    for (slong j = 0; j < q; j++)
        status |= gr_add(t, t, gr_vec_entry_srcptr(b, j, R), R);
    for (slong i = 0; i < p; i++)
        status |= gr_sub(t, t, gr_vec_entry_srcptr(a, i, R), R);
    status |= gr_set_other(s, t, R, CC);

    /* decide whether s is a (small) integer */
    int is_int = gr_get_fmpz(m, t, R);
    if (is_int == GR_SUCCESS)
        is_int = fmpz_fits_si(m) ? 1 : -1;
    else if (is_int == GR_DOMAIN || !acb_contains_int(s))
        is_int = 0;
    else
        is_int = -1;

    if (status != GR_SUCCESS || is_int < 0)
    {
        status |= GR_UNABLE;
    }
    else if (!is_int)
    {
        /* {0, ..., q-1}, {s} */
        if (q > 0)
        {
            acb_zero(grp->leader);
            for (slong i = 0; i < q; i++)
            {
                grp->shifts[i].n = i;
                grp->shifts[i].mult = 1;
            }
            grp->nshifts = q;
            grp++;
        }
        acb_set(grp->leader, s);
        grp->shifts[0].n = 0;
        grp->shifts[0].mult = 1;
        grp->nshifts = 1;
        expos->ngroups = grp - grps + 1;
    }
    else
    {
        /* single group {0, ..., q-1} ∪ {m} */
        slong mm = fmpz_get_si(m);
        slong lead = (q > 0) ? FLINT_MIN(0, mm) : mm;
        slong ns = 0;

        acb_set_si(grp->leader, lead);
        if (mm < 0)
        {
            grp->shifts[ns].n = 0;
            grp->shifts[ns].mult = 1;
            ns++;
        }
        for (slong i = 0; i < q; i++)
        {
            grp->shifts[ns].n = i - lead;
            grp->shifts[ns].mult = (i == mm) ? 2 : 1;
            ns++;
        }
        if (mm >= q)
        {
            grp->shifts[ns].n = mm - lead;
            grp->shifts[ns].mult = 1;
            ns++;
        }
        grp->nshifts = ns;
    }

    fmpz_clear(m);
    acb_clear(s);
    GR_TMP_CLEAR(t, R);

    return status;
}

/* Taylor coefficients of order < len of pFq(a; b; ·) at z, |z| < 1, from the
   series (the series-level functions do not include the k! factor, so we
   append a unit lower parameter like acb_hypgeom_pfq does). */
static void
hypgeom_pfq_jet_series(acb_ptr res, acb_srcptr a, slong p, acb_srcptr b,
                       slong q, const acb_t z, slong len, slong prec)
{
    acb_poly_struct * ab = _acb_poly_vec_init(p + q + 2);
    acb_poly_t jet;
    acb_poly_init(jet);

    for (slong i = 0; i < p; i++)
        acb_poly_set_acb(ab + i, a + i);
    for (slong j = 0; j < q; j++)
        acb_poly_set_acb(ab + p + j, b + j);
    acb_poly_one(ab + p + q);
    acb_poly_set_coeff_acb(ab + p + q + 1, 0, z);
    acb_poly_set_coeff_si(ab + p + q + 1, 1, 1);

    acb_hypgeom_pfq_series_direct(jet, ab, p, ab + p, q + 1, ab + p + q + 1,
                                  0, -1, len, prec);

    for (slong i = 0; i < len; i++)
        acb_poly_get_coeff_acb(res + i, jet, i);

    acb_poly_clear(jet);
    _acb_poly_vec_clear(ab, p + q + 2);
}

/* We match pFq with the solutions of the equation at the pivot 3/4. When
   p = q + 1, the local solutions at 1 are (principal branches of)
   x^λ (log x)^k times analytic functions of x = z - 1, cut along z ∈ (0, 1].
   Since principal values on the cut of log are limits from the upper
   half-plane, the matching yields a combination of them that agrees with pFq
   on the closed upper half-disk. The restriction of pFq to the lower half-disk,
   whose boundary values on [1, 2) are the desired values on the cut, is
   obtained by transporting the coefficient vector through the left half of the
   disk, i.e., by applying the counterclockwise local monodromy matrix. */
int
acb_ode_hypgeom_pfq(acb_ptr res, const gr_vec_t a, const gr_vec_t b,
                    gr_ctx_t R, const acb_t z, slong len, slong prec)
{
    slong p = a->length, q = b->length;
    slong order = q + 1;
    slong wp = prec + 32;

    int status = GR_SUCCESS;

    if (p > q + 1)
        return GR_UNABLE;

    if (len <= 0)
        return GR_SUCCESS;

    gr_ctx_t Pol, Dop, CC;
    gr_ctx_init_gr_poly(Pol, R);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ctx_init_complex_acb(CC, wp);

    gr_ore_poly_t dop;
    gr_ore_poly_init(dop, Dop);

    acb_ode_group_struct grps[2];
    acb_ode_group_init(grps + 0, order);
    acb_ode_group_init(grps + 1, 1);
    acb_ode_exponents_t expos;

    acb_t zp;
    acb_init(zp);
    acb_ptr pts = _acb_vec_init(2);
    acb_ptr jet = _acb_vec_init(order);
    acb_ptr lcrt = _acb_vec_init(order);
    acb_ptr ab = _acb_vec_init(p + q);

    acb_mat_struct mats[2];
    acb_mat_init(mats + 0, order, order);
    acb_mat_init(mats + 1, len, order);
    acb_mat_t B, cup, clow, val, mono;
    acb_mat_init(B, order, 1);
    acb_mat_init(cup, order, 1);
    acb_mat_init(clow, order, 1);
    acb_mat_init(val, len, 1);
    acb_mat_init(mono, order, order);

    status |= hypgeom_operator_at_1(dop, a, b, R, Dop);
    status |= hypgeom_exponents_at_1(expos, grps, a, b, R, CC);
    status |= _gr_vec_set_other(ab, a->entries, R, p, CC);
    status |= _gr_vec_set_other(ab + p, b->entries, R, q, CC);

    if (status != GR_SUCCESS)
        goto cleanup;

    /* fundamental matrices at the pivot and at the evaluation point (the
       leading coefficient of dop is ±x^ε (1 + x)^order, ε ∈ {0, 1}) */

    acb_set_d(zp, 0.75);
    acb_sub_ui(pts + 0, zp, 1, wp);
    acb_sub_ui(pts + 1, z, 1, wp);
    for (slong i = 0; i < order; i++)
        acb_set_si(lcrt + i, -1);

    status |= acb_ode_fundamental_matrix_vec(mats, dop, Dop, expos, lcrt,
                                             pts, 2, ACB_ODE_BASIS_ECHELON, 0,
                                             wp);
    if (status != GR_SUCCESS)
        goto cleanup;

    /* coordinates of pFq in the local basis (rows of the fundamental matrices
       are Taylor coefficients, like those of the series) */

    hypgeom_pfq_jet_series(jet, ab, p, ab + p, q, zp, order, wp);
    for (slong i = 0; i < order; i++)
        acb_swap(acb_mat_entry(B, i, 0), jet + i);
    if (!acb_mat_solve(cup, mats + 0, B, wp))
    {
        _acb_vec_indeterminate(res, len);
        goto cleanup;
    }

    /* value at the evaluation point, on the branch(es) determined by its
       position with respect to the real axis; for input balls straddling the
       cut, we return the union of the values of the two branches */

    int need_up = 1, need_low = 0;
    if (p == q + 1)
    {
        const arb_struct * re = acb_realref(pts + 1);
        const arb_struct * im = acb_imagref(pts + 1);
        if (arb_is_negative(im) || (arb_is_zero(im) && arb_is_nonnegative(re)))
            need_up = 0;
        if (!(arb_is_positive(im) || (arb_is_zero(im) && arb_is_negative(re))))
            need_low = 1;
    }

    if (need_up)
    {
        acb_mat_mul(val, mats + 1, cup, wp);
        for (slong i = 0; i < len; i++)
            acb_set(res + i, acb_mat_entry(val, i, 0));
    }

    if (need_low)
    {
        status |= acb_ode_local_monodromy(mono, dop, Dop, expos, 0, wp);
        if (status != GR_SUCCESS)
            goto cleanup;

        acb_mat_mul(clow, mono, cup, wp);
        acb_mat_mul(val, mats + 1, clow, wp);
        for (slong i = 0; i < len; i++)
        {
            if (need_up)
                acb_union(res + i, res + i, acb_mat_entry(val, i, 0), wp);
            else
                acb_set(res + i, acb_mat_entry(val, i, 0));
        }
    }

cleanup:
    acb_mat_clear(mono);
    acb_mat_clear(val);
    acb_mat_clear(clow);
    acb_mat_clear(cup);
    acb_mat_clear(B);
    acb_mat_clear(mats + 1);
    acb_mat_clear(mats + 0);
    _acb_vec_clear(ab, p + q);
    _acb_vec_clear(lcrt, order);
    _acb_vec_clear(jet, order);
    _acb_vec_clear(pts, 2);
    acb_clear(zp);
    acb_ode_group_clear(grps + 1);
    acb_ode_group_clear(grps + 0);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(CC);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);

    return status;
}
