/*
    Copyright (C) 2026 Marc Mezzarobba

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

// dop in d/dx
// a ordinary
int
_acb_ode_transition_matrix_step(
        acb_mat_t mat,
        const gr_poly_struct * dop, slong dop_len,
        gr_ptr a, gr_ptr b, gr_ctx_t Cst,
        slong prec)
{
    int status = GR_SUCCESS;

    ulong flags = 0; // ?

    gr_ptr delta, shifted;
    acb_t _delta;

    gr_ctx_init_complex_acb(CC, prec);
    gr_ctx_init_gr_poly(Pol, Cst);
    gr_ore_poly_ctx_init(Euler, Pol, 0, ORE_ALGEBRA_EULER_DERIVATIVE);

    GR_TMP_INIT_VEC(shifted, dop_len, Pol);
    GR_TMP_INIT(delta, Cst);
    acb_init(_delta);

    for (slong i = 0; i < dop_len; i++)
        status |= gr_poly_taylor_shift(GR_ENTRY(shifted, i, Pol->sizeof_elem),
                                       GR_ENTRY(dop, i, Pol->sizeof_elem),
                                       a, Cst);

    status |= acb_ode_local_diffop(ldop, &flags, shifted, Dop, Euler);


    status |= gr_sub(delta, a, b, Cst);
    status |= gr_set_other(_delta, delta, Cst, CC);

    for (s

    // XXX passer l'op global ?
    _acb_ode_fundamental_matrix_vec(mat, ldop, dop_len, expos, lcrt, bound,
                                    _delta, 1, ACB_ODE_BASIS_ECHELON,
                                    acb_ode_sum_divconquer, flags, prec);

    GR_TMP_CLEAR(delta, Cst);
    acb_clear(_delta);
    gr_ctx_clear(CC);
    gr_ctx_clear(Pol);
    gr_ctx_clear(Euler);
}

int
_acb_ode_transition_matrix_fixed(
        acb_mat_t mat,
        const gr_poly_struct * gdop, slong dop_len,
        gr_srcptr path, slong path_len, gr_ctx_t Cst,
        slong prec)
{
    int status = GR_SUCCESS;

    FLINT_ASSERT(path_len >= 2);

    if (path_len == 2)
    {

    }
    else
    {
        acb_mat_t mat0, mat1;

        acb_mat_init(mat0, mat->r, dop_len - 1);
        acb_mat_init(mat1, dop_len - 1, dop_len - 1);

        slong m = path_len / 2;

        _acb_ode_transition_matrix_fixed(mat0, gdop, dop_len,
                                         path, m + 1, prec);
        _acb_ode_transition_matrix_fixed(mat1, gdop, dop_len,
                                         path + m, path_len - m, prec);

        // XXX comment on résout dans le cas singulier ? tout en haut de l'arbre
        // de récursion ? cas particulier pour un chemin de longueur 1 ?
        // fonction séparée ?
        acb_mat_mul(mat, mat0, mat1, prec);

        acb_mat_clear(mat0);
        acb_mat_clear(mat1);
    }

    return status;
}

int
acb_ode_transition_matrix_fixed(
        acb_mat_t mat,
        const gr_ore_poly_t dop, gr_ore_poly_ctx_t Dop,
        const acb_ode_exponents_t expos,  // ?
        acb_srcptr lcrt,  // ?
        gr_vec_t path, gr_ctx_t Pts,
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

    if (dop->length < 1)
        return GR_DOMAIN;

    slong order = dop->length - 1;
    gr_poly_struct * lc = (gr_poly_struct *) (dop->coeffs) + order;

    if (lcrt == NULL)  // XXX make it mandatory?
    {
        gr_ctx_t bCC;
        gr_ctx_init_complex_acb(bCC, MAG_BITS);

        _lcrt = _acb_vec_init(lc->length - 1);
        status |= acb_ode_lc_roots(_lcrt, dop, Dop, 0, bCC);
        DEBUG_STATUS;
        lcrt = _lcrt;

        gr_ctx_clear(bCC);
    }

    if (GR_ORE_POLY_CTX(Dop)->which_algebra == ORE_ALGEBRA_DERIVATIVE
        && gr_is_zero(lc->coeffs, Cst) == T_FALSE)
        flags |= ACB_ODE_ORDINARY;

    acb_ode_local_diffop_ctx_init(Pol1, Euler, Dop, Cst);
    GR_TMP_INIT(ldop, Euler);

    status |= acb_ode_local_diffop(ldop, dop, Dop, Euler);
    DEBUG_STATUS;

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
            DEBUG_STATUS;

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
    if (gr_ctx_is_exact(Cst))
    {
        gr_ctx_t bCC, bCCx, bCCxT;
        gr_ctx_init_complex_acb(bCC, MAG_BITS);
        acb_ode_local_diffop_ctx_init(bCCx, bCCxT, Dop, bCC);
        GR_TMP_INIT(acb_ldop, bCCxT);

        status |= gr_poly_set_gr_poly_other((gr_ptr) acb_ldop, (gr_ptr) ldop, Pol1, bCCx);
        DEBUG_STATUS;
        wp = choose_wp(acb_ldop, lcrt, pts, npts, prec, bCCxT);

        GR_TMP_CLEAR(acb_ldop, bCCxT);
        gr_ctx_clear(bCCxT);
        gr_ctx_clear(bCCx);
        gr_ctx_clear(bCC);
    }

    gr_ctx_init_complex_acb(wCC, wp);
    acb_ode_local_diffop_ctx_init(wCCx, wCCxT, Dop, wCC);
    GR_TMP_INIT(acb_ldop, wCCxT);

    status |= gr_poly_set_gr_poly_other((gr_ptr) acb_ldop, (gr_ptr) ldop, Pol1, wCCx);
    DEBUG_STATUS;

    /* Check regularity */
    gr_poly_struct * llc = (gr_poly_struct *) (acb_ldop->coeffs) + order;
    if (status == GR_SUCCESS && gr_is_zero(llc->coeffs, wCC) == T_TRUE)
        status = GR_UNABLE;
    DEBUG_STATUS;

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
