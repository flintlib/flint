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
#include "acb_mat.h"
#include "acb_ode.h"
#include "gr.h"
#include "gr_ore_poly.h"
#include "gr_poly.h"

/* dop in Cst[x][d/dx]; a regular, b ordinary; lcrt including zero */
static int
_acb_ode_transition_matrix_step(acb_mat_t mat,
                                const gr_poly_struct * dop, slong dop_len,
                                gr_ctx_t Cst, const acb_ptr lcrt,
                                gr_srcptr a, gr_ctx_t PtA,
                                gr_srcptr b, gr_ctx_t PtB,
                                const acb_ode_exponents_t expos,
                                ulong flags, slong prec)
{
    gr_ctx_t CC, Pol, Dop;
    gr_ore_poly_t shifted;
    acb_t dx;
    slong mult0 = 0;
    int status = GR_SUCCESS;

    FLINT_ASSERT(dop_len > 0);

    if (PtA != Cst)
    {
        /* preferentially coerce into the operator's constant ring, even when
           gr_ctx_cmp_coercion would select the other way */

        {
            gr_ptr a1;
            GR_TMP_INIT(a1, Cst);

            int cvt = gr_set_other(a1, a, PtA, Cst);

            if (cvt == GR_SUCCESS)
                status |= _acb_ode_transition_matrix_step(
                        mat, dop, dop_len, Cst, lcrt, a1, Cst, b, PtB, expos,
                        flags, prec);

            GR_TMP_CLEAR(a1, Cst);

            if (cvt == GR_SUCCESS)
                return status;
        }

        {
            gr_ctx_t Pol1;
            gr_poly_struct * dop1;

            gr_ctx_init_gr_poly(Pol, Cst);
            gr_ctx_init_gr_poly(Pol1, PtA);
            GR_TMP_INIT_VEC(dop1, dop_len, Pol1);

            int cvt = _gr_vec_set_other(dop1, dop, Pol, dop_len, Pol1);
            /* catch conversions to the zero ring and the like */
            if (dop1[dop_len - 1].length != dop[dop_len-1].length)
                cvt = GR_UNABLE;

            if (cvt == GR_SUCCESS)
                status |= _acb_ode_transition_matrix_step(
                        mat, dop1, dop_len, PtA, lcrt, a, PtA, b, PtB, expos,
                        flags, prec);

            GR_TMP_CLEAR_VEC(dop1, dop_len, Pol1);
            gr_ctx_clear(Pol1);
            gr_ctx_clear(Pol);

            if (cvt == GR_SUCCESS)
                return status;
        }

        return GR_UNABLE;
    }

    gr_ctx_init_complex_acb(CC, prec);
    gr_ctx_init_gr_poly(Pol, Cst);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);

    gr_ore_poly_init(shifted, Dop);
    acb_init(dx);
    slong lcdeg = (dop + dop_len - 1)->length - 1;
    acb_ptr llcrt = _acb_vec_init(lcdeg);

    gr_ore_poly_fit_length(shifted, dop_len, Dop);
    _gr_ore_poly_set_length(shifted, dop_len, Dop);

    for (slong i = 0; i < dop_len; i++)
        status |= gr_poly_taylor_shift(GR_ENTRY(shifted->coeffs, i, Pol->sizeof_elem),
                                       GR_ENTRY(dop, i, Pol->sizeof_elem),
                                       a, Cst);

    if (PtB == Cst)
    {
        gr_ptr _dx;
        GR_TMP_INIT(_dx, Cst);

        status |= gr_sub(_dx, b, a, Cst);
        status |= gr_set_other(dx, _dx, Cst, CC);

        GR_TMP_CLEAR(_dx, Cst);
    }
    else
    {
        gr_ptr _a, _b;
        GR_TMP_INIT2(_a, _b, CC);

        status |= gr_set_other(_a, a, PtA, CC);
        status |= gr_set_other(_b, b, PtB, CC);
        status |= gr_sub(dx, _b, _a, CC);

        GR_TMP_CLEAR2(_a, _b, CC);
    }

    status |= _gr_vec_sub_scalar_other(llcrt, lcrt, lcdeg, a, Cst, CC);

    if (!(flags & ACB_ODE_ORDINARY))
    {
        gr_poly_struct * slc =
            ((gr_poly_struct *) shifted->coeffs) + dop_len - 1;

        slong val = -1;
        for (int is_zero = T_TRUE; is_zero == T_TRUE; )
        {
            is_zero = gr_is_zero(gr_poly_coeff_ptr(slc, ++val, Cst), Cst);
            if (is_zero == T_UNKNOWN)
                status = GR_UNABLE;
        }

        for (slong i = 0; i < lcdeg; i++)
            if (acb_contains_zero(llcrt + i))
                acb_swap(llcrt + i, llcrt + mult0++);

        FLINT_ASSERT(mult0 >= val);
        if (mult0 != val)
            status = GR_UNABLE;
    }

    if (status == GR_SUCCESS)
        status |= acb_ode_fundamental_matrix(mat, shifted, Dop, expos,
                                             llcrt + mult0, dx,
                                             ACB_ODE_BASIS_ECHELON,
                                             flags, prec);

    _acb_vec_clear(llcrt, lcdeg);
    acb_clear(dx);
    gr_ore_poly_clear(shifted, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(CC);

    return status;
}

int
_acb_ode_transition_matrix_fixed(acb_mat_t mat,
                                 const gr_poly_struct * dop, slong dop_len,
                                 gr_ctx_t Cst, const acb_ptr lcrt,
                                 gr_srcptr path, slong path_len, gr_ctx_t Pts,
                                 const acb_ode_exponents_t expos0,
                                 const acb_ode_exponents_t expos1,
                                 ulong flags, slong prec)
{
    gr_ptr lcev;
    int status = GR_SUCCESS;

    FLINT_ASSERT(path_len >= 2);
    /* flags &= ~(ACB_ODE_ORDINARY | ACB_ODE_REGULAR); */
    flags = 0;

    GR_TMP_INIT(lcev, Pts);

    gr_srcptr end = GR_ENTRY(path, path_len - 1, Pts->sizeof_elem);
    status |= gr_poly_evaluate_other(lcev, dop + dop_len - 1, end, Pts, Cst);
    truth_t sing = gr_is_zero(lcev, Pts);

    if (path_len == 2)
    {
        gr_srcptr a = GR_ENTRY(path, 0, Pts->sizeof_elem);
        gr_srcptr b = GR_ENTRY(path, 1, Pts->sizeof_elem);

        if (sing == T_TRUE)
        {
            status |= _acb_ode_transition_matrix_step(mat, dop, dop_len, Cst,
                                                      lcrt, b, Pts, a, Pts,
                                                      expos1, flags, prec);
            if (!acb_mat_inv(mat, mat, prec))
                acb_mat_indeterminate(mat);
        }
        else
        {
            status |= _acb_ode_transition_matrix_step(mat, dop, dop_len, Cst,
                                                      lcrt, a, Pts, b, Pts,
                                                      expos0, flags, prec);
        }
    }
    else
    {
        acb_mat_t mat0, mat1;

        acb_mat_init(mat0, dop_len - 1, dop_len - 1);
        acb_mat_init(mat1, dop_len - 1, dop_len - 1);

        /* In principle we might want to solve with a matrix mat1 corresponding
           to analytic continuation several steps away from the singular point
           (as it may be better conditioned), but with a fixed path, this would
           make the constraints on step lengths too involved. */
        if (sing == T_TRUE)
        {
            gr_srcptr prev = GR_ENTRY(path, path_len - 2, Pts->sizeof_elem);

            status |= _acb_ode_transition_matrix_fixed(mat0,
                                                       dop, dop_len, Cst, lcrt,
                                                       path, path_len - 1, Pts,
                                                       expos0, NULL,
                                                       flags, prec);
            status |= _acb_ode_transition_matrix_step(mat1, dop, dop_len, Cst,
                                                      lcrt, end, Pts, prev, Pts,
                                                      expos1, flags, prec);
            if (!acb_mat_solve(mat, mat1, mat0, prec))
                acb_mat_indeterminate(mat);
        }
        else
        {
            slong m = path_len / 2;

            gr_srcptr path1 = GR_ENTRY(path, m, Pts->sizeof_elem);
            status |= _acb_ode_transition_matrix_fixed(mat0,
                                                       dop, dop_len, Cst, lcrt,
                                                       path, m + 1, Pts,
                                                       expos0, NULL,
                                                       flags, prec);
            status |= _acb_ode_transition_matrix_fixed(mat1,
                                                       dop, dop_len, Cst, lcrt,
                                                       path1, path_len - m, Pts,
                                                       NULL, NULL,
                                                       flags, prec);

            acb_mat_mul(mat, mat1, mat0, prec);
        }

        acb_mat_clear(mat0);
        acb_mat_clear(mat1);
    }

    GR_TMP_CLEAR(lcev, Pts);

    return status;
}

int
acb_ode_transition_matrix_fixed(acb_mat_t mat,
                                const gr_ore_poly_t dop, gr_ore_poly_ctx_t Dop,
                                gr_vec_t path, gr_ctx_t Pts,
                                const acb_ode_exponents_t expos0,
                                const acb_ode_exponents_t expos1,
                                ulong flags, slong prec)
{
    gr_ctx_struct * Cst, * Pol;
    gr_ctx_t CC;
    acb_struct * lcrt;

    int status = GR_SUCCESS;

    status |= gr_ore_poly_ctx_over_gr_poly_base_ptrs(&Cst, &Pol, Dop);
    if (status != GR_SUCCESS)
        return status;

    if (GR_ORE_POLY_CTX(Dop)->which_algebra != ORE_ALGEBRA_DERIVATIVE)
    {
        gr_ctx_t Dop1;
        gr_ore_poly_t dop1;
        slong p;

        gr_ore_poly_ctx_init(Dop1, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
        gr_ore_poly_init(dop1, Dop1);

        status |= gr_ore_poly_convert(dop1, &p, dop, Dop1, Dop);
        status |= acb_ode_transition_matrix_fixed(mat, dop1, Dop1, path, Pts,
                                                  expos0, expos1, flags, prec);

        gr_ore_poly_clear(dop1, Dop1);
        gr_ctx_clear(Dop1);

        return status;
    }

    if (dop->length == 0)
        return GR_DOMAIN;

    if (dop->length - 1 > mat->r || dop->length - 1 > mat->c)
        return GR_UNABLE;  /* support computing only part of the matrix??? */

    if (path->length == 0)
        return GR_DOMAIN;
    else if (path->length == 1)
    {
        acb_mat_one(mat);
        return GR_SUCCESS;
    }

    gr_poly_struct * lc = GR_ENTRY(dop->coeffs, dop->length - 1,
                                   sizeof(gr_poly_t));
    slong lcdeg = lc->length - 1;

    gr_ctx_init_complex_acb(CC, prec);
    lcrt = _acb_vec_init(lcdeg);

    status |= acb_ode_lc_roots(lcrt, dop, Dop, 1, CC);

    if (status == GR_SUCCESS)
        status |= _acb_ode_transition_matrix_fixed(
                mat, dop->coeffs, dop->length, Cst, lcrt, path->entries,
                path->length, Pts, expos0, expos1, flags, prec);

    _acb_vec_clear(lcrt, lcdeg);
    gr_ctx_clear(CC);

    return status;
}
