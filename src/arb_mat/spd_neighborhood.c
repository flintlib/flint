/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

void
arb_mat_spd_neighborhood(arf_t r, const arb_mat_t mat, slong prec)
{
    slong g = acb_mat_nrows(mat);
    arb_t norm;
    arf_t b;
    fmpz_t e;
    arb_mat_t test;
    slong j, k;
    int valid;

    arb_init(norm);
    arf_init(b);
    fmpz_init(e);
    arb_mat_init(test, g, g);

    /* Take a guess at r */
    arb_mat_spd_eig_lbound_arf(b, mat, prec);
    if (arf_cmp_si(b, 0) <= 0)
    {
        fmpz_set_si(e, -prec); /* See stopping condition below */
    }
    else
    {
        arb_mat_max_norm(norm, mat, prec);
        arb_mul_si(norm, norm, g, prec);
        arb_pow_ui(norm, norm, g, prec);
        arb_div_arf(norm, norm, b, prec);
        arb_inv(norm, norm, prec);
        arb_get_lbound_arf(r, norm, prec);
        arf_frexp(r, e, r);
    }
    arf_one(r);
    arf_mul_2exp_fmpz(r, r, e);

    /* Get associated b */
    arb_mat_set(test, mat);
    arb_mat_add_error_arf(test, r);
    arb_mat_spd_eig_lbound_arf(b, test, prec);

    /* Increase r until it gets invalid */
    while (arf_cmp_si(b, 0) > 0)
    {
        arf_mul_2exp_si(r, r, 1);
        fmpz_add_si(e, e, 1);
        arb_mat_set(test, mat);
        arb_mat_add_error_arf(test, r);
        arb_mat_spd_eig_lbound_arf(b, test, prec);
    }

    /* Then reduce r until valid, or we reach prec */
    while (arf_cmp_si(b, 0) <= 0 && fmpz_cmp_si(e, -prec) > 0)
    {
        arf_mul_2exp_si(r, r, -1);
        fmpz_add_si(e, e, -1);
        arb_mat_set(test, mat);
        arb_mat_add_error_arf(test, r);
        arb_mat_spd_eig_lbound_arf(b, test, prec);
    }

    if (arf_cmp_si(b, 0) <= 0) /* Could not find a valid radius */
    {
        arf_zero(r);
    }

    arb_clear(norm);
    arf_clear(b);
    fmpz_clear(e);
    arb_mat_clear(test);
}
