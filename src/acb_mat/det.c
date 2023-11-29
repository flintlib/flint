/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

static void
_acb_mat_det_cofactor_2x2(acb_t t, const acb_mat_t A, slong prec)
{
    acb_mul   (t, acb_mat_entry(A, 0, 0), acb_mat_entry(A, 1, 1), prec);
    acb_submul(t, acb_mat_entry(A, 0, 1), acb_mat_entry(A, 1, 0), prec);
}

static void
_acb_mat_det_cofactor_3x3(acb_t t, const acb_mat_t A, slong prec)
{
    acb_t a;
    acb_init(a);

    acb_mul   (a, acb_mat_entry(A, 1, 0), acb_mat_entry(A, 2, 1), prec);
    acb_submul(a, acb_mat_entry(A, 1, 1), acb_mat_entry(A, 2, 0), prec);
    acb_mul   (t, a, acb_mat_entry(A, 0, 2), prec);

    acb_mul   (a, acb_mat_entry(A, 1, 2), acb_mat_entry(A, 2, 0), prec);
    acb_submul(a, acb_mat_entry(A, 1, 0), acb_mat_entry(A, 2, 2), prec);
    acb_addmul(t, a, acb_mat_entry(A, 0, 1), prec);

    acb_mul   (a, acb_mat_entry(A, 1, 1), acb_mat_entry(A, 2, 2), prec);
    acb_submul(a, acb_mat_entry(A, 1, 2), acb_mat_entry(A, 2, 1), prec);
    acb_addmul(t, a, acb_mat_entry(A, 0, 0), prec);

    acb_clear(a);
}

void
acb_mat_det(acb_t det, const acb_mat_t A, slong prec)
{
    slong n;

    if (!acb_mat_is_square(A))
    {
        flint_throw(FLINT_ERROR, "acb_mat_det: a square matrix is required!\n");
    }

    n = acb_mat_nrows(A);

    if (n == 0)
    {
        acb_one(det);
    }
    else if (n == 1)
    {
        acb_set_round(det, acb_mat_entry(A, 0, 0), prec);
    }
    else if (n == 2)
    {
        _acb_mat_det_cofactor_2x2(det, A, prec);
    }
    else if (!acb_mat_is_finite(A))
    {
        acb_indeterminate(det);
    }
    else if (acb_mat_is_tril(A) || acb_mat_is_triu(A))
    {
        acb_mat_diag_prod(det, A, prec);
    }
    else if (n == 3)
    {
        _acb_mat_det_cofactor_3x3(det, A, prec);
        /* note: 4x4 performs worse than LU */
    }
    else
    {
        if (n <= 14 || prec > 10.0 * n)
            acb_mat_det_lu(det, A, prec);
        else
            acb_mat_det_precond(det, A, prec);
    }
}
