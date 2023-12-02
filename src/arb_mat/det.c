/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

static void
_arb_mat_det_cofactor_2x2(arb_t t, const arb_mat_t A, slong prec)
{
    arb_mul   (t, arb_mat_entry(A, 0, 0), arb_mat_entry(A, 1, 1), prec);
    arb_submul(t, arb_mat_entry(A, 0, 1), arb_mat_entry(A, 1, 0), prec);
}

static void
_arb_mat_det_cofactor_3x3(arb_t t, const arb_mat_t A, slong prec)
{
    arb_t a;
    arb_init(a);

    arb_mul   (a, arb_mat_entry(A, 1, 0), arb_mat_entry(A, 2, 1), prec);
    arb_submul(a, arb_mat_entry(A, 1, 1), arb_mat_entry(A, 2, 0), prec);
    arb_mul   (t, a, arb_mat_entry(A, 0, 2), prec);

    arb_mul   (a, arb_mat_entry(A, 1, 2), arb_mat_entry(A, 2, 0), prec);
    arb_submul(a, arb_mat_entry(A, 1, 0), arb_mat_entry(A, 2, 2), prec);
    arb_addmul(t, a, arb_mat_entry(A, 0, 1), prec);

    arb_mul   (a, arb_mat_entry(A, 1, 1), arb_mat_entry(A, 2, 2), prec);
    arb_submul(a, arb_mat_entry(A, 1, 2), arb_mat_entry(A, 2, 1), prec);
    arb_addmul(t, a, arb_mat_entry(A, 0, 0), prec);

    arb_clear(a);
}

void
arb_mat_det(arb_t det, const arb_mat_t A, slong prec)
{
    slong n;

    if (!arb_mat_is_square(A))
    {
        flint_throw(FLINT_ERROR, "arb_mat_det: a square matrix is required!\n");
    }

    n = arb_mat_nrows(A);

    if (n == 0)
    {
        arb_one(det);
    }
    else if (n == 1)
    {
        arb_set_round(det, arb_mat_entry(A, 0, 0), prec);
    }
    else if (n == 2)
    {
        _arb_mat_det_cofactor_2x2(det, A, prec);
    }
    else if (!arb_mat_is_finite(A))
    {
        arb_indeterminate(det);
    }
    else if (arb_mat_is_tril(A) || arb_mat_is_triu(A))
    {
        arb_mat_diag_prod(det, A, prec);
    }
    else if (n == 3)
    {
        _arb_mat_det_cofactor_3x3(det, A, prec);
        /* note: 4x4 performs worse than LU */
    }
    else
    {
        if (n <= 10 || prec > 10.0 * n)
            arb_mat_det_lu(det, A, prec);
        else
            arb_mat_det_precond(det, A, prec);
    }
}

