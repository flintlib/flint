/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "perm.h"
#include "arb_mat.h"

static void
_apply_permutation(arb_mat_t A, slong * P, slong n)
{
    arb_ptr * Atmp;
    slong i;
    Atmp = flint_malloc(sizeof(arb_ptr) * n);
    for (i = 0; i < n; i++) Atmp[i] = A->rows[P[i]];
    for (i = 0; i < n; i++) A->rows[i] = Atmp[i];
    flint_free(Atmp);
}

/* Enclosure of det(I + eps) using Gershgorin circles.
   Can be improved. */
void
arb_mat_det_one_gershgorin(arb_t det, const arb_mat_t A)
{
    slong n, i, j;
    arb_t t;
    mag_t r, e, f;

    n = arb_mat_nrows(A);

    arb_init(t);
    mag_init(r);
    mag_init(e);
    mag_init(f);

    for (i = 0; i < n; i++)
    {
        mag_zero(e);

        for (j = 0; j < n; j++)
        {
            if (i == j)
            {
                arb_sub_ui(t, arb_mat_entry(A, i, j), 1, MAG_BITS);
                arb_get_mag(f, t);
            }
            else
            {
                arb_get_mag(f, arb_mat_entry(A, i, j));
            }

            mag_add(e, e, f);
        }

        mag_max(r, r, e);
    }

    /* (1 + eps)^n - 1 <= expm1(n*eps) */
    mag_mul_ui(r, r, n);
    mag_expm1(r, r);

    arf_one(arb_midref(det));
    mag_set(arb_radref(det), r);

    arb_clear(t);
    mag_clear(r);
    mag_clear(e);
    mag_clear(f);
}

void
arb_mat_det_precond(arb_t det, const arb_mat_t A, slong prec)
{
    arb_mat_t LU, Linv, Uinv;
    arb_t detU;
    slong n;
    slong *P;

    n = arb_mat_nrows(A);

    if (n == 0)
    {
        arb_one(det);
        return;
    }

    P = _perm_init(n);

    arb_mat_init(LU, n, n);

    if (!arb_mat_approx_lu(P, LU, A, prec))
    {
        /* Fallback. */
        arb_mat_det_lu(det, A, prec);
    }
    else
    {
        arb_mat_init(Linv, n, n);
        arb_mat_init(Uinv, n, n);
        arb_init(detU);

        arb_mat_one(Linv);
        arb_mat_approx_solve_tril(Linv, LU, Linv, 1, prec);
        arb_mat_one(Uinv);
        arb_mat_approx_solve_triu(Uinv, LU, Uinv, 0, prec);

        arb_mat_diag_prod(detU, Uinv, prec);

        arb_mat_mul(LU, A, Uinv, prec);
        _apply_permutation(LU, P, n);
        arb_mat_mul(Uinv, Linv, LU, prec);

        arb_mat_det_one_gershgorin(det, Uinv);
        if (_perm_parity(P, n))
            arb_neg(det, det);

        arb_div(det, det, detU, prec);

        if (arb_contains_zero(det))
        {
            /* Run the interval LU algorithm. This can give a much better
               bound if the Gaussian elimination manages to work through
               several rows, and it is not that expensive. */
            arb_mat_det_lu(detU, A, prec);
            if (mag_cmp(arb_radref(detU), arb_radref(det)) < 0)
                arb_set(det, detU);
        }

        arb_mat_clear(Linv);
        arb_mat_clear(Uinv);
        arb_clear(detU);
    }

    _perm_clear(P);
    arb_mat_clear(LU);
}
