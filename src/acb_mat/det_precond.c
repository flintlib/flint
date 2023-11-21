/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "perm.h"
#include "acb_mat.h"

static void
_apply_permutation(acb_mat_t A, slong * P, slong n)
{
    acb_ptr * Atmp;
    slong i;
    Atmp = flint_malloc(sizeof(acb_ptr) * n);
    for (i = 0; i < n; i++) Atmp[i] = A->rows[P[i]];
    for (i = 0; i < n; i++) A->rows[i] = Atmp[i];
    flint_free(Atmp);
}

/* Enclosure of det(I + eps) using Gershgorin circles.
   Can be improved. */
void
acb_mat_det_one_gershgorin(acb_t det, const acb_mat_t A)
{
    slong n, i, j;
    acb_t t;
    mag_t r, e, f;

    n = acb_mat_nrows(A);

    acb_init(t);
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
                acb_sub_ui(t, acb_mat_entry(A, i, j), 1, MAG_BITS);
                acb_get_mag(f, t);
            }
            else
            {
                acb_get_mag(f, acb_mat_entry(A, i, j));
            }

            mag_add(e, e, f);
        }

        mag_max(r, r, e);
    }

    /* (1 + eps)^n - 1 <= expm1(n*eps) */
    mag_mul_ui(r, r, n);
    mag_expm1(r, r);

    acb_one(det);
    mag_set(arb_radref(acb_realref(det)), r);
    mag_set(arb_radref(acb_imagref(det)), r);

    acb_clear(t);
    mag_clear(r);
    mag_clear(e);
    mag_clear(f);
}

void
acb_mat_det_precond(acb_t det, const acb_mat_t A, slong prec)
{
    acb_mat_t LU, Linv, Uinv;
    acb_t detU;
    slong n;
    slong *P;

    n = acb_mat_nrows(A);

    if (n == 0)
    {
        acb_one(det);
        return;
    }

    P = _perm_init(n);

    acb_mat_init(LU, n, n);

    if (!acb_mat_approx_lu(P, LU, A, prec))
    {
        /* Fallback. */
        acb_mat_det_lu(det, A, prec);
    }
    else
    {
        acb_mat_init(Linv, n, n);
        acb_mat_init(Uinv, n, n);
        acb_init(detU);

        acb_mat_one(Linv);
        acb_mat_approx_solve_tril(Linv, LU, Linv, 1, prec);
        acb_mat_one(Uinv);
        acb_mat_approx_solve_triu(Uinv, LU, Uinv, 0, prec);

        acb_mat_diag_prod(detU, Uinv, prec);

        acb_mat_mul(LU, A, Uinv, prec);
        _apply_permutation(LU, P, n);
        acb_mat_mul(Uinv, Linv, LU, prec);

        acb_mat_det_one_gershgorin(det, Uinv);
        if (acb_mat_is_real(A))
            arb_zero(acb_imagref(det));

        if (_perm_parity(P, n))
            acb_neg(det, det);

        acb_div(det, det, detU, prec);

        if (acb_contains_zero(det))
        {
            mag_t rad1, rad2;

            /* Run the interval LU algorithm. This can give a much better
               bound if the Gaussian elimination manages to work through
               several rows, and it is not that expensive. */
            acb_mat_det_lu(detU, A, prec);

            mag_init(rad1);
            mag_init(rad2);

            mag_hypot(rad1, arb_radref(acb_realref(detU)), arb_radref(acb_imagref(detU)));
            mag_hypot(rad2, arb_radref(acb_realref(det)), arb_radref(acb_imagref(det)));

            if (mag_cmp(rad1, rad2) < 0)
                acb_set(det, detU);

            mag_clear(rad1);
            mag_clear(rad2);
        }

        acb_mat_clear(Linv);
        acb_mat_clear(Uinv);
        acb_clear(detU);
    }

    _perm_clear(P);
    acb_mat_clear(LU);
}
