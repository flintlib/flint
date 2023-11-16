/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

void
acb_mat_randtest_eig(acb_mat_t A, flint_rand_t state, acb_srcptr E, slong prec)
{
    slong n, i, j, ebits;
    acb_mat_t U, Q;

    n = acb_mat_nrows(A);
    ebits = 1 + n_randint(state, 5);

    acb_mat_init(U, n, n);
    acb_mat_init(Q, n, n);

    /* Skew-Hermitian matrix */
    acb_mat_randtest(Q, state, prec, 1);
    if (n_randint(state, 2))
        acb_mat_get_mid(Q, Q);
    for (i = 0; i < n; i++)
    {
        for (j = i + 1; j < n; j++)
        {
            acb_neg(acb_mat_entry(Q, i, j), acb_mat_entry(Q, j, i));
            acb_conj(acb_mat_entry(Q, i, j), acb_mat_entry(Q, i, j));
        }
        arb_zero(acb_realref(acb_mat_entry(Q, i, i)));
    }

    acb_mat_exp(Q, Q, prec);

    acb_mat_randtest(U, state, prec, ebits);
    if (n_randint(state, 2))
        acb_mat_get_mid(U, U);
    for (i = 0; i < n; i++)
        for (j = 0; j < i; j++)
            acb_zero(acb_mat_entry(U, i, j));

    for (i = 0; i < n; i++)
        acb_set(acb_mat_entry(U, i, i), E + i);

    acb_mat_mul(U, Q, U, prec);
    acb_mat_conjugate_transpose(Q, Q);
    acb_mat_mul(A, U, Q, prec);

    acb_mat_clear(U);
    acb_mat_clear(Q);
}
