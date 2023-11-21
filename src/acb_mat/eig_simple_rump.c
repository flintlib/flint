/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

int
acb_mat_eig_simple_rump(acb_ptr E, acb_mat_t L, acb_mat_t R,
    const acb_mat_t A, acb_srcptr E_approx, const acb_mat_t R_approx, slong prec)
{
    slong i, j, n;
    acb_mat_t X, R2;
    int result;

    n = acb_mat_nrows(A);

    if (n == 0)
        return 1;

    if (n == 1)
    {
        acb_set_round(E, acb_mat_entry(A, 0, 0), prec);
        if (L != NULL)
            acb_one(acb_mat_entry(L, 0, 0));
        if (R != NULL)
            acb_one(acb_mat_entry(R, 0, 0));
        return 1;
    }

    acb_mat_init(X, n, 1);
    acb_mat_init(R2, n, n);

    result = 1;

    for (i = 0; i < n && result; i++)
    {
        for (j = 0; j < n; j++)
            acb_set(acb_mat_entry(X, j, 0), acb_mat_entry(R_approx, j, i));

        acb_mat_eig_enclosure_rump(E + i, NULL, X, A, E_approx + i, X, prec);

        if (!acb_is_finite(E + i))
            result = 0;

        for (j = 0; j < i; j++)
            if (acb_overlaps(E + i, E + j))
                result = 0;

        for (j = 0; j < n; j++)
            acb_set(acb_mat_entry(R2, j, i), acb_mat_entry(X, j, 0));
    }

    if (R != NULL)
    {
        if (result)
            acb_mat_set(R, R2);
        else
            acb_mat_indeterminate(R);
    }

    if (L != NULL)
    {
        if (!result || !acb_mat_inv(L, R2, prec))
            acb_mat_indeterminate(L);
    }

    if (!result)
        _acb_vec_indeterminate(E, n);

    acb_mat_clear(X);
    acb_mat_clear(R2);

    return result;
}
