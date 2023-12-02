/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

static truth_t
ca_check_is_zero_fast(const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_QQ(x, ctx))
    {
        return fmpq_is_zero(CA_FMPQ(x)) ? T_TRUE : T_FALSE;
    }
    else
    {
        return T_UNKNOWN;
    }
}

static truth_t ca_check_is_zero_and_simplify(ca_t x, ca_ctx_t ctx)
{
    truth_t res = ca_check_is_zero_fast(x, ctx);

    if (res == T_UNKNOWN)
    {
        res = ca_check_is_zero(x, ctx);

        if (res == T_TRUE)
            ca_zero(x, ctx);
    }

    return res;
}

truth_t
ca_mat_find_pivot(slong * pivot_row, ca_mat_t mat, slong start_row, slong end_row, slong column, ca_ctx_t ctx)
{
    slong best_row, i;
    truth_t is_zero;
    int unknown;

    if (end_row <= start_row)
        flint_throw(FLINT_ERROR, "(%s): end_row <= start_row\n", __func__);

    /* First find the simplest element that is not trivially zero. With high probability
       this will actually be nonzero. */
    best_row = -1;

    for (i = start_row; i < end_row; i++)
    {
        is_zero = ca_check_is_zero_fast(ca_mat_entry(mat, i, column), ctx);

        if (is_zero != T_TRUE)
        {
            if (best_row == -1 || ca_cmp_repr(ca_mat_entry(mat, i, column), ca_mat_entry(mat, best_row, column), ctx) < 0)
            {
                best_row = i;
            }
        }
    }

    if (best_row != -1)
    {
        is_zero = ca_check_is_zero_and_simplify(ca_mat_entry(mat, best_row, column), ctx);

        if (is_zero == T_FALSE)
        {
            *pivot_row = best_row;
            return T_TRUE;
        }
    }

    /* If the above failed, go through all elements again and do more expensive checks.
       Todo: 1) support in-place simplifications.
             2) consider sorting above and traversing all entries in order of
                simplicity (not just the simplest element). */

    best_row = -1;
    unknown = 0;

    for (i = start_row; i < end_row; i++)
    {
        is_zero = ca_check_is_zero_and_simplify(ca_mat_entry(mat, i, column), ctx);

        if (is_zero == T_FALSE)
        {
            if (best_row == -1 || ca_cmp_repr(ca_mat_entry(mat, i, column), ca_mat_entry(mat, best_row, column), ctx) < 0)
            {
                best_row = i;
            }
        }

        if (is_zero == T_UNKNOWN)
            unknown = 1;
    }

    if (best_row == -1)
    {
        *pivot_row = -1;
        if (unknown)
            return T_UNKNOWN;
        else
            return T_FALSE;
    }
    else
    {
        *pivot_row = best_row;
        return T_TRUE;
    }
}
