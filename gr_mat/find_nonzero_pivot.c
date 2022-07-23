/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"

/*
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
*/

int
gr_cmp_repr(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    return 0;
}

#include "arf.h"

int
gr_mat_find_nonzero_pivot(slong * pivot_row, gr_mat_t mat, slong start_row, slong end_row, slong column, gr_ctx_t ctx)
{
    slong best_row, i;
    int unknown;
    truth_t is_zero;
    slong sz;

    if (end_row <= start_row)
        return GR_DOMAIN;

/* arf_mat */
#if 0
    {
        {
            slong best_row, i;

            best_row = -1;

            for (i = start_row; i < end_row; i++)
            {
                if (!arf_is_zero(((arf_srcptr) (mat->rows[i])) + column))
                {
                    if (best_row == -1)
                    {
                        best_row = i;
                    }
                    /* todo: should take the radius into account */
                    else if (arf_cmpabs(((arf_srcptr) (mat->rows[i])) + column,
                                        ((arf_srcptr) (mat->rows[best_row])) + column) > 0)
                    {
                        best_row = i;
                    }
                }
            }

            *pivot_row = best_row;
            return GR_SUCCESS;
        }
    }
#endif


/* todo: for ca-type rings */
#if 0
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
#endif

    /* If the above failed, go through all elements again and do more expensive checks.
       Todo: 1) support in-place simplifications.
             2) consider sorting above and traversing all entries in order of
                simplicity (not just the simplest element). */

    best_row = -1;
    unknown = 0;

    sz = ctx->sizeof_elem;

    for (i = start_row; i < end_row; i++)
    {
        is_zero = gr_is_zero(GR_MAT_ENTRY(mat, i, column, sz), ctx);

        if (is_zero != T_UNKNOWN)
        {
            if (is_zero == T_FALSE)
            {
                if (best_row == -1 || gr_cmp_repr(GR_MAT_ENTRY(mat, i, column, sz), GR_MAT_ENTRY(mat, best_row, column, sz), ctx) < 0)
                {
                    best_row = i;
                }
            }
        }
        else
        {
            unknown = 1;
        }
    }

    if (best_row == -1)
    {
        *pivot_row = -1;
        if (unknown)
            return GR_UNABLE;
        else
            return GR_DOMAIN;
    }
    else
    {
        *pivot_row = best_row;
        return GR_SUCCESS;
    }
}
