/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"
#include "gr_vec.h"

int
gr_mat_set_str(gr_mat_t mat, const char * s, int resize, gr_ctx_t ctx)
{
    const char * lo;
    const char * hi;
    const char * seg;
    slong r, c, i;
    int status;

    /* Validate the outer brackets. */
    status = _gr_str_strip_brackets(&lo, &hi, s, '[', ']');
    if (status != GR_SUCCESS)
        return status;

    /* Empty matrix "[]": zero rows, ambiguous column count. With resize it
       becomes 0 x 0; without resize it is compatible with any 0 x c shape. */
    if (_gr_str_blank(lo, hi))
    {
        if (resize)
        {
            gr_mat_clear(mat, ctx);
            gr_mat_init(mat, 0, 0, ctx);
            return GR_SUCCESS;
        }

        return (gr_mat_nrows(mat, ctx) == 0) ? GR_SUCCESS : GR_DOMAIN;
    }

    /* Number of rows = number of top-level entries of the whole string. */
    status = gr_vec_str_count_entries(&r, s, ctx);
    if (status != GR_SUCCESS)
        return status;

    /* Number of columns = number of entries in the first row. The first-row
       pointer carries trailing characters (the remaining rows), which
       gr_vec_str_count_entries ignores. */
    {
        const char * p = lo;
        while (p < hi && _gr_str_is_space(*p))
            p++;
        status = gr_vec_str_count_entries(&c, p, ctx);
        if (status != GR_SUCCESS)
            return status;
    }

    if (resize)
    {
        if (mat->r != r || mat->c != c)
        {
            gr_mat_clear(mat, ctx);
            gr_mat_init(mat, r, c, ctx);
        }
    }
    else if (gr_mat_nrows(mat, ctx) != r || gr_mat_ncols(mat, ctx) != c)
    {
        return GR_DOMAIN;
    }

    /* Evaluate each row. _gr_vec_set_str enforces exactly c entries. */
    seg = lo;
    for (i = 0; i < r; i++)
    {
        int err = 0;
        const char * sep = _gr_str_find_sep(seg, hi, &err);
        const char * a = seg;
        const char * b = sep;

        _gr_str_trim(&a, &b);   /* a points at this row's '[' */

        status = _gr_vec_set_str(gr_mat_entry_ptr(mat, i, 0, ctx), a, c, ctx);
        if (status != GR_SUCCESS)
            break;

        seg = sep + 1;
    }

    if (status != GR_SUCCESS && resize)
    {
        gr_mat_clear(mat, ctx);
        gr_mat_init(mat, 0, 0, ctx);
    }

    return status;
}
