/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fexpr.h"

/*
   Sets res_view to a view of expr replaced by {xs, ys, len} and returns whether
   any replacement has been made.

   If no replacement has been made, res_view will refer to the data of expr,
   with an alloc field of 0.

   If expr as a whole has been replaced with some entry in ys, res_view will
   refer to the data of that value in ys, also with an alloc field of 0.

   Otherwise, expr will be initialized to a new expression, with a nonzero
   alloc field.
 */
int
_fexpr_replace_vec(fexpr_t res_view, const fexpr_t expr, fexpr_srcptr xs, fexpr_srcptr ys, slong len)
{
    slong i, nargs;
    fexpr_t func, new_func, arg;
    fexpr_struct tmp_args[4];
    fexpr_ptr args;
    int changed;

    for (i = 0; i < len; i++)
    {
        if (fexpr_equal(expr, xs + i))
        {
            res_view->data = ys[i].data;
            res_view->alloc = 0;
            return 1;
        }
    }

    if (fexpr_is_atom(expr))
    {
        res_view->data = expr->data;
        res_view->alloc = 0;
        return 0;
    }

    nargs = fexpr_nargs(expr);

    fexpr_view_func(func, expr);
    changed = _fexpr_replace_vec(new_func, func, xs, ys, len);

    if (nargs >= 5)
        args = flint_malloc(sizeof(fexpr_struct) * nargs);
    else
        args = tmp_args;

    if (nargs >= 1)
    {
        fexpr_view_arg(arg, expr, 0);

        for (i = 0; i < nargs; i++)
        {
            changed |= _fexpr_replace_vec(args + i, arg, xs, ys, len);

            if (i < nargs - 1)
                fexpr_view_next(arg);
        }
    }

    if (changed)
    {
        fexpr_init(res_view);
        fexpr_call_vec(res_view, new_func, args, nargs);

        if (new_func->alloc != 0)
            flint_free(new_func->data);

        for (i = 0; i < nargs; i++)
        {
            if (args[i].alloc != 0)
                flint_free(args[i].data);
        }
    }
    else
    {
        res_view->data = expr->data;
        res_view->alloc = 0;
    }

    if (nargs >= 5)
        flint_free(args);

    return changed;
}

int
fexpr_replace_vec(fexpr_t res, const fexpr_t expr, const fexpr_vec_t xs, const fexpr_vec_t ys)
{
    int changed;
    slong num_rules;
    fexpr_t res_view;

    num_rules = xs->length;

    if (num_rules != ys->length)
    {
        flint_throw(FLINT_ERROR, "fexpr_replace_vec: vectors don't match\n");
    }

    if (num_rules == 0)
    {
        fexpr_set(res, expr);
        return 0;
    }

    changed = _fexpr_replace_vec(res_view, expr, xs->entries, ys->entries, num_rules);

    if (changed)
    {
        if (res_view->alloc != 0)
        {
            /* We have new data */
            fexpr_swap(res, res_view);
            fexpr_clear(res_view);
        }
        else
        {
            /* Must be a view of some ys; we assume that res is not aliased
               with any entry in ys. */
            fexpr_set(res, res_view);
        }
    }
    else
    {
        fexpr_set(res, expr);
    }

    return changed;
}

int
fexpr_replace(fexpr_t res, const fexpr_t expr, const fexpr_t x, const fexpr_t y)
{
    int changed;
    fexpr_t res_view;

    changed = _fexpr_replace_vec(res_view, expr, x, y, 1);

    if (changed)
    {
        if (res_view->alloc != 0)
        {
            /* We have new data */
            fexpr_swap(res, res_view);
            fexpr_clear(res_view);
        }
        else
        {
            /* Must be a view of some ys; we assume that res is not aliased
               with any entry in ys. */
            fexpr_set(res, res_view);
        }
    }
    else
    {
        fexpr_set(res, expr);
    }

    return changed;
}

int
fexpr_replace2(fexpr_t res, const fexpr_t expr, const fexpr_t x1, const fexpr_t y1, const fexpr_t x2, const fexpr_t y2)
{
    int changed;
    fexpr_t res_view;
    fexpr_struct tmp[4];

    tmp[0] = *x1;
    tmp[1] = *x2;
    tmp[2] = *y1;
    tmp[3] = *y2;

    changed = _fexpr_replace_vec(res_view, expr, tmp, tmp + 2, 2);

    if (changed)
    {
        if (res_view->alloc != 0)
        {
            /* We have new data */
            fexpr_swap(res, res_view);
            fexpr_clear(res_view);
        }
        else
        {
            /* Must be a view of some ys; we assume that res is not aliased
               with any entry in ys. */
            fexpr_set(res, res_view);
        }
    }
    else
    {
        fexpr_set(res, expr);
    }

    return changed;
}
