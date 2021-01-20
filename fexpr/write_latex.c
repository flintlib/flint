/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fexpr.h"
#include "fexpr_builtin.h"

static int
fexpr_is_neg_integer(const fexpr_t expr)
{
    ulong head = expr->data[0];

    if (FEXPR_TYPE(head) == FEXPR_TYPE_SMALL_INT)
        return (((slong) head) >> FEXPR_TYPE_BITS) < 0;

    if (FEXPR_TYPE(head) == FEXPR_TYPE_BIG_INT_NEG)
        return 1;

    return 0;
}

static int
fexpr_need_parens_in_mul(const fexpr_t expr, slong arg_index)
{
    if (fexpr_is_atom(expr))
    {
        if (arg_index == 0)
            return 0;

        return fexpr_is_neg_integer(expr);
    }
    else
    {
        fexpr_t func;
        fexpr_view_func(func, expr);

        if (fexpr_is_builtin_symbol(func, FEXPR_Add))
            return 1;

        if (fexpr_is_builtin_symbol(func, FEXPR_Sub))
            return 1;

        if (fexpr_is_builtin_symbol(func, FEXPR_Neg))
            return (arg_index != 0);

        if (fexpr_is_builtin_symbol(func, FEXPR_Pos))
            return (arg_index != 0);

        return 0;
    }
}

static int
fexpr_is_builtin_call(const fexpr_t expr, slong i)
{
    fexpr_t func;

    if (fexpr_is_atom(expr))
        return 0;

    fexpr_view_func(func, expr);

    return fexpr_is_builtin_symbol(func, i);
}

/* todo */
static int
fexpr_need_cdot_before_factor(const fexpr_t expr)
{
    if (fexpr_is_integer(expr))
        return 1;

    if (fexpr_is_builtin_call(expr, FEXPR_Mul) && fexpr_nargs(expr) >= 1)
    {
        fexpr_t first;
        fexpr_view_arg(first, expr, 0);
        return fexpr_need_cdot_before_factor(first);
    }

    if (fexpr_is_builtin_call(expr, FEXPR_Neg) ||
        fexpr_is_builtin_call(expr, FEXPR_Pos))
        return 1;

    return 0;
}

void
fexpr_write_latex_mul(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    fexpr_t arg;
    slong i, len;
    int need_parens;

    len = fexpr_nargs(expr);

    if (len == 0)
        calcium_write(out, "(1)");

    fexpr_view_arg(arg, expr, 0);

    for (i = 0; i < len; i++)
    {
        need_parens = fexpr_need_parens_in_mul(arg, i);

        if (need_parens)
            calcium_write(out, "\\left(");

        fexpr_write_latex(out, arg, flags);

        if (need_parens)
            calcium_write(out, "\\right)");

        if (i < len - 1)
        {
            fexpr_view_next(arg);

            if (fexpr_need_cdot_before_factor(arg))
                calcium_write(out, " \\cdot ");
            else
                calcium_write(out, " ");
        }
    }
}

void
fexpr_write_latex_call(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    fexpr_t view;
    slong i, nargs;

    nargs = fexpr_nargs(expr);
    fexpr_view_func(view, expr);
    fexpr_write_latex(out, view, flags);

    calcium_write(out, "\\!\\left(");

    for (i = 0; i < nargs; i++)
    {
        fexpr_view_next(view);
        fexpr_write_latex(out, view, flags);
        if (i < nargs - 1)
            calcium_write(out, ", ");
    }

    calcium_write(out, "\\right)");
}

void fexpr_write_latex_subscript_call(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    
}

void
fexpr_write_latex_infix(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    fexpr_t func, arg;
    slong i, nargs;

    nargs = fexpr_nargs(expr);

    fexpr_view_func(func, expr);
    fexpr_view_func(arg, expr);

    for (i = 0; i < nargs; i++)
    {
        fexpr_view_next(arg);
        fexpr_write_latex(out, arg, flags);
        if (i < nargs - 1)
        {
            calcium_write(out, " ");
            fexpr_write_latex(out, func, flags);
            calcium_write(out, " ");
        }
    }
}


void
fexpr_write_latex(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    if (fexpr_is_atom(expr))
    {
        if (fexpr_is_integer(expr))
        {
            fexpr_write(out, expr);
        }
        else if (fexpr_is_string(expr))
        {
            flint_abort();
        }
        else if (fexpr_is_any_builtin_symbol(expr))
        {
            slong i;

            i = FEXPR_BUILTIN_ID(expr->data[0]);

            if (strcmp(fexpr_builtin_table[i].latex_string, ""))
            {
                calcium_write(out, fexpr_builtin_table[i].latex_string);
            }
            else
            {
                calcium_write(out, "\\operatorname{");
                calcium_write(out, fexpr_builtin_table[i].string);
                calcium_write(out, "}");
            }
        }
        else
        {
            fexpr_write(out, expr);
        }
    }
    else
    {
        fexpr_t func;
        slong i;

        fexpr_view_func(func, expr);
        i = FEXPR_BUILTIN_ID(func->data[0]);

        if (i != -1 && fexpr_builtin_table[i].latex_writer != NULL)
        {
            (fexpr_builtin_table[i].latex_writer)(out, expr, flags);
        }
        else
        {
            fexpr_write_latex_call(out, expr, flags);
        }
    }
}

void
fexpr_print_latex(const fexpr_t expr, ulong flags)
{
    calcium_stream_t t;
    calcium_stream_init_file(t, stdout);
    fexpr_write_latex(t, expr, flags);
}

char *
fexpr_get_str_latex(const fexpr_t expr, ulong flags)
{
    calcium_stream_t t;
    calcium_stream_init_str(t);
    fexpr_write_latex(t, expr, flags);
    return t->s;
}
