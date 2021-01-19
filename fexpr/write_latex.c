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
        else if (fexpr_is_builtin_symbol(expr))
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
    }
    else
    {
        fexpr_write_latex_call(out, expr, flags);
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
