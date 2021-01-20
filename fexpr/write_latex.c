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
    {
        calcium_write(out, "(1)");
        return;
    }

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

static int
fexpr_need_parens_in_numerator(const fexpr_t expr)
{
    if (fexpr_is_atom(expr))
    {
        return 0;
    }
    else
    {
        fexpr_t func;
        fexpr_view_func(func, expr);

        if (fexpr_is_builtin_symbol(func, FEXPR_Add))
            return 1;

        if (fexpr_is_builtin_symbol(func, FEXPR_Sub))
            return 1;

        return 0;
    }
}

static int
fexpr_need_parens_in_denominator(const fexpr_t expr)
{
    if (fexpr_is_atom(expr))
    {
        return 0;
    }
    else
    {
        fexpr_t func;
        fexpr_view_func(func, expr);

        if (fexpr_is_builtin_symbol(func, FEXPR_Add))
            return 1;

        if (fexpr_is_builtin_symbol(func, FEXPR_Sub))
            return 1;

        if (fexpr_is_builtin_symbol(func, FEXPR_Mul))
            return 1;

        if (fexpr_is_builtin_symbol(func, FEXPR_Div))
            return 1;

        return 0;
    }
}

static int
fexpr_can_extract_leading_sign(const fexpr_t expr)
{
    if (fexpr_is_atom(expr))
    {
        return fexpr_is_neg_integer(expr);
    }
    else
    {
        fexpr_t func;
        fexpr_view_func(func, expr);

        if (fexpr_is_builtin_symbol(func, FEXPR_Neg))
            return 1;

        if (fexpr_is_builtin_symbol(func, FEXPR_Pos))
            return 1;

        if (fexpr_is_builtin_symbol(func, FEXPR_Mul) &&
            fexpr_nargs(expr) >= 1)
        {
            fexpr_view_next(func);
            return fexpr_can_extract_leading_sign(func);
        }

        return 0;
    }
}

void
fexpr_write_latex_div(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    fexpr_t num, den;

    /* Expect exactly 2 arguments */
    if (fexpr_nargs(expr) != 2)
    {
        fexpr_write_latex_call(out, expr, flags);
        return;
    }

    fexpr_view_arg(num, expr, 0);
    fexpr_view_arg(den, expr, 1);

    if (flags & FEXPR_LATEX_SMALL)
    {
        int pnum, pden;

        pnum = fexpr_need_parens_in_numerator(num);
        pden = fexpr_need_parens_in_denominator(den);

        if (pnum)
            calcium_write(out, "\\left(");
        fexpr_write_latex(out, num, flags);
        if (pnum)
            calcium_write(out, "\\right)");

        calcium_write(out, " / ");

        if (pden)
            calcium_write(out, "\\left(");
        fexpr_write_latex(out, den, flags);
        if (pden)
            calcium_write(out, "\\right)");
    }
    else
    {
        if (fexpr_can_extract_leading_sign(num))
        {
            char * s = fexpr_get_str_latex(num, flags);

            if (s[0] == '+' || s[0] == '-')
            {
                char tmp[2];
                tmp[0] = s[0];
                tmp[1] = '\0';
                calcium_write(out, tmp);
                calcium_write(out, "\\frac{");
                calcium_write(out, s + 1);
                calcium_write(out, "}{");
                fexpr_write_latex(out, den, flags);
                calcium_write(out, "}");
            }
            else
            {
                calcium_write(out, "\\frac{");
                fexpr_write_latex(out, num, flags);
                calcium_write(out, "}{");
                fexpr_write_latex(out, den, flags);
                calcium_write(out, "}");
            }

            flint_free(s);
        }
        else
        {
            calcium_write(out, "\\frac{");
            fexpr_write_latex(out, num, flags);
            calcium_write(out, "}{");
            fexpr_write_latex(out, den, flags);
            calcium_write(out, "}");
        }
    }
}

void
fexpr_write_latex_neg_pos(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    fexpr_t arg;

    /* Expect exactly 1 argument */
    if (fexpr_nargs(expr) != 1)
    {
        fexpr_write_latex_call(out, expr, flags);
        return;
    }

    if (fexpr_is_builtin_call(expr, FEXPR_Pos))
        calcium_write(out, "+");
    else
        calcium_write(out, "-");

    fexpr_view_arg(arg, expr, 0);

    if (fexpr_is_builtin_call(arg, FEXPR_Neg) ||
        fexpr_is_builtin_call(arg, FEXPR_Pos) ||
        fexpr_is_builtin_call(arg, FEXPR_Add) ||
        fexpr_is_builtin_call(arg, FEXPR_Sub) ||
        fexpr_is_neg_integer(arg))
    {
        calcium_write(out, "\\left(");
        fexpr_write_latex(out, arg, flags);
        calcium_write(out, "\\right)");
    }
    else
    {
        fexpr_write_latex(out, arg, flags);
    }
}

static int
fexpr_need_parens_in_add(const fexpr_t expr)
{
    if (fexpr_is_atom(expr))
    {
        return 0;
    }
    else
    {
        fexpr_t func;
        fexpr_view_func(func, expr);

        if (fexpr_is_builtin_symbol(func, FEXPR_Sub))
            return 1;

        if (fexpr_is_builtin_symbol(func, FEXPR_Neg))
            return 1;

        return 0;
    }
}

void
fexpr_write_latex_add(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    fexpr_t arg;
    slong i, len;

    len = fexpr_nargs(expr);

    if (len == 0)
    {
        calcium_write(out, "(0)");
        return;
    }

    fexpr_view_arg(arg, expr, 0);

    for (i = 0; i < len; i++)
    {
        if (i == 0)
        {
            fexpr_write_latex(out, arg, flags);
        }
        else
        {
            int need_parens = fexpr_need_parens_in_add(arg);

            if (need_parens)
            {
                calcium_write(out, " + \\left(");
                fexpr_write_latex(out, arg, flags);
                calcium_write(out, "\\right)");
            }
            else
            {
                char * s = fexpr_get_str_latex(arg, flags);

                if (s[0] == '+' || s[0] == '-')
                {
                    calcium_write(out, s);
                }
                else
                {
                    calcium_write(out, " + ");
                    calcium_write(out, s);
                }

                flint_free(s);
            }
        }

        fexpr_view_next(arg);
    }
}

static int
fexpr_need_parens_in_sub(const fexpr_t expr)
{
    if (fexpr_is_atom(expr))
    {
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
            return 1;

        if ((fexpr_is_builtin_symbol(func, FEXPR_Mul) ||
            fexpr_is_builtin_symbol(func, FEXPR_Div)) &&
            fexpr_nargs(expr) >= 1)
        {
            fexpr_t arg;
            fexpr_view_arg(arg, expr, 0);
            return fexpr_can_extract_leading_sign(arg);
        }

        return 0;
    }
}

void
fexpr_write_latex_sub(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    fexpr_t arg;
    slong i, len;

    len = fexpr_nargs(expr);

    if (len == 0)
    {
        calcium_write(out, "(0)");
        return;
    }

    fexpr_view_arg(arg, expr, 0);

    for (i = 0; i < len; i++)
    {
        if (i == 0)
        {
            fexpr_write_latex(out, arg, flags);
        }
        else
        {
            int need_parens = fexpr_need_parens_in_sub(arg);

            if (need_parens)
            {
                calcium_write(out, " - \\left(");
                fexpr_write_latex(out, arg, flags);
                calcium_write(out, "\\right)");
            }
            else
            {
                calcium_write(out, " - ");
                fexpr_write_latex(out, arg, flags);
            }
        }

        fexpr_view_next(arg);
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
