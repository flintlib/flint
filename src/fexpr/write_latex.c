/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "fexpr.h"
#include "fexpr_builtin.h"

#ifdef __GNUC__
# define memcpy __builtin_memcpy
# define strcmp __builtin_strcmp
# define strlen __builtin_strlen
#else
# include <string.h>
#endif

void fexpr_write_latex_symbol(int * subscript, calcium_stream_t out, const fexpr_t expr, ulong flags);
int _fexpr_is_symbol_with_trailing_underscore(const fexpr_t expr);

const char * fexpr_get_symbol_str_pointer(char * tmp, const fexpr_t expr)
{
    slong i;
    ulong head = expr->data[0];

    if (FEXPR_TYPE(head) == FEXPR_TYPE_SMALL_SYMBOL)
    {
        if (((head >> 8) & 0xff) == 0)
        {
            i = head >> 16;
            return fexpr_builtin_table[i].string;
        }

        tmp[FEXPR_SMALL_SYMBOL_LEN] = '\0';

        for (i = 0; i < FEXPR_SMALL_SYMBOL_LEN; i++)
        {
            tmp[i] = (head >> ((i + 1) * 8));
            if (tmp[i] == '\0')
                break;
        }

        return tmp;
    }
    else if (FEXPR_TYPE(head) == FEXPR_TYPE_BIG_SYMBOL)
    {
        return (const char *) (expr->data + 1);
    }
    else
    {
        flint_throw(FLINT_ERROR, "fexpr_get_symbol_str_pointer: a symbol is required\n");
    }
}

int _fexpr_is_symbol_with_trailing_underscore(const fexpr_t expr)
{
    const char *s;
    char tmp[FEXPR_SMALL_SYMBOL_LEN + 1];
    slong len;

    if (!fexpr_is_symbol(expr))
        return 0;

    s = fexpr_get_symbol_str_pointer(tmp, expr);
    len = strlen(s);

    return (len > 1 && s[len - 1] == '_');
}

int _fexpr_is_symbol_with_internal_underscore(const fexpr_t expr)
{
    const char *s;
    char tmp[FEXPR_SMALL_SYMBOL_LEN + 1];
    slong i, len;

    if (!fexpr_is_symbol(expr))
        return 0;

    s = fexpr_get_symbol_str_pointer(tmp, expr);
    len = strlen(s);

    if (len >= 3)
    {
        for (i = 1; i < len - 1; i++)
            if (s[i] == '_')
                return i;
    }

    return 0;
}

/*
static int
fexpr_view_call0(fexpr_t func, const fexpr_t expr)
{
    slong nargs;
    nargs = fexpr_nargs(expr);
    if (nargs != 0)
        return 0;
    fexpr_view_func(func, expr);
    return 1;
}
*/

static int
fexpr_view_call1(fexpr_t func, fexpr_t x1, const fexpr_t expr)
{
    slong nargs;
    nargs = fexpr_nargs(expr);
    if (nargs != 1)
        return 0;
    fexpr_view_func(func, expr);
    *x1 = *func;
    fexpr_view_next(x1);
    return 1;
}

/*
static int
fexpr_view_call2(fexpr_t func, fexpr_t x1, fexpr_t x2, const fexpr_t expr)
{
    slong nargs;
    nargs = fexpr_nargs(expr);
    if (nargs != 2)
        return 0;
    fexpr_view_func(func, expr);
    *x1 = *func;
    fexpr_view_next(x1);
    *x2 = *x1;
    fexpr_view_next(x2);
    return 1;
}
*/

static int
fexpr_view_call3(fexpr_t func, fexpr_t x1, fexpr_t x2, fexpr_t x3, const fexpr_t expr)
{
    slong nargs;
    nargs = fexpr_nargs(expr);
    if (nargs != 3)
        return 0;
    fexpr_view_func(func, expr);
    *x1 = *func;
    fexpr_view_next(x1);
    *x2 = *x1;
    fexpr_view_next(x2);
    *x3 = *x2;
    fexpr_view_next(x3);
    return 1;
}

/*
static int
fexpr_view_call4(fexpr_t func, fexpr_t x1, fexpr_t x2, fexpr_t x3, fexpr_t x4, const fexpr_t expr)
{
    slong nargs;
    nargs = fexpr_nargs(expr);
    if (nargs != 4)
        return 0;
    fexpr_view_func(func, expr);
    *x1 = *func;
    fexpr_view_next(x1);
    *x2 = *x1;
    fexpr_view_next(x2);
    *x3 = *x2;
    fexpr_view_next(x3);
    *x4 = *x3;
    fexpr_view_next(x4);
    return 1;
}
*/

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

/* todo */
static int
fexpr_need_cdot_before_factor(const fexpr_t expr)
{
    if (fexpr_is_integer(expr))
        return 1;

    if (fexpr_is_builtin_symbol(expr, FEXPR_Infinity) ||
        fexpr_is_builtin_symbol(expr, FEXPR_UnsignedInfinity))
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

    if (fexpr_is_builtin_call(expr, FEXPR_Pow) && fexpr_nargs(expr) == 2)
    {
        fexpr_t first;
        fexpr_view_arg(first, expr, 0);
        if (fexpr_is_integer(first))
            return 1;
    }

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
        /* 1 * x * y or -1 * x * y  ->  x * y,  - x * y */
        if (i == 0 && len >= 2 && fexpr_is_integer(arg) &&
            (fexpr_equal_si(arg, 1) || fexpr_equal_si(arg, -1)))
        {
            if (fexpr_equal_si(arg, -1))
                calcium_write(out, "-");

            fexpr_view_next(arg);

            if (fexpr_need_cdot_before_factor(arg))
                calcium_write(out, "1 \\cdot ");
            else
                calcium_write(out, " ");
        }
        else
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

static int
fexpr_power_base_is_safe(const fexpr_t base)
{
    if (fexpr_is_atom(base))
    {
        if (fexpr_is_neg_integer(base))
            return 0;

        return 1;
    }
    else
    {
        slong id;
        fexpr_t func;

        if (fexpr_is_any_builtin_call(base))
        {
            fexpr_view_func(func, base);
            id = FEXPR_BUILTIN_ID(func->data[0]);

            if (id == FEXPR_Abs) return 1;
            if (id == FEXPR_Binomial) return 1;
            if (id == FEXPR_Matrix) return 1;
            if (id == FEXPR_Matrix2x2) return 1;
            if (id == FEXPR_Parentheses) return 1;
            if (id == FEXPR_Brackets) return 1;
            if (id == FEXPR_Braces) return 1;
            if (id == FEXPR_AngleBrackets) return 1;
        }

        return 0;
    }
}

void
_fexpr_write_latex_pow(calcium_stream_t out, const fexpr_t base, const fexpr_t expo, ulong flags)
{
    if (fexpr_is_any_builtin_call(base) && fexpr_nargs(base) == 1)
    {
        fexpr_t func, arg;

        fexpr_view_func(func, base);

        /* (f(x))^n written as f^n(x) for some standard functions */
        switch (FEXPR_BUILTIN_ID(func->data[0]))
        {
            case FEXPR_Sin:
            case FEXPR_Cos:
            case FEXPR_Tan:
            case FEXPR_Cot:
            case FEXPR_Sec:
            case FEXPR_Csc:
            case FEXPR_Sinh:
            case FEXPR_Cosh:
            case FEXPR_Tanh:
            case FEXPR_Coth:
            case FEXPR_Sech:
            case FEXPR_Csch:
            case FEXPR_Log:
            case FEXPR_Sinc:
            case FEXPR_DedekindEta:
                fexpr_write_latex(out, func, flags);
                calcium_write(out, "^{");
                fexpr_write_latex(out, expo, flags | FEXPR_LATEX_SMALL);
                fexpr_view_arg(arg, base, 0);
                calcium_write(out, "}\\!\\left(");
                fexpr_write_latex(out, arg, flags);
                calcium_write(out, "\\right)");
                return;
            default:
                break;
        }
    }

    /* todo: generalize */
    if ((fexpr_is_builtin_call(base, FEXPR_JacobiTheta) || fexpr_is_builtin_call(base, FEXPR_JacobiThetaQ))
            && fexpr_nargs(base) == 3)
    {
        fexpr_t func, x, y, z;

        fexpr_view_func(func, base);
        fexpr_view_arg(x, base, 0);
        fexpr_view_arg(y, base, 1);
        fexpr_view_arg(z, base, 2);
        calcium_write(out, "\\theta_{");
        fexpr_write_latex(out, x, flags | FEXPR_LATEX_SMALL);
        calcium_write(out, "}^{");
        fexpr_write_latex(out, expo, flags | FEXPR_LATEX_SMALL);
        calcium_write(out, "}\\!\\left(");
        fexpr_write_latex(out, y, flags);
        calcium_write(out, ", ");
        fexpr_write_latex(out, z, flags);
        calcium_write(out, "\\right)");
        return;
    }

    if (fexpr_is_symbol(base))
    {
        slong underscore_pos;
        char * s;

        /*  Need to change {x_3}^{y} -> x^{y}_{3} for correct rendering */
        underscore_pos = _fexpr_is_symbol_with_internal_underscore(base);

        if (underscore_pos > 0)
        {
            char *sym, *sub;

            s = fexpr_get_symbol_str(base);
            s[underscore_pos] = '\0';
            sym = s;
            sub = s + underscore_pos + 1;

            if (underscore_pos == 1)
            {
                calcium_write(out, sym);
            }
            else
            {
                calcium_write(out, "\\operatorname{");
                calcium_write(out, sym);
                calcium_write(out, "}");
            }

            calcium_write(out, "^{");
            fexpr_write_latex(out, expo, flags | FEXPR_LATEX_SMALL);
            calcium_write(out, "}_{");

            calcium_write(out, sub);

            calcium_write(out, "}");
            return;
        }
    }

    if (fexpr_is_any_builtin_call(base))
    {
        slong id;
        fexpr_t func, x, y;

        fexpr_view_func(func, base);

        id = FEXPR_BUILTIN_ID(func->data[0]);

        if (fexpr_builtin_table[id].latex_writer == fexpr_write_latex_subscript && fexpr_nargs(base) == 1)
        {
            fexpr_view_arg(x, base, 0);

            fexpr_write_latex(out, func, flags);
            calcium_write(out, "_{");
            fexpr_write_latex(out, x, flags | FEXPR_LATEX_SMALL);
            calcium_write(out, "}^{");
            fexpr_write_latex(out, expo, flags | FEXPR_LATEX_SMALL);
            calcium_write(out, "}");
            return;
        }

        if (fexpr_builtin_table[id].latex_writer == fexpr_write_latex_subscript_call && fexpr_nargs(base) == 2)
        {
            fexpr_view_arg(x, base, 0);
            fexpr_view_arg(y, base, 1);

            fexpr_write_latex(out, func, flags);
            calcium_write(out, "_{");
            fexpr_write_latex(out, x, flags | FEXPR_LATEX_SMALL);
            calcium_write(out, "}^{");
            fexpr_write_latex(out, expo, flags | FEXPR_LATEX_SMALL);
            calcium_write(out, "}\\!\\left(");
            fexpr_write_latex(out, y, flags);
            calcium_write(out, "\\right)");
            return;
        }
    }

    if (fexpr_is_builtin_call(base, FEXPR_Subscript) && fexpr_nargs(base) == 2)
    {
        fexpr_t func, x;

        fexpr_view_arg(func, base, 0);
        fexpr_view_arg(x, base, 1);

        calcium_write(out, "{");
        fexpr_write_latex(out, func, flags);
        calcium_write(out, "}_{");
        fexpr_write_latex(out, x, flags | FEXPR_LATEX_SMALL);
        calcium_write(out, "}^{");
        fexpr_write_latex(out, expo, flags | FEXPR_LATEX_SMALL);
        calcium_write(out, "}");
        return;
    }

    if (fexpr_nargs(base) == 1)
    {
        fexpr_t func, x;
        int subscript;

        fexpr_view_func(func, base);

        if (_fexpr_is_symbol_with_trailing_underscore(func))
        {
            fexpr_view_arg(x, base, 0);
            fexpr_write_latex_symbol(&subscript, out, func, flags);

            calcium_write(out, "_{");
            fexpr_write_latex(out, x, flags | FEXPR_LATEX_SMALL);
            calcium_write(out, "}^{");
            fexpr_write_latex(out, expo, flags | FEXPR_LATEX_SMALL);
            calcium_write(out, "}");

            return;
        }
    }

    if (fexpr_power_base_is_safe(base))
    {
        calcium_write(out, "{");
        fexpr_write_latex(out, base, flags);
        calcium_write(out, "}^{");
        fexpr_write_latex(out, expo, flags | FEXPR_LATEX_SMALL);
        calcium_write(out, "}");
    }
    else
    {
        calcium_write(out, "{\\left(");
        fexpr_write_latex(out, base, flags);
        calcium_write(out, "\\right)}^{");
        fexpr_write_latex(out, expo, flags | FEXPR_LATEX_SMALL);
        calcium_write(out, "}");
    }
}

void
fexpr_write_latex_pow(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    if (fexpr_nargs(expr) == 2)
    {
        fexpr_t base, expo;

        fexpr_view_arg(base, expr, 0);
        fexpr_view_arg(expo, expr, 1);

        _fexpr_write_latex_pow(out, base, expo, flags);
    }
    else
    {
        fexpr_write_latex_call(out, expr, flags);
    }
}

static int
fexpr_show_exp_as_power(const fexpr_t expr)
{
    fexpr_t func, arg;
    slong i, nargs;

    if (fexpr_is_atom(expr))
        return 1;

    fexpr_view_func(func, expr);

    /* todo: more systematic solution */
    if (fexpr_is_builtin_symbol(func, FEXPR_Pos) ||
        fexpr_is_builtin_symbol(func, FEXPR_Neg) ||
        fexpr_is_builtin_symbol(func, FEXPR_Add) ||
        fexpr_is_builtin_symbol(func, FEXPR_Sub) ||
        fexpr_is_builtin_symbol(func, FEXPR_Mul) ||
        fexpr_is_builtin_symbol(func, FEXPR_Div) ||
        fexpr_is_builtin_symbol(func, FEXPR_Pow) ||
        fexpr_is_builtin_symbol(func, FEXPR_Abs) ||
        fexpr_is_builtin_symbol(func, FEXPR_RealAbs) ||
        fexpr_is_builtin_symbol(func, FEXPR_Sqrt) ||
        fexpr_is_builtin_symbol(func, FEXPR_Re) ||
        fexpr_is_builtin_symbol(func, FEXPR_Im) ||
        fexpr_is_builtin_symbol(func, FEXPR_Log))
    {
        nargs = fexpr_nargs(expr);

        if (fexpr_is_builtin_symbol(func, FEXPR_Div) && nargs == 2)
        {
            fexpr_view_arg(arg, expr, 1);
            if (!fexpr_is_atom(arg))
                return 0;
        }

        fexpr_view_arg(arg, expr, 0);

        for (i = 0; i < nargs; i++)
        {
            if (!fexpr_show_exp_as_power(arg))
                return 0;

            fexpr_view_next(arg);
        }

        return 1;
    }

    return 0;
}

void
fexpr_write_latex_exp(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    fexpr_t arg;
    slong nargs;

    nargs = fexpr_nargs(expr);

    if (nargs == 1)
    {
        fexpr_view_arg(arg, expr, 0);

        if (fexpr_show_exp_as_power(arg))
        {
            calcium_write(out, "e^{");
            fexpr_write_latex(out, arg, flags | FEXPR_LATEX_SMALL);
            calcium_write(out, "}");
        }
        else
        {
            calcium_write(out, "\\exp\\!\\left(");
            fexpr_write_latex(out, arg, flags);
            calcium_write(out, "\\right)");
        }

        return;
    }

    fexpr_write_latex_call(out, expr, flags);
}

void
fexpr_write_latex_factorial(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    if (fexpr_nargs(expr) == 1)
    {
        fexpr_t func, arg;

        fexpr_view_func(func, expr);
        fexpr_view_arg(arg, expr, 0);

        if (fexpr_is_symbol(arg) || (fexpr_is_integer(arg) && !fexpr_is_neg_integer(arg)))
        {
            fexpr_write_latex(out, arg, flags);
        }
        else
        {
            calcium_write(out, "\\left(");
            fexpr_write_latex(out, arg, flags);
            calcium_write(out, "\\right)");
        }

        if (fexpr_is_builtin_symbol(func, FEXPR_DoubleFactorial))
            calcium_write(out, "!!");
        else
            calcium_write(out, "!");
    }
    else
    {
        fexpr_write_latex_call(out, expr, flags);
    }
}

void
fexpr_write_latex_sum_product(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    slong nargs, forexpr_nargs;
    fexpr_t f, forexpr, var, low, high, domain, predicate;
    int have_predicate = 0;
    int have_domain = 0;
    int have_low_high = 0;
    int need_parens;

    nargs = fexpr_nargs(expr);

    if (nargs != 2 && nargs != 3)
    {
        fexpr_write_latex_call(out, expr, flags);
        return;
    }

    fexpr_view_arg(f, expr, 0);
    fexpr_view_arg(forexpr, expr, 1);

    if (nargs == 3)
    {
        fexpr_view_arg(predicate, expr, 2);
        have_predicate = 1;
    }

    forexpr_nargs = fexpr_nargs(forexpr);

    if (forexpr_nargs != 2 && forexpr_nargs != 3)
    {
        fexpr_write_latex_call(out, expr, flags);
        return;
    }

    fexpr_view_arg(var, forexpr, 0);

    if (forexpr_nargs == 2)
    {
        fexpr_view_arg(domain, forexpr, 1);
        have_domain = 1;
    }
    else
    {
        fexpr_view_arg(low, forexpr, 1);
        fexpr_view_arg(high, forexpr, 2);
        have_low_high = 1;
    }

    if (fexpr_is_builtin_call(expr, FEXPR_Sum))
        calcium_write(out, "\\sum");
    else
        calcium_write(out, "\\prod");

    if (have_domain)
    {
        if (have_predicate)
        {
            calcium_write(out, "_{\\textstyle{");
            fexpr_write_latex(out, var, flags);
            calcium_write(out, "  \\in ");
            fexpr_write_latex(out, domain, flags);
            calcium_write(out, " \\atop ");
            fexpr_write_latex(out, predicate, flags | FEXPR_LATEX_SMALL);
            calcium_write(out, "}}");
        }
        else
        {
            calcium_write(out, "_{");
            fexpr_write_latex(out, var, flags);
            calcium_write(out, "  \\in ");
            fexpr_write_latex(out, domain, flags | FEXPR_LATEX_SMALL);
            calcium_write(out, "}");
        }
    }
    else if (have_low_high)
    {
        if (have_predicate)
        {
            calcium_write(out, "_{\\textstyle{");
            fexpr_write_latex(out, var, flags);
            calcium_write(out, "=");
            fexpr_write_latex(out, low, flags | FEXPR_LATEX_SMALL);
            calcium_write(out, " \\atop ");
            fexpr_write_latex(out, predicate, flags | FEXPR_LATEX_SMALL);
            calcium_write(out, "}}^{");
            fexpr_write_latex(out, high, flags | FEXPR_LATEX_SMALL);
            calcium_write(out, "}");
        }
        else
        {
            calcium_write(out, "_{");
            fexpr_write_latex(out, var, flags);
            calcium_write(out, "=");
            fexpr_write_latex(out, low, flags | FEXPR_LATEX_SMALL);
            calcium_write(out, "}^{");
            fexpr_write_latex(out, high, flags | FEXPR_LATEX_SMALL);
            calcium_write(out, "}");
        }
    }

    calcium_write(out, " ");

    need_parens = fexpr_is_builtin_call(f, FEXPR_Add) ||
                  fexpr_is_builtin_call(f, FEXPR_Sub);

    if (need_parens)
        calcium_write(out, "\\left(");

    fexpr_write_latex(out, f, flags);

    if (need_parens)
        calcium_write(out, "\\right)");
}

void
fexpr_write_latex_divsum(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    slong nargs, forexpr_nargs;
    fexpr_t f, forexpr, var, high, predicate;
    slong expected_forargs;
    int have_predicate = 0;
    int need_parens;

    nargs = fexpr_nargs(expr);

    if (nargs != 2 && nargs != 3)
    {
        fexpr_write_latex_call(out, expr, flags);
        return;
    }

    fexpr_view_arg(f, expr, 0);
    fexpr_view_arg(forexpr, expr, 1);

    if (nargs == 3)
    {
        fexpr_view_arg(predicate, expr, 2);
        have_predicate = 1;
    }

    forexpr_nargs = fexpr_nargs(forexpr);

    if (fexpr_is_builtin_call(expr, FEXPR_DivisorSum) || fexpr_is_builtin_call(expr, FEXPR_DivisorProduct))
        expected_forargs = 2;
    else
        expected_forargs = 1;

    if (forexpr_nargs != expected_forargs)
    {
        fexpr_write_latex_call(out, expr, flags);
        return;
    }

    fexpr_view_arg(var, forexpr, 0);

    if (forexpr_nargs == 2)
        fexpr_view_arg(high, forexpr, 1);

    if (fexpr_is_builtin_call(expr, FEXPR_DivisorSum) ||
        fexpr_is_builtin_call(expr, FEXPR_PrimeSum))
        calcium_write(out, "\\sum_{");
    else
        calcium_write(out, "\\prod_{");

    if (fexpr_is_builtin_call(expr, FEXPR_DivisorSum) ||
        fexpr_is_builtin_call(expr, FEXPR_DivisorProduct))
    {
        fexpr_write_latex(out, var, flags | FEXPR_LATEX_SMALL);
        calcium_write(out, " \\mid ");
        fexpr_write_latex(out, high, flags | FEXPR_LATEX_SMALL);
        if (have_predicate)
        {
            calcium_write(out, ",\\, ");
            fexpr_write_latex(out, predicate, flags | FEXPR_LATEX_SMALL);
        }
    }
    else
    {
        if (have_predicate)
            fexpr_write_latex(out, predicate, flags | FEXPR_LATEX_SMALL);
        else
            fexpr_write_latex(out, var, flags | FEXPR_LATEX_SMALL);
    }

    calcium_write(out, "} ");

    need_parens = fexpr_is_builtin_call(f, FEXPR_Add) ||
                  fexpr_is_builtin_call(f, FEXPR_Sub);

    if (need_parens)
        calcium_write(out, "\\left(");

    fexpr_write_latex(out, f, flags);

    if (need_parens)
        calcium_write(out, "\\right)");
}

void
fexpr_write_latex_integral(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    if (fexpr_nargs(expr) == 2)
    {
        fexpr_t f, iter, var, low, high, domain;
        int need_parens;

        fexpr_view_arg(f, expr, 0);
        fexpr_view_arg(iter, expr, 1);

        need_parens = fexpr_is_builtin_call(f, FEXPR_Add) ||
                      fexpr_is_builtin_call(f, FEXPR_Sub);

        if (fexpr_is_builtin_call(iter, FEXPR_For))
        {
            if (fexpr_nargs(iter) == 2)
            {
                fexpr_view_arg(var, iter, 0);
                fexpr_view_arg(domain, iter, 1);

                calcium_write(out, "\\int_{");
                fexpr_write_latex(out, var, flags | FEXPR_LATEX_SMALL);
                calcium_write(out, " \\in ");
                fexpr_write_latex(out, domain, flags | FEXPR_LATEX_SMALL);
                calcium_write(out, "} ");
                if (need_parens)
                    calcium_write(out, "\\left(");
                fexpr_write_latex(out, f, flags);
                if (need_parens)
                    calcium_write(out, "\\right)");
                calcium_write(out, " \\, d");
                fexpr_write_latex(out, var, flags);
                return;
            }

            if (fexpr_nargs(iter) == 3)
            {
                fexpr_view_arg(var, iter, 0);
                fexpr_view_arg(low, iter, 1);
                fexpr_view_arg(high, iter, 2);

                calcium_write(out, "\\int_{");
                fexpr_write_latex(out, low, flags | FEXPR_LATEX_SMALL);
                calcium_write(out, "}^{");
                fexpr_write_latex(out, high, flags | FEXPR_LATEX_SMALL);
                calcium_write(out, "} ");
                if (need_parens)
                    calcium_write(out, "\\left(");
                fexpr_write_latex(out, f, flags);
                if (need_parens)
                    calcium_write(out, "\\right)");
                calcium_write(out, " \\, d");
                fexpr_write_latex(out, var, flags);
                return;
            }
        }
    }

    fexpr_write_latex_call(out, expr, flags);
}

void
fexpr_write_latex_limit(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    fexpr_t op, formula, forexpr, var, point, predicate;
    slong nargs, id;
    int have_predicate = 0;
    int parens;

    nargs = fexpr_nargs(expr);

    if (nargs != 2 && nargs != 3)
    {
        fexpr_write_latex_call(out, expr, flags);
        return;
    }

    fexpr_view_func(op, expr);
    fexpr_view_arg(formula, expr, 0);
    fexpr_view_arg(forexpr, expr, 1);

    if (fexpr_nargs(forexpr) != 2)
    {
        fexpr_write_latex_call(out, expr, flags);
        return;
    }

    fexpr_view_arg(var, forexpr, 0);
    fexpr_view_arg(point, forexpr, 1);

    if (nargs == 3)
    {
        fexpr_view_arg(predicate, expr, 2);
        have_predicate = 1;
    }

    id = FEXPR_BUILTIN_ID(op->data[0]);

    if (id == FEXPR_SequenceLimitInferior)
        calcium_write(out, "\\liminf_{");
    else if (id == FEXPR_SequenceLimitSuperior)
        calcium_write(out, "\\limsup_{");
    else
        calcium_write(out, "\\lim_{");

    fexpr_write_latex(out, var, flags);
    calcium_write(out, " \\to ");

    if (id == FEXPR_LeftLimit || id == FEXPR_RightLimit)
        calcium_write(out, "{");

    fexpr_write_latex(out, point, flags | FEXPR_LATEX_SMALL);

    if (id == FEXPR_LeftLimit)
        calcium_write(out, "}^{-}");
    if (id == FEXPR_RightLimit)
        calcium_write(out, "}^{+}");

    if (have_predicate)
    {
        calcium_write(out, ",\\,");
        fexpr_write_latex(out, predicate, flags | FEXPR_LATEX_SMALL);
    }

    calcium_write(out, "} ");

    parens = (fexpr_is_builtin_call(formula, FEXPR_Add) ||
              fexpr_is_builtin_call(formula, FEXPR_Sub));

    if (parens)
        calcium_write(out, "\\left[");

    fexpr_write_latex(out, formula, flags);

    if (parens)
        calcium_write(out, "\\right]");
}

void
fexpr_write_latex_residue(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    fexpr_t op, formula, forexpr, var, point;
    slong nargs, id;
    int parens;

    nargs = fexpr_nargs(expr);

    if (nargs != 2)
    {
        fexpr_write_latex_call(out, expr, flags);
        return;
    }

    fexpr_view_func(op, expr);
    fexpr_view_arg(formula, expr, 0);
    fexpr_view_arg(forexpr, expr, 1);

    if (fexpr_nargs(forexpr) != 2)
    {
        fexpr_write_latex_call(out, expr, flags);
        return;
    }

    fexpr_view_arg(var, forexpr, 0);
    fexpr_view_arg(point, forexpr, 1);

    id = FEXPR_BUILTIN_ID(op->data[0]);

    if (id == FEXPR_ComplexZeroMultiplicity)
        calcium_write(out, "\\mathop{\\operatorname{ord}}\\limits_{");
    else
        calcium_write(out, "\\mathop{\\operatorname{res}}\\limits_{");

    fexpr_write_latex(out, var, flags);
    calcium_write(out, "=");
    fexpr_write_latex(out, point, flags);
    calcium_write(out, "} ");

    parens = (fexpr_is_builtin_call(formula, FEXPR_Add) ||
              fexpr_is_builtin_call(formula, FEXPR_Sub));

    if (parens)
        calcium_write(out, "\\left[");

    fexpr_write_latex(out, formula, flags);

    if (parens)
        calcium_write(out, "\\right]");
}

void
_fexpr_write_latex_derivative(calcium_stream_t out, const fexpr_t f, const fexpr_t subscript, const fexpr_t order, ulong flags)
{
    if (fexpr_equal_ui(order, 1))
    {
        fexpr_write_latex(out, f, flags);
        calcium_write(out, "'");
    }
    else if (fexpr_equal_ui(order, 2))
    {
        fexpr_write_latex(out, f, flags);
        calcium_write(out, "''");
    }
    else if (fexpr_equal_ui(order, 3))
    {
        fexpr_write_latex(out, f, flags);
        calcium_write(out, "'''");
    }
    else
    {
        calcium_write(out, "{");
        fexpr_write_latex(out, f, flags);
        calcium_write(out, "}^{(");
        fexpr_write_latex(out, order, flags);
        calcium_write(out, ")}");
    }

    if (subscript != NULL)
    {
        calcium_write(out, "_{");
        fexpr_write_latex(out, subscript, flags | FEXPR_LATEX_SMALL);
        calcium_write(out, "}");
    }
}

void
fexpr_write_latex_derivative(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    fexpr_t formula, forexpr, var, point, order;
    slong nargs, forexpr_nargs;
    int parens;
    ulong tmp;

    nargs = fexpr_nargs(expr);

    if (nargs != 2)
    {
        fexpr_write_latex_call(out, expr, flags);
        return;
    }

    fexpr_view_arg(formula, expr, 0);
    fexpr_view_arg(forexpr, expr, 1);

    forexpr_nargs = fexpr_nargs(forexpr);

    if (forexpr_nargs != 2 && forexpr_nargs != 3)
    {
        fexpr_write_latex_call(out, expr, flags);
        return;
    }

    fexpr_view_arg(var, forexpr, 0);
    fexpr_view_arg(point, forexpr, 1);

    if (forexpr_nargs == 3)
    {
        fexpr_view_arg(order, forexpr, 2);
    }
    else
    {
        order->data = &tmp;
        order->alloc = 1;
        fexpr_set_ui(order, 1);
    }

    /* d/dx f(x) -> f'(x) */
    if (fexpr_nargs(formula) == 1)
    {
        fexpr_t f, x;
        fexpr_view_func(f, formula);
        fexpr_view_arg(x, formula, 0);

        if (fexpr_equal(x, var) &&
            fexpr_is_symbol(f) && !fexpr_is_builtin_symbol(f, FEXPR_Exp)
                               && !fexpr_is_builtin_symbol(f, FEXPR_Sqrt))
        {
            _fexpr_write_latex_derivative(out, f, NULL, order, flags);

            calcium_write(out, "\\!\\left(");
            fexpr_write_latex(out, point, flags);
            calcium_write(out, "\\right)");
            return;
        }
    }

    /* d/dx f_n(x) -> f'_n(x) */
    if (fexpr_nargs(formula) == 2)
    {
        fexpr_t f, n, x;
        slong id;
        fexpr_view_func(f, formula);
        fexpr_view_arg(n, formula, 0);
        fexpr_view_arg(x, formula, 1);

        if (fexpr_equal(x, var) && fexpr_is_any_builtin_symbol(f))
        {
            id = FEXPR_BUILTIN_ID(f->data[0]);

            if (fexpr_builtin_table[id].latex_writer == fexpr_write_latex_subscript_call)
            {
                _fexpr_write_latex_derivative(out, f, n, order, flags);

                calcium_write(out, "\\!\\left(");
                fexpr_write_latex(out, point, flags);
                calcium_write(out, "\\right)");
                return;
            }
        }
    }

    if (!fexpr_equal(var, point))
        calcium_write(out, "\\left[");

    if (fexpr_equal_ui(order, 1))
    {
        calcium_write(out, "\\frac{d}{d ");
        fexpr_write_latex(out, var, flags);
        calcium_write(out, "}\\, ");
    }
    else
    {
        calcium_write(out, "\\frac{d^{");
        fexpr_write_latex(out, order, flags);
        calcium_write(out, "}}{{d ");
        fexpr_write_latex(out, var, flags);
        calcium_write(out, "}^{");
        fexpr_write_latex(out, order, flags);
        calcium_write(out, "}}\\, ");
    }

    parens = (fexpr_is_builtin_call(formula, FEXPR_Add) ||
              fexpr_is_builtin_call(formula, FEXPR_Sub));

    if (parens)
        calcium_write(out, "\\left[");

    fexpr_write_latex(out, formula, flags);

    if (parens)
        calcium_write(out, "\\right]");

    if (!fexpr_equal(var, point))
    {
        calcium_write(out, " \\right]_{");
        fexpr_write_latex(out, var, flags);
        calcium_write(out, " = ");
        fexpr_write_latex(out, point, flags | FEXPR_LATEX_SMALL);
        calcium_write(out, "}");
    }
}

void
fexpr_write_latex_setop(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    fexpr_t op, formula, forexpr, var, domain, predicate;
    const char * ops;
    slong nargs, id;
    int have_predicate = 0, have_domain = 0;
    int parens;

    nargs = fexpr_nargs(expr);

    fexpr_view_func(op, expr);
    id = FEXPR_BUILTIN_ID(op->data[0]);

    switch (id)
    {
        case FEXPR_Minimum: ops = "\\min"; break;
        case FEXPR_Maximum: ops = "\\max"; break;
        case FEXPR_ArgMin: ops = "\\operatorname{arg\\,min}"; break;
        case FEXPR_ArgMax: ops = "\\operatorname{arg\\,max}"; break;
        case FEXPR_ArgMinUnique: ops = "\\operatorname{arg\\,min*}"; break;
        case FEXPR_ArgMaxUnique: ops = "\\operatorname{arg\\,max*}"; break;
        case FEXPR_Infimum: ops = "\\operatorname{inf}"; break;
        case FEXPR_Supremum: ops = "\\operatorname{sup}"; break;
        case FEXPR_Zeros: ops = "\\operatorname{zeros}"; break;
        case FEXPR_UniqueZero: ops = "\\operatorname{zero*}"; break;
        case FEXPR_Solutions: ops = "\\operatorname{solutions}"; break;
        case FEXPR_UniqueSolution: ops = "\\operatorname{solution*}"; break;
        case FEXPR_Poles: ops = "\\operatorname{poles}"; break;
        default: ops = "";
    }

    if (nargs == 1)
    {
        fexpr_t arg;
        fexpr_view_arg(arg, expr, 0);

        calcium_write(out, ops);

        if (fexpr_is_builtin_call(arg, FEXPR_Set))
        {
            calcium_write(out, " ");
            fexpr_write_latex(out, arg, flags);
        }
        else
        {
            calcium_write(out, "\\left(");
            fexpr_write_latex(out, arg, flags);
            calcium_write(out, "\\right)");
        }
        return;
    }

    if (nargs != 2 && nargs != 3)
    {
        fexpr_write_latex_call(out, expr, flags);
        return;
    }

    fexpr_view_arg(formula, expr, 0);
    fexpr_view_arg(forexpr, expr, 1);

    if (fexpr_nargs(forexpr) != 1 && fexpr_nargs(forexpr) != 2)
    {
        fexpr_write_latex_call(out, expr, flags);
        return;
    }

    fexpr_view_arg(var, forexpr, 0);

    if (fexpr_nargs(forexpr) == 1)
    {
        have_domain = 0;
    }
    else
    {
        fexpr_view_arg(domain, forexpr, 1);
        have_domain = 1;
    }

    if (nargs == 3)
    {
        fexpr_view_arg(predicate, expr, 2);
        have_predicate = 1;
    }


    calcium_write(out, "\\mathop{");
    calcium_write(out, ops);
    calcium_write(out, "\\,}\\limits_{");

    if (have_domain)
    {
        fexpr_write_latex(out, var, flags | FEXPR_LATEX_SMALL);
        calcium_write(out, " \\in ");
        fexpr_write_latex(out, domain, flags | FEXPR_LATEX_SMALL);
    }

    if (have_predicate)
    {
        if (have_domain)
            calcium_write(out, ",\\,");

        fexpr_write_latex(out, predicate, flags | FEXPR_LATEX_SMALL);
    }

    calcium_write(out, "} ");

    parens = (fexpr_is_builtin_call(formula, FEXPR_Add) ||
              fexpr_is_builtin_call(formula, FEXPR_Sub) ||
              fexpr_is_builtin_call(formula, FEXPR_Neg) ||
              fexpr_is_builtin_call(formula, FEXPR_Sum) ||
              fexpr_is_builtin_call(formula, FEXPR_Product) ||
              fexpr_is_builtin_call(formula, FEXPR_Integral));

    if (parens)
        calcium_write(out, "\\left[");

    fexpr_write_latex(out, formula, flags);

    if (parens)
        calcium_write(out, "\\right]");
}

void
fexpr_write_latex_logic(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    fexpr_t arg;
    slong i, nargs;

    nargs = fexpr_nargs(expr);

    if (fexpr_is_builtin_call(expr, FEXPR_Not) && nargs == 1)
    {
        fexpr_view_arg(arg, expr, 0);

        if (flags & FEXPR_LATEX_LOGIC)
            calcium_write(out, "\\neg ");
        else
            calcium_write(out, "\\operatorname{not} ");

        if (fexpr_is_atom(arg))
        {
            fexpr_write_latex(out, arg, flags);
        }
        else
        {
            if (!(flags & FEXPR_LATEX_LOGIC))
                calcium_write(out, "\\,");

            calcium_write(out, "\\left(");
            fexpr_write_latex(out, arg, flags);
            calcium_write(out, "\\right)");
        }

        return;
    }

    if (fexpr_is_builtin_call(expr, FEXPR_Or) && nargs >= 1)
    {
        fexpr_view_arg(arg, expr, 0);

        for (i = 0; i < nargs; i++)
        {
            if (fexpr_is_builtin_call(arg, FEXPR_And) ||
                fexpr_is_builtin_call(arg, FEXPR_Or) ||
                fexpr_is_builtin_call(arg, FEXPR_Not))
            {
                calcium_write(out, "\\left(");
                fexpr_write_latex(out, arg, flags);
                calcium_write(out, "\\right)");
            }
            else
            {
                fexpr_write_latex(out, arg, flags);
            }

            if (i < nargs - 1)
            {
                if (flags & FEXPR_LATEX_LOGIC)
                    calcium_write(out, " \\,\\lor\\, ");
                else
                    calcium_write(out, " \\;\\mathbin{\\operatorname{or}}\\; ");

                fexpr_view_next(arg);
            }
        }

        return;
    }

    if (fexpr_is_builtin_call(expr, FEXPR_And) && nargs >= 1)
    {
        fexpr_view_arg(arg, expr, 0);

        for (i = 0; i < nargs; i++)
        {
            if (fexpr_is_builtin_call(arg, FEXPR_And) ||
                fexpr_is_builtin_call(arg, FEXPR_Or) ||
                fexpr_is_builtin_call(arg, FEXPR_All) ||
                fexpr_is_builtin_call(arg, FEXPR_Exists))
            {
                calcium_write(out, "\\left(");
                fexpr_write_latex(out, arg, flags);
                calcium_write(out, "\\right)");
            }
            else
            {
                fexpr_write_latex(out, arg, flags);
            }

            if (i < nargs - 1)
            {
                if (flags & FEXPR_LATEX_LOGIC)
                    calcium_write(out, " \\,\\land\\, ");
                else if (flags & FEXPR_LATEX_SMALL)  /* see see ff190c */
                    calcium_write(out, " ,\\, ");
                else
                    calcium_write(out, " \\;\\mathbin{\\operatorname{and}}\\; ");

                fexpr_view_next(arg);
            }
        }

        return;
    }

    if (fexpr_is_builtin_call(expr, FEXPR_Implies) && nargs == 2)
    {
        fexpr_t arg1, arg2;
        fexpr_view_arg(arg1, expr, 0);
        fexpr_view_arg(arg2, expr, 1);

        if (!(fexpr_is_atom(arg1) || fexpr_is_builtin_call(arg1, FEXPR_Element)))
        {
            calcium_write(out, "\\left(");
            fexpr_write_latex(out, arg1, flags);
            calcium_write(out, "\\right)");
        }
        else
        {
            fexpr_write_latex(out, arg1, flags);
        }

        calcium_write(out, " \\;\\implies\\; ");

        if (!(fexpr_is_atom(arg2) || fexpr_is_builtin_call(arg2, FEXPR_Element)))
        {
            calcium_write(out, "\\left(");
            fexpr_write_latex(out, arg2, flags);
            calcium_write(out, "\\right)");
        }
        else
        {
            fexpr_write_latex(out, arg2, flags);
        }

        return;
    }

    if (fexpr_is_builtin_call(expr, FEXPR_Equivalent) && nargs >= 1)
    {
        fexpr_view_func(arg, expr);

        for (i = 0; i < nargs; i++)
        {
            fexpr_view_next(arg);

            if (!fexpr_is_atom(arg))
                calcium_write(out, "\\left(");
            fexpr_write_latex(out, arg, flags);
            if (!fexpr_is_atom(arg))
                calcium_write(out, "\\right)");

            if (i < nargs - 1)
            {
                calcium_write(out, " \\iff ");
            }
        }

        return;
    }

    if ((fexpr_is_builtin_call(expr, FEXPR_All) || fexpr_is_builtin_call(expr, FEXPR_Exists))
        && (nargs == 2 || nargs == 3))
    {
        fexpr_t func, forarg, var, domain, condition;

        fexpr_view_arg(func, expr, 0);
        fexpr_view_arg(forarg, expr, 1);
        if (nargs == 3)
            fexpr_view_arg(condition, expr, 2);

        if (fexpr_nargs(forarg) == 1 || fexpr_nargs(forarg) == 2)
        {
            int have_domain = (fexpr_nargs(forarg) == 2);

            fexpr_view_arg(var, forarg, 0);
            if (have_domain)
                fexpr_view_arg(domain, forarg, 1);

            if (flags & FEXPR_LATEX_LOGIC)
            {
                if (fexpr_is_builtin_call(expr, FEXPR_All))
                    calcium_write(out, "\\forall ");
                else
                    calcium_write(out, "\\exists ");

                fexpr_write_latex(out, var, flags);

                if (have_domain)
                {
                    calcium_write(out, " \\in ");
                    fexpr_write_latex(out, domain, flags);
                }

                if (nargs == 3)
                {
                    calcium_write(out, ", \\,");
                    fexpr_write_latex(out, condition, flags);
                }

                calcium_write(out, " : \\, ");
                fexpr_write_latex(out, func, flags);
            }
            else
            {
                fexpr_write_latex(out, func, flags);

                if (fexpr_is_builtin_call(expr, FEXPR_All))
                    calcium_write(out, " \\;\\text{ for all } ");
                else
                    calcium_write(out, " \\;\\text{ for some } ");

                fexpr_write_latex(out, var, flags);

                if (have_domain)
                {
                    calcium_write(out, " \\in ");
                    fexpr_write_latex(out, domain, flags);
                }

                if (nargs == 3)
                {
                    calcium_write(out, " \\text{ with } ");
                    fexpr_write_latex(out, condition, flags);
                }
            }

            return;
        }
    }

    if (fexpr_is_builtin_call(expr, FEXPR_Logic) && nargs == 1)
    {
        fexpr_t arg;
        fexpr_view_arg(arg, expr, 0);
        fexpr_write_latex(out, arg, flags | FEXPR_LATEX_LOGIC);
        return;
    }

    /* todo: move this */
    if (fexpr_is_builtin_call(expr, FEXPR_CongruentMod) && nargs == 3)
    {
        fexpr_t arg;
        fexpr_view_arg(arg, expr, 0);
        fexpr_write_latex(out, arg, flags);
        calcium_write(out, " \\equiv ");
        fexpr_view_next(arg);
        fexpr_write_latex(out, arg, flags);
        calcium_write(out, " \\pmod {");
        fexpr_view_next(arg);
        fexpr_write_latex(out, arg, flags);
        calcium_write(out, " }");
        return;
    }

    fexpr_write_latex_call(out, expr, flags);
}

void
fexpr_write_latex_cases(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    fexpr_t arg, value, condition;
    slong i, nargs;

    calcium_write(out, "\\begin{cases} ");

    nargs = fexpr_nargs(expr);
    fexpr_view_arg(arg, expr, 0);

    for (i = 0; i < nargs; i++)
    {
        if (fexpr_nargs(arg) != 2)
            continue;

        fexpr_view_arg(value, arg, 0);
        fexpr_view_arg(condition, arg, 1);

        fexpr_write_latex(out, value, flags /* | FEXPR_LATEX_SMALL */);
        calcium_write(out, ", & ");

        if (fexpr_is_builtin_symbol(condition, FEXPR_Otherwise))
            calcium_write(out, "\\text{otherwise}");
        else
            fexpr_write_latex(out, condition, flags /* | FEXPR_LATEX_SMALL */);

        calcium_write(out, "\\\\");

        if (i < nargs - 1)
            fexpr_view_next(arg);
    }

    calcium_write(out, " \\end{cases}");
}

void
fexpr_write_latex_simple(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    slong i;
    const char * a;
    const char * b;
    fexpr_t func, arg;

    if (fexpr_nargs(expr) != 1 || !fexpr_is_any_builtin_call(expr))
    {
        fexpr_write_latex_call(out, expr, flags);
        return;
    }

    fexpr_view_func(func, expr);
    fexpr_view_arg(arg, expr, 0);

    i = FEXPR_BUILTIN_ID(func->data[0]);

    switch (i)
    {
        case FEXPR_Sqrt:
            a = "\\sqrt{";
            b = "}";
            break;
        case FEXPR_Conjugate:
            a = "\\overline{";
            b = "}";
            break;
        case FEXPR_Cardinality:
        case FEXPR_Length:
            a = "\\# ";
            b = "";
            break;
        case FEXPR_Abs:
        case FEXPR_RealAbs:
            a = "\\left|";
            b = "\\right|";
            break;
        case FEXPR_Floor:
            a = "\\left\\lfloor ";
            b = " \\right\\rfloor";
            break;
        case FEXPR_Ceil:
            a = "\\left\\lceil ";
            b = "\\right\\rceil";
            break;
        case FEXPR_Parentheses:
            a = "\\left(";
            b = "\\right)";
            break;
        case FEXPR_Brackets:
            a = "\\left[";
            b = "\\right]";
            break;
        case FEXPR_Braces:
            a = "\\left\\{";
            b = "\\right\\}";
            break;
        case FEXPR_AngleBrackets:
            a = "\\left\\langle ";
            b = "\\right\\rangle";
            break;
        case FEXPR_IsEven:
            a = "";
            b = " \\text{ even}";
            break;
        case FEXPR_IsOdd:
            a = "";
            b = " \\text{ odd}";
            break;
        case FEXPR_IsPrime:
            a = "";
            b = " \\text{ prime}";
            break;
        case FEXPR_Enclosure:
            a = "\\left(\\in ";
            b = " \\right)";
            break;
        case FEXPR_Guess:
            a = "\\left(?\\; ";
            b = " \\right)";
            break;
        case FEXPR_Approximation:
            a = "\\left(\\approx ";
            b = " \\right)";
            break;
        default:
            fexpr_write_latex_call(out, expr, flags);
            return;
    }

    calcium_write(out, a);
    fexpr_write_latex(out, arg, flags);
    calcium_write(out, b);
}

void
_fexpr_write_latex_simple2(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    slong i;
    const char * a;
    const char * b;
    const char * c;
    fexpr_t func, arg1, arg2;

    if (fexpr_nargs(expr) != 2 || !fexpr_is_any_builtin_call(expr))
    {
        fexpr_write_latex_call(out, expr, flags);
        return;
    }

    fexpr_view_func(func, expr);
    fexpr_view_arg(arg1, expr, 0);
    fexpr_view_arg(arg2, expr, 1);

    i = FEXPR_BUILTIN_ID(func->data[0]);

    switch (i)
    {
        case FEXPR_RisingFactorial:
            a = "\\left(";
            b = "\\right)_{";
            c = "}";
            break;
        case FEXPR_FallingFactorial:
            a = "\\left(";
            b = "\\right)^{\\underline{";
            c = "}}";
            break;
        case FEXPR_Binomial:
            a = "{";
            b = " \\choose ";
            c = "}";
            break;
        case FEXPR_StirlingCycle:
            a = "\\left[{";
            b = " \\atop ";
            c = "}\\right]";
            break;
        case FEXPR_StirlingS1:
            a = "s\\!\\left(";
            b = ", ";
            c = "\\right)";
            break;
        case FEXPR_StirlingS2:
            a = "\\left\\{{";
            b = " \\atop ";
            c = "}\\right\\}";
            break;
        case FEXPR_LegendreSymbol:
        case FEXPR_JacobiSymbol:
        case FEXPR_KroneckerSymbol:
            a = "\\left(\\frac{";
            b = "}{";
            c = "}\\right)";
            break;
        case FEXPR_Interval:
            a = "\\left[";
            b = ", ";
            c = "\\right]";
            break;
        case FEXPR_OpenInterval:
            a = "\\left(";
            b = ", ";
            c = "\\right)";
            break;
        case FEXPR_ClosedOpenInterval:
            a = "\\left[";
            b = ", ";
            c = "\\right)";
            break;
        case FEXPR_OpenClosedInterval:
            a = "\\left(";
            b = ", ";
            c = "\\right]";
            break;
        case FEXPR_RealBall:
            a = "\\left[";
            b = " \\pm ";
            c = "\\right]";
            break;
        case FEXPR_OpenRealBall:
            a = "\\left(";
            b = " \\pm ";
            c = "\\right)";
            break;
        case FEXPR_KroneckerDelta:
            a = "\\delta_{(";
            b = ",";
            c = ")}";
            break;
        case FEXPR_Subscript:
        case FEXPR_Item:
            a = "{";
            b = "}_{";
            c = "}";
            break;
        default:
            fexpr_write_latex_call(out, expr, flags);
            return;
    }

    calcium_write(out, a);
    fexpr_write_latex(out, arg1, flags);
    calcium_write(out, b);
    fexpr_write_latex(out, arg2, flags);
    calcium_write(out, c);
}

void
fexpr_write_latex_simple2(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    _fexpr_write_latex_simple2(out, expr, flags);
}

void
fexpr_write_latex_simple2_small(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    _fexpr_write_latex_simple2(out, expr, flags | FEXPR_LATEX_SMALL);
}

void
fexpr_write_latex_collection(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    fexpr_t arg;
    slong nargs;

    nargs = fexpr_nargs(expr);

    /* Set comprehension */
    if (fexpr_is_builtin_call(expr, FEXPR_Set) && (nargs == 2 || nargs == 3))
    {
        slong for_nargs;

        fexpr_view_arg(arg, expr, 1);
        for_nargs = fexpr_nargs(arg);

        if (fexpr_is_builtin_call(arg, FEXPR_For) && (for_nargs == 2 || for_nargs == 1))
        {
            fexpr_t func, var, domain, predicate;

            fexpr_view_arg(func, expr, 0);
            fexpr_view_arg(var, arg, 0);
            if (for_nargs == 2)
                fexpr_view_arg(domain, arg, 1);

            calcium_write(out, "\\left\\{ ");

            fexpr_write_latex(out, func, flags);

            if (for_nargs == 2 || nargs == 3)
                calcium_write(out, " : ");

            if (for_nargs == 2)
            {
                fexpr_write_latex(out, var, flags);
                calcium_write(out, " \\in ");
                fexpr_write_latex(out, domain, flags);
            }

            if (nargs == 3)
            {
                fexpr_view_arg(predicate, expr, 2);
                if (for_nargs == 2)
                    calcium_write(out, "\\,\\mathbin{\\operatorname{and}}\\, ");
                fexpr_write_latex(out, predicate, flags);
            }

            calcium_write(out, " \\right\\}");

            return;
        }
    }

    if (fexpr_is_builtin_call(expr, FEXPR_Set))
        calcium_write(out, "\\left\\{");
    else if (fexpr_is_builtin_call(expr, FEXPR_Tuple))
        calcium_write(out, "\\left(");
    else if (fexpr_is_builtin_call(expr, FEXPR_List))
        calcium_write(out, "\\left[");
    else if (fexpr_is_builtin_call(expr, FEXPR_Lattice))
        calcium_write(out, "\\Lambda_{(");

    {
        fexpr_t func, var, a, b;
        slong for_nargs;

        if (nargs == 2)
        {
            fexpr_view_arg(arg, expr, 1);
            for_nargs = fexpr_nargs(arg);
        }

        /* {f(n) for n=a..b} -- todo: special cases for infinity, explicit integer a and b */
        if (nargs == 2 && fexpr_is_builtin_call(arg, FEXPR_For) && for_nargs == 3)
        {
            fexpr_t first, second, last, t, a1;
            fmpz_t n;

            fexpr_view_arg(func, expr, 0);
            fexpr_view_arg(var, arg, 0);
            fexpr_view_arg(a, arg, 1);
            fexpr_view_arg(b, arg, 2);

            fexpr_init(first);
            fexpr_init(second);
            fexpr_init(last);
            fexpr_init(a1);
            fexpr_init(t);
            fmpz_init(n);

            if (fexpr_is_integer(a))
            {
                fexpr_get_fmpz(n, a);
                fmpz_add_ui(n, n, 1);
                fexpr_set_fmpz(a1, n);
            }
            else
            {
                fexpr_set_ui(t, 1);
                fexpr_add(a1, a, t);
            }

            fexpr_replace(first, func, var, a);
            fexpr_replace(second, func, var, a1);
            fexpr_replace(last, func, var, b);

            fexpr_write_latex(out, first, flags);
            calcium_write(out, ", ");
            fexpr_write_latex(out, second, flags);
            calcium_write(out, ", \\ldots, ");
            fexpr_write_latex(out, last, flags);

            fexpr_clear(first);
            fexpr_clear(second);
            fexpr_clear(last);
            fexpr_clear(a1);
            fexpr_clear(t);
            fmpz_clear(n);
        }
        else if (nargs >= 1)
        {
            slong i;

            fexpr_view_arg(arg, expr, 0);

            for (i = 0; i < nargs; i++)
            {
                fexpr_write_latex(out, arg, flags);

                if (i < nargs - 1)
                {
                    calcium_write(out, ", ");
                    fexpr_view_next(arg);
                }
            }
        }
    }

    if (fexpr_is_builtin_call(expr, FEXPR_Set))
        calcium_write(out, "\\right\\}");
    else if (fexpr_is_builtin_call(expr, FEXPR_Tuple))
        calcium_write(out, "\\right)");
    else if (fexpr_is_builtin_call(expr, FEXPR_List))
        calcium_write(out, "\\right]");
    else if (fexpr_is_builtin_call(expr, FEXPR_Lattice))
        calcium_write(out, ")}");
}

void
fexpr_write_latex_range(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    fexpr_t a, b;

    if (fexpr_is_builtin_call(expr, FEXPR_IntegersGreaterEqual) && fexpr_nargs(expr) == 1)
    {
        fexpr_view_arg(a, expr, 0);
        calcium_write(out, "\\mathbb{Z}_{\\ge ");
        fexpr_write_latex(out, a, flags | FEXPR_LATEX_SMALL);
        calcium_write(out, "}");
        return;
    }

    if (fexpr_is_builtin_call(expr, FEXPR_IntegersLessEqual) && fexpr_nargs(expr) == 1)
    {
        fexpr_view_arg(a, expr, 0);
        if (fexpr_is_integer(a))
        {
            fmpz_t n;
            fmpz_init(n);
            fexpr_get_fmpz(n, a);
            calcium_write(out, "\\{");
            calcium_write_fmpz(out, n);
            calcium_write(out, ", ");
            fmpz_sub_ui(n, n, 1);
            calcium_write_fmpz(out, n);
            calcium_write(out, ", \\ldots\\}");
            fmpz_clear(n);
        }
        else
        {
            calcium_write(out, "\\mathbb{Z}_{\\le ");
            fexpr_write_latex(out, a, flags | FEXPR_LATEX_SMALL);
            calcium_write(out, "}");
        }
        return;
    }

    if (fexpr_is_builtin_call(expr, FEXPR_Range) && fexpr_nargs(expr) == 2)
    {
        fexpr_view_arg(a, expr, 0);
        fexpr_view_arg(b, expr, 1);

        if (fexpr_is_integer(a))
        {
            fmpz_t n;
            fmpz_init(n);
            fexpr_get_fmpz(n, a);
            calcium_write(out, "\\{");
            calcium_write_fmpz(out, n);
            calcium_write(out, ", ");
            fmpz_add_ui(n, n, 1);
            calcium_write_fmpz(out, n);
            calcium_write(out, ", \\ldots, ");
            fexpr_write_latex(out, b, flags);
            calcium_write(out, "\\}");
            fmpz_clear(n);
        }
        else
        {
            calcium_write(out, "\\{");
            fexpr_write_latex(out, a, flags);
            calcium_write(out, ", ");
            fexpr_write_latex(out, a, flags);
            calcium_write(out, " + 1, \\ldots, ");
            fexpr_write_latex(out, b, flags);
            calcium_write(out, "\\}");
        }
        return;
    }

    fexpr_write_latex_call(out, expr, flags);
}

void
fexpr_write_latex_matrix(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    fexpr_t arg, row, elem;
    slong i, j, nargs, nrows, ncols;

    nargs = fexpr_nargs(expr);

    if (fexpr_is_builtin_call(expr, FEXPR_RowMatrix) || fexpr_is_builtin_call(expr, FEXPR_ColumnMatrix))
    {
        int isrow = fexpr_is_builtin_call(expr, FEXPR_RowMatrix);

        calcium_write(out, "\\displaystyle{\\begin{pmatrix}");

        if (nargs > 0)
        {
            fexpr_view_arg(elem, expr, 0);

            for (i = 0; i < nargs; i++)
            {
                fexpr_write_latex(out, elem, flags);

                if (i < nargs - 1)
                {
                    if (isrow)
                        calcium_write(out, " & ");
                    else
                        calcium_write(out, " \\\\ ");

                    fexpr_view_next(elem);
                }
            }
        }

        calcium_write(out, "\\end{pmatrix}}");
        return;
    }

    if (fexpr_is_builtin_call(expr, FEXPR_DiagonalMatrix))
    {
        calcium_write(out, "\\displaystyle{\\begin{pmatrix}");

        if (nargs > 0)
        {
            fexpr_view_arg(elem, expr, 0);

            for (i = 0; i < nargs; i++)
            {
                for (j = 0; j < i; j++)
                    calcium_write(out, " & ");

                fexpr_write_latex(out, elem, flags);

                for (j = i + 1; j < nargs; j++)
                    calcium_write(out, " & ");

                if (i < nargs - 1)
                {
                    calcium_write(out, " \\\\ ");
                    fexpr_view_next(elem);
                }
            }
        }

        calcium_write(out, "\\end{pmatrix}}");
        return;
    }

    if (fexpr_is_builtin_call(expr, FEXPR_Matrix2x2) && nargs == 4)
    {
        calcium_write(out, "\\displaystyle{\\begin{pmatrix}");
        fexpr_view_arg(elem, expr, 0);
        fexpr_write_latex(out, elem, flags);
        calcium_write(out, " & ");
        fexpr_view_next(elem);
        fexpr_write_latex(out, elem, flags);
        calcium_write(out, " \\\\ ");
        fexpr_view_next(elem);
        fexpr_write_latex(out, elem, flags);
        calcium_write(out, " & ");
        fexpr_view_next(elem);
        fexpr_write_latex(out, elem, flags);
        calcium_write(out, "\\end{pmatrix}}");
        return;
    }

    if (fexpr_is_builtin_call(expr, FEXPR_Matrix) && nargs == 3)
    {
        fexpr_t for1, for2, f1, f2, i, a, b, j, c, d;

        fexpr_view_arg(for1, expr, 1);
        fexpr_view_arg(for2, expr, 2);

        if (fexpr_view_call3(f1, i, a, b, for1) &&
            fexpr_view_call3(f2, j, c, d, for2) &&
            fexpr_is_builtin_symbol(f1, FEXPR_For) &&
            fexpr_is_builtin_symbol(f2, FEXPR_For))
        {
            fexpr_t a1, c1, x;
            fmpz_t n;

            fmpz_init(n);
            fexpr_init(a1);
            fexpr_init(c1);
            fexpr_init(x);

            fexpr_view_arg(arg, expr, 0);

            /* a1 = a + 1 */
            if (fexpr_is_integer(a))
            {
                fexpr_get_fmpz(n, a);
                fmpz_add_ui(n, n, 1);
                fexpr_set_fmpz(a1, n);
            }
            else
            {
                fexpr_set_ui(x, 1);
                fexpr_add(a1, a, x);
            }

            /* c1 = c + 1 */
            if (fexpr_is_integer(c))
            {
                fexpr_get_fmpz(n, c);
                fmpz_add_ui(n, n, 1);
                fexpr_set_fmpz(c1, n);
            }
            else
            {
                fexpr_set_ui(x, 1);
                fexpr_add(c1, c, x);
            }

            calcium_write(out, "\\displaystyle{\\begin{pmatrix} ");

            fexpr_replace2(x, arg, i, a, j, c);
            fexpr_write_latex(out, x, flags);
            calcium_write(out, " & ");

            fexpr_replace2(x, arg, i, a, j, c1);
            fexpr_write_latex(out, x, flags);
            calcium_write(out, " & \\cdots & ");

            fexpr_replace2(x, arg, i, a, j, d);
            fexpr_write_latex(out, x, flags);
            calcium_write(out, " \\\\ ");

            fexpr_replace2(x, arg, i, a1, j, c);
            fexpr_write_latex(out, x, flags);
            calcium_write(out, " & ");

            fexpr_replace2(x, arg, i, a1, j, c1);
            fexpr_write_latex(out, x, flags);
            calcium_write(out, " & \\cdots & ");

            fexpr_replace2(x, arg, i, a1, j, d);
            fexpr_write_latex(out, x, flags);
            calcium_write(out, " \\\\ ");

            calcium_write(out, "\\vdots & \\vdots & \\ddots & \\vdots \\\\ ");

            fexpr_replace2(x, arg, i, b, j, c);
            fexpr_write_latex(out, x, flags);
            calcium_write(out, " & ");

            fexpr_replace2(x, arg, i, b, j, c1);
            fexpr_write_latex(out, x, flags);
            calcium_write(out, " & \\cdots & ");

            fexpr_replace2(x, arg, i, b, j, d);
            fexpr_write_latex(out, x, flags);
            calcium_write(out, " \\end{pmatrix}}");

            fmpz_clear(n);
            fexpr_clear(a1);
            fexpr_clear(c1);
            fexpr_clear(x);

            return;
        }
    }

    if (fexpr_is_builtin_call(expr, FEXPR_Matrix) && nargs >= 1)
    {
        fexpr_view_arg(arg, expr, 0);

        /* Todo: handle columns too */
        if (fexpr_is_builtin_call(arg, FEXPR_Row))
        {
            nrows = fexpr_nargs(expr);

            calcium_write(out, "\\displaystyle{\\begin{pmatrix}");

            fexpr_view_arg(row, expr, 0);
            for (i = 0; i < nrows; i++)
            {
                ncols = fexpr_nargs(row);

                if (ncols >= 0)
                {
                    fexpr_view_arg(elem, row, 0);

                    for (j = 0; j < ncols; j++)
                    {
                        fexpr_write_latex(out, elem, flags);

                        if (j < ncols - 1)
                        {
                            calcium_write(out, " & ");
                            fexpr_view_next(elem);
                        }
                    }
                }

                if (i < nrows - 1)
                {
                    calcium_write(out, " \\\\");
                    fexpr_view_next(row);
                }
            }

            calcium_write(out, "\\end{pmatrix}}");
            return;
        }
    }

    if (nargs == 1)
    {
        fexpr_view_arg(arg, expr, 0);

        if (fexpr_is_builtin_call(arg, FEXPR_Tuple) ||
            fexpr_is_builtin_call(arg, FEXPR_List))
        {
            nrows = fexpr_nargs(arg);

            calcium_write(out, "\\displaystyle{\\begin{pmatrix}");

            fexpr_view_arg(row, arg, 0);
            for (i = 0; i < nrows; i++)
            {
                ncols = fexpr_nargs(row);

                if (ncols >= 0)
                {
                    fexpr_view_arg(elem, row, 0);

                    for (j = 0; j < ncols; j++)
                    {
                        fexpr_write_latex(out, elem, flags);

                        if (j < ncols - 1)
                        {
                            calcium_write(out, " & ");
                            fexpr_view_next(elem);
                        }
                    }
                }

                if (i < nrows - 1)
                {
                    calcium_write(out, " \\\\");
                    fexpr_view_next(row);
                }
            }

            calcium_write(out, "\\end{pmatrix}}");
            return;
        }
    }

    fexpr_write_latex_call(out, expr, flags);
}

void
fexpr_write_latex_decimal(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    fexpr_t arg;
    char * s;
    slong i, len;

    if (fexpr_nargs(expr) == 1)
    {
        fexpr_view_arg(arg, expr, 0);

        if (fexpr_is_string(arg))
        {
            s = fexpr_get_string(arg);
            len = strlen(s);

            for (i = 0; i < len; i++)
            {
                if (s[i] == 'e' || s[i] == 'E')
                {
                    s[i] = '\0';
                    calcium_write(out, s);
                    calcium_write(out, " \\cdot 10^{");
                    calcium_write(out, s + i + 1);
                    calcium_write(out, "}");
                    flint_free(s);
                    return;
                }
            }

            calcium_write(out, s);
            flint_free(s);
            return;
        }
    }

    fexpr_write_latex_call(out, expr, flags);
}

/* Write (x) */
void
_fexpr_write_latex_call1(calcium_stream_t out, const fexpr_t x, ulong flags)
{
    if (fexpr_is_atom(x))
    {
        calcium_write(out, "(");
        fexpr_write_latex(out, x, flags);
        calcium_write(out, ")");
    }
    else
    {
        calcium_write(out, "\\!\\left(");
        fexpr_write_latex(out, x, flags);
        calcium_write(out, "\\right)");
    }
}

static void
_write_poly(calcium_stream_t out, const fexpr_t pol, ulong flags)
{
    fexpr_t c;
    slong i, d;

    d = fexpr_nargs(pol) - 1;

    if (d < 0)
    {
        calcium_write(out, "0");
        return;
    }

    for (i = d; i >= 0; i--)
    {
        fexpr_view_arg(c, pol, i);

        if (fexpr_equal_si(c, 0))
            continue;

        if (fexpr_equal_si(c, 1))
        {
            if (i == 0)
                calcium_write(out, "+1");
            else if (i != d)
                calcium_write(out, "+");
        }
        else if (fexpr_equal_si(c, -1))
        {
            if (i == 0)
                calcium_write(out, "-1");
            else
                calcium_write(out, "-");
        }
        else
        {
            if (fexpr_need_parens_in_mul(c, 0))
            {
                if (i != d)
                    calcium_write(out, "+ ");
                calcium_write(out, "\\left(");
                fexpr_write_latex(out, c, flags);
                calcium_write(out, "\\right)");
            }
            else
            {
                if (fexpr_is_integer(c))
                {
                    if (!(fexpr_can_extract_leading_sign(c) || i == d))
                        calcium_write(out, "+");

                    fexpr_write_latex(out, c, flags);
                }
                else
                {
                    char * s = fexpr_get_str_latex(c, flags);

                    if (s[0] == '+' || s[0] == '-')
                    {
                        calcium_write(out, s);
                    }
                    else
                    {
                        if (i != d)
                            calcium_write(out, " + ");
                        calcium_write(out, s);
                    }

                    flint_free(s);
                }
            }
        }

        if (i == 1)
        {
            calcium_write(out, " x");
        }
        else if (i >= 2)
        {
            calcium_write(out, " x^{");
            calcium_write_si(out, i);
            calcium_write(out, "}");
        }
    }
}


void
fexpr_write_latex_misc_special(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    if (fexpr_is_builtin_call(expr, FEXPR_Polynomial) && fexpr_nargs(expr) == 1)
    {
        fexpr_t pol;

        fexpr_view_arg(pol, expr, 0);

        _write_poly(out, pol, flags);

        return;
    }

    if (fexpr_is_builtin_call(expr, FEXPR_PolynomialRootNearest) && fexpr_nargs(expr) == 2)
    {
        fexpr_t pol, point;

        fexpr_view_arg(pol, expr, 0);
        fexpr_view_arg(point, expr, 1);

/*
        calcium_write(out, "\\left(x \\approx ");
        fexpr_write_latex(out, point, flags);
        calcium_write(out, ", \\,");
        _write_poly(out, pol, flags);
        calcium_write(out, "= 0\\right)");
        return;
*/

        calcium_write(out, "\\left(\\text{Root }\\, ");
        calcium_write(out, "x \\approx {");
        fexpr_write_latex(out, point, flags);
        calcium_write(out, "} \\;\\text{ of } \\;{");
        _write_poly(out, pol, flags);
        calcium_write(out, "}\\right)");

        return;
    }

    if (fexpr_is_builtin_call(expr, FEXPR_PolynomialRootIndexed) && fexpr_nargs(expr) == 2)
    {
        fexpr_t pol, point;

        fexpr_view_arg(pol, expr, 0);
        fexpr_view_arg(point, expr, 1);

        calcium_write(out, "\\left(\\text{Root \\#}");
        fexpr_write_latex(out, point, flags);
        calcium_write(out, " \\text{ of }\\, {");
        _write_poly(out, pol, flags);
        calcium_write(out, "}\\right)");

        return;
    }

/*
    if (fexpr_is_builtin_call(expr, FEXPR_AlgebraicNumberSerialized) && fexpr_nargs(expr) == 2)
    {
        fexpr_t pol, index;

        fexpr_view_arg(pol, expr, 0);
        fexpr_view_arg(index, expr, 1);

        calcium_write(out, "\\left(\\text{Algebraic } ");
        fexpr_write_latex(out, pol, flags);
        calcium_write(out, "_{");
        fexpr_write_latex(out, index, flags);
        calcium_write(out, "}\\right)");

        return;
    }
*/

    if (fexpr_is_builtin_call(expr, FEXPR_Call) && fexpr_nargs(expr) == 2)
    {
        fexpr_t f, x;

        fexpr_view_arg(f, expr, 0);
        fexpr_view_arg(x, expr, 1);

        fexpr_write_latex(out, f, flags);
        _fexpr_write_latex_call1(out, x, flags);
        return;
    }

    if (fexpr_is_builtin_call(expr, FEXPR_CallIndeterminate) && fexpr_nargs(expr) == 3)
    {
        fexpr_t f, x, v;

        fexpr_view_arg(f, expr, 0);
        fexpr_view_arg(x, expr, 1);
        fexpr_view_arg(v, expr, 2);

        fexpr_write_latex(out, f, flags);
        _fexpr_write_latex_call1(out, v, flags);
        return;
    }

    if (fexpr_is_builtin_call(expr, FEXPR_IndefiniteIntegralEqual) && fexpr_nargs(expr) == 3)
    {
        fexpr_t fx, gx, x;

        fexpr_view_arg(fx, expr, 0);
        fexpr_view_arg(gx, expr, 1);
        fexpr_view_arg(x, expr, 2);

        calcium_write(out, "\\int ");
        fexpr_write_latex(out, fx, flags);
        calcium_write(out, " \\, d");
        fexpr_write_latex(out, x, flags);
        calcium_write(out, " = ");
        fexpr_write_latex(out, gx, flags);
        calcium_write(out, " + \\mathcal{C}");
        return;
    }

    if (fexpr_is_builtin_call(expr, FEXPR_CurvePath) && fexpr_nargs(expr) == 2)
    {
        fexpr_t f, forexpr, x, a, b;
        slong forexpr_nargs;

        fexpr_view_arg(f, expr, 0);
        fexpr_view_arg(forexpr, expr, 1);
        forexpr_nargs = fexpr_nargs(forexpr);

        if (forexpr_nargs == 3)
        {
            fexpr_view_arg(x, forexpr, 0);
            fexpr_view_arg(a, forexpr, 1);
            fexpr_view_arg(b, forexpr, 2);

            calcium_write(out, "\\left(");
            fexpr_write_latex(out, f, flags);
            calcium_write(out, ",\\, ");
            fexpr_write_latex(out, x, flags);
            calcium_write(out, " : ");
            fexpr_write_latex(out, a, flags);
            calcium_write(out, " \\rightsquigarrow ");
            fexpr_write_latex(out, b, flags);
            calcium_write(out, "\\right)");
            return;
        }
    }

    if (fexpr_is_builtin_call(expr, FEXPR_AnalyticContinuation) && fexpr_nargs(expr) == 2)
    {
        fexpr_t f, forexpr, x, a, b;
        slong forexpr_nargs;

        fexpr_view_arg(f, expr, 0);
        fexpr_view_arg(forexpr, expr, 1);
        forexpr_nargs = fexpr_nargs(forexpr);

        if (forexpr_nargs == 2)
        {
            fexpr_view_arg(x, forexpr, 0);
            fexpr_view_arg(a, forexpr, 1);

            calcium_write(out, "\\mathop{\\text{Continuation}}\\limits_{\\displaystyle{");
            fexpr_write_latex(out, x, flags);
            calcium_write(out, ": ");
            fexpr_write_latex(out, a, flags);
            calcium_write(out, "}} \\, ");
            fexpr_write_latex(out, f, flags);
            return;
        }

        if (forexpr_nargs == 3)
        {
            fexpr_view_arg(x, forexpr, 0);
            fexpr_view_arg(a, forexpr, 1);
            fexpr_view_arg(b, forexpr, 2);

            calcium_write(out, "\\mathop{\\text{Continuation}}\\limits_{\\displaystyle{");
            fexpr_write_latex(out, x, flags);
            calcium_write(out, ": ");
            fexpr_write_latex(out, a, flags);
            calcium_write(out, " \\rightsquigarrow ");
            fexpr_write_latex(out, b, flags);
            calcium_write(out, "}} \\, ");
            fexpr_write_latex(out, f, flags);
            return;
        }
    }

    if (fexpr_is_builtin_call(expr, FEXPR_QSeriesCoefficient) && fexpr_nargs(expr) == 5)
    {
        fexpr_t f, tau, q, n, qdef;
        int parens;

        fexpr_view_arg(f, expr, 0);
        fexpr_view_arg(tau, expr, 1);
        fexpr_view_arg(q, expr, 2);
        fexpr_view_arg(n, expr, 3);
        fexpr_view_arg(qdef, expr, 4);
        parens = fexpr_is_builtin_call(f, FEXPR_Add) ||
                 fexpr_is_builtin_call(f, FEXPR_Sub);

        calcium_write(out, "[{");
        fexpr_write_latex(out, q, flags);
        calcium_write(out, "}^{");
        fexpr_write_latex(out, n, flags | FEXPR_LATEX_SMALL);
        calcium_write(out, "}] ");
        if (parens)
            calcium_write(out, "\\left(");
        fexpr_write_latex(out, f, flags);
        if (parens)
            calcium_write(out, "\\right)");
        calcium_write(out, " \\; \\left(");
        fexpr_write_latex(out, qdef, flags);
        calcium_write(out, "\\right)");
        return;
    }

    if (fexpr_is_builtin_call(expr, FEXPR_EqualQSeriesEllipsis) && fexpr_nargs(expr) == 5)
    {
        fexpr_t f, tau, q, ser, qdef;

        fexpr_view_arg(f, expr, 0);
        fexpr_view_arg(tau, expr, 1);
        fexpr_view_arg(q, expr, 2);
        fexpr_view_arg(ser, expr, 3);
        fexpr_view_arg(qdef, expr, 4);

        fexpr_write_latex(out, f, flags);
        calcium_write(out, " = ");
        fexpr_write_latex(out, ser, flags | FEXPR_LATEX_SMALL);
        calcium_write(out, " + \\ldots \\; \\text{ where } ");
        fexpr_write_latex(out, qdef, flags);
        return;
    }

    if (fexpr_is_builtin_call(expr, FEXPR_Coefficient) && fexpr_nargs(expr) == 3)
    {
        fexpr_t f, x, n;
        int parens;

        fexpr_view_arg(f, expr, 0);
        fexpr_view_arg(x, expr, 1);
        fexpr_view_arg(n, expr, 2);
        parens = fexpr_is_builtin_call(f, FEXPR_Add) ||
                 fexpr_is_builtin_call(f, FEXPR_Sub);

        calcium_write(out, "[{");
        fexpr_write_latex(out, x, flags);
        calcium_write(out, "}^{");
        fexpr_write_latex(out, n, flags | FEXPR_LATEX_SMALL);
        calcium_write(out, "}] ");
        if (parens)
            calcium_write(out, "\\left(");
        fexpr_write_latex(out, f, flags);
        if (parens)
            calcium_write(out, "\\right)");
        return;
    }

    if (fexpr_is_builtin_call(expr, FEXPR_DiscreteLog) && fexpr_nargs(expr) == 3)
    {
        fexpr_t n, b, q;

        fexpr_view_arg(n, expr, 0);
        fexpr_view_arg(b, expr, 1);
        fexpr_view_arg(q, expr, 2);

        calcium_write(out, "(\\epsilon : {");
        /* todo: may need parentheses */
        fexpr_write_latex(out, b, flags);
        calcium_write(out, "}^{\\epsilon} \\equiv ");
        fexpr_write_latex(out, n, flags);
        calcium_write(out, " \\text{ mod }");
        fexpr_write_latex(out, q, flags);
        calcium_write(out, ")");
        return;
    }

    if (fexpr_is_builtin_call(expr, FEXPR_AsymptoticTo) && fexpr_nargs(expr) == 4)
    {
        fexpr_t a, b, c, d;

        fexpr_view_arg(a, expr, 0);
        fexpr_view_arg(b, expr, 1);
        fexpr_view_arg(c, expr, 2);
        fexpr_view_arg(d, expr, 3);

        fexpr_write_latex(out, a, flags);
        calcium_write(out, " \\sim ");
        fexpr_write_latex(out, b, flags);
        calcium_write(out, ", \\; ");
        fexpr_write_latex(out, c, flags);
        calcium_write(out, " \\to ");
        fexpr_write_latex(out, d, flags);
        return;
    }

    if ((fexpr_is_builtin_call(expr, FEXPR_CoulombF) || fexpr_is_builtin_call(expr, FEXPR_CoulombG)) && fexpr_nargs(expr) == 3)
    {
        fexpr_t l, eta, z;

        fexpr_view_arg(l, expr, 0);
        fexpr_view_arg(eta, expr, 1);
        fexpr_view_arg(z, expr, 2);

        if (fexpr_is_builtin_call(expr, FEXPR_CoulombF))
            calcium_write(out, "F_{");
        else
            calcium_write(out, "G_{");
        fexpr_write_latex(out, l, flags | FEXPR_LATEX_SMALL);
            calcium_write(out, ",");
        fexpr_write_latex(out, eta, flags | FEXPR_LATEX_SMALL);
        calcium_write(out, "}");
        _fexpr_write_latex_call1(out, z, flags);
        return;
    }

    if (fexpr_is_builtin_call(expr, FEXPR_CoulombH) && fexpr_nargs(expr) == 4)
    {
        fexpr_t omega, l, eta, z;

        fexpr_view_arg(omega, expr, 0);
        fexpr_view_arg(l, expr, 1);
        fexpr_view_arg(eta, expr, 2);
        fexpr_view_arg(z, expr, 3);

        calcium_write(out, "H^{");
        if (fexpr_equal_si(omega, +1))
            calcium_write(out, "+");
        else if (fexpr_equal_si(omega, -1))
            calcium_write(out, "-");
        else
            fexpr_write_latex(out, omega, flags | FEXPR_LATEX_SMALL);

        calcium_write(out, "}_{");
        fexpr_write_latex(out, l, flags | FEXPR_LATEX_SMALL);
            calcium_write(out, ",");
        fexpr_write_latex(out, eta, flags | FEXPR_LATEX_SMALL);
        calcium_write(out, "}");
        _fexpr_write_latex_call1(out, z, flags);
        return;
    }

    if (fexpr_is_builtin_call(expr, FEXPR_Repeat) && fexpr_nargs(expr) >= 2)
    {
        slong i, n;
        fexpr_t x;

        n = fexpr_nargs(expr) - 1;

        calcium_write(out, "\\underbrace{");
        for (i = 0; i < n; i++)
        {
            fexpr_view_arg(x, expr, i);
            fexpr_write_latex(out, x, flags);
            if (i < n - 1)
                calcium_write(out, ", ");
        }

        calcium_write(out, ", \\ldots, ");
        for (i = 0; i < n; i++)
        {
            fexpr_view_arg(x, expr, i);
            fexpr_write_latex(out, x, flags);
            if (i < n - 1)
                calcium_write(out, ", ");
        }

        calcium_write(out, "}_{");

        if (n > 1)
        {
            calcium_write(out, "\\left(");
            for (i = 0; i < n; i++)
            {
                fexpr_view_arg(x, expr, i);
                fexpr_write_latex(out, x, flags);
                if (i < n - 1)
                    calcium_write(out, ", ");
            }
            calcium_write(out, "\\right) \\; ");
        }

        fexpr_view_arg(x, expr, n);
        fexpr_write_latex(out, x, flags);
        calcium_write(out, " \\text{ times}}");

        return;
    }

    if (fexpr_is_builtin_call(expr, FEXPR_Matrices))
    {
        fexpr_t R, n, m;

        if (fexpr_nargs(expr) == 2 || fexpr_nargs(expr) == 3)
        {
            fexpr_view_arg(R, expr, 0);
            fexpr_view_arg(n, expr, 1);

            calcium_write(out, "\\operatorname{M}_{");
            fexpr_write_latex(out, n, flags);
            if (fexpr_nargs(expr) == 3)
            {
                calcium_write(out, " \\times ");
                fexpr_view_arg(m, expr, 2);
                fexpr_write_latex(out, m, flags | FEXPR_LATEX_SMALL);
            }

            calcium_write(out, "}");
            _fexpr_write_latex_call1(out, R, flags);

            return;
        }
    }

    if (fexpr_is_builtin_call(expr, FEXPR_DirichletCharacter))
    {
        fexpr_t a, b, n;

        if (fexpr_nargs(expr) == 2 || fexpr_nargs(expr) == 3)
        {
            fexpr_view_arg(a, expr, 0);
            fexpr_view_arg(b, expr, 1);

            calcium_write(out, "\\chi_{");
            fexpr_write_latex(out, a, flags | FEXPR_LATEX_SMALL);
            calcium_write(out, " \\, . \\, ");
            fexpr_write_latex(out, b, flags | FEXPR_LATEX_SMALL);
            calcium_write(out, "}");

            if (fexpr_nargs(expr) == 3)
            {
                fexpr_view_arg(n, expr, 2);
                _fexpr_write_latex_call1(out, n, flags);
            }

            return;
        }
    }

    if (fexpr_is_builtin_call(expr, FEXPR_LambertW))
    {
        fexpr_t f, x, n, r;

        if (fexpr_nargs(expr) == 1)
        {
            fexpr_view_arg(n, expr, 0);
            calcium_write(out, "W");
            _fexpr_write_latex_call1(out, n, flags);
            return;
        }

        if (fexpr_nargs(expr) == 2)
        {
            fexpr_view_arg(x, expr, 0);
            fexpr_view_arg(n, expr, 1);
            calcium_write(out, "W_{");
            fexpr_write_latex(out, n, flags | FEXPR_LATEX_SMALL);
            calcium_write(out, "}");
            _fexpr_write_latex_call1(out, x, flags);
            return;
        }

        if (fexpr_nargs(expr) == 3)
        {
            fexpr_view_func(f, expr);
            fexpr_view_arg(x, expr, 0);
            fexpr_view_arg(n, expr, 1);
            fexpr_view_arg(r, expr, 2);
            _fexpr_write_latex_derivative(out, f, NULL, r, flags);
            calcium_write(out, "_{");
            fexpr_write_latex(out, n, flags | FEXPR_LATEX_SMALL);
            calcium_write(out, "}");
            _fexpr_write_latex_call1(out, x, flags);
            return;
        }
    }

    if (fexpr_is_builtin_call(expr, FEXPR_SloaneA) && fexpr_nargs(expr) == 2)
    {
        fexpr_t x, n;

        fexpr_view_arg(x, expr, 0);
        fexpr_view_arg(n, expr, 1);

        calcium_write(out, "\\text{");

        if (fexpr_is_integer(x))
        {
            fmpz_t t;
            fmpz_init(t);
            fexpr_get_fmpz(t, x);
            if (fmpz_cmp_ui(t, 10) < 0)
                calcium_write(out, "A00000");
            else if (fmpz_cmp_ui(t, 100) < 0)
                calcium_write(out, "A0000");
            else if (fmpz_cmp_ui(t, 1000) < 0)
                calcium_write(out, "A000");
            else if (fmpz_cmp_ui(t, 10000) < 0)
                calcium_write(out, "A00");
            else if (fmpz_cmp_ui(t, 100000) < 0)
                calcium_write(out, "A0");
            else
                calcium_write(out, "A");
            calcium_write_fmpz(out, t);
            fmpz_clear(t);
        }
        else if (fexpr_is_string(x))
        {
            char * s = fexpr_get_string(x);

            if (s[0] != 'A')
                calcium_write(out, "A");

            calcium_write_free(out, s);
        }
        else
        {
            calcium_write(out, "A00000");
            fexpr_write_latex(out, x, flags);
        }

        calcium_write(out, "}\\!\\left(");
        fexpr_write_latex(out, n, flags);
        calcium_write(out, "\\right)");
        return;
    }

    if (fexpr_is_builtin_call(expr, FEXPR_EqualAndElement))
    {
        fexpr_t a, b, c;

        if (fexpr_nargs(expr) == 3)
        {
            fexpr_view_arg(a, expr, 0);
            fexpr_view_arg(b, expr, 1);
            fexpr_view_arg(c, expr, 2);

            fexpr_write_latex(out, a, flags);
            calcium_write(out, " = ");
            fexpr_write_latex(out, b, flags);
            calcium_write(out, " \\in ");
            fexpr_write_latex(out, c, flags);

            return;
        }
    }

    if (fexpr_is_builtin_call(expr, FEXPR_EqualNearestDecimal))
    {
        fexpr_t a, b, c;

        if (fexpr_nargs(expr) == 3)
        {
            fexpr_view_arg(a, expr, 0);
            fexpr_view_arg(b, expr, 1);
            fexpr_view_arg(c, expr, 2);

            fexpr_write_latex(out, a, flags);
            calcium_write(out, " = ");
            fexpr_write_latex(out, b, flags);
            calcium_write(out, " \\;\\, {\\scriptstyle (\\text{nearest } ");
            fexpr_write_latex(out, c, flags);
            calcium_write(out, " \\text{ digits})}");

            return;
        }
    }

    if (fexpr_is_builtin_call(expr, FEXPR_IsHolomorphicOn) ||
        fexpr_is_builtin_call(expr, FEXPR_IsMeromorphicOn))
    {
        if (fexpr_nargs(expr) == 2)
        {
            fexpr_t func, forexpr;
            fexpr_view_arg(func, expr, 0);
            fexpr_view_arg(forexpr, expr, 1);

            if (fexpr_nargs(forexpr) == 2)
            {
                fexpr_t var, domain;
                fexpr_view_arg(var, forexpr, 0);
                fexpr_view_arg(domain, forexpr, 1);

                fexpr_write_latex(out, func, flags);
                if (fexpr_is_builtin_call(expr, FEXPR_IsHolomorphicOn))
                    calcium_write(out, " \\text{ is holomorphic on } ");
                else
                    calcium_write(out, " \\text{ is meromorphic on } ");
                fexpr_write_latex(out, var, flags);
                calcium_write(out, " \\in ");
                fexpr_write_latex(out, domain, flags);
                return;
            }
        }
    }

    fexpr_write_latex_call(out, expr, flags);
}

void
fexpr_write_latex_show_form(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    fexpr_t f, arg;

    if (fexpr_view_call1(f, arg, expr) && fexpr_is_builtin_symbol(f, FEXPR_ShowExpandedNormalForm))
    {
        fexpr_t v;
        fexpr_init(v);
        fexpr_expanded_normal_form(v, arg, 0);
        fexpr_write_latex(out, v, flags);
        fexpr_clear(v);
        return;
    }

    fexpr_write_latex_call(out, expr, flags);
}

void
fexpr_write_latex_alg_structure(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    fexpr_t func, arg;
    const char *a;
    const char *b;
    slong i, nargs;

    nargs = fexpr_nargs(expr);

    if (nargs <= 1)
    {
        fexpr_write_latex_call(out, expr, flags);
        return;
    }

    fexpr_view_func(func, expr);
    fexpr_view_arg(arg, expr, 0);

    i = FEXPR_BUILTIN_ID(func->data[0]);

    switch (i)
    {
        case FEXPR_Polynomials:
            a = "[";
            b = "]";
            break;
        case FEXPR_PolynomialFractions:
            a = "(";
            b = ")";
            break;
        case FEXPR_FormalPowerSeries:
            a = "[[";
            b = "]]";
            break;
        case FEXPR_FormalLaurentSeries:
            a = "(\\!(";
            b = ")\\!)";
            break;
        case FEXPR_FormalPuiseuxSeries:
            a = "\\!\\left\\langle\\!\\left\\langle ";
            b = " \\right\\rangle\\!\\right\\rangle";
            break;
        default:
            fexpr_write_latex_call(out, expr, flags);
            return;
    }

    fexpr_write_latex(out, arg, flags);
    calcium_write(out, a);

    if (nargs >= 1)
    {
        slong i;

        fexpr_view_next(arg);

        if (fexpr_is_builtin_call(arg, FEXPR_Tuple))
        {
            nargs = fexpr_nargs(arg);
            fexpr_view_arg(arg, arg, 0);
        }
        else
        {
            nargs--;
        }

        for (i = 0; i < nargs; i++)
        {
            fexpr_write_latex(out, arg, flags);

            if (i < nargs - 1)
            {
                calcium_write(out, ", ");
                fexpr_view_next(arg);
            }
        }
    }

    calcium_write(out, b);
}

void
fexpr_write_latex_where(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    fexpr_t f, arg, x, val;
    slong i, nargs;

    nargs = fexpr_nargs(expr);

    if (nargs > 0)
    {
        fexpr_view_arg(f, expr, 0);
        fexpr_write_latex(out, f, flags);
    }

    if (nargs > 1)
    {
        calcium_write(out, "\\; \\text{ where } ");

        fexpr_view_arg(arg, expr, 1);

        for (i = 1; i < nargs; i++)
        {
            if (fexpr_nargs(arg) == 2)
            {
                fexpr_view_arg(x, arg, 0);
                fexpr_view_arg(val, arg, 1);

                fexpr_write_latex(out, x, flags);
                calcium_write(out, " = ");
                fexpr_write_latex(out, val, flags);

                if (i < nargs - 1)
                {
                    calcium_write(out, ",\\;");
                    fexpr_view_next(arg);
                }
            }
        }
    }
}

static int
_fexpr_all_arguments_small(const fexpr_t expr)
{
    fexpr_t arg;
    slong i, nargs;

    nargs = fexpr_nargs(expr);
    fexpr_view_arg(arg, expr, 0);

    for (i = 0; i < nargs; i++)
    {
        if (!fexpr_is_atom(arg))
            return 0;

        fexpr_view_next(arg);
    }

    return 1;
}

void
fexpr_write_latex_symbol(int * subscript, calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    if (fexpr_is_any_builtin_symbol(expr))
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

        *subscript = 0;
    }
    else if (fexpr_is_symbol(expr))
    {
        const char *s;
        char tmp[FEXPR_SMALL_SYMBOL_LEN + 1];
        slong i, len;

        s = fexpr_get_symbol_str_pointer(tmp, expr);
        len = strlen(s);

        if (len > 1 && s[len - 1] == '_')
        {
            char * tmp2;
            tmp2 = flint_malloc(len);
            memcpy(tmp2, tmp, len - 1);
            tmp2[len - 1] = '\0';
            calcium_write(out, tmp2);
            *subscript = 1;
            flint_free(tmp2);
        }
        else if (len == 1)
        {
            calcium_write(out, s);
            *subscript = 0;
        }
        else
        {
            /* Look for internal underscore */
            slong pos = 0;
            for (i = 1; i < len - 1; i++)
            {
                if (s[i] == '_')
                {
                    pos = i;
                    break;
                }
            }

            if (pos == 0)
            {
                calcium_write(out, "\\operatorname{");
                calcium_write(out, s);
                calcium_write(out, "}");
                *subscript = 0;
            }
            else
            {
                char * tmp2;
                tmp2 = flint_malloc(len);
                memcpy(tmp2, tmp, pos);
                tmp2[pos] = '\0';

                if (pos == 1)
                {
                    calcium_write(out, tmp2);
                }
                else
                {
                    calcium_write(out, "\\operatorname{");
                    calcium_write(out, tmp2);
                    calcium_write(out, "}");
                }

                calcium_write(out, "_{");
                calcium_write(out, s + pos + 1);
                calcium_write(out, "}");

                flint_free(tmp2);
                *subscript = 0;
            }
        }
    }
    else
    {
        if (fexpr_is_builtin_call(expr, FEXPR_Add) ||
            fexpr_is_builtin_call(expr, FEXPR_Sub) ||
            fexpr_is_builtin_call(expr, FEXPR_Mul) ||
            fexpr_is_builtin_call(expr, FEXPR_Div) ||
            fexpr_is_builtin_call(expr, FEXPR_Neg) ||
            fexpr_is_builtin_call(expr, FEXPR_Pos))
        {
            calcium_write(out, "\\left(");
            fexpr_write_latex(out, expr, flags);
            calcium_write(out, "\\right)");
            *subscript = 0;
        }
        else
        {
            fexpr_write_latex(out, expr, flags);
            *subscript = 0;
        }
    }
}

void
fexpr_write_latex_call(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    fexpr_t view;
    slong i, nargs;
    int small, subscript;

    nargs = fexpr_nargs(expr);

    fexpr_view_func(view, expr);
    fexpr_write_latex_symbol(&subscript, out, view, flags);

    if (subscript)
    {
        calcium_write(out, "_{");
        for (i = 0; i < nargs; i++)
        {
            fexpr_view_next(view);
            fexpr_write_latex(out, view, flags | FEXPR_LATEX_SMALL);
            if (i < nargs - 1)
                calcium_write(out, ", ");
        }
        calcium_write(out, "}");
    }
    else
    {
        fexpr_view_next(view);

        /* Function of matrix already has parentheses */
        if (nargs == 1)
        {
            if (fexpr_is_builtin_call(view, FEXPR_Matrix) || fexpr_is_builtin_call(view, FEXPR_Matrix2x2))
            {
                calcium_write(out, " ");
                fexpr_write_latex(out, view, flags);
                return;
            }
        }

        small = _fexpr_all_arguments_small(expr);

        if (small)
            calcium_write(out, "(");
        else
            calcium_write(out, "\\!\\left(");

        for (i = 0; i < nargs; i++)
        {
            fexpr_write_latex(out, view, flags);
            if (i < nargs - 1)
            {
                calcium_write(out, ", ");
                fexpr_view_next(view);
            }
        }

        if (small)
            calcium_write(out, ")");
        else
            calcium_write(out, "\\right)");
    }
}

void fexpr_write_latex_subscript_call(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    fexpr_t view;
    slong i, nargs;
    int subscript;

    nargs = fexpr_nargs(expr);

    fexpr_view_func(view, expr);
    fexpr_write_latex_symbol(&subscript, out, view, flags);

    if (nargs >= 1)
    {
        calcium_write(out, "_{");
        fexpr_view_next(view);
        fexpr_write_latex(out, view, flags | FEXPR_LATEX_SMALL);
        calcium_write(out, "}");
    }

    if (nargs >= 2)
    {
        calcium_write(out, "\\!\\left(");

        for (i = 1; i < nargs; i++)
        {
            fexpr_view_next(view);
            fexpr_write_latex(out, view, flags);
            if (i < nargs - 1)
                calcium_write(out, ", ");
        }

        calcium_write(out, "\\right)");
    }
}

void
fexpr_write_latex_subscript(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    fexpr_t view;
    slong i, nargs;
    int subscript;

    nargs = fexpr_nargs(expr);

    fexpr_view_func(view, expr);
    fexpr_write_latex_symbol(&subscript, out, view, flags);

    calcium_write(out, "_{");
    for (i = 0; i < nargs; i++)
    {
        fexpr_view_next(view);
        fexpr_write_latex(out, view, flags | FEXPR_LATEX_SMALL);
        if (i < nargs - 1)
            calcium_write(out, ", ");
    }
    calcium_write(out, "}");
}

/* f(x) */
/* f(x,n) -> nth derivative */
void
fexpr_write_latex_call1_optional_derivative(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    fexpr_t func, x, order;
    slong nargs;

    nargs = fexpr_nargs(expr);

    if (nargs == 2)
    {
        fexpr_view_func(func, expr);
        fexpr_view_arg(x, expr, 0);
        fexpr_view_arg(order, expr, 1);
        _fexpr_write_latex_derivative(out, func, NULL, order, flags);
        calcium_write(out, "\\!\\left(");
        fexpr_write_latex(out, x, flags);
        calcium_write(out, "\\right)");
    }
    else
    {
        fexpr_write_latex_call(out, expr, flags);
    }
}

/* f(x,y) */
/* f(x,y,n) -> nth derivative */
void
fexpr_write_latex_call2_optional_derivative(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    fexpr_t func, x, y, order;
    slong nargs;

    nargs = fexpr_nargs(expr);

    if (nargs == 3)
    {
        fexpr_view_func(func, expr);
        fexpr_view_arg(x, expr, 0);
        fexpr_view_arg(y, expr, 1);
        fexpr_view_arg(order, expr, 2);
        _fexpr_write_latex_derivative(out, func, NULL, order, flags);
        calcium_write(out, "\\!\\left(");
        fexpr_write_latex(out, x, flags);
        calcium_write(out, ", ");
        fexpr_write_latex(out, y, flags);
        calcium_write(out, "\\right)");
    }
    else
    {
        fexpr_write_latex_call(out, expr, flags);
    }
}

/* f_x(y) */
/* f_x(y,n) -> nth derivative */
void
fexpr_write_latex_sub1_call1_optional_derivative(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    fexpr_t func, x, y, order;
    slong nargs;

    nargs = fexpr_nargs(expr);

    if (nargs == 3)
    {
        fexpr_view_func(func, expr);
        fexpr_view_arg(x, expr, 0);
        fexpr_view_arg(y, expr, 1);
        fexpr_view_arg(order, expr, 2);
        _fexpr_write_latex_derivative(out, func, NULL, order, flags);
        calcium_write(out, "_{");
        fexpr_write_latex(out, x, flags | FEXPR_LATEX_SMALL);
        calcium_write(out, "}");
        calcium_write(out, "\\!\\left(");
        fexpr_write_latex(out, y, flags);
        calcium_write(out, "\\right)");
    }
    else
    {
        fexpr_write_latex_subscript_call(out, expr, flags);
    }
}

/* f_x(y,z) */
/* f_x(y,z,n) -> nth derivative */
void
fexpr_write_latex_sub1_call2_optional_derivative(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    fexpr_t func, x, y, z, order;
    slong nargs;

    nargs = fexpr_nargs(expr);

    if (nargs == 4)
    {
        fexpr_view_func(func, expr);
        fexpr_view_arg(x, expr, 0);
        fexpr_view_arg(y, expr, 1);
        fexpr_view_arg(z, expr, 2);
        fexpr_view_arg(order, expr, 3);
        _fexpr_write_latex_derivative(out, func, NULL, order, flags);
        calcium_write(out, "_{");
        fexpr_write_latex(out, x, flags | FEXPR_LATEX_SMALL);
        calcium_write(out, "}");
        calcium_write(out, "\\!\\left(");
        fexpr_write_latex(out, y, flags);
        calcium_write(out, ", ");
        fexpr_write_latex(out, z, flags);
        calcium_write(out, "\\right)");
    }
    else
    {
        fexpr_write_latex_subscript_call(out, expr, flags);
    }
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

        if (fexpr_is_builtin_call(arg, FEXPR_Step) && fexpr_nargs(arg) == 2)
        {
            fexpr_t x, forexpr, n, a, b, na, na1, nb;
            fmpz_t t;

            fexpr_view_arg(x, arg, 0);
            fexpr_view_arg(forexpr, arg, 1);

            if (fexpr_nargs(forexpr) == 3)
            {
                fexpr_view_arg(n, forexpr, 0);
                fexpr_view_arg(a, forexpr, 1);
                fexpr_view_arg(b, forexpr, 2);

                fexpr_init(na);
                fexpr_init(na1);
                fexpr_init(nb);
                fmpz_init(t);

                fexpr_replace(na, x, n, a);
                if (fexpr_is_integer(a))
                {
                    fexpr_get_fmpz(t, a);
                    fmpz_add_ui(t, t, 1);
                    fexpr_set_fmpz(na1, t);
                    fexpr_swap(nb, na1);
                    fexpr_replace(na1, x, n, nb);
                }
                else
                {
                    fexpr_set_ui(nb, 1);
                    fexpr_add(na1, a, nb);
                    fexpr_swap(nb, na1);
                    fexpr_replace(na1, x, n, nb);
                }
                fexpr_replace(nb, x, n, b);

                fexpr_write_latex(out, na, flags);
                calcium_write(out, " ");
                fexpr_write_latex(out, func, flags);
                calcium_write(out, " ");
                fexpr_write_latex(out, na1, flags);
                calcium_write(out, " ");
                fexpr_write_latex(out, func, flags);
                calcium_write(out, " \\ldots ");
                fexpr_write_latex(out, func, flags);
                calcium_write(out, " ");
                fexpr_write_latex(out, nb, flags);

                fexpr_clear(na);
                fexpr_clear(na1);
                fexpr_clear(nb);
                fmpz_clear(t);
            }
            else
            {
                fexpr_write_latex(out, arg, flags);
            }
        }
        else
        {
            fexpr_write_latex(out, arg, flags);
        }

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
            char * s;

            /* todo: escape strings */
            s = fexpr_get_string(expr);

            calcium_write(out, "\\text{``");
            calcium_write(out, s);
            calcium_write(out, "''}");
            flint_free(s);
        }
        else
        {
            int subscript;
            fexpr_write_latex_symbol(&subscript, out, expr, flags);
        }
    }
    else
    {
        fexpr_t func;
        slong i;

        fexpr_view_func(func, expr);

        if (fexpr_is_any_builtin_symbol(func))
        {
            i = FEXPR_BUILTIN_ID(func->data[0]);

            if (fexpr_builtin_table[i].latex_writer != NULL)
            {
                (fexpr_builtin_table[i].latex_writer)(out, expr, flags);
                return;
            }
        }

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
