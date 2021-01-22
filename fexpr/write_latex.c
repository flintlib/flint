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

static int
fexpr_is_any_builtin_call(const fexpr_t expr)
{
    fexpr_t func;

    if (fexpr_is_atom(expr))
        return 0;

    fexpr_view_func(func, expr);

    return fexpr_is_any_builtin_symbol(func);
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
        /* todo: Parentheses, Braces, Brackets, ... */
        if (fexpr_is_builtin_call(base, FEXPR_Abs))
            return 1;
        if (fexpr_is_builtin_call(base, FEXPR_Binomial))
            return 1;
        if (fexpr_is_builtin_call(base, FEXPR_Matrix2x2))
            return 1;

        return 0;
    }
}

/*
Special cases todo:
- Subscripted functions
- ...
*/
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

    fexpr_write_latex_call(out, expr, flags);
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
        flint_printf("fexpr_get_symbol_str_pointer: a symbol is required\n");
        flint_abort();
    }
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
        slong len;

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
            free(tmp2);
        }
        else
        {
            calcium_write(out, s);
            *subscript = 0;
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
        small = _fexpr_all_arguments_small(expr);

        if (small)
            calcium_write(out, "(");
        else
            calcium_write(out, "\\!\\left(");

        for (i = 0; i < nargs; i++)
        {
            fexpr_view_next(view);
            fexpr_write_latex(out, view, flags);
            if (i < nargs - 1)
                calcium_write(out, ", ");
        }

        if (small)
            calcium_write(out, ")");
        else
            calcium_write(out, "\\right)");
    }
}

void fexpr_write_latex_subscript_call(calcium_stream_t out, const fexpr_t expr, ulong flags)
{
    
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
