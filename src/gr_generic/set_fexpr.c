/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fexpr.h"
#include "fexpr_builtin.h"
#include "qqbar.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_special.h"

#define BINARY_OP(gr_func) \
    if (nargs == 2) \
    { \
        GR_TMP_INIT(t, ctx); \
        fexpr_view_arg(arg, expr, 0); \
        status = gr_set_fexpr(res, inputs, outputs, arg, ctx); \
        if (status == GR_SUCCESS) \
        { \
            fexpr_view_next(arg); \
            status = gr_set_fexpr(t, inputs, outputs, arg, ctx); \
            if (status == GR_SUCCESS) \
                status = gr_func(res, res, t, ctx); \
        } \
        GR_TMP_CLEAR(t, ctx); \
        return status; \
    } \
    return GR_DOMAIN; \

#define BINARY_OP_WITH_FMPZ(gr_func, gr_func_fmpz) \
    if (nargs == 2) \
    { \
        GR_TMP_INIT(t, ctx); \
        fexpr_view_arg(arg, expr, 0); \
        status = gr_set_fexpr(res, inputs, outputs, arg, ctx); \
        if (status == GR_SUCCESS) \
        { \
            fexpr_view_next(arg); \
            if (fexpr_is_integer(arg)) \
            { \
                fmpz_t n; \
                fmpz_init(n); \
                fexpr_get_fmpz(n, arg); \
                status = gr_func_fmpz(res, res, n, ctx); \
                fmpz_clear(n); \
            } \
            else \
            { \
                status = gr_set_fexpr(t, inputs, outputs, arg, ctx); \
                if (status == GR_SUCCESS) \
                    status = gr_func(res, res, t, ctx); \
            } \
        } \
        GR_TMP_CLEAR(t, ctx); \
        return status; \
    } \
    return GR_DOMAIN; \


#define UNARY_OP(gr_func) \
    if (nargs == 1) \
    { \
        fexpr_view_arg(arg, expr, 0); \
        status = gr_set_fexpr(res, inputs, outputs, arg, ctx); \
        if (status == GR_SUCCESS) \
            status = gr_func(res, res, ctx); \
        return status; \
    } \
    return GR_DOMAIN; \

#define NARY_OP(gr_func, gr_empty_func) \
    if (nargs == 0) \
        return gr_empty_func(res, ctx); \
    fexpr_view_arg(arg, expr, 0); \
    status = gr_set_fexpr(res, inputs, outputs, arg, ctx); \
    if (status == GR_SUCCESS && nargs > 1) \
    { \
        /* todo: divide and conquer */ \
        GR_TMP_INIT(t, ctx); \
        for (i = 1; i < nargs && status == GR_SUCCESS; i++) \
        { \
            fexpr_view_next(arg); \
            status = gr_set_fexpr(t, inputs, outputs, arg, ctx); \
            if (status == GR_SUCCESS) \
                status = gr_func(res, res, t, ctx); \
        } \
        GR_TMP_CLEAR(t, ctx); \
    } \
    return status; \

int
gr_generic_set_fexpr(gr_ptr res, fexpr_vec_t inputs, gr_vec_t outputs, const fexpr_t expr, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (fexpr_is_integer(expr))
    {
        fmpz_t t;
        fmpz_init(t);
        fexpr_get_fmpz(t, expr);
        status = gr_set_fmpz(res, t, ctx);
        fmpz_clear(t);
        return status;
    }

    if (fexpr_is_any_builtin_symbol(expr))
    {
        slong op = FEXPR_BUILTIN_ID(expr->data[0]);

        switch (op)
        {
            case FEXPR_Pi:
                return gr_pi(res, ctx);
            case FEXPR_NumberI:
                return gr_i(res, ctx);
            case FEXPR_NumberE:
                status |= gr_one(res, ctx);
                status |= gr_exp(res, res, ctx);
                return status;
            case FEXPR_Euler:
                return gr_euler(res, ctx);
            case FEXPR_GoldenRatio:
                /* todo: have a builtin */
                status |= gr_set_ui(res, 5, ctx);
                status |= gr_sqrt(res, res, ctx);
                status |= gr_add_ui(res, res, 1, ctx);
                status |= gr_div_ui(res, res, 2, ctx);
                return status;
            case FEXPR_CatalanConstant:
                return gr_catalan(res, ctx);
            case FEXPR_KhinchinConstant:
                return gr_khinchin(res, ctx);
            case FEXPR_GlaisherConstant:
                return gr_glaisher(res, ctx);
            case FEXPR_Infinity:
                return gr_pos_inf(res, ctx);
            case FEXPR_UnsignedInfinity:
                return gr_uinf(res, ctx);
            case FEXPR_Undefined:
                return gr_undefined(res, ctx);
            case FEXPR_Unknown:
                return gr_unknown(res, ctx);
        }

        return GR_UNABLE;
    }

    if (fexpr_is_symbol(expr))
    {
        slong i, num_defs;

        num_defs = inputs->length;

        /* Treat local definitions as a stack, more recent ones
           overriding older ones */
        for (i = num_defs - 1; i >= 0; i--)
        {
            if (fexpr_equal(expr, fexpr_vec_entry(inputs, i)))
            {
                return gr_set(res, gr_vec_entry_srcptr(outputs, i, ctx), ctx);
            }
        }

        return GR_UNABLE;
    }

    if (fexpr_is_any_builtin_call(expr))
    {
        fexpr_t func, arg;
        slong op, i, nargs;
        int status;
        gr_ptr t;

        nargs = fexpr_nargs(expr);
        fexpr_view_func(func, expr);
        op = FEXPR_BUILTIN_ID(func->data[0]);

        /* Parse local definitions */
        if (op == FEXPR_Where)
        {
            slong num_previous_defs;

            num_previous_defs = inputs->length;
            status = GR_SUCCESS;

            /* Parse in reverse order; this assumes definitions are in the form
               x = f(y,z), y = g(z), z = h  which is what gr_get_fexpr currently
               generates. Should this work in both directions (or any order)?
            */

            for (i = nargs - 1; i >= 1; i--)
            {
                fexpr_t defn, symbol, value;

                fexpr_view_arg(defn, expr, i);

                if (!fexpr_is_builtin_call(defn, FEXPR_Def) || fexpr_nargs(defn) != 2)
                {
                    status = GR_DOMAIN;
                    break;
                }

                fexpr_view_arg(symbol, defn, 0);
                fexpr_view_arg(value, defn, 1);

                status = gr_set_fexpr(res, inputs, outputs, value, ctx);
                if (status != GR_SUCCESS)
                    break;

                fexpr_vec_append(inputs, symbol);
                status = gr_vec_append(outputs, res, ctx);
                if (status != GR_SUCCESS)
                    break;
            }

            if (status == GR_SUCCESS)
            {
                fexpr_view_arg(arg, expr, 0);
                status = gr_set_fexpr(res, inputs, outputs, arg, ctx);
            }

            /* We are done with the local definitions, so erase anything new. */
            fexpr_vec_set_length(inputs, num_previous_defs);
            gr_vec_set_length(outputs, num_previous_defs, ctx);

            return status;
        }

        switch (op)
        {
            /* todo: generalize to non-algebraics; handle large decimals efficiently */
            case FEXPR_Decimal:
            case FEXPR_PolynomialRootIndexed:
            case FEXPR_PolynomialRootNearest:
            case FEXPR_AlgebraicNumberSerialized:
                {
                    qqbar_t a;
                    qqbar_init(a);
                    status = qqbar_set_fexpr(a, expr);
                    if (status == GR_SUCCESS)
                    {
                        gr_ctx_t QQbar;
                        gr_ctx_init_complex_qqbar(QQbar);  /* no need to free */
                        status = gr_set_other(res, a, QQbar, ctx);
                    }
                    qqbar_clear(a);
                    return status;
                }

            case FEXPR_Pos: UNARY_OP(gr_set)
            case FEXPR_Neg: UNARY_OP(gr_neg)
            case FEXPR_Sub: BINARY_OP(gr_sub)
            case FEXPR_Div: BINARY_OP(gr_div)
            case FEXPR_Pow: BINARY_OP_WITH_FMPZ(gr_pow, gr_pow_fmpz)
            case FEXPR_Sqrt: UNARY_OP(gr_sqrt)
            case FEXPR_Exp: UNARY_OP(gr_exp)
            case FEXPR_Log: UNARY_OP(gr_log)
            case FEXPR_Sin: UNARY_OP(gr_sin)
            case FEXPR_Cos: UNARY_OP(gr_cos)
            case FEXPR_Tan: UNARY_OP(gr_tan)
            case FEXPR_Cot: UNARY_OP(gr_cot)
            case FEXPR_Atan: UNARY_OP(gr_atan)
            case FEXPR_Acos: UNARY_OP(gr_acos)
            case FEXPR_Asin: UNARY_OP(gr_asin)
            case FEXPR_Sign: UNARY_OP(gr_sgn)
            case FEXPR_Csgn: UNARY_OP(gr_csgn)
            case FEXPR_Arg: UNARY_OP(gr_arg)
            case FEXPR_Abs: UNARY_OP(gr_abs)
            case FEXPR_Re: UNARY_OP(gr_re)
            case FEXPR_Im: UNARY_OP(gr_im)
            case FEXPR_Conjugate: UNARY_OP(gr_conj)
            case FEXPR_Floor: UNARY_OP(gr_floor)
            case FEXPR_Ceil: UNARY_OP(gr_ceil)
            case FEXPR_Gamma: UNARY_OP(gr_gamma)
            case FEXPR_Erf: UNARY_OP(gr_erf)
            case FEXPR_Erfc: UNARY_OP(gr_erfc)
            case FEXPR_Erfi: UNARY_OP(gr_erfi)

            case FEXPR_Add: NARY_OP(gr_add, gr_zero)
            case FEXPR_Mul: NARY_OP(gr_mul, gr_one)
            case FEXPR_GCD: NARY_OP(gr_gcd, gr_zero)
        }
    }

    return 0;
}

/*
int
gr_set_fexpr(gr_ptr res, const fexpr_t expr, gr_ctx_t ctx)
{
    int status;

    fexpr_vec_t inputs;
    gr_vec_t outputs;

    fexpr_vec_init(inputs, 0);
    gr_vec_init(outputs, 0, ctx);

    status = _gr_set_fexpr(res, inputs, outputs, expr, ctx);

    fexpr_vec_clear(inputs);
    gr_vec_clear(outputs, ctx);

    return status;
}
*/
