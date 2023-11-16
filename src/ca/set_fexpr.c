/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fexpr.h"
#include "fexpr_builtin.h"
#include "ca.h"
#include "ca_ext.h"
#include "ca_vec.h"

#define BINARY_OP(ca_func) \
    if (nargs == 2) \
    { \
        ca_init(t, ctx); \
        fexpr_view_arg(arg, expr, 0); \
        success = _ca_set_fexpr(res, inputs, outputs, arg, ctx); \
        if (success) \
        { \
            fexpr_view_next(arg); \
            success = _ca_set_fexpr(t, inputs, outputs, arg, ctx); \
            if (success) \
                ca_func(res, res, t, ctx); \
        } \
        ca_clear(t, ctx); \
        return success; \
    } \
    return 0; \

#define UNARY_OP(ca_func) \
    if (nargs == 1) \
    { \
        fexpr_view_arg(arg, expr, 0); \
        success = _ca_set_fexpr(res, inputs, outputs, arg, ctx); \
        if (success) \
            ca_func(res, res, ctx); \
        return success; \
    } \
    return 0; \

int
_ca_set_fexpr(ca_t res, fexpr_vec_t inputs, ca_vec_t outputs, const fexpr_t expr, ca_ctx_t ctx)
{
    if (fexpr_is_integer(expr))
    {
        _ca_make_fmpq(res, ctx);
        fexpr_get_fmpz(CA_FMPQ_NUMREF(res), expr);
        fmpz_one(CA_FMPQ_DENREF(res));
        return 1;
    }

    if (fexpr_is_any_builtin_symbol(expr))
    {
        slong op = FEXPR_BUILTIN_ID(expr->data[0]);

        switch (op)
        {
            case FEXPR_Pi:
                ca_pi(res, ctx);
                return 1;
            case FEXPR_NumberI:
                ca_i(res, ctx);
                return 1;
            case FEXPR_NumberE:
                ca_one(res, ctx);
                ca_exp(res, res, ctx);
                return 1;
            case FEXPR_Euler:
                ca_euler(res, ctx);
                return 1;
            case FEXPR_GoldenRatio:
                ca_sqrt_ui(res, 5, ctx);
                ca_add_ui(res, res, 1, ctx);
                ca_div_ui(res, res, 2, ctx);
                return 1;
            case FEXPR_Infinity:
                ca_pos_inf(res, ctx);
                return 1;
            case FEXPR_UnsignedInfinity:
                ca_uinf(res, ctx);
                return 1;
            case FEXPR_Undefined:
                ca_undefined(res, ctx);
                return 1;
            case FEXPR_Unknown:
                ca_unknown(res, ctx);
                return 1;
        }

        return 0;
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
                ca_set(res, ca_vec_entry(outputs, i), ctx);
                return 1;
            }
        }

        return 0;
    }

    if (fexpr_is_any_builtin_call(expr))
    {
        fexpr_t func, arg;
        slong op, i, nargs;
        int success;
        ca_t t;

        nargs = fexpr_nargs(expr);
        fexpr_view_func(func, expr);
        op = FEXPR_BUILTIN_ID(func->data[0]);

        /* Parse local definitions */
        if (op == FEXPR_Where)
        {
            slong num_previous_defs;

            num_previous_defs = inputs->length;
            success = 1;

            /* Parse in reverse order; this assumes definitions are in the form
               x = f(y,z), y = g(z), z = h  which is what ca_get_fexpr currently
               generates. Should this work in both directions (or any order)?
            */

            for (i = nargs - 1; i >= 1; i--)
            {
                fexpr_t defn, symbol, value;

                fexpr_view_arg(defn, expr, i);

                if (!fexpr_is_builtin_call(defn, FEXPR_Def) || fexpr_nargs(defn) != 2)
                {
                    success = 0;
                    break;
                }

                fexpr_view_arg(symbol, defn, 0);
                fexpr_view_arg(value, defn, 1);

                success = _ca_set_fexpr(res, inputs, outputs, value, ctx);
                if (!success)
                    break;

                fexpr_vec_append(inputs, symbol);
                ca_vec_append(outputs, res, ctx);
            }

            if (success)
            {
                fexpr_view_arg(arg, expr, 0);
                success = _ca_set_fexpr(res, inputs, outputs, arg, ctx);
            }

            /* We are done with the local definitions, so erase anything new. */
            fexpr_vec_set_length(inputs, num_previous_defs);
            ca_vec_set_length(outputs, num_previous_defs, ctx);

            return success;
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
                    success = qqbar_set_fexpr(a, expr);
                    if (success)
                        ca_set_qqbar(res, a, ctx);
                    qqbar_clear(a);
                    return success;
                }

            case FEXPR_Pos: UNARY_OP(ca_set)
            case FEXPR_Neg: UNARY_OP(ca_neg)
            case FEXPR_Sub: BINARY_OP(ca_sub)
            case FEXPR_Div: BINARY_OP(ca_div)
            case FEXPR_Pow: BINARY_OP(ca_pow)
            case FEXPR_Sqrt: UNARY_OP(ca_sqrt)
            case FEXPR_Exp: UNARY_OP(ca_exp)
            case FEXPR_Log: UNARY_OP(ca_log)
            case FEXPR_Sin: UNARY_OP(ca_sin)
            case FEXPR_Cos: UNARY_OP(ca_cos)
            case FEXPR_Tan: UNARY_OP(ca_tan)
            case FEXPR_Cot: UNARY_OP(ca_cot)
            case FEXPR_Atan: UNARY_OP(ca_atan)
            case FEXPR_Acos: UNARY_OP(ca_acos)
            case FEXPR_Asin: UNARY_OP(ca_asin)
            case FEXPR_Sign: UNARY_OP(ca_sgn)
            case FEXPR_Csgn: UNARY_OP(ca_csgn)
            case FEXPR_Arg: UNARY_OP(ca_arg)
            case FEXPR_Abs: UNARY_OP(ca_abs)
            case FEXPR_Re: UNARY_OP(ca_re)
            case FEXPR_Im: UNARY_OP(ca_im)
            case FEXPR_Conjugate: UNARY_OP(ca_conj)
            case FEXPR_Floor: UNARY_OP(ca_floor)
            case FEXPR_Ceil: UNARY_OP(ca_ceil)
            case FEXPR_Gamma: UNARY_OP(ca_gamma)
            case FEXPR_Erf: UNARY_OP(ca_erf)
            case FEXPR_Erfc: UNARY_OP(ca_erfc)
            case FEXPR_Erfi: UNARY_OP(ca_erfi)

            case FEXPR_Add:
                if (nargs == 0)
                {
                    ca_zero(res, ctx);
                    return 1;
                }

                fexpr_view_arg(arg, expr, 0);
                success = _ca_set_fexpr(res, inputs, outputs, arg, ctx);

                if (success && nargs > 1)
                {
                    /* todo: divide and conquer for large nargs? */
                    ca_init(t, ctx);
                    for (i = 1; i < nargs && success; i++)
                    {
                        fexpr_view_next(arg);
                        success = _ca_set_fexpr(t, inputs, outputs, arg, ctx);
                        if (success)
                            ca_add(res, res, t, ctx);
                    }
                    ca_clear(t, ctx);
                }
                return success;

            case FEXPR_Mul:
                if (nargs == 0)
                {
                    ca_one(res, ctx);
                    return 1;
                }

                fexpr_view_arg(arg, expr, 0);
                success = _ca_set_fexpr(res, inputs, outputs, arg, ctx);

                if (success && nargs > 1)
                {
                    /* todo: divide and conquer for large nargs? */
                    ca_init(t, ctx);
                    for (i = 1; i < nargs && success; i++)
                    {
                        fexpr_view_next(arg);
                        success = _ca_set_fexpr(t, inputs, outputs, arg, ctx);
                        if (success)
                            ca_mul(res, res, t, ctx);
                    }
                    ca_clear(t, ctx);
                }
                return success;
        }
    }

    return 0;
}

int
ca_set_fexpr(ca_t res, const fexpr_t expr, ca_ctx_t ctx)
{
    int success;

    fexpr_vec_t inputs;
    ca_vec_t outputs;

    fexpr_vec_init(inputs, 0);
    ca_vec_init(outputs, 0, ctx);

    success = _ca_set_fexpr(res, inputs, outputs, expr, ctx);

    fexpr_vec_clear(inputs);
    ca_vec_clear(outputs, ctx);

    return success;
}
