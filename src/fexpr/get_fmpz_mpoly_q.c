/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_q.h"
#include "fexpr.h"
#include "fexpr_builtin.h"

int
fexpr_get_fmpz_mpoly_q(fmpz_mpoly_q_t res, const fexpr_t expr, const fexpr_vec_t vars, const fmpz_mpoly_ctx_t ctx)
{
    if (fexpr_is_integer(expr))
    {
        fmpz_t c;
        fmpz_init(c);
        fexpr_get_fmpz(c, expr);
        fmpz_mpoly_q_set_fmpz(res, c, ctx);
        fmpz_clear(c);
        return 1;
    }
    else
    {
        slong i, nargs;
        ulong op_head;
        fexpr_t func, arg;
        int success;

        if (fexpr_is_arithmetic_operation(expr))
        {
            fmpz_mpoly_q_t A, B;

            nargs = fexpr_nargs(expr);

            if (nargs == 0)
            {
                /* not implemented */
                return 0;
            }

            success = 1;

            fmpz_mpoly_q_init(A, ctx);
            fmpz_mpoly_q_init(B, ctx);

            fexpr_view_func(func, expr);
            op_head = func->data[0];

            fexpr_view_arg(arg, expr, 0);
            success = fexpr_get_fmpz_mpoly_q(res, arg, vars, ctx);
            if (!success)
                goto cleanup1;

            /* todo: verify nargs == 1 */
            if (op_head == FEXPR_SYMBOL_Neg)
            {
                fmpz_mpoly_q_neg(res, res, ctx);
                goto cleanup1;
            }

            for (i = 1; i < nargs; i++)
            {
                fexpr_view_next(arg);

                success = fexpr_get_fmpz_mpoly_q(A, arg, vars, ctx);
                if (!success)
                    goto cleanup1;

                if (op_head == FEXPR_SYMBOL_Add)
                    fmpz_mpoly_q_add(B, res, A, ctx);
                else if (op_head == FEXPR_SYMBOL_Sub)
                    fmpz_mpoly_q_sub(B, res, A, ctx);
                else if (op_head == FEXPR_SYMBOL_Mul)
                    fmpz_mpoly_q_mul(B, res, A, ctx);
                else if (op_head == FEXPR_SYMBOL_Div)
                {
                    /* formal division by zero */
                    if (fmpz_mpoly_q_is_zero(A, ctx))
                        success = 0;
                    else
                        fmpz_mpoly_q_div(B, res, A, ctx);
                }

                fmpz_mpoly_q_swap(res, B, ctx);
            }

cleanup1:
            fmpz_mpoly_q_clear(A, ctx);
            fmpz_mpoly_q_clear(B, ctx);

            return success;
        }

        if (fexpr_is_builtin_call(expr, FEXPR_Pow) && (fexpr_nargs(expr) == 2))
        {
            fexpr_t base, exp;
            int success;

            fexpr_view_arg(base, expr, 0);
            fexpr_view_arg(exp, expr, 1);

            if (fexpr_is_integer(exp))
            {
                fmpz_t c;

                success = fexpr_get_fmpz_mpoly_q(res, base, vars, ctx);
                if (!success)
                    return 0;

                fmpz_init(c);
                fexpr_get_fmpz(c, exp);
                if (fmpz_sgn(c) < 0)
                {
                    if (fmpz_mpoly_q_is_zero(res, ctx))
                    {
                        success = 0;
                        goto cleanup2;
                    }

                    fmpz_neg(c, c);
                    fmpz_mpoly_q_inv(res, res, ctx);
                }

                success = (fmpz_mpoly_pow_fmpz(fmpz_mpoly_q_numref(res), fmpz_mpoly_q_numref(res), c, ctx) &&
                           fmpz_mpoly_pow_fmpz(fmpz_mpoly_q_denref(res), fmpz_mpoly_q_denref(res), c, ctx));

cleanup2:
                fmpz_clear(c);
                return success;
            }
        }

        for (i = 0; i < vars->length; i++)
        {
            if (fexpr_equal(expr, fexpr_vec_entry(vars, i)))
            {
                fmpz_mpoly_q_gen(res, i, ctx);
                return 1;
            }
        }

        return 0;
    }
}
