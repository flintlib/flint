/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"
#include "ca_ext.h"
#include "ca_vec.h"

void
ca_sin_cos_special(ca_t res1, ca_t res2, const ca_t x, ca_ctx_t ctx)
{
    truth_t t1, t2;

    if (ca_check_is_signed_inf(x, ctx) == T_TRUE)
    {
        t1 = ca_check_is_pos_i_inf(x, ctx);

        if (t1 == T_TRUE)
        {
            if (res1 != NULL) ca_pos_i_inf(res1, ctx);
            if (res2 != NULL) ca_pos_inf(res2, ctx);
            return;
        }

        t2 = ca_check_is_neg_i_inf(x, ctx);

        if (t2 == T_TRUE)
        {
            if (res1 != NULL) ca_neg_i_inf(res1, ctx);
            if (res2 != NULL) ca_pos_inf(res2, ctx);
            return;
        }

        if (t1 == T_FALSE && t2 == T_FALSE)
        {
            if (res1 != NULL) ca_undefined(res1, ctx);
            if (res2 != NULL) ca_undefined(res2, ctx);
            return;
        }
    }

    if (ca_check_is_undefined(x, ctx) == T_TRUE || ca_check_is_uinf(x, ctx) == T_TRUE)
    {
        if (res1 != NULL) ca_undefined(res1, ctx);
        if (res2 != NULL) ca_undefined(res2, ctx);
        return;
    }

    if (res1 != NULL) ca_unknown(res1, ctx);
    if (res2 != NULL) ca_unknown(res2, ctx);
}

void
ca_sin_cos_exponential(ca_t res1, ca_t res2, const ca_t x, ca_ctx_t ctx)
{
    ca_t ix, y, t;

    if (CA_IS_SPECIAL(x))
    {
        ca_sin_cos_special(res1, res2, x, ctx);
        return;
    }

    ca_init(ix, ctx);
    ca_init(y, ctx);
    ca_init(t, ctx);

    ca_i(ix, ctx);
    ca_mul(ix, x, ix, ctx);
    ca_exp(y, ix, ctx);
    ca_inv(t, y, ctx);

    if (res2 != NULL)
    {
        ca_add(res2, y, t, ctx);
        ca_div_ui(res2, res2, 2, ctx);
    }

    if (res1 != NULL)
    {
        ca_sub(res1, y, t, ctx);
        ca_div_ui(res1, res1, 2, ctx);
        ca_neg_i(t, ctx);
        ca_mul(res1, res1, t, ctx);
    }

    ca_clear(ix, ctx);
    ca_clear(y, ctx);
    ca_clear(t, ctx);
}

void
ca_sin_cos_direct(ca_t res1, ca_t res2, const ca_t x, ca_ctx_t ctx)
{
    ca_t t, pi;
    fmpq_t v;

    if (CA_IS_SPECIAL(x))
    {
        ca_sin_cos_special(res1, res2, x, ctx);
        return;
    }

    ca_init(t, ctx);
    ca_init(pi, ctx);
    fmpq_init(v);

    ca_pi(pi, ctx);
    ca_div(t, x, pi, ctx);

    if (ca_get_fmpq(v, t, ctx))
    {
        /* For now, only convert trivial values */
        if (fmpz_cmp_ui(fmpq_denref(v), 6) <= 0 && !fmpz_equal_ui(fmpq_denref(v), 5))
        {
            slong p, q;
            qqbar_t a;

            q = fmpz_get_ui(fmpq_denref(v));
            p = fmpz_fdiv_ui(fmpq_numref(v), 2 * q);

            qqbar_init(a);

            if (res1 != NULL)
            {
                qqbar_sin_pi(a, p, q);
                ca_set_qqbar(res1, a, ctx);
            }

            if (res2 != NULL)
            {
                qqbar_cos_pi(a, p, q);
                ca_set_qqbar(res2, a, ctx);
            }

            qqbar_clear(a);
        }
        else
        {
            ca_mul_fmpq(t, pi, v, ctx);

            if (fmpq_sgn(v) > 0)
            {
                if (res1 != NULL) _ca_function_fx(res1, CA_Sin, t, ctx);
                if (res2 != NULL) _ca_function_fx(res2, CA_Cos, t, ctx);
            }
            else
            {
                ca_neg(t, t, ctx);
                if (res1 != NULL) _ca_function_fx(res1, CA_Sin, t, ctx);
                if (res2 != NULL) _ca_function_fx(res2, CA_Cos, t, ctx);
                if (res1 != NULL) ca_neg(res1, res1, ctx);
            }
        }
    }
    else
    {
        if (res1 != NULL) _ca_function_fx(res1, CA_Sin, x, ctx);
        if (res2 != NULL) _ca_function_fx(res2, CA_Cos, x, ctx);
    }

    ca_clear(pi, ctx);
    ca_clear(t, ctx);
    fmpq_clear(v);
}

#if 0
void
ca_sin_cos_direct_exp_hack(ca_t res1, ca_t res2, const ca_t x, ca_ctx_t ctx)
{
    ca_t i, ix, exp, sin, cos, t, u;

    ca_init(i, ctx);
    ca_init(ix, ctx);
    ca_init(exp, ctx);
    ca_init(sin, ctx);
    ca_init(cos, ctx);
    ca_init(t, ctx);
    ca_init(u, ctx);

    ca_i(i, ctx);
    ca_mul(ix, i, x, ctx);

    _ca_function_fx(sin, CA_Sin, x, ctx);
    _ca_function_fx(cos, CA_Cos, x, ctx);
    _ca_function_fx(exp, CA_Exp, ix, ctx);

    /* todo: create field in one step! */
    ca_merge_fields(t, u, sin, exp, ctx);
    ca_swap(sin, t, ctx);
    ca_merge_fields(t, u, sin, i, ctx);
    ca_swap(sin, t, ctx);
    ca_merge_fields(t, u, sin, cos, ctx);
    ca_swap(sin, t, ctx);
    ca_swap(cos, u, ctx);

    if (res1 != NULL) ca_swap(res1, sin, ctx);
    if (res2 != NULL) ca_swap(res2, cos, ctx);

    ca_clear(i, ctx);
    ca_clear(ix, ctx);
    ca_clear(exp, ctx);
    ca_clear(sin, ctx);
    ca_clear(cos, ctx);
    ca_clear(t, ctx);
    ca_clear(u, ctx);
}
#endif


void
ca_sin_cos_tangent(ca_t res1, ca_t res2, const ca_t x, ca_ctx_t ctx)
{
    ca_t t, u, v;

    if (CA_IS_SPECIAL(x))
    {
        ca_sin_cos_special(res1, res2, x, ctx);
        return;
    }

    ca_init(t, ctx);
    ca_init(u, ctx);
    ca_init(v, ctx);

    ca_div_ui(t, x, 2, ctx);
    ca_tan_direct(t, t, ctx);

    if (CA_IS_SPECIAL(t))
    {
        ca_sin_cos_direct(res1, res2, x, ctx);
    }
    else
    {
        ca_sqr(u, t, ctx);
        ca_add_ui(v, u, 1, ctx);
        ca_inv(v, v, ctx);

        if (res1 != NULL)
        {
            ca_mul(res1, t, v, ctx);
            ca_mul_ui(res1, res1, 2, ctx);
        }

        if (res2 != NULL)
        {
            ca_ui_sub(u, 1, u, ctx);
            ca_mul(res2, u, v, ctx);
        }
    }

    ca_clear(t, ctx);
    ca_clear(u, ctx);
    ca_clear(v, ctx);
}

void
ca_sin_cos(ca_t res1, ca_t res2, const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        ca_sin_cos_special(res1, res2, x, ctx);
        return;
    }

    if (CA_IS_QQ(x, ctx) && fmpq_is_zero(CA_FMPQ(x)))
    {
        if (res1 != NULL) ca_zero(res1, ctx);
        if (res2 != NULL) ca_one(res2, ctx);
        return;
    }

    if (ctx->options[CA_OPT_TRIG_FORM] == CA_TRIG_EXPONENTIAL)
    {
        ca_sin_cos_exponential(res1, res2, x, ctx);
        return;
    }

    if (ctx->options[CA_OPT_TRIG_FORM] == CA_TRIG_TANGENT)
    {
        ca_sin_cos_tangent(res1, res2, x, ctx);
        return;
    }

    ca_sin_cos_direct(res1, res2, x, ctx);
}

void
ca_sin(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    ca_sin_cos(res, NULL, x, ctx);
}

void
ca_cos(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    ca_sin_cos(NULL, res, x, ctx);
}

void
ca_tan_special(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    if (ca_check_is_signed_inf(x, ctx) == T_TRUE)
    {
        ca_t s;

        ca_init(s, ctx);
        ca_sgn(s, x, ctx);
        ca_im(s, s, ctx);
        ca_sgn(s, s, ctx);

        if (ca_check_is_one(s, ctx) == T_TRUE)
        {
            ca_i(res, ctx);
        }
        else if (ca_check_is_neg_one(s, ctx) == T_TRUE)
        {
            ca_neg_i(res, ctx);
        }
        else if (ca_check_is_zero(s, ctx) == T_TRUE)
        {
            ca_undefined(res, ctx);
        }
        else
        {
            ca_unknown(res, ctx);
        }

        ca_clear(s, ctx);
        return;
    }

    if (ca_is_unknown(res, ctx))
        ca_unknown(res, ctx);
    else
        ca_undefined(res, ctx);
}

void
ca_tan_direct(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    ca_t t, u;
    truth_t pole;

    if (CA_IS_SPECIAL(x))
    {
        ca_tan_special(res, x, ctx);
        return;
    }

    ca_init(t, ctx);
    ca_init(u, ctx);

    ca_pi(t, ctx);
    ca_div(t, x, t, ctx);

    if (ca_check_is_integer(t, ctx) == T_TRUE)
    {
        ca_zero(res, ctx);
    }
    else
    {
        ca_set_d(u, 0.5, ctx);
        ca_add(u, u, t, ctx);

        pole = ca_check_is_integer(u, ctx);

        if (pole == T_TRUE)
        {
            ca_uinf(res, ctx);
        }
        else if (pole == T_UNKNOWN)
        {
            ca_unknown(res, ctx);
        }
        else
        {
            fmpq_t v;
            fmpq_init(v);

            /* For now, only convert trivial values */
            if (ca_get_fmpq(v, t, ctx) &&
                (fmpz_equal_ui(fmpq_denref(v), 3) ||
                 fmpz_equal_ui(fmpq_denref(v), 4) ||
                 fmpz_equal_ui(fmpq_denref(v), 6) ||
                 fmpz_equal_ui(fmpq_denref(v), 8) ||
                 fmpz_equal_ui(fmpq_denref(v), 12)))
            {
                slong p, q;
                qqbar_t a;

                q = fmpz_get_ui(fmpq_denref(v));
                p = fmpz_fdiv_ui(fmpq_numref(v), q);

                qqbar_init(a);
                qqbar_tan_pi(a, p, q);
                ca_set_qqbar(res, a, ctx);
                qqbar_clear(a);
            }
            else
            {
                /* todo: sign symmetry, ... */
                _ca_function_fx(res, CA_Tan, x, ctx);
            }

            fmpq_clear(v);
        }
    }

    ca_clear(t, ctx);
    ca_clear(u, ctx);
}

void
ca_tan_exponential(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        ca_tan_special(res, x, ctx);
    }
    else
    {
        ca_t s, c;
        ca_init(s, ctx);
        ca_init(c, ctx);
        ca_sin_cos_exponential(s, c, x, ctx);
        ca_div(res, s, c, ctx);
        ca_clear(s, ctx);
        ca_clear(c, ctx);
        return;
    }
}

void
ca_tan_sine_cosine(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        ca_tan_special(res, x, ctx);
    }
    else
    {
        ca_t s, c;
        ca_init(s, ctx);
        ca_init(c, ctx);
        ca_sin_cos_direct(s, c, x, ctx);
        ca_div(res, s, c, ctx);
        ca_clear(s, ctx);
        ca_clear(c, ctx);
        return;
    }
}

void
ca_tan(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        ca_tan_special(res, x, ctx);
    }
    else if (ctx->options[CA_OPT_TRIG_FORM] == CA_TRIG_EXPONENTIAL)
    {
        ca_tan_exponential(res, x, ctx);
    }
    else if (ctx->options[CA_OPT_TRIG_FORM] == CA_TRIG_SINE_COSINE)
    {
        ca_tan_sine_cosine(res, x, ctx);
    }
    else
    {
        ca_tan_direct(res, x, ctx);
    }
}

void
ca_cot(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    ca_tan(res, x, ctx);
    ca_inv(res, res, ctx);
}
