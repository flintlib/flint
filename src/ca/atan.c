/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"

void
ca_atan_special(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    if (ca_check_is_signed_inf(x, ctx) == T_TRUE)
    {
        ca_t s;
        ca_init(s, ctx);
        ca_csgn(s, x, ctx);

        if (ca_check_is_one(s, ctx) == T_TRUE)
        {
            ca_pi(res, ctx);
            ca_div_ui(res, res, 2, ctx);
        }
        else if (ca_check_is_neg_one(s, ctx) == T_TRUE)
        {
            ca_pi(res, ctx);
            ca_div_ui(res, res, 2, ctx);
            ca_neg(res, res, ctx);
        }
        else
        {
            ca_unknown(res, ctx);
        }

        ca_clear(s, ctx);
        return;
    }

    if (ca_check_is_uinf(x, ctx) == T_TRUE || ca_check_is_undefined(x, ctx) == T_TRUE)
    {
        ca_undefined(res, ctx);
        return;
    }

    ca_unknown(res, ctx);
    return;
}

static int
_ca_atan_rational(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    qqbar_t v;
    slong p;
    ulong q;
    int success;

    qqbar_init(v);

    /* todo: rule out non-tangents more quickly */
    if (ca_get_qqbar(v, x, ctx) && qqbar_atan_pi(&p, &q, v))
    {
        ca_pi(res, ctx);
        ca_mul_si(res, res, p, ctx);
        ca_div_ui(res, res, q, ctx);
        success = 1;
    }
    else
    {
        success = 0;
    }

    qqbar_clear(v);
    return success;
}

void
ca_atan_logarithm(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    ca_t t, u, v;
    acb_t z;
    arb_t one, minus_one;

    if (CA_IS_SPECIAL(x))
    {
        ca_atan_special(res, x, ctx);
        return;
    }

    if (_ca_atan_rational(res, x, ctx))
        return;

    acb_init(z);
    arb_init(one);
    arb_init(minus_one);

    ca_init(t, ctx);
    ca_init(u, ctx);
    ca_init(v, ctx);

    ca_i(t, ctx);
    ca_mul(u, x, t, ctx);
    /* v = 1 + i x */
    ca_add_ui(v, u, 1, ctx);
    /* res = 1 - i x */
    ca_sub_ui(res, u, 1, ctx);
    ca_neg(res, res, ctx);

    ca_get_acb(z, x, ctx->options[CA_OPT_LOW_PREC], ctx);
    arb_set_si(one, 1);
    arb_set_si(minus_one, -1);

    if (arb_lt(acb_imagref(z), one))
    {
        /* atan(x) = i/2 log((1-ix)/(1+ix)) */
        ca_div(res, res, v, ctx);
        ca_log(res, res, ctx);
        ca_mul(res, res, t, ctx);
        ca_div_ui(res, res, 2, ctx);
    }
    else if (arb_gt(acb_imagref(z), minus_one))
    {
        /* atan(x) = -i/2 log((1+ix)/(1-ix)) */
        ca_div(res, v, res, ctx);
        ca_log(res, res, ctx);
        ca_mul(res, res, t, ctx);
        ca_div_ui(res, res, 2, ctx);
        ca_neg(res, res, ctx);
    }
    else
    {
        /* atan(x) = i/2 (log(1-ix) - log(1+ix)) */
        ca_log(res, res, ctx);
        ca_log(v, v, ctx);
        ca_sub(res, res, v, ctx);
        ca_mul(res, res, t, ctx);
        ca_div_ui(res, res, 2, ctx);
    }

    ca_clear(t, ctx);
    ca_clear(u, ctx);
    ca_clear(v, ctx);

    acb_clear(z);
    arb_clear(one);
    arb_clear(minus_one);
}

void
ca_atan_direct(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    truth_t pole;

    if (CA_IS_SPECIAL(x))
    {
        ca_atan_special(res, x, ctx);
        return;
    }

    if (_ca_atan_rational(res, x, ctx))
        return;

    pole = ca_check_is_i(x, ctx);

    if (pole == T_TRUE)
    {
        ca_pos_i_inf(res, ctx);
        return;
    }

    if (pole == T_UNKNOWN)
    {
        ca_unknown(res, ctx);
        return;
    }

    pole = ca_check_is_neg_i(x, ctx);

    if (pole == T_TRUE)
    {
        ca_neg_i_inf(res, ctx);
        return;
    }

    if (pole == T_UNKNOWN)
    {
        ca_unknown(res, ctx);
        return;
    }

    /* todo: csgn normalization, reflection...? */
    _ca_function_fx(res, CA_Atan, x, ctx);
}

void
ca_atan(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    if (ctx->options[CA_OPT_TRIG_FORM] == CA_TRIG_EXPONENTIAL)
    {
        ca_atan_logarithm(res, x, ctx);
    }
    else
    {
        ca_atan_direct(res, x, ctx);
    }
}
