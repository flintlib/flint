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
ca_asin_special(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    if (ca_check_is_signed_inf(x, ctx) == T_TRUE)
    {
        ca_t s;
        ca_init(s, ctx);

        /* -i*csgn(i*z)*inf */
        ca_i(s, ctx);
        ca_mul(res, x, s, ctx);
        ca_csgn(res, res, ctx);
        ca_mul(res, res, s, ctx);
        ca_neg(res, res, ctx);
        ca_pos_inf(s, ctx);
        ca_mul(res, res, s, ctx);

        ca_clear(s, ctx);
        return;
    }

    if (ca_check_is_uinf(x, ctx) == T_TRUE || ca_check_is_undefined(x, ctx) == T_TRUE)
    {
        ca_set(res, x, ctx);
        return;
    }

    ca_unknown(res, ctx);
    return;
}

static int
_ca_asin_rational(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    qqbar_t v;
    slong p;
    ulong q;
    int success;

    qqbar_init(v);

    /* todo: rule out non-sines more quickly */
    if (ca_get_qqbar(v, x, ctx) && qqbar_asin_pi(&p, &q, v))
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
ca_asin_logarithm(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    ca_t t, u;

    if (CA_IS_SPECIAL(x))
    {
        ca_asin_special(res, x, ctx);
        return;
    }

    if (_ca_asin_rational(res, x, ctx))
        return;

    /* asin(x) = -i log(ix + sqrt(1-x^2)) */

    ca_init(t, ctx);
    ca_init(u, ctx);

    ca_sqr(t, x, ctx);
    ca_ui_sub(t, 1, t, ctx);
    ca_sqrt(t, t, ctx);
    ca_i(u, ctx);
    ca_mul(u, u, x, ctx);
    ca_add(t, t, u, ctx);
    ca_log(t, t, ctx);
    ca_i(u, ctx);
    ca_mul(res, t, u, ctx);
    ca_neg(res, res, ctx);

    ca_clear(t, ctx);
    ca_clear(u, ctx);
}

void
ca_asin_direct(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        ca_asin_special(res, x, ctx);
        return;
    }

    if (_ca_asin_rational(res, x, ctx))
        return;

    /* todo: csgn normalization, reflection...? */
    _ca_function_fx(res, CA_Asin, x, ctx);
}

void
ca_asin(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    if (ctx->options[CA_OPT_TRIG_FORM] == CA_TRIG_EXPONENTIAL)
    {
        ca_asin_logarithm(res, x, ctx);
    }
/* todo:
    else if (ctx->options[CA_OPT_TRIG_FORM] == CA_TRIG_TANGENT)
    {
        ca_asin_arctangent(res, x, ctx);
    }
*/
    else
    {
        ca_asin_direct(res, x, ctx);
    }
}

void
ca_acos_logarithm(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    ca_t t;
    ca_init(t, ctx);
    ca_pi(t, ctx);
    ca_div_ui(t, t, 2, ctx);
    ca_asin_logarithm(res, x, ctx);
    ca_sub(res, t, res, ctx);
    ca_clear(t, ctx);
}

void
ca_acos_direct(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    ca_t t;
    ca_init(t, ctx);
    ca_pi(t, ctx);
    ca_div_ui(t, t, 2, ctx);
    ca_asin_direct(res, x, ctx);
    ca_sub(res, t, res, ctx);
    ca_clear(t, ctx);
}

void
ca_acos(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    ca_t t;
    ca_init(t, ctx);
    ca_pi(t, ctx);
    ca_div_ui(t, t, 2, ctx);
    ca_asin(res, x, ctx);
    ca_sub(res, t, res, ctx);
    ca_clear(t, ctx);
}
