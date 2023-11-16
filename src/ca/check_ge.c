/*
    Copyright (C) 2020, 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"

#define CMP_UNDEFINED -2
#define CMP_UNKNOWN -3

int
_ca_cmp(const ca_t x, const ca_t y, ca_ctx_t ctx)
{
    acb_t v, w;
    truth_t x_real, y_real;
    slong prec, prec_limit;
    int result;

    if (CA_IS_QQ(x, ctx) && CA_IS_QQ(y, ctx))
    {
        result = fmpq_cmp(CA_FMPQ(x), CA_FMPQ(y));
        if (result < 0) result = -1;
        if (result > 0) result = 1;
        return result;
    }

    if (CA_IS_SPECIAL(x) || CA_IS_SPECIAL(y))
    {
        if (ca_check_is_pos_inf(x, ctx) == T_TRUE)
        {
            if (ca_check_is_pos_inf(y, ctx) == T_TRUE)
                return 0;
            if (ca_check_is_neg_inf(y, ctx) == T_TRUE)
                return 1;
            y_real = ca_check_is_real(y, ctx);
            if (y_real == T_TRUE)
                return 1;
            if (y_real == T_FALSE)
                return CMP_UNDEFINED;
            return CMP_UNKNOWN;
        }

        if (ca_check_is_neg_inf(x, ctx) == T_TRUE)
        {
            if (ca_check_is_neg_inf(y, ctx) == T_TRUE)
                return 0;
            if (ca_check_is_pos_inf(y, ctx) == T_TRUE)
                return -1;
            y_real = ca_check_is_real(y, ctx);
            if (y_real == T_TRUE)
                return -1;
            if (y_real == T_FALSE)
                return CMP_UNDEFINED;
            return CMP_UNKNOWN;
        }

        if (ca_check_is_pos_inf(y, ctx) == T_TRUE)
        {
            if (ca_check_is_pos_inf(x, ctx) == T_TRUE)
                return 0;
            if (ca_check_is_neg_inf(x, ctx) == T_TRUE)
                return -1;
            x_real = ca_check_is_real(x, ctx);
            if (x_real == T_TRUE)
                return -1;
            if (x_real == T_FALSE)
                return CMP_UNDEFINED;
            return CMP_UNKNOWN;
        }

        if (ca_check_is_neg_inf(y, ctx) == T_TRUE)
        {
            if (ca_check_is_neg_inf(x, ctx) == T_TRUE)
                return 0;
            if (ca_check_is_pos_inf(x, ctx) == T_TRUE)
                return 1;
            x_real = ca_check_is_real(x, ctx);
            if (x_real == T_TRUE)
                return 1;
            if (x_real == T_FALSE)
                return CMP_UNDEFINED;
            return CMP_UNKNOWN;
        }

        /* Not comparable */
        if (ca_check_is_undefined(x, ctx) == T_TRUE ||
            ca_check_is_undefined(y, ctx) == T_TRUE ||
            ca_check_is_uinf(x, ctx) == T_TRUE ||
            ca_check_is_uinf(y, ctx) == T_TRUE ||
            (ca_check_is_signed_inf(x, ctx) == T_TRUE && ca_check_is_pos_inf(x, ctx) == T_FALSE && ca_check_is_neg_inf(x, ctx) == T_FALSE) ||
            (ca_check_is_signed_inf(y, ctx) == T_TRUE && ca_check_is_pos_inf(y, ctx) == T_FALSE && ca_check_is_neg_inf(y, ctx) == T_FALSE) ||
            (ca_check_is_number(x, ctx) == T_TRUE && ca_check_is_real(x, ctx) == T_FALSE) ||
            (ca_check_is_number(y, ctx) == T_TRUE && ca_check_is_real(y, ctx) == T_FALSE))
        {
            return CMP_UNDEFINED;
        }

        return CMP_UNKNOWN;
    }

    result = CMP_UNKNOWN;
    x_real = y_real = T_UNKNOWN;

    acb_init(v);
    acb_init(w);

    prec_limit = ctx->options[CA_OPT_PREC_LIMIT];
    prec_limit = FLINT_MAX(prec_limit, 64);

    for (prec = 64; (prec <= prec_limit) && (result == CMP_UNKNOWN); prec *= 2)
    {
        ca_get_acb_raw(v, x, prec, ctx);
        ca_get_acb_raw(w, y, prec, ctx);

        if (arb_is_zero(acb_imagref(v)))
            x_real = T_TRUE;
        else if (!arb_contains_zero(acb_imagref(v)))
            x_real = T_FALSE;

        if (arb_is_zero(acb_imagref(w)))
            y_real = T_TRUE;
        else if (!arb_contains_zero(acb_imagref(w)))
            y_real = T_FALSE;

        if (x_real == T_FALSE || y_real == T_FALSE)
        {
            /* Not comparable */
            result = CMP_UNDEFINED;
            break;
        }

        /* Force a verification that we have comparable numbers. */
        if (x_real == T_UNKNOWN && prec == 64)
            x_real = ca_check_is_real(x, ctx);
        if (y_real == T_UNKNOWN && prec == 64)
            y_real = ca_check_is_real(y, ctx);

        if (x_real == T_FALSE || y_real == T_FALSE)
        {
            /* Not comparable */
            result = CMP_UNDEFINED;
            break;
        }

        if (x_real == T_TRUE && y_real == T_TRUE)
        {
            if (arb_gt(acb_realref(v), acb_realref(w)))
            {
                result = 1;
                break;
            }
            else if (arb_lt(acb_realref(v), acb_realref(w)))
            {
                result = -1;
                break;
            }
        }

        /* try qqbar computation */
        /* todo: precision to do this should depend on complexity of the polynomials, degree of the elements... */
        if (prec == 64)
        {
            if (ca_can_evaluate_qqbar(x, ctx) && ca_can_evaluate_qqbar(y, ctx))
            {
                qqbar_t t, u;

                qqbar_init(t);
                qqbar_init(u);

                if (ca_get_qqbar(t, x, ctx))
                {
                    if (!qqbar_is_real(t))
                    {
                        /* Not comparable */
                        result = CMP_UNDEFINED;
                    }
                    else if (ca_get_qqbar(u, y, ctx))
                    {
                        if (!qqbar_is_real(u))
                        {
                            /* Not comparable */
                            result = CMP_UNDEFINED;
                        }
                        else
                        {
                            result = qqbar_cmp_re(t, u);
                            if (result < 0) result = -1;
                            if (result > 0) result = 1;
                        }
                    }
                }

                qqbar_clear(t);
                qqbar_clear(u);
            }
        }
    }

    /* Todo: subtract, compute sign? Symbolic cancellations may help. */

    if (result == CMP_UNKNOWN && x_real == T_TRUE && y_real == T_TRUE)
    {
        if (ca_check_equal(x, y, ctx) == T_TRUE)
            result = 0;
    }

    acb_clear(v);
    acb_clear(w);

    return result;
}

truth_t
ca_check_ge(const ca_t x, const ca_t y, ca_ctx_t ctx)
{
    int c = _ca_cmp(x, y, ctx);

    if (c == CMP_UNKNOWN)
        return T_UNKNOWN;

    if (c == CMP_UNDEFINED)
        return T_FALSE;

    return (c >= 0) ? T_TRUE : T_FALSE;
}

truth_t
ca_check_gt(const ca_t x, const ca_t y, ca_ctx_t ctx)
{
    int c = _ca_cmp(x, y, ctx);

    if (c == CMP_UNKNOWN)
        return T_UNKNOWN;

    if (c == CMP_UNDEFINED)
        return T_FALSE;

    return (c > 0) ? T_TRUE : T_FALSE;
}

truth_t
ca_check_lt(const ca_t x, const ca_t y, ca_ctx_t ctx)
{
    int c = _ca_cmp(x, y, ctx);

    if (c == CMP_UNKNOWN)
        return T_UNKNOWN;

    if (c == CMP_UNDEFINED)
        return T_FALSE;

    return (c < 0) ? T_TRUE : T_FALSE;
}

truth_t
ca_check_le(const ca_t x, const ca_t y, ca_ctx_t ctx)
{
    int c = _ca_cmp(x, y, ctx);

    if (c == CMP_UNKNOWN)
        return T_UNKNOWN;

    if (c == CMP_UNDEFINED)
        return T_FALSE;

    return (c <= 0) ? T_TRUE : T_FALSE;
}
