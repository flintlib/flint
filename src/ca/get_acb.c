/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_fmpz_poly.h"
#include "ca.h"

void
ca_get_acb(acb_t res, const ca_t x, slong prec, ca_ctx_t ctx)
{
    slong wp, initial, maxprec, exact_check_prec;

    acb_indeterminate(res);

    initial = prec * 1.05 + 30;
    maxprec = FLINT_MAX(2 * initial, ctx->options[CA_OPT_PREC_LIMIT]);
    exact_check_prec = FLINT_MIN(maxprec, initial * 16);

    for (wp = initial; wp <= maxprec; wp *= 2)
    {
        ca_get_acb_raw(res, x, wp, ctx);

        if (acb_rel_accuracy_bits(res) >= prec)
            break;

        if (wp == exact_check_prec)
        {
            if (ca_check_is_zero(x, ctx) == T_TRUE)
            {
                acb_zero(res);
                break;
            }
        }
    }
}

void
ca_get_acb_accurate_parts(acb_t res, const ca_t x, slong prec, ca_ctx_t ctx)
{
    slong wp, initial, maxprec, exact_check_prec;
    int re_ok, im_ok;

    acb_indeterminate(res);

    initial = prec * 1.05 + 30;
    maxprec = FLINT_MAX(2 * initial, ctx->options[CA_OPT_PREC_LIMIT]);
    exact_check_prec = FLINT_MIN(maxprec, initial * 16);

    for (wp = initial; wp <= maxprec; wp *= 2)
    {
        ca_get_acb_raw(res, x, wp, ctx);

        re_ok = arb_rel_accuracy_bits(acb_realref(res)) >= prec;
        im_ok = arb_rel_accuracy_bits(acb_imagref(res)) >= prec;

        if (re_ok && im_ok)
            break;

        if (wp == exact_check_prec)
        {
            if (acb_rel_accuracy_bits(res) < prec && ca_check_is_zero(x, ctx) == T_TRUE)
            {
                acb_zero(res);
                break;
            }

            if (re_ok && ca_check_is_real(x, ctx) == T_TRUE)
            {
                arb_zero(acb_imagref(res));
                break;
            }

            if (im_ok && ca_check_is_imaginary(x, ctx) == T_TRUE)
            {
                arb_zero(acb_realref(res));
                break;
            }
        }
    }
}

char *
ca_get_decimal_str(const ca_t x, slong digits, ulong flags, ca_ctx_t ctx)
{
    calcium_stream_t t;
    acb_t v;

    digits = FLINT_MAX(digits, 1);

    acb_init(v);
    calcium_stream_init_str(t);

    if (flags & 1)
        ca_get_acb_accurate_parts(v, x, digits * 3.333 + 1, ctx);
    else
        ca_get_acb(v, x, digits * 3.333 + 1, ctx);

    if (acb_is_finite(v))
        calcium_write_acb(t, v, digits, ARB_STR_NO_RADIUS);
    else
        calcium_write(t, "?");

    acb_clear(v);
    return t->s;
}
