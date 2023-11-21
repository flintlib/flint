/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"

void
ca_div_fmpq(ca_t res, const ca_t x, const fmpq_t y, ca_ctx_t ctx)
{
    ca_field_srcptr field;

    if (CA_IS_SPECIAL(x))
    {
        if (CA_IS_SIGNED_INF(x))
        {
            if (fmpq_is_zero(y))
                ca_uinf(res, ctx);
            else if (fmpq_sgn(y) > 0)
                ca_set(res, x, ctx);
            else
                ca_neg(res, x, ctx);
        }
        else
        {
            ca_set(res, x, ctx);
        }

        return;
    }

    if (fmpq_is_zero(y))
    {
        truth_t x_zero = ca_check_is_zero(x, ctx);

        if (x_zero == T_TRUE)
            ca_undefined(res, ctx);
        else if (x_zero == T_FALSE)
            ca_uinf(res, ctx);
        else
            ca_unknown(res, ctx);
        return;
    }

    if (CA_IS_QQ(x, ctx))
    {
        _ca_make_fmpq(res, ctx);
        fmpq_div(CA_FMPQ(res), CA_FMPQ(x), y);
        return;
    }
    else
    {
        field = CA_FIELD(x, ctx);
        _ca_make_field_element(res, field, ctx);

        if (CA_FIELD_IS_NF(field))
        {
            nf_elem_scalar_div_fmpq(CA_NF_ELEM(res), CA_NF_ELEM(x), y, CA_FIELD_NF(field));
        }
        else
        {
            fmpz_mpoly_q_div_fmpq(CA_MPOLY_Q(res), CA_MPOLY_Q(x), y, CA_FIELD_MCTX(field, ctx));
        }
    }
}

void
ca_div_fmpz(ca_t res, const ca_t x, const fmpz_t y, ca_ctx_t ctx)
{
    fmpq_t t;
    *fmpq_numref(t) = *y;
    *fmpq_denref(t) = 1;
    ca_div_fmpq(res, x, t, ctx);
}

void
ca_div_ui(ca_t res, const ca_t x, ulong y, ca_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init_set_ui(t, y);
    ca_div_fmpz(res, x, t, ctx);
    fmpz_clear(t);
}

void
ca_div_si(ca_t res, const ca_t x, slong y, ca_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init_set_si(t, y);
    ca_div_fmpz(res, x, t, ctx);
    fmpz_clear(t);
}

void
ca_fmpq_div(ca_t res, const fmpq_t x, const ca_t y, ca_ctx_t ctx)
{
    ca_t t;
    ca_init(t, ctx);
    ca_set_fmpq(t, x, ctx);
    ca_div(res, t, y, ctx);
    ca_clear(t, ctx);
}

void
ca_fmpz_div(ca_t res, const fmpz_t x, const ca_t y, ca_ctx_t ctx)
{
    ca_t t;
    ca_init(t, ctx);
    ca_set_fmpz(t, x, ctx);
    ca_div(res, t, y, ctx);
    ca_clear(t, ctx);
}

void
ca_si_div(ca_t res, slong x, const ca_t y, ca_ctx_t ctx)
{
    ca_t t;
    ca_init(t, ctx);
    ca_set_si(t, x, ctx);
    ca_div(res, t, y, ctx);
    ca_clear(t, ctx);
}

void
ca_ui_div(ca_t res, ulong x, const ca_t y, ca_ctx_t ctx)
{
    ca_t t;
    ca_init(t, ctx);
    ca_set_ui(t, x, ctx);
    ca_div(res, t, y, ctx);
    ca_clear(t, ctx);
}

void
ca_div(ca_t res, const ca_t x, const ca_t y, ca_ctx_t ctx)
{
    ca_field_srcptr xfield, yfield, zfield;
    truth_t x_is_zero, y_is_zero;

    xfield = CA_FIELD(x, ctx);
    yfield = CA_FIELD(y, ctx);

    if (CA_IS_QQ(x, ctx) && (xfield == yfield))
    {
        if (fmpq_is_zero(CA_FMPQ(y)))
        {
            if (fmpq_is_zero(CA_FMPQ(x)))
                ca_undefined(res, ctx);
            else
                ca_uinf(res, ctx);
        }
        else
        {
            _ca_make_fmpq(res, ctx);
            fmpq_div(CA_FMPQ(res), CA_FMPQ(x), CA_FMPQ(y));
        }
        return;
    }

    if (CA_IS_QQ(y, ctx))
    {
        if (res == y)
        {
            fmpq_t t;
            fmpq_init(t);
            fmpq_set(t, CA_FMPQ(y));
            ca_div_fmpq(res, x, t, ctx);
            fmpq_clear(t);
        }
        else
        {
            ca_div_fmpq(res, x, CA_FMPQ(y), ctx);
        }
        return;
    }

    if (CA_IS_SPECIAL(x) || CA_IS_SPECIAL(y))
    {
        ca_t t;
        ca_init(t, ctx);
        ca_inv(t, y, ctx);
        ca_mul(res, x, t, ctx);
        ca_clear(t, ctx);
        return;
    }

    y_is_zero = ca_check_is_zero(y, ctx);

    if (y_is_zero == T_TRUE)
    {
        x_is_zero = ca_check_is_zero(x, ctx);

        if (x_is_zero == T_FALSE)
            ca_uinf(res, ctx);
        else if (x_is_zero == T_TRUE)
            ca_undefined(res, ctx);
        else
            ca_unknown(res, ctx);
        return;
    }
    else if (y_is_zero == T_UNKNOWN)
    {
        ca_unknown(res, ctx);
        return;
    }

    if (xfield == yfield)
    {
        zfield = xfield;
        _ca_make_field_element(res, zfield, ctx);

        if (CA_FIELD_IS_NF(zfield))
        {
            nf_elem_div(CA_NF_ELEM(res), CA_NF_ELEM(x), CA_NF_ELEM(y), CA_FIELD_NF(zfield));
        }
        else
        {
            fmpz_mpoly_q_div(CA_MPOLY_Q(res), CA_MPOLY_Q(x), CA_MPOLY_Q(y), CA_FIELD_MCTX(zfield, ctx));
            _ca_mpoly_q_reduce_ideal(CA_MPOLY_Q(res), zfield, ctx);
            _ca_mpoly_q_simplify_fraction_ideal(CA_MPOLY_Q(res), zfield, ctx);
        }

        ca_condense_field(res, ctx);
        return;
    }

    {
        ca_t t;
        ca_init(t, ctx);
        ca_inv(t, y, ctx);
        ca_mul(res, x, t, ctx);
        ca_clear(t, ctx);
    }
}
