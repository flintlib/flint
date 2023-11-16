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
ca_mul_fmpq(ca_t res, const ca_t x, const fmpq_t y, ca_ctx_t ctx)
{
    ca_field_srcptr field;

    if (CA_IS_SPECIAL(x))
    {
        if (CA_IS_SIGNED_INF(x))
        {
            if (fmpq_is_zero(y))
                ca_undefined(res, ctx);
            else if (fmpq_sgn(y) > 0)
                ca_set(res, x, ctx);
            else
                ca_neg(res, x, ctx);
        }
        else if (CA_IS_UNSIGNED_INF(x))
        {
            if (fmpq_is_zero(y))
                ca_undefined(res, ctx);
            else
                ca_set(res, x, ctx);
        }
        else
        {
            ca_set(res, x, ctx);
        }

        return;
    }

    if (fmpq_is_zero(y))
    {
        ca_zero(res, ctx);
        return;
    }

    if (CA_IS_QQ(x, ctx))
    {
        _ca_make_fmpq(res, ctx);
        fmpq_mul(CA_FMPQ(res), CA_FMPQ(x), y);
    }
    else
    {
        field = CA_FIELD(x, ctx);
        _ca_make_field_element(res, field, ctx);

        if (CA_FIELD_IS_NF(field))
        {
            nf_elem_scalar_mul_fmpq(CA_NF_ELEM(res), CA_NF_ELEM(x), y, CA_FIELD_NF(field));
        }
        else
        {
            fmpz_mpoly_q_mul_fmpq(CA_MPOLY_Q(res), CA_MPOLY_Q(x), y, CA_FIELD_MCTX(field, ctx));
        }
    }
}

void
ca_mul_fmpz(ca_t res, const ca_t x, const fmpz_t y, ca_ctx_t ctx)
{
    fmpq_t t;
    *fmpq_numref(t) = *y;
    *fmpq_denref(t) = 1;
    ca_mul_fmpq(res, x, t, ctx);
}

void
ca_mul_ui(ca_t res, const ca_t x, ulong y, ca_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init_set_ui(t, y);
    ca_mul_fmpz(res, x, t, ctx);
    fmpz_clear(t);
}

void
ca_mul_si(ca_t res, const ca_t x, slong y, ca_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init_set_si(t, y);
    ca_mul_fmpz(res, x, t, ctx);
    fmpz_clear(t);
}

void
ca_mul(ca_t res, const ca_t x, const ca_t y, ca_ctx_t ctx)
{
    ca_field_srcptr xfield, yfield, zfield;

    xfield = CA_FIELD(x, ctx);
    yfield = CA_FIELD(y, ctx);

    if (CA_IS_QQ(x, ctx) && (xfield == yfield))
    {
        _ca_make_fmpq(res, ctx);
        fmpq_mul(CA_FMPQ(res), CA_FMPQ(x), CA_FMPQ(y));
        return;
    }

    if (CA_IS_QQ(y, ctx))
    {
        if (res == y)
        {
            fmpq_t t;
            fmpq_init(t);
            fmpq_set(t, CA_FMPQ(y));
            ca_mul_fmpq(res, x, t, ctx);
            fmpq_clear(t);
        }
        else
        {
            ca_mul_fmpq(res, x, CA_FMPQ(y), ctx);
        }
        return;
    }

    if (CA_IS_QQ(x, ctx))
    {
        if (res == x)
        {
            fmpq_t t;
            fmpq_init(t);
            fmpq_set(t, CA_FMPQ(x));
            ca_mul_fmpq(res, y, t, ctx);
            fmpq_clear(t);
        }
        else
        {
            ca_mul_fmpq(res, y, CA_FMPQ(x), ctx);
        }
        return;
    }

    if (CA_IS_SPECIAL(x) || CA_IS_SPECIAL(y))
    {
        if (CA_IS_UNDEFINED(x) || CA_IS_UNDEFINED(y))
        {
            ca_undefined(res, ctx);
            return;
        }

        if (CA_IS_UNKNOWN(x) || CA_IS_UNKNOWN(y))
        {
            ca_unknown(res, ctx);
            return;
        }

        if (CA_IS_UNSIGNED_INF(x) && CA_IS_INF(y))
        {
            ca_uinf(res, ctx);
            return;
        }

        if (CA_IS_UNSIGNED_INF(y) && CA_IS_INF(x))
        {
            ca_uinf(res, ctx);
            return;
        }

        if (CA_IS_UNSIGNED_INF(x) && !CA_IS_SPECIAL(y))
        {
            truth_t zero = ca_check_is_zero(y, ctx);

            if (zero == T_TRUE)
                ca_undefined(res, ctx);
            else if (zero == T_FALSE)
                ca_uinf(res, ctx);
            else
                ca_unknown(res, ctx);
            return;
        }

        if (CA_IS_UNSIGNED_INF(y) && !CA_IS_SPECIAL(x))
        {
            truth_t zero = ca_check_is_zero(x, ctx);

            if (zero == T_TRUE)
                ca_undefined(res, ctx);
            else if (zero == T_FALSE)
                ca_uinf(res, ctx);
            else
                ca_unknown(res, ctx);
            return;
        }

        if (CA_IS_SIGNED_INF(x) || CA_IS_SIGNED_INF(y))
        {
            ca_t t, u;
            truth_t xzero, yzero;

            xzero = ca_check_is_zero(x, ctx);
            yzero = ca_check_is_zero(y, ctx);

            if (xzero == T_TRUE || yzero == T_TRUE)
            {
                ca_undefined(res, ctx);
                return;
            }

            if (xzero == T_UNKNOWN || yzero == T_UNKNOWN)
            {
                ca_unknown(res, ctx);
                return;
            }

            ca_init(t, ctx);
            ca_init(u, ctx);
            ca_sgn(t, x, ctx);
            ca_sgn(u, y, ctx);

            ca_mul(res, t, u, ctx);
            if (ca_check_is_number(res, ctx) == T_TRUE)
                res->field |= CA_INF;

            ca_clear(t, ctx);
            ca_clear(u, ctx);

            return;
        }

        ca_unknown(res, ctx);
        return;
    }

    if (xfield == yfield)
    {
        zfield = xfield;
        _ca_make_field_element(res, zfield, ctx);

        if (CA_FIELD_IS_NF(zfield))
        {
            nf_elem_mul(CA_NF_ELEM(res), CA_NF_ELEM(x), CA_NF_ELEM(y), CA_FIELD_NF(zfield));
        }
        else
        {
            fmpz_mpoly_q_mul(CA_MPOLY_Q(res), CA_MPOLY_Q(x), CA_MPOLY_Q(y), CA_FIELD_MCTX(zfield, ctx));
            _ca_mpoly_q_reduce_ideal(CA_MPOLY_Q(res), zfield, ctx);
            _ca_mpoly_q_simplify_fraction_ideal(CA_MPOLY_Q(res), zfield, ctx);
        }

        ca_condense_field(res, ctx);

        return;
    }

    {
        ca_t t, u;

        ca_init(t, ctx);
        ca_init(u, ctx);

        ca_merge_fields(t, u, x, y, ctx);
        ca_mul(res, t, u, ctx);
        ca_condense_field(res, ctx);

        ca_clear(t, ctx);
        ca_clear(u, ctx);
    }
}

