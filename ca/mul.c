/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"

void
ca_mul_fmpq(ca_t res, const ca_t x, const fmpq_t y, ca_ctx_t ctx)
{
    ulong xfield;
    slong zfield;
    ca_field_type_t type;

    xfield = x->field;

    if (CA_IS_SPECIAL(x))
    {
        if (xfield & CA_SIGNED_INF)
        {
            if (fmpq_is_zero(y))
                ca_undefined(res, ctx);
            else if (fmpq_sgn(y) > 0)
                ca_set(res, x, ctx);
            else
                ca_neg(res, x, ctx);
        }
        else if (xfield & CA_UNSIGNED_INF)
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

    if (xfield == CA_FIELD_ID_QQ)
    {
        _ca_make_fmpq(res, ctx);
        fmpq_mul(CA_FMPQ(res), CA_FMPQ(x), y);
        return;
    }

    zfield = xfield;
    type = ctx->fields[zfield].type;
    _ca_make_field_element(res, zfield, ctx);

    if (type == CA_FIELD_TYPE_NF)
    {
        nf_elem_scalar_mul_fmpq(CA_NF_ELEM(res), CA_NF_ELEM(x), y, CA_FIELD_NF(ctx->fields + zfield));
    }
    else if (type == CA_FIELD_TYPE_FUNC)
    {
        fmpz_mpoly_q_mul_fmpq(CA_MPOLY_Q(res), CA_MPOLY_Q(x), y, ctx->mctx + 0);
    }
    else if (type == CA_FIELD_TYPE_MULTI)
    {
        fmpz_mpoly_q_mul_fmpq(CA_MPOLY_Q(res), CA_MPOLY_Q(x), y, CA_FIELD_MCTX(ctx->fields + zfield, ctx));
    }
    else
    {
        flint_printf("ca_mul_fmpq: unknown field type\n");
        flint_abort();
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
    ulong xfield, yfield;
    slong zfield;
    ca_field_type_t type;

    xfield = x->field;
    yfield = y->field;

    if (xfield == CA_FIELD_ID_QQ && yfield == CA_FIELD_ID_QQ)
    {
        _ca_make_fmpq(res, ctx);
        fmpq_mul(CA_FMPQ(res), CA_FMPQ(x), CA_FMPQ(y));
        return;
    }

    if (yfield == CA_FIELD_ID_QQ)
    {
        ca_mul_fmpq(res, x, CA_FMPQ(y), ctx);
        return;
    }

    if (xfield == CA_FIELD_ID_QQ)
    {
        ca_mul_fmpq(res, y, CA_FMPQ(x), ctx);
        return;
    }

    if (CA_IS_SPECIAL(x) || CA_IS_SPECIAL(y))
    {
        if ((xfield & CA_UNDEFINED) || (yfield & CA_UNDEFINED))
        {
            ca_undefined(res, ctx);
            return;
        }

        if ((xfield & CA_UNSIGNED_INF) && ((yfield & CA_UNSIGNED_INF) || (yfield & CA_SIGNED_INF)))
        {
            ca_uinf(res, ctx);
            return;
        }

        if ((yfield & CA_UNSIGNED_INF) && ((xfield & CA_UNSIGNED_INF) || (xfield & CA_SIGNED_INF)))
        {
            ca_uinf(res, ctx);
            return;
        }

        if ((xfield & CA_UNSIGNED_INF) && !CA_IS_SPECIAL(y))
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

        if ((yfield & CA_UNSIGNED_INF) && !CA_IS_SPECIAL(x))
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

        if ((xfield & CA_SIGNED_INF) && (yfield & CA_SIGNED_INF))
        {
            /* ca_sign + mul */
        }

        if ((xfield & CA_SIGNED_INF) && !CA_IS_SPECIAL(y))
        {
            /* ca_sign + mul */
        }

        if ((yfield & CA_SIGNED_INF) && !CA_IS_SPECIAL(x))
        {
            /* ca_sign + mul */
        }

        ca_unknown(res, ctx);
        return;
    }

    if (xfield == yfield)
    {
        zfield = xfield;

        type = ctx->fields[zfield].type;
        _ca_make_field_element(res, zfield, ctx);

        if (type == CA_FIELD_TYPE_NF)
        {
            nf_elem_mul(CA_NF_ELEM(res), CA_NF_ELEM(x), CA_NF_ELEM(y), CA_FIELD_NF(ctx->fields + zfield));
        }
        else if (type == CA_FIELD_TYPE_FUNC)
        {
            fmpz_mpoly_q_mul(CA_MPOLY_Q(res), CA_MPOLY_Q(x), CA_MPOLY_Q(y), ctx->mctx + 0);
        }
        else if (type == CA_FIELD_TYPE_MULTI)
        {
            fmpz_mpoly_q_mul(CA_MPOLY_Q(res), CA_MPOLY_Q(x), CA_MPOLY_Q(y), CA_FIELD_MCTX(ctx->fields + zfield, ctx));
        }
        else
        {
            flint_printf("ca_mul: unknown field type\n");
            flint_abort();
        }

        ca_condense_field(res, ctx);

        return;
    }

    /* todo: subfields, merge fields */

    ca_unknown(res, ctx);
    return;
}

