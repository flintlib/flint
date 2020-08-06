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
ca_div_fmpq(ca_t res, const ca_t x, const fmpq_t y, ca_ctx_t ctx)
{
    ulong xfield;
    slong zfield;

    xfield = x->field;

    if (CA_IS_SPECIAL(x))
    {
        if (xfield & CA_SIGNED_INF)
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

    if (xfield == CA_FIELD_ID_QQ)
    {
        _ca_make_fmpq(res, ctx);
        fmpq_div(CA_FMPQ(res), CA_FMPQ(x), y);
        return;
    }

    zfield = xfield;
    _ca_make_field_element(res, zfield, ctx);

    if (CA_FIELD_IS_NF(ctx->fields + zfield))
    {
        nf_elem_scalar_div_fmpq(CA_NF_ELEM(res), CA_NF_ELEM(x), y, CA_FIELD_NF(ctx->fields + zfield));
    }
    else
    {
        fmpz_mpoly_q_div_fmpq(CA_MPOLY_Q(res), CA_MPOLY_Q(x), y, CA_FIELD_MCTX(ctx->fields + zfield, ctx));
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
    ulong xfield, yfield;
    slong zfield;
    truth_t x_is_zero, y_is_zero;

    xfield = x->field;
    yfield = y->field;

    if (xfield == CA_FIELD_ID_QQ && yfield == CA_FIELD_ID_QQ)
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

    if (yfield == CA_FIELD_ID_QQ)
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

        if (CA_FIELD_IS_NF(ctx->fields + zfield))
        {
            nf_elem_div(CA_NF_ELEM(res), CA_NF_ELEM(x), CA_NF_ELEM(y), CA_FIELD_NF(ctx->fields + zfield));
        }
        else
        {
            slong i, n;

            fmpz_mpoly_q_div(CA_MPOLY_Q(res), CA_MPOLY_Q(x), CA_MPOLY_Q(y), CA_FIELD_MCTX(ctx->fields + zfield, ctx));

            /* todo: improve, deduplicate this code */

            n = CA_FIELD_IDEAL_LENGTH(ctx->fields + zfield);

            if (n != 0)
            {
                fmpz_mpoly_struct ** I;
                fmpz_mpoly_struct ** Q;
                fmpq_t scale;

                I = flint_malloc(sizeof(fmpz_mpoly_struct *) * n);
                for (i = 0; i < n; i++)
                    I[i] = CA_FIELD_IDEAL_POLY(ctx->fields + zfield, i);

                Q = flint_malloc(sizeof(fmpz_mpoly_struct *) * n);
                for (i = 0; i < n; i++)
                {
                    Q[i] = flint_malloc(sizeof(fmpz_mpoly_struct));
                    fmpz_mpoly_init(Q[i], CA_FIELD_MCTX(ctx->fields + zfield, ctx));
                }

                fmpq_init(scale);

                fmpz_mpoly_quasidivrem_ideal(fmpq_denref(scale), Q, fmpz_mpoly_q_numref(CA_MPOLY_Q(res)), fmpz_mpoly_q_numref(CA_MPOLY_Q(res)), I, n, CA_FIELD_MCTX(ctx->fields + zfield, ctx));
                fmpz_mpoly_quasidivrem_ideal(fmpq_numref(scale), Q, fmpz_mpoly_q_denref(CA_MPOLY_Q(res)), fmpz_mpoly_q_denref(CA_MPOLY_Q(res)), I, n, CA_FIELD_MCTX(ctx->fields + zfield, ctx));

                fmpq_canonicalise(scale);
                fmpz_mpoly_q_canonicalise(CA_MPOLY_Q(res), CA_FIELD_MCTX(ctx->fields + zfield, ctx));
                fmpz_mpoly_q_mul_fmpq(CA_MPOLY_Q(res), CA_MPOLY_Q(res), scale, CA_FIELD_MCTX(ctx->fields + zfield, ctx));

                for (i = 0; i < n; i++)
                {
                    fmpz_mpoly_clear(Q[i], CA_FIELD_MCTX(ctx->fields + zfield, ctx));
                    flint_free(Q[i]);
                }

                flint_free(Q);
                flint_free(I);

                fmpq_clear(scale);
            }
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

