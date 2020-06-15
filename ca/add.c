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
ca_add_fmpq(ca_t res, const ca_t x, const fmpq_t y, ca_ctx_t ctx)
{
    ulong xfield;
    slong zfield;
    ca_field_type_t type;

    if (fmpq_is_zero(y) || CA_IS_SPECIAL(x))
    {
        ca_set(res, x, ctx);
        return;
    }

    xfield = x->field;

    if (xfield == CA_FIELD_ID_QQ)
    {
        _ca_make_fmpq(res, ctx);
        fmpq_add(CA_FMPQ(res), CA_FMPQ(x), y);
        return;
    }

    zfield = xfield;
    type = ctx->fields[zfield].type;
    _ca_make_field_element(res, zfield, ctx);

    if (type == CA_FIELD_TYPE_NF)
    {
        nf_elem_add_fmpq(CA_NF_ELEM(res), CA_NF_ELEM(x), y, CA_FIELD_NF(ctx->fields + zfield));
    }
    else if (type == CA_FIELD_TYPE_FUNC)
    {
        fmpz_mpoly_q_add_fmpq(CA_MPOLY_Q(res), CA_MPOLY_Q(x), y, ctx->mctx + 0);
    }
    else if (type == CA_FIELD_TYPE_MULTI)
    {
        fmpz_mpoly_q_add_fmpq(CA_MPOLY_Q(res), CA_MPOLY_Q(x), y, CA_FIELD_MCTX(ctx->fields + zfield, ctx));
    }
    else
    {
        flint_printf("ca_add_fmpq: unknown field type\n");
        flint_abort();
    }
}

void
ca_add_fmpz(ca_t res, const ca_t x, const fmpz_t y, ca_ctx_t ctx)
{
    fmpq_t t;
    *fmpq_numref(t) = *y;
    *fmpq_denref(t) = 1;
    ca_add_fmpq(res, x, t, ctx);
}

void
ca_add_ui(ca_t res, const ca_t x, ulong y, ca_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init_set_ui(t, y);
    ca_add_fmpz(res, x, t, ctx);
    fmpz_clear(t);
}

void
ca_add_si(ca_t res, const ca_t x, slong y, ca_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init_set_si(t, y);
    ca_add_fmpz(res, x, t, ctx);
    fmpz_clear(t);
}

void
ca_add(ca_t res, const ca_t x, const ca_t y, ca_ctx_t ctx)
{
    ulong xfield, yfield;
    slong zfield;
    ca_field_type_t type;

    xfield = x->field;
    yfield = y->field;

    if (xfield == CA_FIELD_ID_QQ && yfield == CA_FIELD_ID_QQ)
    {
        _ca_make_fmpq(res, ctx);
        fmpq_add(CA_FMPQ(res), CA_FMPQ(x), CA_FMPQ(y));
        return;
    }

    if (yfield == CA_FIELD_ID_QQ)
    {
        if (res == y)
        {
            fmpq_t t;
            fmpq_init(t);
            fmpq_set(t, CA_FMPQ(y));
            ca_add_fmpq(res, x, t, ctx);
            fmpq_clear(t);
        }
        else
        {
            ca_add_fmpq(res, x, CA_FMPQ(y), ctx);
        }
        return;
    }

    if (xfield == CA_FIELD_ID_QQ)
    {
        if (res == x)
        {
            fmpq_t t;
            fmpq_init(t);
            fmpq_set(t, CA_FMPQ(x));
            ca_add_fmpq(res, y, t, ctx);
            fmpq_clear(t);
        }
        else
        {
            ca_add_fmpq(res, y, CA_FMPQ(x), ctx);
        }
        return;
    }

    if (CA_IS_SPECIAL(x) || CA_IS_SPECIAL(y))
    {
        if ((xfield & CA_UNDEFINED) || (yfield & CA_UNDEFINED))
        {
            ca_undefined(res, ctx);
            return;
        }

        if (!CA_IS_SPECIAL(y))
        {
            ca_set(res, x, ctx);
            return;
        }

        if (!CA_IS_SPECIAL(x))
        {
            ca_set(res, y, ctx);
            return;
        }

        if ((xfield & CA_SIGNED_INF) && (yfield & CA_SIGNED_INF))
        {
            truth_t eq = ca_check_equal(x, y, ctx);

            if (eq == T_TRUE)         /* sum of infs with same sign */
                ca_set(res, x, ctx);
            else if (eq == T_FALSE)   /* sum of infs with different sign */
                ca_undefined(res, ctx);
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
            nf_elem_add(CA_NF_ELEM(res), CA_NF_ELEM(x), CA_NF_ELEM(y), CA_FIELD_NF(ctx->fields + zfield));
        }
        else if (type == CA_FIELD_TYPE_FUNC)
        {
            fmpz_mpoly_q_add(CA_MPOLY_Q(res), CA_MPOLY_Q(x), CA_MPOLY_Q(y), ctx->mctx + 0);
        }
        else if (type == CA_FIELD_TYPE_MULTI)
        {
            slong i, n;

            fmpz_mpoly_q_add(CA_MPOLY_Q(res), CA_MPOLY_Q(x), CA_MPOLY_Q(y), CA_FIELD_MCTX(ctx->fields + zfield, ctx));

            /* todo: improve, deduplicate this code */

            n = ctx->fields[zfield].data.multi.ideal_len;

            if (n != 0)
            {
                fmpz_mpoly_struct ** I;
                fmpz_mpoly_struct ** Q;
                fmpq_t scale;

                I = flint_malloc(sizeof(fmpz_mpoly_struct *) * n);
                for (i = 0; i < n; i++)
                    I[i] = &ctx->fields[zfield].data.multi.ideal[i];

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
        else
        {
            flint_printf("ca_add: unknown field type\n");
            flint_abort();
        }

        ca_condense_field(res, ctx);

        return;
    }

    {
        ca_t t, u;

        ca_init(t, ctx);
        ca_init(u, ctx);

        ca_merge_fields(t, u, x, y, ctx);

/*
        printf("merged\n");
        ca_print(t, ctx); printf("\n");
        ca_print(u, ctx); printf("\n\n");
*/

        ca_add(res, t, u, ctx);

/*
        printf("added\n");
        ca_print(res, ctx); printf("\n");
*/

        ca_condense_field(res, ctx);

/*
        printf("condensed\n");
        ca_print(res, ctx); printf("\n");
*/

        ca_clear(t, ctx);
        ca_clear(u, ctx);
    }
}

void
ca_sub(ca_t res, const ca_t x, const ca_t y, ca_ctx_t ctx)
{
    ca_t t;
    ca_init(t, ctx);
    ca_neg(t, y, ctx);
    ca_add(res, x, t, ctx);
    ca_clear(t, ctx);
}

void
ca_sub_fmpq(ca_t res, const ca_t x, const fmpq_t y, ca_ctx_t ctx)
{
    fmpq_t t;
    fmpq_init(t);
    fmpq_neg(t, y);
    ca_add_fmpq(res, x, t, ctx);
    fmpq_clear(t);
}

void
ca_sub_fmpz(ca_t res, const ca_t x, const fmpz_t y, ca_ctx_t ctx)
{
    fmpq_t t;
    *fmpq_numref(t) = *y;
    *fmpq_denref(t) = 1;
    ca_sub_fmpq(res, x, t, ctx);
}

void
ca_sub_ui(ca_t res, const ca_t x, ulong y, ca_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init_set_ui(t, y);
    ca_sub_fmpz(res, x, t, ctx);
    fmpz_clear(t);
}

void
ca_sub_si(ca_t res, const ca_t x, slong y, ca_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init_set_si(t, y);
    ca_sub_fmpz(res, x, t, ctx);
    fmpz_clear(t);
}

