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
ca_print(const ca_t x, const ca_ctx_t ctx)
{
    slong field;

    if (CA_IS_SPECIAL(x))
    {
        if (x->field & CA_UNDEFINED)
        {
            flint_printf("Undefined");
        }
        else if (x->field & CA_UNKNOWN)
        {
            flint_printf("Unknown");
        }
        else if (x->field & CA_UNSIGNED_INF)
        {
            flint_printf("UnsignedInfinity");
        }
        else  if (x->field & CA_SIGNED_INF)
        {
            ca_t sgn;
            sgn->field = x->field & ~CA_SPECIAL;
            sgn->elem = x->elem;

            if (CA_FIELD_IS_QQ(sgn, ctx))
            {
                if (fmpq_sgn(CA_FMPQ(sgn)) > 0)
                    flint_printf("+Infinity");
                else
                    flint_printf("-Infinity");
            }
            else if (CA_FIELD_IS_QQ_I(sgn, ctx))
            {
                if (fmpz_sgn(QNF_ELEM_NUMREF(CA_NF_ELEM(sgn)) + 1) > 0)
                    flint_printf("+I * Infinity");
                else
                    flint_printf("-I * Infinity");
            }
            else
            {
                flint_printf("Infinity * (");
                ca_print(sgn, ctx);
                flint_printf(")");
            }
        }

        return;
    }
    else
    {
        field = x->field;
    }

    if (field == CA_FIELD_ID_QQ)
    {
        fmpq_print(CA_FMPQ(x));
    }
    else if (CA_FIELD_IS_NF(ctx->fields + field))
    {
        nf_elem_print_pretty(CA_NF_ELEM(x), CA_FIELD_NF(ctx->fields + field), "x");
    }
    else
    {
        /* todo: could use depth to select different characters */
        if (CA_FIELD_LENGTH(ctx->fields + field) == 1)
        {
            const char * xx = "x1";

            fmpz_mpoly_q_print_pretty(CA_MPOLY_Q(x), &xx, CA_FIELD_MCTX(ctx->fields + field, ctx));
        }
        else
        {
            fmpz_mpoly_q_print_pretty(CA_MPOLY_Q(x), NULL, CA_FIELD_MCTX(ctx->fields + field, ctx));
        }
    }

    flint_printf("  in  ");
    ca_field_print(ctx->fields + field, ctx);
}

