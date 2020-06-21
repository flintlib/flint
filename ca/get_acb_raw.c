/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_fmpz_poly.h"
#include "ca.h"

/*
todo: use cached enclosure  &K->data.func.enclosure
*/

void
ca_field_func_get_acb_raw(acb_t res, ca_field_t K, slong prec, ca_ctx_t ctx)
{
    /* todo: verify args_len for functions... */

    switch (K->data.func.func)
    {
        case CA_Pi:
            acb_const_pi(res, prec);
            break;
        case CA_Exp:
            ca_get_acb_raw(res, K->data.func.args, prec, ctx);
            acb_exp(res, res, prec);
            break;
        case CA_Log:
            ca_get_acb_raw(res, K->data.func.args, prec, ctx);
            acb_log(res, res, prec);
            break;

        default:
            flint_printf("ca_field_func_get_acb_raw: unknown function\n");
            flint_abort();
    }
}

void
ca_get_acb_raw(acb_t res, const ca_t x, slong prec, ca_ctx_t ctx)
{
    slong field_index;
    ca_field_type_t type;
    ulong xfield;

    xfield = x->field;

    if (xfield == CA_FIELD_ID_QQ)
    {
        acb_set_fmpq(res, CA_FMPQ(x), prec);
        return;
    }

    if (xfield == CA_FIELD_ID_QQ_I)
    {
        const fmpz *n, *d;

        n = QNF_ELEM_NUMREF(CA_NF_ELEM(x));
        d = QNF_ELEM_DENREF(CA_NF_ELEM(x));

        if (fmpz_is_one(d))
        {
            arb_set_round_fmpz(acb_realref(res), n, prec);
            arb_set_round_fmpz(acb_imagref(res), n + 1, prec);
        }
        else
        {
            arb_fmpz_div_fmpz(acb_realref(res), n, d, prec);
            arb_fmpz_div_fmpz(acb_imagref(res), n + 1, d, prec);
        }
        return;
    }

    if (CA_IS_SPECIAL(x))
    {
        acb_indeterminate(res);
        return;
    }

    field_index = xfield;
    type = ctx->fields[field_index].type;

    if (type == CA_FIELD_TYPE_QQ)
    {
        flint_abort();  /* QQ should be unique and caught above */
    }
    else if (type == CA_FIELD_TYPE_NF)
    {
        if (CA_FIELD_NF(ctx->fields + field_index)->flag & NF_LINEAR)
            flint_abort();

        qqbar_cache_enclosure(CA_FIELD_NF_QQBAR(ctx->fields + field_index), prec);
        qqbar_get_acb(res, CA_FIELD_NF_QQBAR(ctx->fields + field_index), prec);

        if (CA_FIELD_NF(ctx->fields + field_index)->flag & NF_QUADRATIC)
        {
            _arb_fmpz_poly_evaluate_acb(res, QNF_ELEM_NUMREF(CA_NF_ELEM(x)), 2, res, prec);
            acb_div_fmpz(res, res, QNF_ELEM_DENREF(CA_NF_ELEM(x)), prec);
        }
        else
        {
            _arb_fmpz_poly_evaluate_acb(res, NF_ELEM_NUMREF(CA_NF_ELEM(x)), NF_ELEM(CA_NF_ELEM(x))->length, res, prec);
            acb_div_fmpz(res, res, NF_ELEM_DENREF(CA_NF_ELEM(x)), prec);
        }
    }
    else if (type == CA_FIELD_TYPE_FUNC)
    {
        ca_field_func_get_acb_raw(res, ctx->fields + field_index, prec, ctx);
        fmpz_mpoly_q_evaluate_acb(res, CA_MPOLY_Q(x), res, prec, ctx->mctx + 0);
    }
    else if (type == CA_FIELD_TYPE_MULTI)
    {
        acb_ptr v;
        slong i, field_i, n;

        n = (ctx->fields + field_index)-> data.multi.len;
        v = _acb_vec_init(n);

        for (i = 0; i < n; i++)
        {
            field_i =  (ctx->fields + field_index)-> data.multi.ext[i];

            if ((ctx->fields + field_i)->type == CA_FIELD_TYPE_NF)
            {
                qqbar_cache_enclosure(CA_FIELD_NF_QQBAR(ctx->fields + field_i), prec);
                qqbar_get_acb(v + i, CA_FIELD_NF_QQBAR(ctx->fields + field_i), prec);
            }
            else if ((ctx->fields + field_i)->type == CA_FIELD_TYPE_FUNC)
            {
                ca_field_func_get_acb_raw(v + i, ctx->fields + field_i, prec, ctx);
            }
            else
            {
                flint_abort();
            }
        }

        fmpz_mpoly_q_evaluate_acb(res, CA_MPOLY_Q(x), v, prec, CA_FIELD_MCTX(ctx->fields + field_index, ctx));

        _acb_vec_clear(v, n);
    }
}

