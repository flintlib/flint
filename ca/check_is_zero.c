/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"

/* stupid algorithm, just to have something working... */
int
fmpz_mpoly_evaluate_qqbar(qqbar_t res, const fmpz_mpoly_t pol, qqbar_srcptr x, slong deg_limit, slong bits_limit, const fmpz_mpoly_ctx_t ctx)
{
    slong i, j, len, nvars;
    qqbar_t s, t, u;
    ulong * exp;
    int success;

    len = fmpz_mpoly_length(pol, ctx);

    if (len == 0)
    {
        qqbar_zero(res);
        return 1;
    }

    if (len == 1 && fmpz_mpoly_is_fmpz(pol, ctx))
    {
        qqbar_set_fmpz(res, pol->coeffs);
        return 1;
    }

    success = 0;

    nvars = ctx->minfo->nvars;
    exp = flint_malloc(sizeof(ulong) * nvars);

    qqbar_init(s);
    qqbar_init(t);
    qqbar_init(u);

    for (i = 0; i < len; i++)
    {
        fmpz_mpoly_get_term_exp_ui(exp, pol, i, ctx);

        qqbar_one(t);

        for (j = 0; j < nvars; j++)
        {
            if (exp[j] == 1)
            {
                if ((double) qqbar_degree(t) * (double) qqbar_degree(x + j) > (double) deg_limit ||
                    (double) qqbar_height_bits(t) + (double) qqbar_height_bits(x + j) > (double) bits_limit)
                    goto cleanup;

                qqbar_mul(t, t, x + j);
            }
            else if (exp[j] >= 2)
            {
                if ((double) qqbar_height_bits(x + j) * (double) exp[j] > bits_limit)
                    goto cleanup;

                qqbar_pow_ui(u, x + j, exp[j]);

                if ((double) qqbar_degree(t) * (double) qqbar_degree(u) > (double) deg_limit ||
                    (double) qqbar_height_bits(t) + (double) qqbar_height_bits(u) > (double) bits_limit)
                    goto cleanup;

                qqbar_mul(t, t, u);
            }
        }

        qqbar_mul_fmpz(t, t, pol->coeffs + i);

        if ((double) qqbar_degree(s) * (double) qqbar_degree(t) > (double) deg_limit ||
            (double) qqbar_height_bits(s) + (double) qqbar_height_bits(t) > (double) bits_limit)
            goto cleanup;

        qqbar_add(s, s, t);
    }

    success = 1;
    qqbar_swap(res, s);

cleanup:
    flint_free(exp);

    qqbar_clear(s);
    qqbar_clear(t);
    qqbar_clear(u);

    return success;
}

truth_t
_ca_check_is_zero_qqbar(const ca_t x, ca_ctx_t ctx)
{
    slong i, len, deg_limit, bits_limit;
    qqbar_ptr xs;
    qqbar_t y;
    truth_t res;

    deg_limit = ctx->options[CA_OPT_QQBAR_DEG_LIMIT];
    bits_limit = 10 * ctx->options[CA_OPT_PREC_LIMIT]; /* xxx */

    len = CA_FIELD_LENGTH(ctx->fields + x->field);

    for (i = 0; i < len; i++)
    {
        /* todo: could allow symbolic functions that allow evaluation to qqbar */
        if (!CA_EXT_IS_QQBAR(CA_FIELD_GET_EXT(ctx->fields + x->field, i)))
            return T_UNKNOWN;
    }

    res = T_UNKNOWN;
    xs = (qqbar_struct *) flint_malloc(sizeof(qqbar_struct) * len);
    qqbar_init(y);

    for (i = 0; i < len; i++)
        xs[i] = *CA_EXT_QQBAR(CA_FIELD_GET_EXT(ctx->fields + x->field, i));

    if (fmpz_mpoly_evaluate_qqbar(y, fmpz_mpoly_q_numref(CA_MPOLY_Q(x)), xs, deg_limit, bits_limit, CA_FIELD_MCTX(ctx->fields + x->field, ctx)))
    {
        res = qqbar_is_zero(y) ? T_TRUE : T_FALSE;
    }

    flint_free(xs);
    qqbar_clear(y);

    return res;
}

truth_t
ca_check_is_zero(const ca_t x, ca_ctx_t ctx)
{
    acb_t v;
    truth_t res;
    slong prec, prec_limit;

    if (CA_IS_SPECIAL(x))
    {
        if (ca_is_unknown(x, ctx))
            return T_UNKNOWN;

        return T_FALSE;
    }

    if (x->field == CA_FIELD_ID_QQ)
    {
        if (fmpq_is_zero(CA_FMPQ(x)))
            return T_TRUE;
        else
            return T_FALSE;
    }

    if (x->field == CA_FIELD_ID_QQ_I)
    {
        const fmpz *n;

        n = QNF_ELEM_NUMREF(CA_NF_ELEM(x));

        if (fmpz_is_zero(n) && fmpz_is_zero(n + 1))
            return T_TRUE;

        return T_FALSE;
    }

    if (CA_FIELD_IS_NF(ctx->fields + x->field))
    {
        if (nf_elem_is_zero(CA_NF_ELEM(x), CA_FIELD_NF(ctx->fields + x->field)))
            return T_TRUE;
        else
            return T_FALSE;
    }

    res = T_UNKNOWN;

    /* todo: in the following, we should extract the numerator; the denominator is irrelevant */

    acb_init(v);

    prec_limit = ctx->options[CA_OPT_PREC_LIMIT];
    prec_limit = FLINT_MAX(prec_limit, 64);

    for (prec = 64; (prec <= prec_limit) && (res == T_UNKNOWN); prec *= 2)
    {
        ca_get_acb_raw(v, x, prec, ctx);

        if (!acb_contains_zero(v))
        {
            res = T_FALSE;
        }

        /* try qqbar computation */
        /* todo: precision to do this should depend on complexity of the polynomials, degree of the elements... */
        if (prec == 64)
        {
            res = _ca_check_is_zero_qqbar(x, ctx);
        }
    }

    acb_clear(v);

    return res;
}

