/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"

int fmpz_mpoly_evaluate_qqbar(qqbar_t res, const fmpz_mpoly_t pol, qqbar_srcptr x, slong deg_limit, slong bits_limit, const fmpz_mpoly_ctx_t ctx);

int
ca_get_qqbar(qqbar_t res, const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        return 0;
    }
    else if (CA_FIELD_IS_QQ(x, ctx))
    {
        qqbar_set_fmpq(res, CA_FMPQ(x));
        return 1;
    }
    else
    {
        if (CA_FIELD_IS_NF(CA_FIELD(x, ctx)))
        {
            const fmpz * num;
            const fmpz * den;
            slong len;

            if (CA_FIELD_NF(CA_FIELD(x, ctx))->flag & NF_LINEAR)
            {
                num = (fmpz *) LNF_ELEM_NUMREF(CA_NF_ELEM(x));
                den = LNF_ELEM_DENREF(CA_NF_ELEM(x));
                len = 1;
            }
            else if (CA_FIELD_NF(CA_FIELD(x, ctx))->flag & NF_QUADRATIC)
            {
                num = (fmpz *) QNF_ELEM_NUMREF(CA_NF_ELEM(x));
                den = QNF_ELEM_DENREF(CA_NF_ELEM(x));
                len = 2;
            }
            else
            {
                num = (fmpz *) NF_ELEM_NUMREF(CA_NF_ELEM(x));
                den = NF_ELEM_DENREF(CA_NF_ELEM(x));
                len = NF_ELEM(CA_NF_ELEM(x))->length;
            }

            _qqbar_evaluate_fmpq_poly(res, num, den, len, CA_FIELD_NF_QQBAR(CA_FIELD(x, ctx)));

            return 1;
        }
        else
        {
            slong i, len, deg_limit, bits_limit;
            qqbar_ptr xs;
            qqbar_t y;
            int success;

            deg_limit = ctx->options[CA_OPT_QQBAR_DEG_LIMIT];
            bits_limit = 10 * ctx->options[CA_OPT_PREC_LIMIT]; /* xxx */

            len = CA_FIELD_LENGTH(CA_FIELD(x, ctx));

            for (i = 0; i < len; i++)
            {
                if (!CA_EXT_IS_QQBAR(CA_FIELD_GET_EXT(CA_FIELD(x, ctx), i)))
                    return 0;
            }

            success = 0;
            xs = (qqbar_struct *) flint_malloc(sizeof(qqbar_struct) * len);
            qqbar_init(y);

            for (i = 0; i < len; i++)
                xs[i] = *CA_EXT_QQBAR(CA_FIELD_GET_EXT(CA_FIELD(x, ctx), i));

            if (fmpz_mpoly_evaluate_qqbar(y, fmpz_mpoly_q_numref(CA_MPOLY_Q(x)), xs, deg_limit, bits_limit, CA_FIELD_MCTX(CA_FIELD(x, ctx), ctx)))
            {
                if (fmpz_mpoly_evaluate_qqbar(res, fmpz_mpoly_q_denref(CA_MPOLY_Q(x)), xs, deg_limit, bits_limit, CA_FIELD_MCTX(CA_FIELD(x, ctx), ctx)))
                {
                    if (qqbar_binop_within_limits(y, res, deg_limit, bits_limit))
                    {
                        qqbar_div(res, y, res);
                        success = 1;
                    }
                }
            }

            flint_free(xs);
            qqbar_clear(y);

            return success;
        }
    }
}

int
ca_get_fmpq(fmpq_t res, const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        return 0;
    }
    else if (CA_FIELD_IS_QQ(x, ctx))
    {
        fmpq_set(res, CA_FMPQ(x));
        return 1;
    }
    else if (CA_FIELD_IS_NF(CA_FIELD(x, ctx)))
    {
        if (nf_elem_is_rational(CA_NF_ELEM(x), CA_FIELD_NF(CA_FIELD(x, ctx))))
        {
            const fmpz * num;
            const fmpz * den;
            slong len;

            if (CA_FIELD_NF(CA_FIELD(x, ctx))->flag & NF_LINEAR)
            {
                num = (fmpz *) LNF_ELEM_NUMREF(CA_NF_ELEM(x));
                den = LNF_ELEM_DENREF(CA_NF_ELEM(x));
                fmpz_set(fmpq_numref(res), num);
                fmpz_set(fmpq_denref(res), den);
            }
            else if (CA_FIELD_NF(CA_FIELD(x, ctx))->flag & NF_QUADRATIC)
            {
                num = (fmpz *) QNF_ELEM_NUMREF(CA_NF_ELEM(x));
                den = QNF_ELEM_DENREF(CA_NF_ELEM(x));
                fmpz_set(fmpq_numref(res), num);
                fmpz_set(fmpq_denref(res), den);
            }
            else
            {
                num = (fmpz *) NF_ELEM_NUMREF(CA_NF_ELEM(x));
                den = NF_ELEM_DENREF(CA_NF_ELEM(x));
                len = NF_ELEM(CA_NF_ELEM(x))->length;
                if (len == 0)
                {
                    fmpq_zero(res);
                }
                else
                {
                    fmpz_set(fmpq_numref(res), num);
                    fmpz_set(fmpq_denref(res), den);
                }
            }

            return 1;
        }
        else
        {
            return 0;
        }
    }
    else
    {
        /* todo: exclude complex numbers, obviously irrational numbers
           before evaluating */
        int success;
        qqbar_t t;
        qqbar_init(t);
        success = ca_get_qqbar(t, x, ctx);
        if (success && qqbar_is_rational(t))
        {
            fmpz_neg(fmpq_numref(res), QQBAR_COEFFS(t));
            fmpz_set(fmpq_denref(res), QQBAR_COEFFS(t) + 1);
            success = 1;
        }
        qqbar_clear(t);
        return success;
    }
}

int
ca_get_fmpz(fmpz_t res, const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        return 0;
    }
    else if (CA_FIELD_IS_QQ(x, ctx))
    {
        if (fmpz_is_one(fmpq_denref(CA_FMPQ(x))))
        {
            fmpz_set(res, fmpq_numref(CA_FMPQ(x)));
            return 1;
        }

        return 0;
    }
    else if (CA_FIELD_IS_NF(CA_FIELD(x, ctx)))
    {
        if (nf_elem_is_integer(CA_NF_ELEM(x), CA_FIELD_NF(CA_FIELD(x, ctx))))
        {
            const fmpz * num;
            slong len;

            if (CA_FIELD_NF(CA_FIELD(x, ctx))->flag & NF_LINEAR)
            {
                num = (fmpz *) LNF_ELEM_NUMREF(CA_NF_ELEM(x));
                fmpz_set(res, num);
            }
            else if (CA_FIELD_NF(CA_FIELD(x, ctx))->flag & NF_QUADRATIC)
            {
                num = (fmpz *) QNF_ELEM_NUMREF(CA_NF_ELEM(x));
                fmpz_set(res, num);
            }
            else
            {
                num = (fmpz *) NF_ELEM_NUMREF(CA_NF_ELEM(x));
                len = NF_ELEM(CA_NF_ELEM(x))->length;
                if (len == 0)
                    fmpz_zero(res);
                else
                    fmpz_set(res, num);
            }

            return 1;
        }
        else
        {
            return 0;
        }
    }
    else
    {
        /* todo: exclude (numerically) obvious non-integers before evaluating */
        int success;
        qqbar_t t;
        qqbar_init(t);
        success = ca_get_qqbar(t, x, ctx);
        if (success && qqbar_is_integer(t))
        {
            fmpz_neg(res, QQBAR_COEFFS(t));
            success = 1;
        }
        qqbar_clear(t);
        return success;
    }
}

